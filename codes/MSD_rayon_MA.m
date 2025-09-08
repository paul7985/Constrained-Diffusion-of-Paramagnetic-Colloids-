Interval = 15;
MaxDisp = 10;
Radii = linspace(40, 80, 10);

% Sélection des fichiers CSV
[filenames, pathname] = uigetfile('*.csv', 'Select Results.csv files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames};
end

for f = 1:length(filenames)
    filename = filenames{f}; 
    fullpath = fullfile(pathname, filename);
    
    % Lecture du fichier CSV
    A = importdata(fullpath, ',', 1);
    i_Frame = find(strcmp('Frame', A.colheaders));
    i_X = find(strcmp('X', A.colheaders));
    i_Y = find(strcmp('Y', A.colheaders));
    
    time = (A.data(:, i_Frame) - 1) * Interval;
    InputForTrack = [A.data(:, i_X), A.data(:, i_Y), time];
    InputForTrack(any(isnan(InputForTrack), 2), :) = [];  % Suppression des NaN
    
    % Tracking
    AllTraj = track(InputForTrack, MaxDisp);
    NrOfTrajs = max(AllTraj(:, 4));
    
    % Extraction des trajectoires originales
    original_trajectories = cell(NrOfTrajs, 1);
    for n = 1:NrOfTrajs
        idx = AllTraj(:, 4) == n;
        original_trajectories{n} = AllTraj(idx, [1 2 3]);
    end
    
    % Préparation de la figure
    figure;
    hold on;
    
    % Détermination du max_lag
    max_frames = max(AllTraj(:, 3)) / Interval;
    max_lag = floor(max_frames);
    
    % Boucle sur les différents rayons
    for r = 1:length(Radii)
        Radius = Radii(r)
        corrected_trajectories = cell(NrOfTrajs, 1);
        
        for n = 1:NrOfTrajs
            corrected_trajectories{n} = correct_local_drift_individual(original_trajectories{n}, AllTraj, Radius);
        end
        
        % Calcul et tracé du MSD
        [msd_local, time_lags] = compute_msd(corrected_trajectories, max_lag, Interval);
        plot(time_lags, msd_local, 'DisplayName', sprintf('Radius = %.0f', Radius));
    end
    
    % Finalisation du tracé
    title(['MSD - Local Correction: ' filename], 'Interpreter', 'none');
    xlabel('Time Lag (s)'); ylabel('MSD'); grid on;
    legend show;
    set(gcf, 'Position', [100, 100, 1200, 700]);
    hold off;
end

function corrected_traj = correct_local_drift_individual(traj, AllTraj, Radius)
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    
    for i = 1:n_points
        x_i = traj(i, 1);
        y_i = traj(i, 2);
        t_i = traj(i, 3);
        
        distances = sqrt((AllTraj(:, 1) - x_i).^2 + (AllTraj(:, 2) - y_i).^2);
        neighbors_idx = find(distances < Radius);
        
        vx = []; vy = [];
        for j = 1:length(neighbors_idx)-1
            idx1 = neighbors_idx(j);
            idx2 = neighbors_idx(j+1);
            dt = AllTraj(idx2, 3) - AllTraj(idx1, 3);
            if dt > 0
                vx(end+1) = (AllTraj(idx2, 1) - AllTraj(idx1, 1)) / dt;
                vy(end+1) = (AllTraj(idx2, 2) - AllTraj(idx1, 2)) / dt;
            end
        end
        
        if ~isempty(vx)
            ux_local = mean(vx);
            uy_local = mean(vy);
            x_corrected(i) = x_i - ux_local * t_i;
            y_corrected(i) = y_i - uy_local * t_i;
        else
            x_corrected(i) = x_i;
            y_corrected(i) = y_i;
        end
    end
    corrected_traj = [x_corrected, y_corrected, traj(:, 3)];
end

function [msd, time_lags] = compute_msd(trajectories, max_lag, interval)
    % Calcule le MSD (Mean Square Displacement) pour un ensemble de trajectoires
    % trajectories: cellule contenant des matrices [x, y, t] pour chaque trajectoire
    % max_lag: nombre maximal de décalages temporels à considérer
    % interval: intervalle de temps entre les frames
    
    n_trajectories = length(trajectories);
    time_lags = (1:max_lag) * interval;
    msd = zeros(1, max_lag);
    n_counts = zeros(1, max_lag);
    
    % Pour chaque trajectoire
    for traj_idx = 1:n_trajectories
        traj = trajectories{traj_idx};
        n_points = size(traj, 1);
        
        if n_points < 2
            continue;  % Ignorer les trajectoires trop courtes
        end
        
        % Vérifier que les coordonnées sont valides
        if any(any(~isfinite(traj(:, 1:2))))
            fprintf('Trajectoire %d contient des valeurs non finies\n', traj_idx);
            continue;
        end
        
        % Pour chaque décalage temporel
        for lag = 1:min(max_lag, n_points-1)
            % Pour chaque point dans la trajectoire
            for i = 1:(n_points - lag)
                % Calcul du déplacement carré
                dx = traj(i+lag, 1) - traj(i, 1);
                dy = traj(i+lag, 2) - traj(i, 2);
                
                % Vérifier que le déplacement est valide
                if isfinite(dx) && isfinite(dy)
                    squared_displacement = dx^2 + dy^2;
                    
                    % Accumulation pour le calcul de la moyenne
                    msd(lag) = msd(lag) + squared_displacement;
                    n_counts(lag) = n_counts(lag) + 1;
                end
            end
        end
    end
    
    % Calcul de la moyenne
    for lag = 1:max_lag
        if n_counts(lag) > 0
            msd(lag) = msd(lag) / n_counts(lag);
        else
            msd(lag) = NaN;
            fprintf('Avertissement: Aucun point valide trouvé pour le lag %d\n', lag);
        end
    end
end
