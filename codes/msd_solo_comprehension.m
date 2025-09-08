% Paramètres
Interval = 15;
MaxDisp = 10;
Radii = linspace(0, 1500, 11);  

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

    % Détermination du max_lag
    max_frames = max(AllTraj(:, 3)) / Interval;
    max_lag = floor(max_frames);

    % Choix d'une particule aléatoire parmi les trajectoires complètes
    rng('shuffle');
    selected_idx = randi(length(original_trajectories));
    % selected_idx = 1206;
    selected_original_traj = original_trajectories{selected_idx};

    % Affichage des trajectoires corrigées
    figure(1); clf;
    % plot(selected_original_traj(:,1), selected_original_traj(:,2), 'k--', 'DisplayName', 'Originale');
    hold on;
    colors = lines(length(Radii));

    figure(2); clf; hold on;  % Pour MSD
    for r = 1:length(Radii)
        Radius = Radii(r);

        % Correction avec affichage d'infos
        [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = ...
            correct_local_drift_individual_verbose(selected_original_traj, AllTraj, Radius);

        % Affichage des infos dans le terminal
        fprintf('--- Rayon = %.0f ---\n', Radius);
        fprintf('Voisins détectés (total) : %d\n', total_neighbors);
        fprintf('Voisins moyens par image : %.2f\n', avg_neighbors_per_frame);
        fprintf('Vitesse moyenne locale détectée : %.4f px/s\n\n', avg_speed);

        % Trajectoires corrigées
        figure(1);
        plot(corrected_traj(:,1), corrected_traj(:,2), '-', 'Color', colors(r,:), ...
             'DisplayName', sprintf('Corrigée R=%.0f', Radius));

        % Calcul et tracé du MSD avec étiquette sur la courbe
        [msd_single, time_lags] = compute_msd({corrected_traj}, max_lag, Interval);
        figure(2);
        plot(time_lags, msd_single, 'Color', colors(r,:), 'DisplayName', sprintf('R=%.0f', Radius));
        
        % Ajout d'une étiquette sur la courbe MSD
        valid_indices = find(~isnan(msd_single)); % Tous les indices non-NaN
        if ~isempty(valid_indices)
            target_idx = valid_indices(max(1, round(0.90 * length(valid_indices)))); 
            text(time_lags(target_idx), msd_single(target_idx), ...
                sprintf('R=%.0f, %.2f voisins/image', Radius, avg_neighbors_per_frame), ...
                'Color', colors(r,:), 'FontSize', 10, 'VerticalAlignment', 'bottom', ...
                'HorizontalAlignment', 'left');
        end
    end

    % Finalisation figures
    figure(1);
    title(sprintf('Trajectoires corrigées - Particule #%d', selected_idx));
    xlabel('X (px)'); ylabel('Y (px)'); legend show; axis equal; grid on;

    figure(2);
    title(sprintf('MSD - Particule #%d', selected_idx));
    xlabel('Time Lag (s)'); ylabel('MSD (px²)'); legend show; grid on;
    set(gcf, 'Position', [100, 100, 1200, 700]);
end

% --- Fonction de correction locale avec retour d'infos
function [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = correct_local_drift_individual_verbose(traj, AllTraj, Radius)
    % Initialisation
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    vx_all = []; vy_all = [];
    total_neighbors = 0;

    if Radius == 0
        % Cas spécial : calculer la vitesse moyenne de la particule seule
        vx = []; vy = [];
        for i = 1:n_points-1
            dt = traj(i+1, 3) - traj(i, 3);
            if dt > 0
                vx(end+1) = (traj(i+1, 1) - traj(i, 1)) / dt;
                vy(end+1) = (traj(i+1, 2) - traj(i, 2)) / dt;
            end
        end
        
        if ~isempty(vx)
            ux_mean = mean(vx);
            uy_mean = mean(vy);
            for i = 1:n_points
                t_i = traj(i, 3);
                x_corrected(i) = traj(i, 1) - ux_mean * t_i;
                y_corrected(i) = traj(i, 2) - uy_mean * t_i;
            end
            vx_all = ux_mean;
            vy_all = uy_mean;
            total_neighbors = n_points^2 - 1; % Nombre de segments utilisés
        else
            x_corrected = traj(:, 1);
            y_corrected = traj(:, 2);
        end
    else
        % Cas standard : utiliser les voisins dans AllTraj
        for i = 1:n_points
            x_i = traj(i, 1);
            y_i = traj(i, 2);
            t_i = traj(i, 3);

            % Recherche des voisins
            distances = sqrt((AllTraj(:, 1) - x_i).^2 + (AllTraj(:, 2) - y_i).^2);           
            neighbors_idx = find(distances <= Radius);
            total_neighbors = total_neighbors + length(neighbors_idx);

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
                vx_all(end+1) = ux_local;
                vy_all(end+1) = uy_local;
            else
                x_corrected(i) = x_i;
                y_corrected(i) = y_i;
            end
        end
        avg_neighbors_per_frame = (total_neighbors / n_points^2);
    end

    % Construction de la trajectoire corrigée
    corrected_traj = [x_corrected, y_corrected, traj(:, 3)];

    % Calcul de la vitesse moyenne globale
    if ~isempty(vx_all)
        speed_norms = sqrt(vx_all.^2 + vy_all.^2);
        avg_speed = mean(speed_norms);
    else
        avg_speed = 0;
    end
    % Calcul du nombre moyen de voisins par image
    avg_neighbors_per_frame = (total_neighbors / n_points^2);
end

% --- Fonction MSD
function [msd, time_lags] = compute_msd(trajectories, max_lag, interval)
    n_trajectories = length(trajectories);
    time_lags = (1:max_lag) * interval;
    msd = zeros(1, max_lag);
    n_counts = zeros(1, max_lag);

    for traj_idx = 1:n_trajectories
        traj = trajectories{traj_idx};
        n_points = size(traj, 1);
        if n_points < 2 || any(any(~isfinite(traj(:, 1:2))))
            continue;
        end
        for lag = 1:min(max_lag, n_points - 1)
            for i = 1:(n_points - lag)
                dx = traj(i+lag, 1) - traj(i, 1);
                dy = traj(i+lag, 2) - traj(i, 2);
                if isfinite(dx) && isfinite(dy)
                    squared_displacement = dx^2 + dy^2;
                    msd(lag) = msd(lag) + squared_displacement;
                    n_counts(lag) = n_counts(lag) + 1;
                end
            end
        end
    end
    for lag = 1:max_lag
        if n_counts(lag) > 0
            msd(lag) = msd(lag) / n_counts(lag);
        else
            msd(lag) = NaN;
            fprintf('Avertissement: Aucun point valide trouvé pour le lag %d\n', lag);
        end
    end
end