% Paramètres
Interval = 15;
MaxDisp = 10;
pixel_to_um = 1 / 1.8; % Conversion factor: 1 px = 0.5556 µm
pixel_to_um2 = pixel_to_um^2; % Conversion factor for area: 1 px² = 0.3086 µm²
Radii = linspace(70, 580, 5) * pixel_to_um; % Convert radii to µm

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
    InputForTrack = [A.data(:, i_X) * pixel_to_um, A.data(:, i_Y) * pixel_to_um, time]; % Convert X, Y to µm
    InputForTrack(any(isnan(InputForTrack), 2), :) = [];  % Suppression des NaN
    
    % Tracking
    AllTraj = track(InputForTrack, MaxDisp * pixel_to_um); % Convert MaxDisp to µm
    
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
    
    % Initialisation pour stocker les MSD pour toutes les particules et tous les rayons
    all_msd_values = cell(NrOfTrajs, length(Radii));
    
    % Boucle sur toutes les particules
    fprintf('Analyse du fichier %s\n', filename);
    fprintf('Nombre total de particules : %d\n', NrOfTrajs);
    
    tic; 
    for n = 1:NrOfTrajs
        selected_original_traj = original_trajectories{n};
        
        if mod(n, 10) == 1 || n == NrOfTrajs
            fprintf('Traitement de la particule %d/%d (%.1f%%)\n', ...
                n, NrOfTrajs, (n/NrOfTrajs)*100);
        end
        
        for r = 1:length(Radii)
            Radius = Radii(r);
            
            % Correction avec affichage d'infos
            [corrected_traj, ~, ~, ~] = ...
                correct_local_drift_individual_verbose(selected_original_traj, original_trajectories, AllTraj, Radius, Interval);
            
            % Calcul du MSD
            [msd_single, ~] = compute_msd({corrected_traj}, max_lag, Interval);
            
            % Stockage du MSD complet, converti en µm²
            all_msd_values{n, r} = msd_single * pixel_to_um2;
        end
    end
    toc; 
    fprintf('Temps d''exécution pour %s : %.2f secondes\n', filename, toc);
    
    % Calcul de la moyenne des MSD pour chaque rayon
    mean_msd = zeros(max_lag, length(Radii));
    for r = 1:length(Radii)
        msd_sum = zeros(max_lag, 1);
        n_counts = zeros(max_lag, 1);
        for n = 1:NrOfTrajs
            msd_single = all_msd_values{n, r};
            for lag = 1:max_lag
                if ~isnan(msd_single(lag))
                    msd_sum(lag) = msd_sum(lag) + msd_single(lag);
                    n_counts(lag) = n_counts(lag) + 1;
                end
            end
        end
        for lag = 1:max_lag
            if n_counts(lag) > 0
                mean_msd(lag, r) = msd_sum(lag) / n_counts(lag);
            else
                mean_msd(lag, r) = NaN;
            end
        end
    end
    
    % Affichage des courbes MSD en fonction de tau pour chaque rayon
    figure; clf;
    hold on;
    colors = lines(length(Radii));
    time_lags = (1:max_lag) * Interval;
    
    for r = 1:length(Radii)
        valid_msd = mean_msd(:, r);
        valid_lags = time_lags(~isnan(valid_msd));
        valid_msd = valid_msd(~isnan(valid_msd));
        if ~isempty(valid_msd)
            plot(valid_lags, valid_msd, '-o', 'Color', colors(r, :), ...
                'LineWidth', 1.5, 'MarkerSize', 4, ...
                'DisplayName', sprintf('Rayon = %.0f µm', Radii(r)));
        end
    end
    
    % Mise en forme du graphique
    xlabel('Time lag \tau (s)');
    ylabel('MSD moyen (µm²)');
    title(sprintf('MSD moyen vs Time Lag - %s', filename));
    grid on;
    legend('show', 'Location', 'best');
    set(gca, 'YScale', 'log');
    hold off;
end

% --- Fonction de correction locale avec retour d'infos
function [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = ...
    correct_local_drift_individual_verbose(traj, original_trajectories, AllTraj, Radius, Interval)
    
    pixel_to_um = 1 / 1.8; % Conversion factor for velocity
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    vx_all = []; vy_all = [];
    neighbors_per_frame = zeros(n_points, 1);

    ntraj = max(AllTraj(:, 4));
    corr_x = 0;
    corr_y = 0;

    x_corrected(1) = traj(1,1);
    y_corrected(1) = traj(1,2);

    for i = 2:n_points
        x_i = traj(i, 1);
        y_i = traj(i, 2);
        t_i = traj(i, 3);

        iv = 0;
        vx = 0;
        vy = 0;

        for itraj = 1:ntraj
            traj2 = original_trajectories{itraj};
            i2 = find(traj2(:,3)==t_i);
            if ~isempty(i2) && i2 > 1
                if sqrt((traj2(i2, 1) - x_i).^2 + (traj2(i2, 2) - y_i).^2) <= Radius
                    iv = iv + 1;
                    vx = vx + (traj2(i2,1) - traj2(i2-1,1));
                    vy = vy + (traj2(i2,2) - traj2(i2-1,2));
                end
            end
        end
        neighbors_per_frame(i) = iv;

        if iv > 0
            vx = vx / (iv * Interval) * pixel_to_um; % Convert velocity to µm/s
            vy = vy / (iv * Interval) * pixel_to_um;
            vx_all = [vx_all; vx];
            vy_all = [vy_all; vy];
        end

        corr_x = corr_x - vx * Interval;
        corr_y = corr_y - vy * Interval;

        x_corrected(i) = x_i + corr_x;
        y_corrected(i) = y_i + corr_y;
    end

    corrected_traj = [x_corrected, y_corrected, traj(:, 3)];

    if ~isempty(vx_all)
        speed_norms = sqrt(vx_all.^2 + vy_all.^2);
        avg_speed = mean(speed_norms);
    else
        avg_speed = 0;
    end
    
    total_neighbors = sum(neighbors_per_frame);
    if n_points > 1
        avg_neighbors_per_frame = mean(neighbors_per_frame(2:end));
    else
        avg_neighbors_per_frame = 0;
    end
end

% --- Fonction MSD
function [msd, time_lags] = compute_msd(trajectories, max_lag, interval)
    pixel_to_um2 = (1 / 1.8)^2; % Conversion factor for MSD
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
            msd(lag) = msd(lag) / n_counts(lag) * pixel_to_um2; % Convert MSD to µm²
        else
            msd(lag) = NaN;
        end
    end
end