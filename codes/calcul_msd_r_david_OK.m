% Paramètres
Interval = 15;
MaxDisp = 10;
pixel_to_um = 1 / 1.8; % Conversion factor: 1 px = 0.5556 µm
pixel_to_um2 = pixel_to_um^2; % Conversion factor for area: 1 px² = 0.3086 µm²
Radii = linspace(0, 580, 100) * pixel_to_um; % Convert radii to µm

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
    NrOfFrames = max(AllTraj(:, 3)) / Interval; % Nombre total d'images
    
    % Initialisation pour stocker les MSD et voisins pour toutes les particules
    all_last_msd_values = zeros(NrOfTrajs, length(Radii));
    all_avg_neighbors_values = zeros(NrOfTrajs, length(Radii));
    
    % Boucle sur toutes les particules
    fprintf('Analyse du fichier %s\n', filename);
    fprintf('Nombre total de particules : %d\n', NrOfTrajs);
    
    tic; 
    for n = 1:NrOfTrajs
        selected_original_traj = original_trajectories{n};
        
        % Affichage de l'avancement toutes les 10 particules
        if mod(n, 10) == 1 || n == NrOfTrajs
            fprintf('Traitement de la particule %d/%d (%.1f%%)\n', ...
                n, NrOfTrajs, (n/NrOfTrajs)*100);
        end
        
        for r = 1:length(Radii)
            Radius = Radii(r);
            
            % Correction avec affichage d'infos
            [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = ...
                correct_local_drift_individual_verbose(selected_original_traj, original_trajectories, AllTraj, Radius, Interval);
            
            % Calcul du MSD
            [msd_single, ~] = compute_msd({corrected_traj}, max_lag, Interval);
            
            % Extraction de la dernière valeur non-NaN du MSD, converti en µm²
            valid_msd = msd_single(~isnan(msd_single));
            if ~isempty(valid_msd)
                all_last_msd_values(n, r) = valid_msd(end) * pixel_to_um2; % Convert to µm²
            else
                all_last_msd_values(n, r) = NaN;
            end
            
            all_avg_neighbors_values(n, r) = avg_neighbors_per_frame;
        end
    end
    toc; 
    fprintf('Temps d''exécution pour %s : %.2f secondes\n', filename, toc);
    
    % Calcul des moyennes et écarts-types
    mean_msd = zeros(1, length(Radii));
    std_msd = zeros(1, length(Radii));
    mean_neighbors = zeros(1, length(Radii));
    std_neighbors = zeros(1, length(Radii));
    
    for r = 1:length(Radii)
        valid_msd = all_last_msd_values(:, r);
        valid_msd = valid_msd(~isnan(valid_msd));
        if ~isempty(valid_msd)
            mean_msd(r) = mean(valid_msd);
            if length(valid_msd) > 1
                std_msd(r) = sqrt(sum((valid_msd - mean_msd(r)).^2) / (length(valid_msd) - 1));
            else
                std_msd(r) = 0;
            end
        else
            mean_msd(r) = NaN;
            std_msd(r) = NaN;
        end
        
        % Moyenne pour voisins
        valid_neighbors = all_avg_neighbors_values(:, r);
        valid_neighbors = valid_neighbors(~isnan(valid_neighbors));
        if ~isempty(valid_neighbors)
            mean_neighbors(r) = mean(valid_neighbors);
            if length(valid_neighbors) > 1
                std_neighbors(r) = sqrt(sum((valid_neighbors - mean_neighbors(r)).^2) / (length(valid_neighbors) - 1));
            else
                std_neighbors(r) = 0;
            end
        else
            mean_neighbors(r) = NaN;
            std_neighbors(r) = NaN;
        end
    end
    
    % Affichage des informations pertinentes dans la console
    fprintf('Résultats pour le fichier %s :\n', filename);
    fprintf('MSD moyen (dernier lag) pour Rayon = %.0f µm : %.2f ± %.2f µm²\n', ...
        Radii(end), mean_msd(end), std_msd(end));
    fprintf('Voisins moyens pour Rayon = %.0f µm : %.2f ± %.2f\n\n', ...
        Radii(end), mean_neighbors(end), std_neighbors(end));
    
    % Affichage des résultats avec barres d'erreur
    figure; clf;

    % Couleurs personnalisées
    msd_color = [0.2 0.6 0.8];       % Bleu
    neighbor_color = [0.8 0.4 0.2];  % Orange

    % Tracé MSD avec zone d'erreur (yyaxis gauche)
    yyaxis left;
    hold on;
    fill([Radii fliplr(Radii)], ...
        [mean_msd + std_msd/sqrt(NrOfTrajs), fliplr(mean_msd - std_msd/sqrt(NrOfTrajs))], ...
        msd_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(Radii, mean_msd, '-o', 'Color', msd_color, 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('MSD moyen (µm²)');

    % Tracé voisins avec zone d'erreur (yyaxis droit)
    yyaxis right;
    hold on;
    fill([Radii fliplr(Radii)], ...
        [mean_neighbors + std_neighbors/sqrt(NrOfFrames), fliplr(mean_neighbors - std_neighbors/sqrt(NrOfFrames))], ...
        neighbor_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(Radii, mean_neighbors, '--s', 'Color', neighbor_color, 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('Voisins moyens');

    % Titres et axes
    title(sprintf('MSD moyen et Voisins vs Rayon - %s', filename));
    xlabel('Rayon (µm)');
    grid on;
    legend({'Zone erreur MSD', 'MSD moyen', 'Zone erreur voisins', 'Voisins moyens'}, 'Location', 'best');
end

% --- Fonction de correction locale avec retour d'infos
function [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = ...
    correct_local_drift_individual_verbose(traj, original_trajectories, AllTraj, Radius, Interval)
    
    pixel_to_um = 1 / 1.8; % Conversion factor for velocity
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    vx_all = []; vy_all = [];
    neighbors_per_frame = zeros(n_points, 1); % Stocke le nombre de voisins par image

    ntraj = max(AllTraj(:, 4));
    corr_x = 0;
    corr_y = 0;

    x_corrected(1) = traj(1,1);
    y_corrected(1) = traj(1,2);

    for i = 2:n_points
        x_i = traj(i, 1);
        y_i = traj(i, 2);
        t_i = traj(i, 3);

        iv = 0; % Nombre de voisins
        vx = 0;
        vy = 0;

        % Recherche des voisins
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
        neighbors_per_frame(i) = iv; % Stocker le nombre de voisins pour cette image

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

    % Construction de la trajectoire corrigée
    corrected_traj = [x_corrected, y_corrected, traj(:, 3)];

    % Calcul de la vitesse moyenne globale
    if ~isempty(vx_all)
        speed_norms = sqrt(vx_all.^2 + vy_all.^2);
        avg_speed = mean(speed_norms);
    else
        avg_speed = 0;
    end
    
    % Calcul du nombre total et moyen de voisins
    total_neighbors = sum(neighbors_per_frame);
    if n_points > 1
        avg_neighbors_per_frame = mean(neighbors_per_frame(2:end)); % Exclure le premier point
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