% Paramètres
Interval = 15;
MaxDisp = 10;
Radii = linspace(0, 1000, 5);  % Test avec 5 rayons, à ajuster à 200 après validation

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
    
    % Initialisation pour stocker les MSD et voisins pour toutes les particules
    all_last_msd_values = zeros(NrOfTrajs, length(Radii));
    all_avg_neighbors_values = zeros(NrOfTrajs, length(Radii));
    
    % Boucle sur toutes les particules
    fprintf('Analyse du fichier %s\n', filename);
    fprintf('Nombre total de particules : %d\n', NrOfTrajs);
    
    tic; % Début du chronométrage
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
                correct_local_drift_individual_verbose(selected_original_traj, AllTraj, Radius, Interval);
            
            % Calcul du MSD
            [msd_single, ~] = compute_msd({corrected_traj}, max_lag, Interval);
            
            % Extraction de la dernière valeur non-NaN du MSD
            valid_msd = msd_single(~isnan(msd_single));
            if ~isempty(valid_msd)
                all_last_msd_values(n, r) = valid_msd(end); % Dernière valeur non-NaN
            else
                all_last_msd_values(n, r) = NaN; % Si aucun MSD valide, utiliser NaN
            end
            
            all_avg_neighbors_values(n, r) = avg_neighbors_per_frame;
        end
    end
    toc; % Fin du chronométrage
    fprintf('Temps d''exécution pour %s : %.2f secondes\n', filename, toc);
    
    % Calcul des moyennes et écarts-types en ignorant les NaN
    mean_msd = zeros(1, length(Radii));
    std_msd = zeros(1, length(Radii));
    mean_neighbors = zeros(1, length(Radii));
    std_neighbors = zeros(1, length(Radii));
    
    for r = 1:length(Radii)
        % Moyenne pour MSD
        valid_msd = all_last_msd_values(:, r);
        valid_msd = valid_msd(~isnan(valid_msd));
        if ~isempty(valid_msd)
            mean_msd(r) = mean(valid_msd);
            % Écart-type pour MSD
            if length(valid_msd) > 1
                std_msd(r) = sqrt(sum((valid_msd - mean_msd(r)).^2) / (length(valid_msd) - 1));
            else
                std_msd(r) = 0; % Écart-type nul si une seule valeur
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
            % Écart-type pour voisins
            if length(valid_neighbors) > 1
                std_neighbors(r) = sqrt(sum((valid_neighbors - mean_neighbors(r)).^2) / (length(valid_neighbors) - 1));
            else
                std_neighbors(r) = 0; % Écart-type nul si une seule valeur
            end
        else
            mean_neighbors(r) = NaN;
            std_neighbors(r) = NaN;
        end
    end
    
    % Affichage des informations pertinentes dans la console
    fprintf('Résultats pour le fichier %s :\n', filename);
    fprintf('MSD moyen (dernier lag) pour Rayon = %.0f px : %.2f ± %.2f px²\n', ...
        Radii(end), mean_msd(end), std_msd(end));
    fprintf('Voisins moyens pour Rayon = %.0f px : %.2f ± %.2f\n\n', ...
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
        [mean_msd + std_msd, fliplr(mean_msd - std_msd)], ...
        msd_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(Radii, mean_msd, '-o', 'Color', msd_color, 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('MSD moyen (px²)');

    % Tracé voisins avec zone d'erreur (yyaxis droit)
    yyaxis right;
    hold on;
    fill([Radii fliplr(Radii)], ...
         [mean_neighbors + std_neighbors, fliplr(mean_neighbors - std_neighbors)], ...
        neighbor_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(Radii, mean_neighbors, '--s', 'Color', neighbor_color, 'LineWidth', 1.5, 'MarkerSize', 6);
    ylabel('Voisins moyens par image');

    %   Titres et axes
    title(sprintf('MSD moyen et Voisins vs Rayon - %s', filename));
    xlabel('Rayon (px)');
    grid on;
    legend({'Zone erreur MSD', 'MSD moyen', 'Zone erreur voisins', 'Voisins moyens'}, 'Location', 'best');

end

% --- Fonction de correction locale avec retour d'infos
function [corrected_traj, avg_speed, total_neighbors, avg_neighbors_per_frame] = ...
    correct_local_drift_individual_verbose(traj, AllTraj, Radius, Interval)
    
    % Initialisation
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    vx_all = []; vy_all = [];
    total_neighbors = 0;
    delta_t = 3;
    % Fenêtre temporelle pour la recherche des voisins (±3 frames)
    temporal_window = delta_t * Interval; % ±3 frames en temps (ms)

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
            total_neighbors = n_points - 1; % Nombre de segments utilisés
        else
            x_corrected = traj(:, 1);
            y_corrected = traj(:, 2);
        end
    else
        % Cas standard : utiliser les voisins dans AllTraj dans une fenêtre temporelle
        for i = 1:n_points
            x_i = traj(i, 1);
            y_i = traj(i, 2);
            t_i = traj(i, 3);

            % Recherche des voisins dans la fenêtre temporelle (±3 frames)
            time_diff = abs(AllTraj(:, 3) - t_i);
            dist2 = (AllTraj(:, 1) - x_i).^2 + (AllTraj(:, 2) - y_i).^2;
            neighbors_idx = find(time_diff <= temporal_window & dist2 <= Radius^2);
            total_neighbors = total_neighbors + length(neighbors_idx);

            vx = []; vy = [];
            % Calcul des vitesses à partir des voisins consécutifs dans la fenêtre
            if length(neighbors_idx) > 1
                % Trier les voisins par temps pour assurer des paires consécutives valides
                [~, sort_order] = sort(AllTraj(neighbors_idx, 3));
                sorted_neighbors_idx = neighbors_idx(sort_order);
                
                for j = 1:length(sorted_neighbors_idx)-1
                    idx1 = sorted_neighbors_idx(j);
                    idx2 = sorted_neighbors_idx(j+1);
                    dt = AllTraj(idx2, 3) - AllTraj(idx1, 3);
                    if dt > 0
                        vx(end+1) = (AllTraj(idx2, 1) - AllTraj(idx1, 1)) / dt;
                        vy(end+1) = (AllTraj(idx2, 2) - AllTraj(idx1, 2)) / dt;
                    end
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
    avg_neighbors_per_frame = total_neighbors / (n_points*delta_t);
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
        if n_points < 10 || any(any(~isfinite(traj(:, 1:2))))
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
        end
    end
end