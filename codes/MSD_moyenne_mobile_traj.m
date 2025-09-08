Interval = 1;
MaxDisp = 10;
Radius = 100; % Rayon pour la moyenne mobile

% Sélection multiple des fichiers CSV
[filenames, pathname] = uigetfile('*.csv', 'Select Results.csv files', 'MultiSelect', 'on');

% Vérification si un seul fichier ou plusieurs ont été sélectionnés
if ischar(filenames)  % Cas d'un seul fichier
    filenames = {filenames};
end

% Boucle sur chaque fichier sélectionné
for f = 1:length(filenames)
    filename = filenames{f};
    plot_trajectories(filename, pathname, Interval, MaxDisp, Radius);
end

function plot_trajectories(filename, pathname, Interval, MaxDisp, Radius)
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata([pathname filename], delimiterIn, headerlinesIn);
    
    i_Frame = find(strcmp('Frame', A.colheaders));
    i_X = find(strcmp('X', A.colheaders));
    i_Y = find(strcmp('Y', A.colheaders));
    
    time = (A.data(:, i_Frame) - 1) * Interval;
    InputForTrack = [A.data(:, i_X), A.data(:, i_Y), time];
    
    % Nettoyage des NaN
    InputForTrack(any(isnan(InputForTrack), 2), :) = [];
    
    % Suivi des particules
    AllTraj = track(InputForTrack, MaxDisp);
    NrOfTrajs = max(AllTraj(:, 4));
    
    % Calcul de la vitesse de drift moyenne
    [ux, uy] = compute_drift(AllTraj);
    
    % Préparation des structures pour stocker les trajectoires
    original_trajectories = cell(NrOfTrajs, 1);
    corrected_trajectories = cell(NrOfTrajs, 1);
    
    % Extraction des trajectoires et application de la correction globale
    for n = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == n);
        SingleTraj = AllTraj(indx, :);
        
        % Trajectoire originale
        original_trajectories{n} = [SingleTraj(:, 1), SingleTraj(:, 2), SingleTraj(:, 3)];
        
        % Trajectoire corrigée globalement
        x_corrected = SingleTraj(:, 1) - ux * SingleTraj(:, 3);
        y_corrected = SingleTraj(:, 2) - uy * SingleTraj(:, 3);
        corrected_trajectories{n} = [x_corrected, y_corrected, SingleTraj(:, 3)];
    end
    
    % Application de la correction locale
    local_corrected_trajectories = cell(NrOfTrajs, 1);
    
    % Pour chaque trajectoire individuelle
    for n = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == n);
        SingleTraj = AllTraj(indx, :);
        
        % Appliquer la correction locale à cette trajectoire spécifique
        [x_local_corrected, y_local_corrected] = correct_local_drift_individual(SingleTraj, AllTraj, Radius);
        
        % Stocker la trajectoire corrigée localement
        local_corrected_trajectories{n} = [x_local_corrected, y_local_corrected, SingleTraj(:, 3)];
    end
    
    % Calcul du MSD pour chaque type de trajectoire
    % Utiliser la longueur maximale des trajectoires au lieu d'une valeur fixe
    max_frames = max(AllTraj(:, 3)) / Interval;
    max_lag = min(60, floor(max_frames)); % Utiliser jusqu'à 60 points ou le maximum disponible
    
    [msd_original, time_lags] = compute_msd(original_trajectories, max_lag, Interval);
    [msd_corrected, ~] = compute_msd(corrected_trajectories, max_lag, Interval);
    [msd_local_corrected, ~] = compute_msd(local_corrected_trajectories, max_lag, Interval);
    
    % Création d'une figure avec six sous-figures (3 pour trajectoires, 3 pour MSD)
    figure;
    
    % Sous-figure 1 : Trajectoires originales
    subplot(2, 3, 1);
    hold on;
    title(sprintf('Original Trajectories\n%s', filename), 'Interpreter', 'none');
    xlabel('X Position');
    ylabel('Y Position');
    
    for n = 1:NrOfTrajs
        plot(original_trajectories{n}(:, 1), original_trajectories{n}(:, 2), '.-');
    end
    hold off;
    
    % Sous-figure 2 : Trajectoires corrigées globalement
    subplot(2, 3, 2);
    hold on;
    title(sprintf('Global Drift-Corrected Trajectories\n%s', filename), 'Interpreter', 'none');
    xlabel('X Position');
    ylabel('Y Position');
    
    for n = 1:NrOfTrajs
        plot(corrected_trajectories{n}(:, 1), corrected_trajectories{n}(:, 2), '.-');
    end
    hold off;
    
    % Sous-figure 3 : Trajectoires corrigées localement
    subplot(2, 3, 3);
    hold on;
    title(sprintf('Local Drift-Corrected Trajectories\n%s', filename), 'Interpreter', 'none');
    xlabel('X Position');
    ylabel('Y Position');
    
    for n = 1:NrOfTrajs
        plot(local_corrected_trajectories{n}(:, 1), local_corrected_trajectories{n}(:, 2), '.-');
    end
    hold off;
    
    % Sous-figure 4 : MSD des trajectoires originales
    subplot(2, 3, 4);
    hold on;
    title('MSD - Original', 'Interpreter', 'none');
    xlabel('Time Lag (s)');
    ylabel('MSD');
    
    % Tracé du MSD 
    plot(time_lags, msd_original, 'b.-', 'LineWidth', 1.5);
    
    % Ajout d'une ligne de référence pour la diffusion normale (pente 1)
    if ~isnan(msd_original(2))
        ref_line = time_lags .* (msd_original(2) / time_lags(2));
        plot(time_lags, ref_line, 'k--', 'DisplayName', 'Slope = 1');
        legend('MSD', 'Slope = 1', 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Sous-figure 5 : MSD des trajectoires corrigées globalement
    subplot(2, 3, 5);
    hold on;
    title('MSD - Global Correction', 'Interpreter', 'none');
    xlabel('Time Lag (s)');
    ylabel('MSD');
    
    % Tracé du MSD 
    plot(time_lags, msd_corrected, 'r.-', 'LineWidth', 1.5);
    
    % Ajout d'une ligne de référence pour la diffusion normale (pente 1)
    if ~isnan(msd_corrected(2))
        ref_line = time_lags .* (msd_corrected(2) / time_lags(2));
        plot(time_lags, ref_line, 'k--', 'DisplayName', 'Slope = 1');
        legend('MSD', 'Slope = 1', 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Sous-figure 6 : MSD des trajectoires corrigées localement
    subplot(2, 3, 6);
    hold on;
    title('MSD - Local Correction', 'Interpreter', 'none');
    xlabel('Time Lag (s)');
    ylabel('MSD');
    
    % Tracé du MSD
    plot(time_lags, msd_local_corrected, 'g.-', 'LineWidth', 1.5);
    
    % Ajout d'une ligne de référence pour la diffusion normale (pente 1)
    if ~isnan(msd_local_corrected(2))
        ref_line = time_lags .* (msd_local_corrected(2) / time_lags(2));
        plot(time_lags, ref_line, 'k--', 'DisplayName', 'Slope = 1');
        legend('MSD', 'Slope = 1', 'Location', 'best');
    end
    
    grid on;
    hold off;
    
    % Ajustement de la taille de la figure
    set(gcf, 'Position', [50, 50, 1400, 800]);
    
    % Ajouter du débogage pour vérifier les valeurs NaN
    fprintf('MSD original first 5 values: %s\n', mat2str(msd_original(1:min(5,length(msd_original)))));
    fprintf('MSD global correction first 5 values: %s\n', mat2str(msd_corrected(1:min(5,length(msd_corrected)))));
    fprintf('MSD local correction first 5 values: %s\n', mat2str(msd_local_corrected(1:min(5,length(msd_local_corrected)))));
end

function [ux, uy] = compute_drift(AllTraj)
    NrOfTrajs = max(AllTraj(:, 4));
    
    total_displacement_x = 0;
    total_displacement_y = 0;
    total_time = 0;
    
    for i = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == i);
        traj = AllTraj(indx, :);
        
        dx = traj(end, 1) - traj(1, 1);
        dy = traj(end, 2) - traj(1, 2);
        dt = traj(end, 3) - traj(1, 3);
        
        if dt > 0
            total_displacement_x = total_displacement_x + dx;
            total_displacement_y = total_displacement_y + dy;
            total_time = total_time + dt;
        end
    end
    
    if total_time > 0
        ux = total_displacement_x / total_time;
        uy = total_displacement_y / total_time;
    else
        ux = 0;
        uy = 0;
    end
end

function [x_corrected, y_corrected] = correct_local_drift_individual(traj, AllTraj, Radius)
    % Correction locale appliquée à une trajectoire spécifique
    n_points = size(traj, 1);
    x_corrected = zeros(n_points, 1);
    y_corrected = zeros(n_points, 1);
    
    for i = 1:n_points
        % Position et temps actuels
        x_i = traj(i, 1);
        y_i = traj(i, 2);
        t_i = traj(i, 3);
        
        % Trouver les voisins spatiaux dans toutes les trajectoires
        distances = sqrt((AllTraj(:, 1) - x_i).^2 + (AllTraj(:, 2) - y_i).^2);
        neighbors_idx = find(distances < Radius);
        
        if length(neighbors_idx) > 1
            % Extraire les coordonnées des voisins
            x_neighbors = AllTraj(neighbors_idx, 1);
            y_neighbors = AllTraj(neighbors_idx, 2);
            t_neighbors = AllTraj(neighbors_idx, 3);
            
            % Calculer les vecteurs vitesse pour chaque paire de points consécutifs
            vx = zeros(length(neighbors_idx)-1, 1);
            vy = zeros(length(neighbors_idx)-1, 1);
            
            for j = 1:length(neighbors_idx)-1
                if t_neighbors(j+1) > t_neighbors(j)
                    dt = t_neighbors(j+1) - t_neighbors(j);
                    vx(j) = (x_neighbors(j+1) - x_neighbors(j)) / dt;
                    vy(j) = (y_neighbors(j+1) - y_neighbors(j)) / dt;
                end
            end
            
            % Filtrer les valeurs non finies
            valid_idx = isfinite(vx) & isfinite(vy);
            
            if sum(valid_idx) > 0
                % Calculer la moyenne des vitesses
                ux_local = mean(vx(valid_idx));
                uy_local = mean(vy(valid_idx));
                
                % Appliquer la correction
                x_corrected(i) = x_i - ux_local * t_i;
                y_corrected(i) = y_i - uy_local * t_i;
            else
                % Si pas de vitesse valide, conserver les coordonnées originales
                x_corrected(i) = x_i;
                y_corrected(i) = y_i;
            end
        else
            % Si pas assez de voisins, conserver les coordonnées originales
            x_corrected(i) = x_i;
            y_corrected(i) = y_i;
        end
    end
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