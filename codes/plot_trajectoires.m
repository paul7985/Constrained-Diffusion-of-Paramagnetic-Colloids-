Interval = 1;
MaxDisp = 10;

% Sélection multiple des fichiers CSV
[filenames, pathname] = uigetfile('*.csv', 'Select Results.csv files', 'MultiSelect', 'on');

% Vérification si un seul fichier ou plusieurs ont été sélectionnés
if ischar(filenames)  % Cas d'un seul fichier
    filenames = {filenames};
end

% Boucle sur chaque fichier sélectionné
for f = 1:length(filenames)
    filename = filenames{f};
    plot_trajectories(filename, pathname, Interval, MaxDisp);
end

function plot_trajectories(filename, pathname, Interval, MaxDisp)
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
    
    % Création d'une figure avec deux sous-figures
    figure;
    
    % Sous-figure 1 : Trajectoires originales
    hold on;
    title(sprintf('Original Trajectories\n%s', filename), 'Interpreter', 'none');
    xlabel('X Position');
    ylabel('Y Position');
    
    for n = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == n);
        SingleTraj = AllTraj(indx, :);
        plot(SingleTraj(:, 1), SingleTraj(:, 2), '.-');
    end
    hold off;
    
    % Sous-figure 2 : Trajectoires corrigées
    % subplot(1, 2, 2);
    % hold on;
    % title(sprintf('Drift-Corrected Trajectories\n%s', filename), 'Interpreter', 'none');
    % xlabel('X Position');
    % ylabel('Y Position');
    
    % for n = 1:NrOfTrajs
    %     indx = find(AllTraj(:, 4) == n);
    %     SingleTraj = AllTraj(indx, :);
    %     x_corrected = SingleTraj(:, 1) - ux * SingleTraj(:, 3);
    %     y_corrected = SingleTraj(:, 2) - uy * SingleTraj(:, 3);
    %     plot(x_corrected, y_corrected, '.-');
    % end
    % hold off;
    
    % % Ajustement de la taille de la figure
    % set(gcf, 'Position', [100, 100, 800, 400]);
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