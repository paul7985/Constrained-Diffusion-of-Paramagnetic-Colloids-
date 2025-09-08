Interval = 15;
MaxDisp = 10;
XLimits = [1400, 1700]; % Limites en X
YLimits = [0, 200]; % Limites en Y

% Sélection de plusieurs fichiers
[filenames, pathname] = uigetfile('*.csv', 'Select files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames}; % Convertir en cellule si un seul fichier
end

numFiles = length(filenames);
msd_means_raw = cell(numFiles, 1);
msd_means_corrected = cell(numFiles, 1);
max_msd_length = 0;

% Traitement des fichiers
for i = 1:numFiles
    [msd_means_raw{i}, msd_means_corrected{i}, rawTrajectories, correctedTrajectories] = compute_MSD(filenames{i}, pathname, Interval, MaxDisp, XLimits, YLimits);
    max_msd_length = max(max_msd_length, length(msd_means_raw{i}));
    
    % Affichage des trajectoires brutes et corrigées
    plot_trajectories(rawTrajectories, correctedTrajectories, filenames{i});
end

time_axis = (1:max_msd_length) * Interval;

% Tracer les MSD bruts et corrigés
figure;
hold on;
for i = 1:numFiles
    plot(time_axis(1:length(msd_means_raw{i})), msd_means_raw{i}, '--o', 'LineWidth', 1.5, 'DisplayName', sprintf('Raw MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Mean Squared Displacements (Raw)');
grid on;

figure;
hold on;
for i = 1:numFiles
    plot(time_axis(1:length(msd_means_corrected{i})), msd_means_corrected{i}, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('Corrected MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Mean Squared Displacements (Corrected)');
grid on;

function [msd_raw, msd_corrected, rawTrajectories, correctedTrajectories] = compute_MSD(filename, pathname, Interval, MaxDisp, XLimits, YLimits)
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata(fullfile(pathname, filename), delimiterIn, headerlinesIn);
    
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
    Trajectory = arrayfun(@(n) AllTraj(AllTraj(:, 4) == n, :), 1:NrOfTrajs, 'UniformOutput', false);
    
    % Filtrage des trajectoires courtes
    Trajectory = Trajectory(cellfun(@(traj) size(traj, 1) > 10, Trajectory));
    nTraj = length(Trajectory);
    
    % Calcul de la vitesse moyenne (ux, uy)
    [ux, uy] = compute_drift(Trajectory);
    
    % Calcul du MSD brut et corrigé
    [msd_raw, rawTrajectories] = cellfun(@(traj) compute_msd(traj, XLimits, YLimits), Trajectory, 'UniformOutput', false);
    [msd_corrected, correctedTrajectories] = cellfun(@(traj) correct_and_compute_msd(traj, ux, uy, XLimits, YLimits), Trajectory, 'UniformOutput', false);
    
    % Retirer les trajectoires qui n'ont pas passé le filtre
    validRaw = ~cellfun(@isempty, msd_raw);
    msd_raw = msd_raw(validRaw);
    rawTrajectories = rawTrajectories(validRaw);
    
    validCorrected = ~cellfun(@isempty, msd_corrected);
    msd_corrected = msd_corrected(validCorrected);
    correctedTrajectories = correctedTrajectories(validCorrected);
    
    % Moyenne des MSDs
    msd_raw = mean_msd(msd_raw);
    msd_corrected = mean_msd(msd_corrected);
end

function [ux, uy] = compute_drift(Trajectory)
    total_displacement_x = 0;
    total_displacement_y = 0;
    total_time = 0;
    
    for i = 1:length(Trajectory)
        mat = Trajectory{i};
        dx = mat(end, 1) - mat(1, 1);
        dy = mat(end, 2) - mat(1, 2);
        dt = mat(end, 3) - mat(1, 3);
        
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

function [msd, traj_filtered] = compute_msd(traj, XLimits, YLimits)
    % Filtrage des positions
    valid_indices = (traj(:, 1) >= XLimits(1) & traj(:, 1) <= XLimits(2)) & ...
                    (traj(:, 2) >= YLimits(1) & traj(:, 2) <= YLimits(2));
                
    x_filtered = traj(valid_indices, 1);
    y_filtered = traj(valid_indices, 2);
    time_filtered = traj(valid_indices, 3);
    
    traj_filtered = [x_filtered, y_filtered, time_filtered];
    
    if length(time_filtered) < 2
        msd = [];
        return;
    end

    % Calcul du MSD
    T = length(time_filtered);
    msd = arrayfun(@(j) mean((x_filtered(1:T-j) - x_filtered(1+j:T)).^2 + ...
                             (y_filtered(1:T-j) - y_filtered(1+j:T)).^2), 1:T-1);
end

function [msd, traj_corrected] = correct_and_compute_msd(traj, ux, uy, XLimits, YLimits)
    x_corrected = traj(:, 1) - ux * traj(:, 3);
    y_corrected = traj(:, 2) - uy * traj(:, 3);

    valid_indices = (x_corrected >= XLimits(1) & x_corrected <= XLimits(2)) & ...
                    (y_corrected >= YLimits(1) & y_corrected <= YLimits(2));
                
    x_filtered = x_corrected(valid_indices);
    y_filtered = y_corrected(valid_indices);
    time_filtered = traj(valid_indices, 3);
    
    traj_corrected = [x_filtered, y_filtered, time_filtered];

    if length(time_filtered) < 2
        msd = [];
        return;
    end

    % Calcul du MSD
    T = length(time_filtered);
    msd = arrayfun(@(j) mean((x_filtered(1:T-j) - x_filtered(1+j:T)).^2 + ...
                             (y_filtered(1:T-j) - y_filtered(1+j:T)).^2), 1:T-1);
end

function msd_mean = mean_msd(msd_list)
    if isempty(msd_list)
        msd_mean = NaN;
        return;
    end
    max_length = max(cellfun(@length, msd_list));
    msd_matrix = NaN(length(msd_list), max_length);
    
    for i = 1:length(msd_list)
        len = length(msd_list{i});
        msd_matrix(i, 1:len) = msd_list{i};
    end
    
    msd_mean = mean(msd_matrix, 1, 'omitnan');
end

function plot_trajectories(rawTrajectories, correctedTrajectories, filename)
    figure;
    subplot(1, 2, 1);
    hold on;
    title(sprintf('Raw Trajectories - %s', filename), 'Interpreter', 'none');
    for i = 1:length(rawTrajectories)
        plot(rawTrajectories{i}(:, 1), rawTrajectories{i}(:, 2), '.-');
    end
    hold off;

    subplot(1, 2, 2);
    hold on;
    title(sprintf('Corrected Trajectories - %s', filename), 'Interpreter', 'none');
    for i = 1:length(correctedTrajectories)
        plot(correctedTrajectories{i}(:, 1), correctedTrajectories{i}(:, 2), '.-');
    end
    hold off;
end
