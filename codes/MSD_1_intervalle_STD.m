Interval = 15;
MaxDisp = 20;

% Sélection de plusieurs fichiers
[filenames, pathname] = uigetfile('*.csv', 'Select files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames}; % Convertir en cellule si un seul fichier
end

numFiles = length(filenames);
msd_means = cell(numFiles, 1);
msd_stds = cell(numFiles, 1);
max_msd_length = 0;

% Calcul du MSD et de l'écart-type pour chaque fichier
for i = 1:numFiles
    [msd_means{i}, msd_stds{i}] = compute_MSD_t0(filenames{i}, pathname, Interval, MaxDisp);
    max_msd_length = max(max_msd_length, length(msd_means{i}));
end

time_axis = (1:max_msd_length) * Interval;

% Tracer les MSD avec bande d'erreur
figure;
hold on;
colors = lines(numFiles);
for i = 1:numFiles
    len = length(msd_means{i});
    fill([time_axis(1:len), fliplr(time_axis(1:len))], ...
         [msd_means{i} + msd_stds{i}, fliplr(msd_means{i} - msd_stds{i})], ...
         colors(i, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_axis(1:len), msd_means{i}, '-o', 'Color', colors(i, :), 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Mean Squared Displacement with Standard Deviation (relative to t_0)');
grid on;

function [msd_mean, msd_std] = compute_MSD_t0(filename, pathname, Interval, MaxDisp)
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
    
    % Calcul du MSD basé sur un t0 commun
    msd = cellfun(@compute_msd_t0_single, Trajectory, 'UniformOutput', false);
    max_msd_length = max(cellfun(@length, msd));
    msd_matrix = NaN(nTraj, max_msd_length);
    
    for i = 1:nTraj
        len = length(msd{i});
        msd_matrix(i, 1:len) = msd{i};
    end
    
    msd_mean = mean(msd_matrix, 1, 'omitnan');
    msd_std = std(msd_matrix, 0, 1, 'omitnan');
end

function msd = compute_msd_t0_single(traj)
    t0_x = traj(1, 1);
    t0_y = traj(1, 2);
    T = size(traj, 1);
    msd = arrayfun(@(j) mean((traj(j, 1) - t0_x).^2 + (traj(j, 2) - t0_y).^2), 1:T);
end
