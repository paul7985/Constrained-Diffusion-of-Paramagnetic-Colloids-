Interval = 15;
MaxDisp = 20;

% Sélection de plusieurs fichiers
[filenames, pathname] = uigetfile('*.csv', 'Select files', 'MultiSelect', 'on');
if ischar(filenames)
    filenames = {filenames}; % Convertir en cellule si un seul fichier
end

numFiles = length(filenames);
msd_means = cell(numFiles, 1);
max_msd_length = 0;

% Calcul du MSD pour chaque fichier
for i = 1:numFiles
    msd_means{i} = compute_MSD(filenames{i}, pathname, Interval, MaxDisp);
    max_msd_length = max(max_msd_length, length(msd_means{i}));
end

time_axis = (1:max_msd_length) * Interval;

% Tracer les MSD sur la même figure
figure;
hold on;
for i = 1:numFiles
    plot(time_axis(1:length(msd_means{i})), msd_means{i}, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Mean Squared Displacement');
grid on;

function msd_mean = compute_MSD(filename, pathname, Interval, MaxDisp)
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
    total_displacement_x = 0;
    total_displacement_y = 0;
    total_time = 0;
    
    for i = 1:nTraj
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
    
    % Correction des trajectoires et calcul du MSD
    msd = cellfun(@(traj) compute_msd_single(traj, ux, uy), Trajectory, 'UniformOutput', false);
    max_msd_length = max(cellfun(@length, msd));
    msd_matrix = NaN(nTraj, max_msd_length);
    
    for i = 1:nTraj
        len = length(msd{i});
        msd_matrix(i, 1:len) = msd{i};
    end
    
    msd_mean = mean(msd_matrix, 1, 'omitnan');
end

function msd = compute_msd_single(traj, ux, uy)
    x_corrected = traj(:, 1) - ux * traj(:, 3);
    y_corrected = traj(:, 2) - uy * traj(:, 3);
    
    T = size(traj, 1);
    msd = arrayfun(@(j) mean((x_corrected(1:T-j) - x_corrected(1+j:T)).^2 + (y_corrected(1:T-j) - y_corrected(1+j:T)).^2), 1:T-1);
end
