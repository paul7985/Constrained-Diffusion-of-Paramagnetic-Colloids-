Interval = 15;
MaxDisp = 10;

% Sélection multiple de fichiers
[filenames, pathname] = uigetfile('*.csv', 'Select Results.csv files', 'MultiSelect', 'on');

% Si un seul fichier est sélectionné, conversion en cell array
if ischar(filenames)
    filenames = {filenames};
end

nFiles = length(filenames);

% Calcul du MSD pour chaque fichier (brut et corrigé)
msd_raw = cell(nFiles, 1);
msd_corrected = cell(nFiles, 1);

for i = 1:nFiles
    [msd_raw{i}, msd_corrected{i}] = compute_MSD(fullfile(pathname, filenames{i}), Interval, MaxDisp);
end

% Déterminer la longueur maximale pour l'axe temporel
max_msd_length = max(cellfun(@length, msd_raw));
time_axis = (1:max_msd_length) * Interval;

% Création de la figure avec deux sous-figures
figure;

% Sous-figure 1 : MSD bruts
subplot(2,1,1);
hold on;
for i = 1:nFiles
    plot(time_axis, msd_raw{i}, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Raw Mean Squared Displacements');
grid on;

% Sous-figure 2 : MSD corrigés
subplot(2,1,2);
hold on;
for i = 1:nFiles
    plot(time_axis, msd_corrected{i}, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filenames{i}));
end
hold off;
legend('Interpreter', 'none');
xlabel('Time lag \tau (s)');
ylabel('MSD');
title('Corrected Mean Squared Displacements');
grid on;

% Ajustement de la disposition
sgtitle('Comparison of Mean Squared Displacements');

function [msd_raw_mean, msd_corrected_mean] = compute_MSD(filepath, Interval, MaxDisp)
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata(filepath, delimiterIn, headerlinesIn);
    
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
    Trajectory = cell(NrOfTrajs, 1);
    
    for n = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == n);
        Trajectory{n} = AllTraj(indx, :);
    end
    
    % Filtrage des trajectoires courtes
    nTraj = 0;
    for i = 1:NrOfTrajs
        if length(Trajectory{i}(:, 1)) > 10
            nTraj = nTraj + 1;
            Traj{nTraj} = Trajectory{i};
        end
    end
    
    % Calcul de la vitesse moyenne (ux, uy) sur toutes les trajectoires
    total_displacement_x = 0;
    total_displacement_y = 0;
    total_time = 0;
    
    for i = 1:nTraj
        mat = Traj{i};
        x = mat(:, 1);
        y = mat(:, 2);
        t = mat(:, 3);
        
        dx = x(end) - x(1);
        dy = y(end) - y(1);
        dt = t(end) - t(1);
        
        if dt > 0
            total_displacement_x = total_displacement_x + dx;
            total_displacement_y = total_displacement_y + dy;
            total_time = total_time + dt;
        end
    end
    
    ux = total_displacement_x / total_time;
    uy = total_displacement_y / total_time;
    
    % Calcul des MSD brut et corrigé
    msd_raw = cell(nTraj, 1);
    msd_corrected = cell(nTraj, 1);
    
    for i = 1:nTraj
        mat = Traj{i};
        x = mat(:, 1);
        y = mat(:, 2);
        t = mat(:, 3);
        
        % Positions corrigées
        x_corrected = x - ux * t;
        y_corrected = y - uy * t;
        
        T = length(t);
        msd_raw{i} = zeros(T-1, 1);
        msd_corrected{i} = zeros(T-1, 1);
        
        for j = 1:T-1
            ind = 1:(T-j);
            % MSD brut
            msd_raw{i}(j) = mean((x(ind+j) - x(ind)).^2 + (y(ind+j) - y(ind)).^2);
            % MSD corrigé
            msd_corrected{i}(j) = mean((x_corrected(ind+j) - x_corrected(ind)).^2 + ...
                                     (y_corrected(ind+j) - y_corrected(ind)).^2);
        end
    end
    
    % Calcul des moyennes
    max_msd_length = max(cellfun(@length, msd_raw));
    msd_raw_matrix = NaN(nTraj, max_msd_length);
    msd_corrected_matrix = NaN(nTraj, max_msd_length);
    
    for i = 1:nTraj
        len = length(msd_raw{i});
        msd_raw_matrix(i, 1:len) = msd_raw{i};
        msd_corrected_matrix(i, 1:len) = msd_corrected{i};
    end
    
    msd_raw_mean = mean(msd_raw_matrix, 1, 'omitnan');
    msd_corrected_mean = mean(msd_corrected_matrix, 1, 'omitnan');
end