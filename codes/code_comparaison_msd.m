Interval = 1;
MaxDisp = 10;

% Sélection des deux fichiers
[filename1, pathname1] = uigetfile('*.csv', 'Select the first file Results.csv');
[filename2, pathname2] = uigetfile('*.csv', 'Select the second file Results.csv');

% Calcul du MSD pour chaque fichier
msd_mean1 = compute_MSD(filename1, pathname1, Interval, MaxDisp);
msd_mean2 = compute_MSD(filename2, pathname2, Interval, MaxDisp);

% Déterminer la longueur maximale des MSD pour normaliser l'affichage
max_msd_length = max(length(msd_mean1), length(msd_mean2));
time_axis = (1:max_msd_length) * Interval;

% Adapter les tailles des vecteurs MSD en les complétant avec NaN si nécessaire
%msd_mean1 = padarray(msd_mean1, [0, max_msd_length - length(msd_mean1)], NaN, 'post');
%msd_mean2 = padarray(msd_mean2, [0, max_msd_length - length(msd_mean2)], NaN, 'post');

% Tracer les deux MSD sur la même figure
figure;
plot(time_axis, msd_mean1, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filename1));
hold on;
plot(time_axis, msd_mean2, '-s', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filename2));
hold off;
legend;
xlabel('Time Lag');
ylabel('Mean MSD');
title('Comparison of Mean Squared Displacements');
grid on;

function msd_mean = compute_MSD(filename, pathname, Interval, MaxDisp)
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
    Trajectory = cell(NrOfTrajs, 1);
    
    for n = 1:NrOfTrajs
        indx = find(AllTraj(:, 4) == n);
        Trajectory{n} = AllTraj(indx, :);
    end
    
    % Filtrage des trajectoires courtes
    nTraj = 0;
    for i = 1:NrOfTrajs
        if length(Trajectory{i}(:, 1)) > 500
            nTraj = nTraj + 1;
            Traj{nTraj} = Trajectory{i};
        end
    end
    
    % Calcul du MSD moyen
    msd = cell(nTraj, 1);
    for i = 1:nTraj
        mat = Traj{i};
        x = mat(:, 1);
        y = mat(:, 2);
        t = mat(:, 3);
        
        m = length(t);
        msd{i} = zeros(m - 1, 1);
        for j = 1:m - 1
            ind = 1:(m - j);
            msd{i}(j) = mean((x(ind + j) - x(ind)).^2 + (y(ind + j) - y(ind)).^2);
        end
    end
    
    max_msd_length = max(cellfun(@length, msd));
    msd_matrix = NaN(nTraj, max_msd_length);
    
    for i = 1:nTraj
        len = length(msd{i});
        msd_matrix(i, 1:len) = msd{i};
    end
    
    msd_mean = mean(msd_matrix, 1,'omitnan');
end