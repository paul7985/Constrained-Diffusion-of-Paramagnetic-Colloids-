Interval = 1;
MaxDisp = 10;


[filename1, pathname1] = uigetfile('*.csv', 'Select the first file Results.csv');
msd_mean1 = compute_MSD(filename1, pathname1, Interval, MaxDisp);
max_msd_length = length(msd_mean1);
time_axis = (1:max_msd_length) * Interval;

% Tracer les MSD sur la même figure
figure;
plot(time_axis, msd_mean1, '-o', 'LineWidth', 1.5, 'DisplayName', sprintf('MSD %s', filename1));
hold on;
legend;
xlabel('Time (s)');
ylabel('MSD');
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
        if length(Trajectory{i}(:, 1)) > 300
            nTraj = nTraj + 1;
            Traj{nTraj} = Trajectory{i};
        end
    end
    
    % Calcul de la vitesse moyenne (ux, uy) sur toutes les trajectoires
    total_displacement_x = 0;
    total_displacement_y = 0;
    total_time = 0;
    n_points = 0;
    
    for i = 1:nTraj
        mat = Traj{i};
        x = mat(:, 1);
        y = mat(:, 2);
        t = mat(:, 3);
        
        % Calcul du déplacement total et du temps total pour cette trajectoire
        dx = x(end) - x(1); % Déplacement total en x
        dy = y(end) - y(1); % Déplacement total en y
        dt = t(end) - t(1); % Durée totale
        
        if dt > 0 % Éviter division par zéro
            total_displacement_x = total_displacement_x + dx;
            total_displacement_y = total_displacement_y + dy;
            total_time = total_time + dt;
            n_points = n_points + 1;
        end
    end
    
    % Vitesse moyenne
    ux = total_displacement_x / total_time; % Vitesse moyenne en x
    uy = total_displacement_y / total_time; % Vitesse moyenne en y
    
    % Correction des trajectoires et calcul du MSD
    msd = cell(nTraj, 1);
    for i = 1:nTraj
        mat = Traj{i};
        x = mat(:, 1);
        y = mat(:, 2);
        t = mat(:, 3);
        
        % Correction des positions : x = x - ux*t, y = y - uy*t
        x_corrected = x - ux * t;
        y_corrected = y - uy * t;
        
        m = length(t);
        msd{i} = zeros(m - 1, 1);
        for j = 1:m - 1
            ind = 1:(m - j);
            msd{i}(j) = mean((x_corrected(ind + j) - x_corrected(ind)).^2 + ...
                            (y_corrected(ind + j) - y_corrected(ind)).^2);
        end
        
        % Affichage du MSD individuel
        figure;
        %plot((1:length(msd{i})) * Interval, msd{i}, '-o', 'LineWidth', 1.5);
        plot((1:length(msd{nTraj-i})) * Interval, msd{nTraj-i}, '-o', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('MSD');
        title(sprintf('MSD for Particle %d', nTraj-i));
        grid on;
        keyboard;
    end
    
    max_msd_length = max(cellfun(@length, msd));
    msd_matrix = NaN(nTraj, max_msd_length);
    
    for i = 1:nTraj
        len = length(msd{i});
        msd_matrix(i, 1:len) = msd{i};
    end
    
    msd_mean = mean(msd_matrix, 1, 'omitnan');
end
