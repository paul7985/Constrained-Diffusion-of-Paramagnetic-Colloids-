Interval = 0.5;
MaxDisp = 10;
dr = 1; % Pas de l'histogramme pour g(r)
r_max = 120; % Distance maximale pour g(r)

% Sélection du fichier
[filename, pathname] = uigetfile('*.csv', 'Select the file Results.csv');

% Calcul de g(r) pour chaque intervalle de temps
[g_r_matrix, r_vals, t_vals] = compute_GR(filename, pathname, Interval, MaxDisp, dr, r_max);

% Tracer g(r) en fonction du temps
figure;
imagesc(t_vals, r_vals, g_r_matrix);
colormap("turbo");
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Time (s)');
ylabel('r');
title('Radial Distribution Function g(r) over Time');
grid on;

function [g_r_matrix, r_vals, t_vals] = compute_GR(filename, pathname, Interval, MaxDisp, dr, r_max)
    delimiterIn = ',';
    headerlinesIn = 1;
    A = importdata([pathname filename], delimiterIn, headerlinesIn);
    
    i_Frame = find(strcmp('Frame', A.colheaders));
    i_X = find(strcmp('X', A.colheaders));
    i_Y = find(strcmp('Y', A.colheaders));
    
    time = (A.data(:, i_Frame) - 1) * Interval;
    InputForTrack = [A.data(:, i_X), A.data(:, i_Y), time];
    
    InputForTrack(any(isnan(InputForTrack), 2), :) = [];
    AllTraj = track(InputForTrack, MaxDisp);
    
    t_vals = unique(time);
    r_vals = dr:dr:r_max;
    g_r_matrix = zeros(length(r_vals), length(t_vals));
    
    for t_idx = 1:length(t_vals)
        t = t_vals(t_idx);
        frame_idx = find(time == t);
        if isempty(frame_idx), continue; end
        
        pos = [A.data(frame_idx, i_X), A.data(frame_idx, i_Y)];
        nPart = size(pos, 1);
        rho = nPart / (max(pos(:,1)) * max(pos(:,2)));
        
        g_r = zeros(size(r_vals));
        
        for j = 1:nPart-1
            for k = j+1:nPart
                r = norm(pos(j, :) - pos(k, :));
                bin = floor(r / dr) + 1;
                if bin <= length(g_r)
                    g_r(bin) = g_r(bin) + 2; % Comptage symétrique
                end
            end
        end
        
        % Normalisation
        for i = 1:length(r_vals)
            r = r_vals(i);
            shell_volume = pi * ((r + dr)^2 - r^2); % Surface en 2D
            g_r(i) = g_r(i) / (nPart * rho * shell_volume);
        end
        
        g_r_matrix(:, t_idx) = g_r;
    end
end