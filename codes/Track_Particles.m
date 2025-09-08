Interval = 10;
MaxDips = 15;
TrackParticles(Interval, MaxDips);
function [Dmean,Traj,uxmean,uymean,umean]=Track_Particles(Interval, MaxDisp);

% Particle Tracking from Macro_movie.ijm ImageJ
% Connects particle trajectoires des particules using the
% algorithm Crocker-Grier-Weeks in track.m 
% En sortie les trajectoire y(t) avec points finaux marquï¿½s '*'
%
%%NB Lyd Mai 2019: Trajectoires considérés >100s; si fichiers<100s,
%%modifier paramètre 100 ligne 85


% INPUTS: 
% The function can be used with two inputs, then :
%    It will pop up a dialogue to retrieve data from Results.csv, containing the
%    particle trajectories X and Y
%    The file Results.csv should contain columns with headers 'Frame', 'X', and 'Y'.
%    Interval:  Time lapse between two consecutive images    
%    MaxDisp:   Maximum particle displacement between two images, in pixels
% If the function are used without input, particle trajectories will be
% read from Trajectories.mat and not calculated
% OUTPUTS  : 
% Particle trajectories (cell array)
% Average drift velocity (ux, uy, u=sqrt(ux^2+uy^2)
% Mean and median diffusion coefficient (corrected for mean drift)

calib=1.8; %%%conversion pix->um
if nargin==2
    [filename,pathname] = uigetfile('*.csv','Select the file Results.csv to be analyzed');
    delimiterIn = ',';
    headerlinesIn = 1;
    %A = importdata([pathname filename],delimiterIn,headerlinesIn);
    A = importdata([pathname filename]);
    i_Frame = find(strcmp('Frame',A.colheaders));
    i_X = find(strcmp('X',A.colheaders));
    i_Y = find(strcmp('Y',A.colheaders));
    
    time=(A.data(:,i_Frame)-1)*Interval;
    InputForTrack=[A.data(:,i_X),A.data(:,i_Y),time];  % define Input (X,Y,t) for track.m
    %Clean NaN data
    i1 = find(isnan(InputForTrack(:,1))==1);
    i2 = find(isnan(InputForTrack(:,2))==1);
    i3 = find(isnan(InputForTrack(:,3))==1);
    inan = union(i1,i2);
    inan = union(inan,i3);
    InputForTrack(inan,:) = [];
    
    AllTraj=track  (InputForTrack,MaxDisp); % TRACKING using track.m
    %%pause;
    % AllTraj is a matrix whose columns are [X, Y, t, traj number]
    NrOfTrajs=max(AllTraj(:,4));
    ['Total Nr of Trajectories = ', num2str(NrOfTrajs)];
    %Plot trajectories
    figure(1)
    clf
    hold on
% 
%    
    for n=1:NrOfTrajs
        indx=find(AllTraj(:,4)==n);
        SingleTraj=AllTraj(indx,:);
        Trajectory{n}=SingleTraj;
%         
plot(SingleTraj(:,1),SingleTraj(:,2),'.-');
       %pause; 
     end
%     
     hold off
    
    save([pathname 'Trajectories.mat'])
    
elseif nargin==0
    [filename,pathname] = uigetfile('*.mat','Select the file Trajectories.mat storing the trajectories');
    load([pathname filename]);
else
    error('Incorrect number of input arguments');
end

%Analyze trajectories
nTrajraw = length(Trajectory);
sz = size(Trajectory);
Traj = cell(sz);

%Eliminate trajectories whose duration is less than 10 seconds
nTraj=0;
for i=1:nTrajraw
    if length(Trajectory{i}(:,1))>round(100/Interval)
        nTraj = nTraj+1;
        Traj{nTraj} = Trajectory{i};
    end
end

% %Plot drifts
% figure
% clf
% hold on
% for i=1:nTraj
%     mat = Traj{i};
%     x = mat(:,1);
%     y = mat(:,2);
%     t = mat(:,3);
%     m = length(t);
%     %plot(x(1),y(1),'o')
%     %plot(x([1,end]),y([1,end]),'-')
%     plot(x,y,'-')
% end
% hold off


%Analyze trajectories For each trajectory, estimate drift and subtract it
ux = zeros(nTraj,1);
uy = zeros(nTraj,1);
D = zeros(nTraj,1);
itotaltime = 0;
for i=1:nTraj
    
    mat = Traj{i};
    x = mat(:,1);
    y = mat(:,2);
    t = mat(:,3);
    itotaltime = max(itotaltime,max(t));
    [ux(i),uy(i),D(i),msd{i}] = analyze_trajectory(x,y,t);
    %for j=1:length(msd{i})
        
    %pause;
end

%Obtain means and standard deviations
uxmean = mean(ux);
uxstd = std(ux);
uymean = mean(uy);
uystd = std(uy);
Dmean = mean(D);
Dstd = std(D);
umean=sqrt(uxmean^2+uymean^2);
ustd=sqrt(uxstd^2+uystd^2);
    
pix_per_um=calib;
display(['Calculations based on ',num2str(nTraj),' trajectories, ',num2str(itotaltime),' s:'])
display(['uxmean=',num2str(uxmean/pix_per_um),' std=, ',num2str(uxstd/pix_per_um),' um/s'])
display(['uymean=',num2str(uymean/pix_per_um),' std=, ',num2str(uystd/pix_per_um),' um/s'])
display(['umean=',num2str(umean/pix_per_um),' std=, ',num2str(ustd/pix_per_um),' um/s'])
display(['Dmean=',num2str(Dmean/pix_per_um^2),' std=, ',num2str(Dstd/pix_per_um^2),' um2/s'])

figure(1)
clf
subplot(3,1,1)
histogram(ux)
title('ux')
subplot(3,1,2)
histogram(uy)
title('uy')
subplot(3,1,3)
histogram(D)
title('D')

return
end

function [ux,uy,D,msd] = analyze_trajectory(x,y,t)
%Computes mean drift and diffusion coefficient

%Time interval
Interval = t(2)-t(1);

%Estimate drift
ux = (x(end)-x(1)) / (t(end)-t(1));
uy = (y(end)-y(1)) / (t(end)-t(1));

%Subtract drift from trajectory
x = x - ux*(t-t(1)); 
y = y - uy*(t-t(1));

%Compute msd
m = length(t);
msd = zeros(m-1,1);
for i=1:m-1
    ind = 1:(m-i);
    msd(i) = mean((x(ind+i)-x(ind)).^2 + (y(ind+i)-y(ind)).^2);
end

%Estimate diffusion coefficient for msd corresponding to the 1/10 shortest
%time lags
imax = ceil((m-1)/10);
D = ([1:imax]'*Interval)\msd(1:imax);
D = D/4; %2D diffusion coefficient

figure(2)
clf
hold on
plot([1:(m-1)]*Interval,msd,'o')
plot([1:(m-1)]*Interval,[1:(m-1)]*Interval*4*D,'-k')
title(['ux=',num2str(ux),'uy=',num2str(uy)])
hold off
% keyboard

return
end