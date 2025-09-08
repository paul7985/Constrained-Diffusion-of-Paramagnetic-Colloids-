%%%Pgm pour charger une image et superposer les particules d�tect�es par
%%%macros ImageJ Macro_order_v2
clear all
close all

choixImg=1;
%%%Plot image
 [FileName,PathName] = uigetfile('Select location of the image');
[X,map] = imread([PathName,FileName],'tif');
figure(11),clf
image(X,'CDataMapping','scaled');
colormap gray
caxis auto
hold on

%%%Chargement data
%[FileName,PathName] = uigetfile('Results_Order.csv','Select location of the Results_Order.csv file');

[FileName,PathName] = uigetfile('Results.csv','Select location of the Results.csv file');

data = csvread([PathName,FileName],1,0);
%jcol = 6; %Use this line if images were analyzed with Macro_order_v3.ijm
jcol = 2; %Use this line if images were analyzed with Macro_densities_tris.ijm
%jcol = 28; %Use this line if images were analyzed with Macro_order_v2.ijm

    N = length(data(:,1));%Total number of particles
    iImgPart = data(:,jcol);%Image number for each particle
    numImg = max(iImgPart);%Total number of images
%     xpart = data(:,3);%Horizontal position of a particle
%     ypart = data(:,4);%Vertical position of a particle
    xpart = data(:,4);%Horizontal position of a particle
    ypart = data(:,5);%Vertical position of a particle
    iparts = find(iImgPart==choixImg);
   
     x = xpart(iparts);
     y = ypart(iparts);
     plot(x,y,'bs')