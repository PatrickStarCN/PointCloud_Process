clc 
clear
close all

%Instructions :

%just uncomment one of the following line and run


%% Points cloud :
clear all
% load Cactus.mat
% load Skull.mat 
% load Standford_Bunny.mat
% load Horse.mat
% load hippo.mat 
% load Elephant.mat
% load chair.mat
% load 7.mat
% load Skull.mat
% load Block.mat
% load gargo50k.mat
b = importdata('reconstructionAABB.txt')
%% plot of the current point cloud
b=[b(:,1),b(:,2),b(:,3)];
figure(1);
hold on
axis equal
title('Points Cloud','fontsize',14)
plot3(b(:,1),b(:,2),b(:,3),'g.')
view(-37.5,30)

%% Run  program
[t]=reconstruction(b);

%% plot of the oyput triangulation
figure(2)
        hold on
        title('Output Triangulation','fontsize',14)
        axis equal
        trisurf(t,b(:,1),b(:,2),b(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
        view(-37.5,30)