% Seeking the average point cloud density
% To find the average number of data points contained in the unit volume to characterize the density of the point cloud.
% The smallest cuboid bounding box that determines the point cloud data

clear all  \
clc;
in = importdata('Horse.txt');%Load point cloud data
%% plot of the current point cloud
figure(1);
hold on
axis equal  %The units of length of axes are set equal
title('Points Cloud','fontsize',14)
plot3(in(:,1),in(:,2),in(:,3),'g.')
grid
view(-37.5,30)
xlabel('X');ylabel('Y');zlabel('Z');

[N,b]=size(in);e=0.8;K=1;  %N is the number of rows£¬b is3
max_x=max(in(:,1));min_x=min(in(:,1));
max_y=max(in(:,2));min_y=min(in(:,2));
max_z=max(in(:,3));min_z=min(in(:,3));
Lx=max_x-min_x;
Ly=max_y-min_y;
Lz=max_z-min_z;

%% Determines the number of data points in a small cube raster
V=Lx*Ly*Lz;   %The maximum size of the bounding box

%the number of data points in a small cube  is£º
n=N/(V+eps);

%Determine the side length Ls of the sub-cube grid, ¦Ë is a scale
%factor, which is used to adjust the side length of the
%sub-cube grid, K is the number of adjacent points

Ls=(e*K/(n+eps))^(1/3);
mm=floor(Lx/(Ls+eps))+1;
nn=floor(Ly/(Ls+eps))+1;
ll=floor(Lz/(Ls+eps))+1;%The point cloud data is divided into m ¡Á n ¡Á l small cube grids
hh=floor((in(:,1)-min_x)/(Ls+eps))+1;
jj=floor((in(:,2)-min_y)/(Ls+eps))+1;
kk=floor((in(:,3)-min_z)/(Ls+eps))+1;   %Make sure that each point cloud is in that square
%Find the bounding box center coordinates
center=zeros(mm*nn*ll,3);
for i=1:mm
    center_x=min_x+Ls*(i-0.5);
    for j=1:nn
        center_y=min_y+Ls*(j-0.5);
        for k=1:ll
            center_z=min_z+Ls*(k-0.5);
            center((i-1)*nn*ll+(j-1)*ll+k,1)=center_x;
            center((i-1)*nn*ll+(j-1)*ll+k,2)=center_y;
            center((i-1)*nn*ll+(j-1)*ll+k,3)=center_z;
        end
    end
end
%Find the center point and the distance between each point cloud
distance=zeros(N,4);
for i=1:N
    distance(i,1)=hh(i);
    distance(i,2)=jj(i);
    distance(i,3)=kk(i);
    distance(i,4)=sqrt((center((hh(i)-1)*nn*ll+(jj(i)-1)*ll+kk(i),1)-in(i,1))^2+(center((hh(i)-1)*nn*ll+(jj(i)-1)*ll+kk(i),2)-in(i,2))^2+(center((hh(i)-1)*nn*ll+(jj(i)-1)*ll+kk(i),3)-in(i,3))^2);
end

[Y mm nn]=unique(distance(:,1:3),'rows')
num=length(mm);
X=zeros(num,1);
for i=1:num
    X(i)=max(distance(nn==i,4));
end
Y=[Y X];

%% plot of the current point cloud
figure(2);
hold off
axis equal
title('Points Cloud','fontsize',14)
plot3(Y(:,1),Y(:,2),Y(:,3),'g.')
fid=fopen('reconstructionAABB.txt','a')
for g=1:length(Y)
    restructionP=Y(g,1:3)
    [m,n]=size( restructionP);
    
    for i=1:1:m
        for j=1:1:n
            if j==n
                fprintf(fid,'%g\n',restructionP(i,j));
            else fprintf(fid,'%g\t',restructionP(i,j));
            end
        end
    end
end
grid
view(-37.5,30)
xlabel('X');ylabel('Y');zlabel('Z');
%
%% Run  program
p=Y(:,1:3);
[t]=MyCrust(p);
[w]=MyCrust(in);
%% plot of the oyput triangulation
figure(3)
hold on
title('Output Triangulation','fontsize',14)
axis equal
trisurf(w,in(:,1),in(:,2),in(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
grid
view(-37.5,30)
xlabel('X');ylabel('Y');zlabel('Z');
%

%% plot of the oyput triangulation
figure(4)
hold on
title('Output Triangulation','fontsize',14)
axis equal
trisurf(t,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
grid
view(-37.5,30)
xlabel('X');ylabel('Y');zlabel('Z');


