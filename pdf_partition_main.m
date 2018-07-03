clear all;
%This program is used for point clouds representation
%
%function used:
%myFindPeaks.m----used to calculate the pdf information
%calculateDis.m----used to calculate the distortion
%findSubSpace.m----used  to find the subspace and form the tree
%
%Input:
%               -- training_martix is a Nx3 array containing the 3D set of points
%              
%Output:
%               -- restructionP(Reconstruction point),distorSum(Sum of distortion),
%               -- countP(counter for the leaf node)
%               -- occVol(Occupied Volume)
%Created by YuYAN and NannanXU

%% load original point clouds data;
training_martix = importdata('ball.txt');

%load the data of different axis
training_martix_x = training_martix(:,1);
training_martix_y = training_martix(:,2);
training_martix_z = training_martix(:,3);

%draw the original point clouds
figure(1)

plot3(training_martix_x,training_martix_y,training_martix_z,'.');


%Create the structure for nodes
space = struct('depth',{0},'leaf',{1},'data',{training_martix},...
    'distor',{inf},'index',{0},'fatherIndex','null','reP',{[0,0,0]},...
    'volume',{inf})
space.data = [training_martix_x,training_martix_y,training_martix_z]
count = 1;
index =1 ;


%% find the sub space and
[subSpace{1},index,count,space(1)]=findSubSpace(space(1),index,count);
%if you want to control the depth 
%please change the length for the loop
% depth=length(loop)+1
% for depth 1 ,comment this part

for ii=1:2
    count=index-1;
    subSpace{ii+1}=[];
    index =1 ;   
    for i=1:count
        if ((length(subSpace{ii}(i).data)>10))
            [k,index,count,subSpace{ii}(i)]=findSubSpace(subSpace{ii}(i),index,count);
            subSpace{ii+1}=[subSpace{ii+1},k];
        end
    end   
end

%% judge if is the leaf nod plot the reconstruction point 

%initialize the -- distorSum(Sum of distortion)
%               -- countP(counter for the leaf node)
%               -- occVol(Occupied Volume)
 
distorSum=0;
countP =0;
occVol =0;

fid=fopen('reconstructionPDF.txt','a');%save the data for reconstruction program

figure(2)
%Traverse all nodes and determine whether it is a leaf node
for i=1:length(subSpace)
    for j=1:length(subSpace{i})
        if(subSpace{i}(j).leaf==1)
            plot3(subSpace{i}(j).reP(1),subSpace{i}(j).reP(2),subSpace{i}(j).reP(3),'.');%plot the reconstruction point
            
            %find the boundary of the sub space
            minX=min(subSpace{i}(j).data(:,1));
            maxX=max(subSpace{i}(j).data(:,1));
            minY=min(subSpace{i}(j).data(:,2));
            maxY=max(subSpace{i}(j).data(:,2));
            minZ=min(subSpace{i}(j).data(:,3));
            maxZ=max(subSpace{i}(j).data(:,3));
            
            %calculate the distortion and occypied volume
            [dis]=calculateDis(subSpace{i}(j).data,minX,minY,minZ,maxX,maxY,maxZ);
            distorSum = dis + distorSum;
            occVol = occVol + subSpace{i}(j).volume;

            %output the reconstrucion point for reconstruction program
            reconstructionP=[subSpace{i}(j).reP(1),subSpace{i}(j).reP(2),subSpace{i}(j).reP(3)];
            [m,n]=size( reconstructionP);    
            for l=1:1:m
                for p=1:1:n
                    if p==n
                         fprintf(fid,'%g\n',reconstructionP(l,p));
                    else fprintf(fid,'%g\t',reconstructionP(l,p));
                    end
                end
            end           
            hold on
            countP = countP+1;
        end
    end
end
fclose(fid);
% average distortion
distorSum=distorSum/countP;






