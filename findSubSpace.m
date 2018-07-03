%this fonction is used to find the subspace
%and update the node information


function [temSubSpace,index,count,space]=findSubSpace(space,index,count)
index1 =1
space.leaf=0
count1=1;
if(~isempty(space.data))
    %caculate the pdf and find trough
    [fx,~] = ksdensity(space.data(:,1));
    x='x';
    [IndMinX,IndMaxX]=myFindPeaks(fx,x);
    y='y';
    [fy,~] = ksdensity(space.data(:,2));
    [IndMinY,IndMaxY]=myFindPeaks(fy,y);
    [fz,~] = ksdensity(space.data(:,3));
    z='z';
    [IndMinZ,IndMaxZ]=myFindPeaks(fz,z);
    
    %find the boundary of subspace
    minX=min(space.data(:,1));
    maxX=max(space.data(:,1));
    minY=min(space.data(:,2));
    maxY=max(space.data(:,2));
    minZ=min(space.data(:,3));
    maxZ=max(space.data(:,3));
    
    [dis]=calculateDis( space.data,minX,minY,minZ,maxX,maxY,maxZ)
    %merge the boundary for temSubSpace
    IndMinX = [minX,minX+IndMinX*(abs(minX)+abs(maxX))/100,maxX];
    IndMinY = [minY,minY+IndMinY*(abs(minY)+abs(maxY))/100,maxY];
    IndMinZ = [minZ,minZ+IndMinZ*(abs(minZ)+abs(maxZ))/100,maxZ];
    
    if(length(IndMinX)==2)
        IndMinX= [minX,(minX+maxX)/2,maxX];
    end
    if(length(IndMinY)==2)
        IndMinY= [minY,(minY+maxY)/2,maxY];
    end
    if(length(IndMinZ)==2)
        IndMinZ= [minZ,(minZ+maxZ)/2,maxZ];
    end
    plot3(space.data(:,1),space.data(:,2),space.data(:,3),'.')
    hold on
    
    
    for l=1:length(IndMinX)-1
        for j=1:length(IndMinY)-1
            for k=1:length(IndMinZ)-1
                sub =find((IndMinX(1,l)<space.data(:,1))&(space.data(:,1)<IndMinX(1,l+1))&...
                    (IndMinY(1,j)<space.data(:,2))&(space.data(:,2)<IndMinY(1,j+1))&...
                    (IndMinZ(1,k)<space.data(:,3))&(space.data(:,3)<IndMinZ(1,k+1)));
                
                if(~sub==0)
                    %update the node information
                    temSubSpace(index1).data = space.data(sub,[1,2,3]);%update the data for subspace
                    temSubSpace(index1).index = index;%update the index for subspace
                    temSubSpace(index1).depth = space.depth+1;%update the depth for subspace
                    temSubSpace(index1).leaf = 1;%update the information of leaf node for subspace
                    
                    temSubSpace(index1).fatherIndex = space.index;
                    space.leaf=0;
                    minX=min(temSubSpace(index1).data(:,1));
                    maxX=max(temSubSpace(index1).data(:,1));
                    minY=min(temSubSpace(index1).data(:,2));
                    maxY=max(temSubSpace(index1).data(:,2));
                    minZ=min(temSubSpace(index1).data(:,3));
                    maxZ=max(temSubSpace(index1).data(:,3));
                    temSubSpace(index1).reP = [mean(temSubSpace(index1).data(:,1)),...
                    mean(temSubSpace(index1).data(:,2)),mean(temSubSpace(index1).data(:,3))];
                    temSubSpace(index1).volume =(maxX-minX)*(maxY-minY)*(maxZ-minZ);
                    [dis]=calculateDis( temSubSpace(index1).data,minX,minY,minZ,maxX,maxY,maxZ);
                    temSubSpace(index1).distor = dis;
                    
                    
%                this part is used for show the block for each subspace
%                also used to debug and to see the default of this algorithme

                    
                   
                    
                    
%                     if ~(isempty(temSubSpace(index1).data))
%                         plot3(temSubSpace(index1).data(:,1),temSubSpace(index1).data(:,2),temSubSpace(index1).data(:,3),'.')
%                         hold on
%                         if((maxY-minY~=0)&&(maxX-minX~=0)&&(maxZ-minZ~=0))
%                             [X,Z]=meshgrid(minX:maxX-minX:maxX,minZ:maxZ-minZ:maxZ);
%                             N=floor((maxY-minY)/(maxY-minY)+1);
%                             Y=minY*ones(N,N);
%                             mesh(X,Y,Z);
%                             alpha(.5)
%                             hold on
%                             
%                             [Y,Z]=meshgrid(minY:maxY-minY:maxY,minZ:maxZ-minZ:maxZ);
%                             N=floor((maxX-minX)/(maxX-minX)+1);
%                             X=minX*ones(N,N);
%                             mesh(X,Y,Z);%画出分割平面
%                             alpha(.5) %设置透明度
%                             hold on
%                             
%                             [X,Y]=meshgrid(minX:maxX-minX:maxX,minY:maxY-minY:maxY);
%                             N=floor((maxZ-minZ)/(maxZ-minZ)+1);
%                             Z=minZ*ones(N,N);
%                             mesh(X,Y,Z);
%                             alpha(.5)
%                             hold on
%                             
%                             [X,Z]=meshgrid(minX:maxX-minX:maxX,minZ:maxZ-minZ:maxZ);
%                             N=floor((maxY-minY)/(maxY-minY)+1);
%                             Y=maxY*ones(N,N);
%                             mesh(X,Y,Z);
%                             alpha(.5)
%                             hold on
%                             
%                             [Y,Z]=meshgrid(minY:maxY-minY:maxY,minZ:maxZ-minZ:maxZ);
%                             N=floor((maxX-minX)/(maxX-minX)+1);
%                             X=maxX*ones(N,N);
%                             mesh(X,Y,Z);%画出分割平面
%                             alpha(.5) %设置透明度
%                             hold on
%                             
%                             [X,Y]=meshgrid(minX:maxX-minX:maxX,minY:maxY-minY:maxY);
%                             N=floor((maxZ-minZ)/(maxZ-minZ)+1);
%                             Z=maxZ*ones(N,N);
%                             mesh(X,Y,Z);
%                             alpha(.5)
%                             hold on
%                         end
%                     end

                count = count+1;
                index = index +1;
                index1=index1+1;
                count1= count1+1
                
                
            end
        end
    end
    
    
    count1=count1-1;
    
end
end
end