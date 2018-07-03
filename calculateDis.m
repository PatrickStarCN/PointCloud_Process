%this foction is used to calculate the distortion
%information for each subspace

function [dis]=calculateDis(subSpaceData,minX,minY,minZ,maxX,maxY,maxZ)
dis=0;
TI = [mean(subSpaceData(:,1)),...
      mean(subSpaceData(:,2)),...
      mean(subSpaceData(:,3))];
  
for i=1:size(subSpaceData,1)
    dis=dis+(sum((subSpaceData(i,:)-TI).^2,2));
end
end