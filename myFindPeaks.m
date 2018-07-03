%this fonction is used for calculating the pdf information
%Output:
%             IndMin:the location of the local minimum
%             IndMax:the location of the local maximum
%
%

function [IndMin,IndMax]=myFindPeaks(data,aixs)
IndMin=find(diff(sign(diff(data)))>0)+1;   %Get the location of the local minimum
IndMax=find(diff(sign(diff(data)))<0)+1;   %Get the location of the local maximum

%this part is used to plot the pdf information
%uncomment only in debug mode to see the infoamtion
%if not in debug mode ,there will be too mant pictures


% figure; hold on; box on;
% plot(1:length(data),data);
% plot(IndMin,data(IndMin),'r^')
% plot(IndMax,data(IndMax),'k*')
% legend('curve','Trough point','Peak point')
% title(['Calculate the peak and valley information aixs=',aixs], 'FontWeight', 'Bold');
% end
