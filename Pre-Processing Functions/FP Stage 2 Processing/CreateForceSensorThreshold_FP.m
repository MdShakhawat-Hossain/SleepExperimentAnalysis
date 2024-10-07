function [thresh] = CreateForceSensorThreshold_FP(forceSensor,Fs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: View the force sensor data and determine what a good value is to binarize movement.
%________________________________________________________________________________________________________________________

y = hilbert(diff(forceSensor));
force = abs(y);
forceThresh = figure;
% figure;
isok = 'n';
ForceTime = (1:length(forceSensor))/Fs;
TotalTime = 15*60*Fs;
% plot(ForceTime(1:TotalTime),force(1:TotalTime),'k');

% for NN = 1:1:10
%     plot(ForceTime(1:TotalTime),force(1:TotalTime),'k');
%     thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
%     binForceSensor = BinarizeForceSensor_FP(forceSensor,thresh);
%     subplot(3,1,1)
%     plot(ForceTime(1:TotalTime),forceSensor(1:TotalTime),'k') 
%     axis tight
%     subplot(3,1,2) 
%     plot(ForceTime(1:TotalTime),force(1:TotalTime),'k')
%     axis tight
%     subplot(3,1,3)
%     plot(ForceTime(1:TotalTime),binForceSensor(1:TotalTime),'k')
%     axis tight
%     isok = input('Is this threshold okay? (y/n): ','s'); disp(' ')
%     if strcmp(isok,'y') == 1
%         break; 
%     end
% end

while strcmp(isok,'y') == 0
    plot(ForceTime(1:TotalTime),force(1:TotalTime),'k');
    drawnow
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
    binForceSensor = BinarizeForceSensor_FP(forceSensor,thresh);
    subplot(3,1,1)
    plot(ForceTime(1:TotalTime),forceSensor(1:TotalTime),'k') 
    axis tight
    subplot(3,1,2) 
    plot(ForceTime(1:TotalTime),force(1:TotalTime),'k')
    axis tight
    subplot(3,1,3)
    plot(ForceTime(1:TotalTime),binForceSensor(1:TotalTime),'k')
    axis tight
    isok = input('Is this threshold okay? (y/n): ','s'); disp(' ')
end
close 
close(forceThresh);

end
