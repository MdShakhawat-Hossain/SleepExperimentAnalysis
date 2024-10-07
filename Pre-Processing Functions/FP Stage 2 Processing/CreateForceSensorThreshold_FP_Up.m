function [thresh] = CreateForceSensorThreshold_FP_Up(forceSensor)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: View the force sensor data and determine what a good value is to binarize movement.
%________________________________________________________________________________________________________________________

y = hilbert(diff(forceSensor));
force = abs(y);
forceThresh = figure;
isok = 'n';
while strcmp(isok,'y') == 0
    plot(force,'k');
    drawnow
    thresh = input('No Threshold to binarize pressure sensor found. Please enter a threshold: '); disp(' ')
    binForceSensor = BinarizeForceSensor_FP(forceSensor,thresh);
    subplot(3,1,1)
    plot(forceSensor,'k'); hold on; yline(thresh,'r');
    axis tight
    subplot(3,1,2)
    plot(force,'k')
    axis tight
    subplot(3,1,3)
    plot(binForceSensor,'k')
    axis tight
    drawnow
    isok = input('Is this threshold okay? (y/n): ','s'); disp(' ')
end
close(forceThresh);

end
