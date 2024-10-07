function [binForceSensor] = BinarizeForceSensor_FP(forceSensor,thresh)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Binarize the Force sensor with a given threshold.
%________________________________________________________________________________________________________________________

y = hilbert(diff(forceSensor));
env = abs(y);
binForceSensor = gt(env,thresh);

end
