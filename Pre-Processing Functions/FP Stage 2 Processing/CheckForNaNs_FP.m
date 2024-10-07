function [] = CheckForNaNs_FP(ProcData)
%________________________________________________________________________________________________________________________
% Edited by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%________________________________________________________________________________________________________________________
%
%   Purpose: Check the most importance CBV ROIs for any NaN values where the IOS camera glitched.
%________________________________________________________________________________________________________________________

% check ROI fields corresponding to the imaging type (single, bilateral hem) 
ROInames = fieldnames(ProcData.data.CBV);
% go through each ROI and check each individual value
for b = 1:length(ROInames)
    nanCheck{b,1} = sum(isnan(ProcData.data.CBV.(ROInames{b,1}))); %#ok<AGROW>
end
% pause the program if an NaN is found. Will need to add an interpolation method here if NaN events are found and
% the specific file needs to be kept
for b = 1:length(nanCheck)
    if nanCheck{b,1} ~= 0 
        disp('WARNING - NaNs found in CBV array'); disp(' ')
        keyboard
    end
end

end