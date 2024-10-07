function [] = AnalyzeBilateralCoherence_Handler(rootFolder,delim,runFromStart)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the spectral coherence between bilateral hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

% create or load results structure
if runFromStart == true
    Results_BilatCoher = [];
elseif runFromStart == false
    % load existing results structure, if it exists
    if exist('Results_BilatCoher.mat','file') == 2
        load('Results_BilatCoher.mat','-mat')
    else
        Results_BilatCoher = [];
    end
end
setName = 'IOS Set A';
expGroups = {'Naive','SSP-SAP','Blank-SAP'};
% determine waitbar length
waitBarLength = 0;
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name}, '.'));
    folderAnimalIDs = {folderList.name};
    waitBarLength = waitBarLength + length(folderAnimalIDs);
end
% run analysis for each animal in the group
cc = 1;
for aa = 1:length(expGroups)
    folderList = dir([expGroups{1,aa} delim setName]);
    folderList = folderList(~startsWith({folderList.name},'.'));
    animalIDs = {folderList.name};
    multiWaitbar('Analyzing bilateral coherence',0,'Color','P'); pause(0.25);
    for bb = 1:length(animalIDs)
        if isfield(Results_BilatCoher,(animalIDs{1,bb})) == false
            [Results_BilatCoher] = AnalyzeBilateralCoherence(animalIDs{1,bb},[expGroups{1,aa} delim 'IOS Set A'],rootFolder,delim,Results_BilatCoher);
        end
        multiWaitbar('Analyzing bilateral coherence','Value',cc/waitBarLength);
        cc = cc + 1;
    end
end
multiWaitbar('CloseAll');

end
