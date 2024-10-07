function [AnalysisResults] = AnalyzeGCaMP7sGammaRelationship(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the relationship between gamma-band power and hemodynamics [GCaMP7s]
%________________________________________________________________________________________________________________________

%% function parameters
dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
cd(dataLocation)
% find and load manual baseline event information
load([ animalID '_Forest_ScoringResults.mat'])
% find and load EventData.mat struct
procDataFileStruct = dir('*_ProcData.mat');
procDataFile = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFile);
% extract/concatenate data from each file
catLH_Gamma = [];
catRH_Gamma = [];
catLH_GCaMP7s = [];
catRH_GCaMP7s = [];
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    for bb = 1:length(ProcData.sleep.parameters.cortical_LH.gammaBandPower)
        LH_Gamma = mean(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{bb,1}{1,1});
        RH_Gamma = mean(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{bb,1}{1,1});
        LH_GCaMP7s = mean(ProcData.sleep.parameters.GCaMP7s.LH{bb,1});
        RH_GCaMP7s = mean(ProcData.sleep.parameters.GCaMP7s.RH{bb,1});
        % group date based on arousal-state classification
        state = ScoringResults.labels{aa,1}{bb,1};
        if strcmp(state,'Not Sleep') == true
            if isfield(catLH_Gamma,'Awake') == false
                catLH_Gamma.Awake = [];
                catRH_Gamma.Awake = [];
                catLH_GCaMP7s.Awake = [];
                catRH_GCaMP7s.Awake = [];
            end
            catLH_Gamma.Awake = cat(1,catLH_Gamma.Awake,LH_Gamma);
            catRH_Gamma.Awake = cat(1,catRH_Gamma.Awake,RH_Gamma);
            catLH_GCaMP7s.Awake = cat(1,catLH_GCaMP7s.Awake,LH_GCaMP7s);
            catRH_GCaMP7s.Awake = cat(1,catRH_GCaMP7s.Awake,RH_GCaMP7s);
        elseif strcmp(state,'NREM Sleep') == true
            if isfield(catLH_Gamma,'NREM') == false
                catLH_Gamma.NREM = [];
                catRH_Gamma.NREM = [];
                catLH_GCaMP7s.NREM = [];
                catRH_GCaMP7s.NREM = [];
            end
            catLH_Gamma.NREM = cat(1,catLH_Gamma.NREM,LH_Gamma);
            catRH_Gamma.NREM = cat(1,catRH_Gamma.NREM,RH_Gamma);
            catLH_GCaMP7s.NREM = cat(1,catLH_GCaMP7s.NREM,LH_GCaMP7s);
            catRH_GCaMP7s.NREM = cat(1,catRH_GCaMP7s.NREM,RH_GCaMP7s);
        elseif strcmp(state,'REM Sleep') == true
            if isfield(catLH_Gamma,'REM') == false
                catLH_Gamma.REM = [];
                catRH_Gamma.REM = [];
                catLH_GCaMP7s.REM = [];
                catRH_GCaMP7s.REM = [];
            end
            catLH_Gamma.REM = cat(1,catLH_Gamma.REM,LH_Gamma);
            catRH_Gamma.REM = cat(1,catRH_Gamma.REM,RH_Gamma);
            catLH_GCaMP7s.REM = cat(1,catLH_GCaMP7s.REM,LH_GCaMP7s);
            catRH_GCaMP7s.REM = cat(1,catRH_GCaMP7s.REM,RH_GCaMP7s);
        end
    end
end
% save results
AnalysisResults.(animalID).GCaMP7svsGamma = [];
AnalysisResults.(animalID).GCaMP7svsGamma.Awake.LH.Gamma = catLH_Gamma.Awake;
AnalysisResults.(animalID).GCaMP7svsGamma.NREM.LH.Gamma = catLH_Gamma.NREM;
AnalysisResults.(animalID).GCaMP7svsGamma.REM.LH.Gamma = catLH_Gamma.REM;
AnalysisResults.(animalID).GCaMP7svsGamma.Awake.RH.Gamma = catRH_Gamma.Awake;
AnalysisResults.(animalID).GCaMP7svsGamma.NREM.RH.Gamma = catRH_Gamma.NREM;
AnalysisResults.(animalID).GCaMP7svsGamma.REM.RH.Gamma = catRH_Gamma.REM;
AnalysisResults.(animalID).GCaMP7svsGamma.Awake.LH.GCaMP7s = catLH_GCaMP7s.Awake;
AnalysisResults.(animalID).GCaMP7svsGamma.NREM.LH.GCaMP7s = catLH_GCaMP7s.NREM;
AnalysisResults.(animalID).GCaMP7svsGamma.REM.LH.GCaMP7s = catLH_GCaMP7s.REM;
AnalysisResults.(animalID).GCaMP7svsGamma.Awake.RH.GCaMP7s = catRH_GCaMP7s.Awake;
AnalysisResults.(animalID).GCaMP7svsGamma.NREM.RH.GCaMP7s = catRH_GCaMP7s.NREM;
AnalysisResults.(animalID).GCaMP7svsGamma.REM.RH.GCaMP7s = catRH_GCaMP7s.REM;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
