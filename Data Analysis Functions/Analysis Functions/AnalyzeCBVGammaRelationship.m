function [AnalysisResults] = AnalyzeCBVGammaRelationship(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the relationship between gamma-band power and hemodynamics [TRITC]
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
catLH_TRITC = [];
catRH_TRITC = [];
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    for bb = 1:length(ProcData.sleep.parameters.cortical_LH.gammaBandPower)
        LH_Gamma = mean(ProcData.sleep.parameters.cortical_LH.specGammaBandPower{bb,1}{1,1});
        RH_Gamma = mean(ProcData.sleep.parameters.cortical_RH.specGammaBandPower{bb,1}{1,1});
        LH_TRITC = mean(ProcData.sleep.parameters.TRITC.LH{bb,1});
        RH_TRITC = mean(ProcData.sleep.parameters.TRITC.RH{bb,1});
        % group date based on arousal-state classification
        state = ScoringResults.labels{aa,1}{bb,1};
        if strcmp(state,'Not Sleep') == true
            if isfield(catLH_Gamma,'Awake') == false
                catLH_Gamma.Awake = [];
                catRH_Gamma.Awake = [];
                catLH_TRITC.Awake = [];
                catRH_TRITC.Awake = [];
            end
            catLH_Gamma.Awake = cat(1,catLH_Gamma.Awake,LH_Gamma);
            catRH_Gamma.Awake = cat(1,catRH_Gamma.Awake,RH_Gamma);
            catLH_TRITC.Awake = cat(1,catLH_TRITC.Awake,LH_TRITC);
            catRH_TRITC.Awake = cat(1,catRH_TRITC.Awake,RH_TRITC);
        elseif strcmp(state,'NREM Sleep') == true
            if isfield(catLH_Gamma,'NREM') == false
                catLH_Gamma.NREM = [];
                catRH_Gamma.NREM = [];
                catLH_TRITC.NREM = [];
                catRH_TRITC.NREM = [];
            end
            catLH_Gamma.NREM = cat(1,catLH_Gamma.NREM,LH_Gamma);
            catRH_Gamma.NREM = cat(1,catRH_Gamma.NREM,RH_Gamma);
            catLH_TRITC.NREM = cat(1,catLH_TRITC.NREM,LH_TRITC);
            catRH_TRITC.NREM = cat(1,catRH_TRITC.NREM,RH_TRITC);
        elseif strcmp(state,'REM Sleep') == true
            if isfield(catLH_Gamma,'REM') == false
                catLH_Gamma.REM = [];
                catRH_Gamma.REM = [];
                catLH_TRITC.REM = [];
                catRH_TRITC.REM = [];
            end
            catLH_Gamma.REM = cat(1,catLH_Gamma.REM,LH_Gamma);
            catRH_Gamma.REM = cat(1,catRH_Gamma.REM,RH_Gamma);
            catLH_TRITC.REM = cat(1,catLH_TRITC.REM,LH_TRITC);
            catRH_TRITC.REM = cat(1,catRH_TRITC.REM,RH_TRITC);
        end
    end
end
% save results
AnalysisResults.(animalID).TRITCvsGamma = [];
AnalysisResults.(animalID).TRITCvsGamma.Awake.LH.Gamma = catLH_Gamma.Awake;
AnalysisResults.(animalID).TRITCvsGamma.NREM.LH.Gamma = catLH_Gamma.NREM;
AnalysisResults.(animalID).TRITCvsGamma.REM.LH.Gamma = catLH_Gamma.REM;
AnalysisResults.(animalID).TRITCvsGamma.Awake.RH.Gamma = catRH_Gamma.Awake;
AnalysisResults.(animalID).TRITCvsGamma.NREM.RH.Gamma = catRH_Gamma.NREM;
AnalysisResults.(animalID).TRITCvsGamma.REM.RH.Gamma = catRH_Gamma.REM;
AnalysisResults.(animalID).TRITCvsGamma.Awake.LH.TRITC = catLH_TRITC.Awake;
AnalysisResults.(animalID).TRITCvsGamma.NREM.LH.TRITC = catLH_TRITC.NREM;
AnalysisResults.(animalID).TRITCvsGamma.REM.LH.TRITC = catLH_TRITC.REM;
AnalysisResults.(animalID).TRITCvsGamma.Awake.RH.TRITC = catRH_TRITC.Awake;
AnalysisResults.(animalID).TRITCvsGamma.NREM.RH.TRITC = catRH_TRITC.NREM;
AnalysisResults.(animalID).TRITCvsGamma.REM.RH.TRITC = catRH_TRITC.REM;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
