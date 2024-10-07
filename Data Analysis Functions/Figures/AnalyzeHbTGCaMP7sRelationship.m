function [AnalysisResults] = AnalyzeHbTGCaMP7sRelationship(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the relationship between Ca2+ activity and hemodynamics 
%________________________________________________________________________________________________________________________
%% function parameters
dataLocation = [rootFolder  '\' animalID '\CombinedImaging\'];
cd(dataLocation)
% find and load manual baseline event information
load([ animalID '_Forest_ScoringResults.mat'])
% find and load EventData.mat struct
procDataFileStruct = dir('*_ProcData.mat');
procDataFile = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFile);
% extract/concatenate data from each file
catLH_TRITC = [];
catRH_TRITC = [];
catLH_GCaMP7s = [];
catRH_GCaMP7s = [];
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID)
    for bb = 1:length(ProcData.sleep.parameters.GCaMP7s.LH)
        LH_TRITC = mean(ProcData.sleep.parameters.TRITC.LH{bb,1});
        RH_TRITC = mean(ProcData.sleep.parameters.TRITC.RH{bb,1});
        LH_GCaMP7s = mean(ProcData.sleep.parameters.GCaMP7s.LH{bb,1});
        RH_GCaMP7s = mean(ProcData.sleep.parameters.GCaMP7s.RH{bb,1});
        % group date based on arousal-state classification
        state = ScoringResults.labels{aa,1}{bb,1};
        if strcmp(state,'Not Sleep') == true
            if isfield(catLH_TRITC,'Awake') == false
                catLH_TRITC.Awake = [];
                catRH_TRITC.Awake = [];
                catLH_GCaMP7s.Awake = [];
                catRH_GCaMP7s.Awake = [];
            end
            catLH_TRITC.Awake = cat(1,catLH_TRITC.Awake,LH_TRITC);
            catRH_TRITC.Awake = cat(1,catRH_TRITC.Awake,RH_TRITC);
            catLH_GCaMP7s.Awake = cat(1,catLH_GCaMP7s.Awake,LH_GCaMP7s);
            catRH_GCaMP7s.Awake = cat(1,catRH_GCaMP7s.Awake,RH_GCaMP7s);
        elseif strcmp(state,'NREM Sleep') == true
            if isfield(catLH_TRITC,'NREM') == false
                catLH_TRITC.NREM = [];
                catRH_TRITC.NREM = [];
                catLH_GCaMP7s.NREM = [];
                catRH_GCaMP7s.NREM = [];
            end
            catLH_TRITC.NREM = cat(1,catLH_TRITC.NREM,LH_TRITC);
            catRH_TRITC.NREM = cat(1,catRH_TRITC.NREM,RH_TRITC);
            catLH_GCaMP7s.NREM = cat(1,catLH_GCaMP7s.NREM,LH_GCaMP7s);
            catRH_GCaMP7s.NREM = cat(1,catRH_GCaMP7s.NREM,RH_GCaMP7s);
        elseif strcmp(state,'REM Sleep') == true
            if isfield(catLH_TRITC,'REM') == false
                catLH_TRITC.REM = [];
                catRH_TRITC.REM = [];
                catLH_GCaMP7s.REM = [];
                catRH_GCaMP7s.REM = [];
            end
            catLH_TRITC.REM = cat(1,catLH_TRITC.REM,LH_TRITC);
            catRH_TRITC.REM = cat(1,catRH_TRITC.REM,RH_TRITC);
            catLH_GCaMP7s.REM = cat(1,catLH_GCaMP7s.REM,LH_GCaMP7s);
            catRH_GCaMP7s.REM = cat(1,catRH_GCaMP7s.REM,RH_GCaMP7s);
        end
    end
end
% save results
AnalysisResults.(animalID).TRITCvsGCaMP7s = [];
AnalysisResults.(animalID).TRITCvsGCaMP7s.Awake.LH.TRITC = catLH_TRITC.Awake;
AnalysisResults.(animalID).TRITCvsGCaMP7s.NREM.LH.TRITC = catLH_TRITC.NREM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.REM.LH.TRITC = catLH_TRITC.REM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.Awake.RH.TRITC = catRH_TRITC.Awake;
AnalysisResults.(animalID).TRITCvsGCaMP7s.NREM.RH.TRITC = catRH_TRITC.NREM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.REM.RH.TRITC = catRH_TRITC.REM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.Awake.LH.GCaMP7s = catLH_GCaMP7s.Awake;
AnalysisResults.(animalID).TRITCvsGCaMP7s.NREM.LH.GCaMP7s = catLH_GCaMP7s.NREM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.REM.LH.GCaMP7s = catLH_GCaMP7s.REM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.Awake.RH.GCaMP7s = catRH_GCaMP7s.Awake;
AnalysisResults.(animalID).TRITCvsGCaMP7s.NREM.RH.GCaMP7s = catRH_GCaMP7s.NREM;
AnalysisResults.(animalID).TRITCvsGCaMP7s.REM.RH.GCaMP7s = catRH_GCaMP7s.REM;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
