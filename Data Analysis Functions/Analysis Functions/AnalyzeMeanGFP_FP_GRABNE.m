function [AnalysisResults] = AnalyzeMeanGFP_FP_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [GFP] during different arousal states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
modelType = 'Manual';
params.minTime.Rest = 10;
params.Offset = 15; 
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 10;
params.minTime.NREM = 30;%30;
params.minTime.REM = 60; %30;
%% only run analysis for valid animal IDs
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end
cd(dataLocation)
% find and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% lowpass filter
samplingRate = RestData.GFP.Z_Ach.RhodamineCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for the whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskPuffCriteria.Fieldname = {'puffDistance'};
WhiskPuffCriteria.Comparison = {'gt'};
WhiskPuffCriteria.Value = {5};
% criteria for the resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% criteria for the stimulation
StimCriteriaA.Value = {'RPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'LPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
%% analyze [GFP] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.GFP.Z_Ach,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.GFP.Z_Ach,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.GFP.Z_Ach.fileIDs(combRestLogical,:);
restEventTimes = RestData.GFP.Z_Ach.eventTimes(combRestLogical,:);
restDurations = RestData.GFP.Z_Ach.durations(combRestLogical,:);
Ach_RestingData = RestData.GFP.Z_Ach.data(combRestLogical,:);
NE_RestingData = RestData.GFP.Z_NE.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[Ach_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(Ach_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
[NE_finalRestData,~,~,~] = RemoveInvalidData_IOS(NE_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter [GFP]
for gg = 1:length(Ach_finalRestData)
    Ach_ProcRestData{gg,1} = filtfilt(sos,g,Ach_finalRestData{gg,1}); %#ok<*AGROW>
    NE_ProcRestData{gg,1} = filtfilt(sos,g,NE_finalRestData{gg,1});
end
% take mean [GFP] during resting epochs
for nn = 1:length(Ach_ProcRestData)
    Ach_restGFPMean(nn,1) = mean(Ach_ProcRestData{nn,1}(1:end));
    NE_restGFPMean(nn,1) = mean(NE_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanGFP.Rest.GFP.MeanAch = Ach_restGFPMean;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.MeanNE = NE_restGFPMean;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.IndAch = Ach_ProcRestData;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.IndNE = NE_ProcRestData;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.FileIDs = finalRestFileIDs;
%% analyze [GFP] during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.GFP.Z_Ach.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.GFP.Z_Ach.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.GFP.Z_Ach.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.GFP.Z_Ach.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.GFP.Z_Ach.whisk.duration(combWhiskLogical,:);
Ach_whiskData = EventData.GFP.Z_Ach.whisk.data(combWhiskLogical,:);
NE_whiskData = EventData.GFP.Z_NE.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[Ach_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(Ach_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
[NE_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(NE_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter [GFP] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(Ach_finalWhiskData,1)
    Ach_ProcWhiskData_temp = filtfilt(sos,g,Ach_finalWhiskData(gg,:));
    Ach_ProcWhiskData(gg,:) = Ach_ProcWhiskData_temp;% - mean(Ach_ProcWhiskData_temp(1:params.Offset*samplingRate));
    NE_ProcWhiskData_temp = filtfilt(sos,g,NE_finalWhiskData(gg,:));
    NE_ProcWhiskData(gg,:) = NE_ProcWhiskData_temp;% - mean(NE_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean [GFP] during whisking epochs from onset through 5 seconds
for nn = 1:size(Ach_ProcWhiskData,1)
    Ach_whiskGFPMean(nn,1) = mean(Ach_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    NE_whiskGFPMean(nn,1) = mean(NE_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    Ach_whiskGFP{nn,1} = Ach_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    NE_whiskGFP{nn,1} = NE_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.MeanAch = (Ach_whiskGFPMean);
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.MeanNE = (NE_whiskGFPMean);
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.IndAch = Ach_whiskGFP;
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.IndNE = NE_whiskGFP;
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.FileIDs = finalWhiskFileIDs;
%% analyze [GFP] during periods of stimulation
% pull data from EventData.mat structure
    if firstHrs == "true"
        Ach_stimFilter = FilterEvents_IOS(EventData.GFP.Z_Ach.stim,StimCriteriaA);
        NE_stimFilter = FilterEvents_IOS(EventData.GFP.Z_NE.stim,StimCriteriaB);
        [Ach_stimGFPData] = EventData.GFP.Z_Ach.stim.data(Ach_stimFilter,:);
        [NE_stimGFPData] = EventData.GFP.Z_NE.stim.data(NE_stimFilter,:);
        [Ach_stimFileIDs] = EventData.GFP.Z_Ach.stim.fileIDs(Ach_stimFilter,:);
        [NE_stimFileIDs] = EventData.GFP.Z_NE.stim.fileIDs(NE_stimFilter,:);
        [Ach_stimEventTimes] = EventData.GFP.Z_Ach.stim.eventTime(Ach_stimFilter,:);
        [NE_stimEventTimes] = EventData.GFP.Z_NE.stim.eventTime(NE_stimFilter,:);
        Ach_stimDurations = zeros(length(Ach_stimEventTimes),1);
        NE_stimDurations = zeros(length(NE_stimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [Ach_finalStimData,Ach_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(Ach_stimGFPData,Ach_stimFileIDs,Ach_stimDurations,Ach_stimEventTimes,ManualDecisions);
        [NE_finalStimData,NE_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(NE_stimGFPData,NE_stimFileIDs,NE_stimDurations,NE_stimEventTimes,ManualDecisions);
        % filter [GFP] and mean-subtract 2 seconds prior to stimulus (left hem)
        for gg = 1:size(Ach_finalStimData,1)
            Ach_ProcStimData_temp = filtfilt(sos,g,Ach_finalStimData(gg,:));
            Ach_ProcStimData(gg,:) = Ach_ProcStimData_temp;% - mean(Ach_ProcStimData_temp(1:params.Offset*samplingRate));
        end
        % filter [GFP] and mean-subtract 2 seconds prior to stimulus (right hem)
        for gg = 1:size(NE_finalStimData,1)
            NE_ProcStimData_temp = filtfilt(sos,g,NE_finalStimData(gg,:));
            NE_ProcStimData(gg,:) = NE_ProcStimData_temp;% - mean(NE_ProcStimData_temp(1:params.Offset*samplingRate));
        end
        % take mean [GFP] 1-2 seconds after stimulation (left hem)
        for nn = 1:size(Ach_ProcStimData,1)
            Ach_stimGFPMean{nn,1} = mean(Ach_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate),2);
            Ach_stimGFP{nn,1} = Ach_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate);
        end
        % take mean [GFP] 1-2 seconds after stimulation (right hem)
        for nn = 1:size(NE_ProcStimData,1)
            NE_stimGFPMean{nn,1} = mean(NE_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate),2);
            NE_stimGFP{nn,1} = NE_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate);
        end
        % save results
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.MeanAch = cell2mat(Ach_stimGFPMean);
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.MeanNE = cell2mat(NE_stimGFPMean);
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.IndAch = Ach_stimGFP;
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.IndNE = NE_stimGFP;
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.Ach_FileIDs = Ach_finalStimFileIDs;
        AnalysisResults.(animalID).MeanGFP.Stim.GFP.NE_FileIDs = NE_finalStimFileIDs;
    end
%% analyze [GFP] during periods of NREM sleep
% pull data from SleepData.mat structure
    if firstHrs == "true"
    [AchremData,nremFileIDs,~] = NRemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NEremData,~,~] = NRemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif firstHrs == "false"
    [AchremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NEremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
% filter and take mean [GFP] during NREM epochs
for nn = 1:length(AchremData)
    AchremGFPMean(nn,1) = mean(filtfilt(sos,g,AchremData{nn,1}(1:end)));
    NEremGFPMean(nn,1) = mean(filtfilt(sos,g,NEremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanGFP.NREM.GFP.MeanAch = AchremGFPMean;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.MeanNE = NEremGFPMean;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.IndAch = AchremData;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.IndNE = NEremData;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.FileIDs = nremFileIDs;
%% analyze [GFP] during periods of REM sleep
% pull data from SleepData.mat structure
if firstHrs == "false"
    [Ach_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Z_Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [NE_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Z_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter and take mean [GFP] during REM epochs
    for nn = 1:length(Ach_remData)
        Ach_remGFPMean(nn,1) = mean(filtfilt(sos,g,Ach_remData{nn,1}(1:end)));
        NE_remGFPMean(nn,1) = mean(filtfilt(sos,g,NE_remData{nn,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).MeanGFP.REM.GFP.MeanAch = Ach_remGFPMean;
    AnalysisResults.(animalID).MeanGFP.REM.GFP.MeanNE = NE_remGFPMean;
    AnalysisResults.(animalID).MeanGFP.REM.GFP.IndAch = Ach_remData;
    AnalysisResults.(animalID).MeanGFP.REM.GFP.IndNE = NE_remData;
    AnalysisResults.(animalID).MeanGFP.REM.GFP.FileIDs = remFileIDs;
end
 %% save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end

end
