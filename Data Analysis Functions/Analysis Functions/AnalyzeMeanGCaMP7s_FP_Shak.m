function [AnalysisResults] = AnalyzeMeanGCaMP7s_FP_Shak(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [GCaMP7s] during different arousal states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 2;
params.minTime.Whisk = params.Offset + 3;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' animalID '\CombinedImaging'];
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
samplingRate = RestData.GCaMP7s.LH.TRITCCamSamplingRate;
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
%% analyze [GCaMP7s] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.GCaMP7s.LH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.GCaMP7s.LH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.GCaMP7s.LH.fileIDs(combRestLogical,:);
restEventTimes = RestData.GCaMP7s.LH.eventTimes(combRestLogical,:);
restDurations = RestData.GCaMP7s.LH.durations(combRestLogical,:);
LH_RestingData = RestData.GCaMP7s.LH.data(combRestLogical,:);
RH_RestingData = RestData.GCaMP7s.RH.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[LH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
[RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter [GCaMP7s]
for gg = 1:length(LH_finalRestData)
    LH_ProcRestData{gg,1} = filtfilt(sos,g,LH_finalRestData{gg,1}); %#ok<*AGROW>
    RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1});
end
% take mean [GCaMP7s] during resting epochs
for nn = 1:length(LH_ProcRestData)
    LH_restGCaMP7sMean(nn,1) = mean(LH_ProcRestData{nn,1}(1:end));
    RH_restGCaMP7sMean(nn,1) = mean(RH_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanGCaMP7s.Rest.GCaMP7s.MeanLH = LH_restGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.Rest.GCaMP7s.MeanRH = RH_restGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.Rest.GCaMP7s.IndLH = LH_ProcRestData;
AnalysisResults.(animalID).MeanGCaMP7s.Rest.GCaMP7s.IndRH = RH_ProcRestData;
AnalysisResults.(animalID).MeanGCaMP7s.Rest.GCaMP7s.FileIDs = finalRestFileIDs;
%% analyze [GCaMP7s] during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.GCaMP7s.LH.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.GCaMP7s.LH.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.GCaMP7s.LH.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.GCaMP7s.LH.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.GCaMP7s.LH.whisk.duration(combWhiskLogical,:);
LH_whiskData = EventData.GCaMP7s.LH.whisk.data(combWhiskLogical,:);
RH_whiskData = EventData.GCaMP7s.RH.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[LH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
[RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter [GCaMP7s] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(LH_finalWhiskData,1)
    LH_ProcWhiskData_temp = filtfilt(sos,g,LH_finalWhiskData(gg,:));
    LH_ProcWhiskData(gg,:) = LH_ProcWhiskData_temp - mean(LH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
    RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean [GCaMP7s] during whisking epochs from onset through 5 seconds
for nn = 1:size(LH_ProcWhiskData,1)
    LH_whiskGCaMP7sMean{nn,1} = mean(LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    RH_whiskGCaMP7sMean{nn,1} = mean(RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    LH_whiskGCaMP7s{nn,1} = LH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    RH_whiskGCaMP7s{nn,1} = RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.MeanLH = cell2mat(LH_whiskGCaMP7sMean);
AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.MeanRH = cell2mat(RH_whiskGCaMP7sMean);
AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.IndLH = LH_whiskGCaMP7s;
AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.IndRH = RH_whiskGCaMP7s;
AnalysisResults.(animalID).MeanGCaMP7s.Whisk.GCaMP7s.FileIDs = finalWhiskFileIDs;
%% analyze [GCaMP7s] during periods of stimulation
% pull data from EventData.mat structure
LH_stimFilter = FilterEvents_IOS(EventData.GCaMP7s.LH.stim,StimCriteriaA);
RH_stimFilter = FilterEvents_IOS(EventData.GCaMP7s.RH.stim,StimCriteriaB);
[LH_stimGCaMP7sData] = EventData.GCaMP7s.LH.stim.data(LH_stimFilter,:);
[RH_stimGCaMP7sData] = EventData.GCaMP7s.RH.stim.data(RH_stimFilter,:);
[LH_stimFileIDs] = EventData.GCaMP7s.LH.stim.fileIDs(LH_stimFilter,:);
[RH_stimFileIDs] = EventData.GCaMP7s.RH.stim.fileIDs(RH_stimFilter,:);
[LH_stimEventTimes] = EventData.GCaMP7s.LH.stim.eventTime(LH_stimFilter,:);
[RH_stimEventTimes] = EventData.GCaMP7s.RH.stim.eventTime(RH_stimFilter,:);
LH_stimDurations = zeros(length(LH_stimEventTimes),1);
RH_stimDurations = zeros(length(RH_stimEventTimes),1);
% keep only the data that occurs within the manually-approved awake regions
[LH_finalStimData,LH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(LH_stimGCaMP7sData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
[RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(RH_stimGCaMP7sData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
% filter [GCaMP7s] and mean-subtract 2 seconds prior to stimulus (left hem)
for gg = 1:size(LH_finalStimData,1)
    LH_ProcStimData_temp = filtfilt(sos,g,LH_finalStimData(gg,:));
    LH_ProcStimData(gg,:) = LH_ProcStimData_temp - mean(LH_ProcStimData_temp(1:params.Offset*samplingRate));
end
% filter [GCaMP7s] and mean-subtract 2 seconds prior to stimulus (right hem)
for gg = 1:size(RH_finalStimData,1)
    RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
    RH_ProcStimData(gg,:) = RH_ProcStimData_temp - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
end
% take mean [GCaMP7s] 1-2 seconds after stimulation (left hem)
for nn = 1:size(LH_ProcStimData,1)
    LH_stimGCaMP7sMean{nn,1} = mean(LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    LH_stimGCaMP7s{nn,1} = LH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
end
% take mean [GCaMP7s] 1-2 seconds after stimulation (right hem)
for nn = 1:size(RH_ProcStimData,1)
    RH_stimGCaMP7sMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    RH_stimGCaMP7s{nn,1} = RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.MeanLH = cell2mat(LH_stimGCaMP7sMean);
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.MeanRH = cell2mat(RH_stimGCaMP7sMean);
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.IndLH = LH_stimGCaMP7s;
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.IndRH = RH_stimGCaMP7s;
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.LH_FileIDs = LH_finalStimFileIDs;
AnalysisResults.(animalID).MeanGCaMP7s.Stim.GCaMP7s.RH_FileIDs = RH_finalStimFileIDs;
%% analyze [GCaMP7s] during periods of NREM sleep
% pull data from SleepData.mat structure
[LHremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GCaMP7s.LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
[RHremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GCaMP7s.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
% filter and take mean [GCaMP7s] during NREM epochs
for nn = 1:length(LHremData)
    LHremGCaMP7sMean(nn,1) = mean(filtfilt(sos,g,LHremData{nn,1}(1:end)));
    RHremGCaMP7sMean(nn,1) = mean(filtfilt(sos,g,RHremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanGCaMP7s.NREM.GCaMP7s.MeanLH = LHremGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.NREM.GCaMP7s.MeanRH = RHremGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.NREM.GCaMP7s.IndLH = LHremData;
AnalysisResults.(animalID).MeanGCaMP7s.NREM.GCaMP7s.IndRH = RHremData;
AnalysisResults.(animalID).MeanGCaMP7s.NREM.GCaMP7s.FileIDs = nremFileIDs;
%% analyze [GCaMP7s] during periods of REM sleep
% pull data from SleepData.mat structure
[LH_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GCaMP7s.LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
[RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GCaMP7s.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
% filter and take mean [GCaMP7s] during REM epochs
for nn = 1:length(LH_remData)
    LH_remGCaMP7sMean(nn,1) = mean(filtfilt(sos,g,LH_remData{nn,1}(1:end)));
    RH_remGCaMP7sMean(nn,1) = mean(filtfilt(sos,g,RH_remData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanGCaMP7s.REM.GCaMP7s.MeanLH = LH_remGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.REM.GCaMP7s.MeanRH = RH_remGCaMP7sMean;
AnalysisResults.(animalID).MeanGCaMP7s.REM.GCaMP7s.IndLH = LH_remData;
AnalysisResults.(animalID).MeanGCaMP7s.REM.GCaMP7s.IndRH = RH_remData;
AnalysisResults.(animalID).MeanGCaMP7s.REM.GCaMP7s.FileIDs = remFileIDs;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
