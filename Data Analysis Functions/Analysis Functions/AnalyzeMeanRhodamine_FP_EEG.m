function [AnalysisResults] = AnalyzeMeanRhodamine_FP_EEG(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [Rhodamine] during different arousal states
%________________________________________________________________________________________________________________________
%% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 15;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;
params.minTime.REM = 30;
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
samplingRate = RestData.Rhodamine.RH.RhodamineCamSamplingRate;
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
%% analyze [Rhodamine] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.Rhodamine.RH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.Rhodamine.RH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.Rhodamine.RH.fileIDs(combRestLogical,:);
restEventTimes = RestData.Rhodamine.RH.eventTimes(combRestLogical,:);
restDurations = RestData.Rhodamine.RH.durations(combRestLogical,:);
RH_RestingData = RestData.Rhodamine.RH.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter [Rhodamine]
for gg = 1:length(RH_finalRestData)
    RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1}); %#ok<*AGROW>
end
% take mean [Rhodamine] during resting epochs
for nn = 1:length(RH_ProcRestData)
    RH_restRhodamineMean(nn,1) = mean(RH_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanRhodamine.Rest.Rhodamine.MeanRH = RH_restRhodamineMean;
AnalysisResults.(animalID).MeanRhodamine.Rest.Rhodamine.IndRH = RH_ProcRestData;
AnalysisResults.(animalID).MeanRhodamine.Rest.Rhodamine.FileIDs = finalRestFileIDs;
%% analyze [Rhodamine] during periods of extended whisking (>5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.Rhodamine.RH.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.Rhodamine.RH.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.Rhodamine.RH.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.Rhodamine.RH.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.Rhodamine.RH.whisk.duration(combWhiskLogical,:);
RH_whiskData = EventData.Rhodamine.RH.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter [Rhodamine] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(RH_finalWhiskData,1)
    RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
    RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp;%  - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean [Rhodamine] during whisking epochs from onset through 5 seconds
for nn = 1:size(RH_ProcWhiskData,1)
    RH_whiskRhodamineMean{nn,1} = mean(RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2); % mean
    RH_whiskRhodamine{nn,1} = RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanRhodamine.Whisk.Rhodamine.MeanRH = cell2mat(RH_whiskRhodamineMean);
AnalysisResults.(animalID).MeanRhodamine.Whisk.Rhodamine.IndRH = RH_whiskRhodamine;
AnalysisResults.(animalID).MeanRhodamine.Whisk.Rhodamine.FileIDs = finalWhiskFileIDs;
%% analyze [Rhodamine] during periods of stimulation
% pull data from EventData.mat structure
RH_stimFilter = FilterEvents_IOS(EventData.Rhodamine.RH.stim,StimCriteriaA);
[RH_stimRhodamineData] = EventData.Rhodamine.RH.stim.data(RH_stimFilter,:);
[RH_stimFileIDs] = EventData.Rhodamine.RH.stim.fileIDs(RH_stimFilter,:);
[RH_stimEventTimes] = EventData.Rhodamine.RH.stim.eventTime(RH_stimFilter,:);
RH_stimDurations = zeros(length(RH_stimEventTimes),1);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(RH_stimRhodamineData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
% filter [Rhodamine] and mean-subtract 2 seconds prior to stimulus (left hem)
for gg = 1:size(RH_finalStimData,1)
    RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
    RH_ProcStimData(gg,:) = RH_ProcStimData_temp;% - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
end

% take mean [Rhodamine] 1-2 seconds after stimulation (left hem)
for nn = 1:size(RH_ProcStimData,1)
    RH_stimRhodamineMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
    RH_stimRhodamine{nn,1} = RH_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
end

% save results
AnalysisResults.(animalID).MeanRhodamine.Stim.Rhodamine.MeanRH = cell2mat(RH_stimRhodamineMean);
AnalysisResults.(animalID).MeanRhodamine.Stim.Rhodamine.IndRH = RH_stimRhodamine;
AnalysisResults.(animalID).MeanRhodamine.Stim.Rhodamine.RH_FileIDs = RH_finalStimFileIDs;
%% analyze [Rhodamine] during periods of NREM sleep
% pull data from SleepData.mat structure
[RHremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
% filter and take mean [Rhodamine] during NREM epochs
for nn = 1:length(RHremData)
    RHremRhodamineMean(nn,1) = mean(filtfilt(sos,g,RHremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanRhodamine.NREM.Rhodamine.MeanRH = RHremRhodamineMean;
AnalysisResults.(animalID).MeanRhodamine.NREM.Rhodamine.IndRH = RHremData;
AnalysisResults.(animalID).MeanRhodamine.NREM.Rhodamine.FileIDs = nremFileIDs;
%% analyze [Rhodamine] during periods of REM sleep
% pull data from SleepData.mat structure
[RH_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
% filter and take mean [Rhodamine] during REM epochs
for nn = 1:length(RH_remData)
    RH_remRhodamineMean(nn,1) = mean(filtfilt(sos,g,RH_remData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanRhodamine.REM.Rhodamine.MeanRH = RH_remRhodamineMean;
AnalysisResults.(animalID).MeanRhodamine.REM.Rhodamine.IndRH = RH_remData;
AnalysisResults.(animalID).MeanRhodamine.REM.Rhodamine.FileIDs = remFileIDs;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
