function [AnalysisResults] = AnalyzeMeanGFP_FP_EEG(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [GFP] during different arousal states (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
modelType = 'Forest';
params.minTime.Rest = 10;
params.Offset = 15; 
params.minTime.Whisk = params.Offset + 2;
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
samplingRate = RestData.GFP.RH.RhodamineCamSamplingRate;
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
[restLogical] = FilterEvents_IOS(RestData.GFP.RH,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.GFP.RH,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.GFP.RH.fileIDs(combRestLogical,:);
restEventTimes = RestData.GFP.RH.eventTimes(combRestLogical,:);
restDurations = RestData.GFP.RH.durations(combRestLogical,:);
RH_RestingData = RestData.GFP.RH.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter [GFP]
for gg = 1:length(RH_finalRestData)
    RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1}); %#ok<*AGROW>
end
% take mean [GFP] during resting epochs
for nn = 1:length(RH_ProcRestData)
    RH_restGFPMean(nn,1) = mean(RH_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanGFP.Rest.GFP.MeanRH = RH_restGFPMean;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.IndRH = RH_ProcRestData;
AnalysisResults.(animalID).MeanGFP.Rest.GFP.FileIDs = finalRestFileIDs;
%% analyze [GFP] during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.GFP.RH.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.GFP.RH.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.GFP.RH.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.GFP.RH.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.GFP.RH.whisk.duration(combWhiskLogical,:);
RH_whiskData = EventData.GFP.RH.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter [GFP] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(RH_finalWhiskData,1)
    RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
    RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp;% - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean [GFP] during whisking epochs from onset through 5 seconds
for nn = 1:size(RH_ProcWhiskData,1)
    RH_whiskGFPMean{nn,1} = mean(RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    RH_whiskGFP{nn,1} = RH_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.MeanRH = cell2mat(RH_whiskGFPMean);
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.IndRH = RH_whiskGFP;
AnalysisResults.(animalID).MeanGFP.Whisk.GFP.FileIDs = finalWhiskFileIDs;
%% analyze [GFP] during periods of stimulation
% pull data from EventData.mat structure
RH_stimFilter = FilterEvents_IOS(EventData.GFP.RH.stim,StimCriteriaA);
[RH_stimGFPData] = EventData.GFP.RH.stim.data(RH_stimFilter,:);
[RH_stimFileIDs] = EventData.GFP.RH.stim.fileIDs(RH_stimFilter,:);
[RH_stimEventTimes] = EventData.GFP.RH.stim.eventTime(RH_stimFilter,:);
RH_stimDurations = zeros(length(RH_stimEventTimes),1);
% keep only the data that occurs within the manually-approved awake regions
[RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(RH_stimGFPData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
% filter [GFP] and mean-subtract 2 seconds prior to stimulus (left hem)
for gg = 1:size(RH_finalStimData,1)
    RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
    RH_ProcStimData(gg,:) = RH_ProcStimData_temp;% - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
end

% take mean [GFP] 1-2 seconds after stimulation (left hem)
for nn = 1:size(RH_ProcStimData,1)
    RH_stimGFPMean{nn,1} = mean(RH_ProcStimData(nn,(params.Offset + 2)*samplingRate:params.minTime.Stim*samplingRate),2);
    RH_stimGFP{nn,1} = RH_ProcStimData(nn,(params.Offset + 2)*samplingRate:params.minTime.Stim*samplingRate);
end

% save results
AnalysisResults.(animalID).MeanGFP.Stim.GFP.MeanRH = cell2mat(RH_stimGFPMean);
AnalysisResults.(animalID).MeanGFP.Stim.GFP.IndRH = RH_stimGFP;
AnalysisResults.(animalID).MeanGFP.Stim.GFP.RH_FileIDs = RH_finalStimFileIDs;
%% analyze [GFP] during periods of NREM sleep
% pull data from SleepData.mat structure
[RHremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
% filter and take mean [GFP] during NREM epochs
for nn = 1:length(RHremData)
    RHremGFPMean(nn,1) = mean(filtfilt(sos,g,RHremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanGFP.NREM.GFP.MeanRH = RHremGFPMean;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.IndRH = RHremData;
AnalysisResults.(animalID).MeanGFP.NREM.GFP.FileIDs = nremFileIDs;
%% analyze [GFP] during periods of REM sleep
% pull data from SleepData.mat structure
[RH_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
% filter and take mean [GFP] during REM epochs
for nn = 1:length(RH_remData)
    RH_remGFPMean(nn,1) = mean(filtfilt(sos,g,RH_remData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanGFP.REM.GFP.MeanRH = RH_remGFPMean;
AnalysisResults.(animalID).MeanGFP.REM.GFP.IndRH = RH_remData;
AnalysisResults.(animalID).MeanGFP.REM.GFP.FileIDs = remFileIDs;
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
