function [AnalysisResults] = AnalyzeMeanCBV_FP_GRABNE_Sleep(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal [CBV] during different arousal states
%________________________________________________________________________________________________________________________
%% function parameters
modelType = 'Manual';
params.minTime.Rest = 10;
params.Offset = 15;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 15;
params.minTime.NREM = 30; %30
params.minTime.REM = 60; %30
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

try
    samplingRate = RestData.CBV.P_ACh.CBVCamSamplingRate;
catch
    samplingRate = RestData.CBV.P_ACh.RhodamineCamSamplingRate;
end

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
%% analyze [CBV] during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.CBV.P_ACh,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.CBV.P_ACh,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.CBV.P_ACh.fileIDs(combRestLogical,:);
restEventTimes = RestData.CBV.P_ACh.eventTimes(combRestLogical,:);
restDurations = RestData.CBV.P_ACh.durations(combRestLogical,:);
ACh_RestingData = RestData.CBV.P_ACh.data(combRestLogical,:);
NE_RestingData = RestData.CBV.P_NE.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[ACh_finalRestData_OLD,finalRestFileIDs_OLD,finalRestDuration_OLD,finalRestEventTimes_OLD] = RemoveInvalidData_IOS(ACh_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
[NE_finalRestData_OLD,~,~,~] = RemoveInvalidData_IOS(NE_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% keep only the data that doesn't have sleep
SleepFileIDs = [];
SleepBinTimes = [];
[ACh_finalRestData,finalRestFileIDs,~,~] = RemoveSleepData(animalID,ACh_finalRestData_OLD,finalRestFileIDs_OLD,finalRestDuration_OLD,finalRestEventTimes_OLD,SleepFileIDs,SleepBinTimes);
[NE_finalRestData,~,~,~] = RemoveSleepData(animalID,NE_finalRestData_OLD,finalRestFileIDs_OLD,finalRestDuration_OLD,finalRestEventTimes_OLD,SleepFileIDs,SleepBinTimes);
% filter [CBV]
for gg = 1:length(ACh_finalRestData)
    ACh_ProcRestData{gg,1} = filtfilt(sos,g,ACh_finalRestData{gg,1}); %#ok<*AGROW>
    NE_ProcRestData{gg,1} = filtfilt(sos,g,NE_finalRestData{gg,1});
end
% take mean [CBV] during resting epochs
for nn = 1:length(ACh_ProcRestData)
    ACh_restCBVMean(nn,1) = mean(ACh_ProcRestData{nn,1}(1:end));
    NE_restCBVMean(nn,1) = mean(NE_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanCBV.Rest.CBV.MeanACh = ACh_restCBVMean;
AnalysisResults.(animalID).MeanCBV.Rest.CBV.MeanNE = NE_restCBVMean;
AnalysisResults.(animalID).MeanCBV.Rest.CBV.IndACh = ACh_ProcRestData;
AnalysisResults.(animalID).MeanCBV.Rest.CBV.IndNE = NE_ProcRestData;
AnalysisResults.(animalID).MeanCBV.Rest.CBV.FileIDs = finalRestFileIDs;
%% analyze [CBV] during periods of brief whisking (<2.5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.CBV.P_ACh.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.CBV.P_ACh.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.CBV.P_ACh.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.CBV.P_ACh.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.CBV.P_ACh.whisk.duration(combWhiskLogical,:);
ACh_whiskData = EventData.CBV.P_ACh.whisk.data(combWhiskLogical,:);
NE_whiskData = EventData.CBV.P_NE.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[ACh_finalWhiskData_OLD,finalWhiskFileIDs_OLD,finalWhiskDuration_OLD,finalWhiskEventTimes_OLD] = RemoveInvalidData_IOS(ACh_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
[NE_finalWhiskData_OLD,~,~,~] = RemoveInvalidData_IOS(NE_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% keep only the data that doesn't have sleep
SleepFileIDs = [];
SleepBinTimes = [];
[ACh_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveSleepData(animalID,ACh_finalWhiskData_OLD,finalWhiskFileIDs_OLD,finalWhiskDuration_OLD,finalWhiskEventTimes_OLD,SleepFileIDs,SleepBinTimes);
[NE_finalWhiskData,~,~,~] = RemoveSleepData(animalID,NE_finalWhiskData_OLD,finalWhiskFileIDs_OLD,finalWhiskDuration_OLD,finalWhiskEventTimes_OLD,SleepFileIDs,SleepBinTimes);


% filter [CBV] and mean-subtract 2 seconds prior to whisk
for gg = 1:size(ACh_finalWhiskData,1)
    ACh_ProcWhiskData_temp = filtfilt(sos,g,ACh_finalWhiskData(gg,:));
    ACh_ProcWhiskData(gg,:) = ACh_ProcWhiskData_temp;%  - mean(ACh_ProcWhiskData_temp(1:params.Offset*samplingRate));
    NE_ProcWhiskData_temp = filtfilt(sos,g,NE_finalWhiskData(gg,:));
    NE_ProcWhiskData(gg,:) = NE_ProcWhiskData_temp;%  - mean(NE_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean [CBV] during whisking epochs from onset through 5 seconds
for nn = 1:size(ACh_ProcWhiskData,1)
    ACh_whiskCBVMean(nn,1) = mean(ACh_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2); % mean
    NE_whiskCBVMean(nn,1) = mean(NE_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);  % mean
end

for nn = 1:size(ACh_ProcWhiskData,1)
    ACh_whiskCBV{nn,1} = ACh_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    NE_whiskCBV{nn,1} = NE_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanCBV.Whisk.CBV.MeanACh = ACh_whiskCBVMean;
AnalysisResults.(animalID).MeanCBV.Whisk.CBV.MeanNE = NE_whiskCBVMean;
AnalysisResults.(animalID).MeanCBV.Whisk.CBV.IndACh = ACh_whiskCBV;
AnalysisResults.(animalID).MeanCBV.Whisk.CBV.IndNE = NE_whiskCBV;
AnalysisResults.(animalID).MeanCBV.Whisk.CBV.FileIDs = finalWhiskFileIDs;
%% analyze [CBV] during periods of stimulation
% if firstHrs == "true"
    % pull data from EventData.mat structure
    ACh_stimFilter = FilterEvents_IOS(EventData.CBV.P_ACh.stim,StimCriteriaA);
    NE_stimFilter = FilterEvents_IOS(EventData.CBV.P_NE.stim,StimCriteriaB);
    [ACh_stimData] = EventData.CBV.P_ACh.stim.data(ACh_stimFilter,:);
    [NE_stimData] = EventData.CBV.P_NE.stim.data(NE_stimFilter,:);
    [ACh_stimFileIDs] = EventData.CBV.P_ACh.stim.fileIDs(ACh_stimFilter,:);
    [NE_stimFileIDs] = EventData.CBV.P_NE.stim.fileIDs(NE_stimFilter,:);
    [ACh_stimEventTimes] = EventData.CBV.P_ACh.stim.eventTime(ACh_stimFilter,:);
    [NE_stimEventTimes] = EventData.CBV.P_NE.stim.eventTime(NE_stimFilter,:);
    ACh_stimDurations = zeros(length(ACh_stimEventTimes),1);
    NE_stimDurations = zeros(length(NE_stimEventTimes),1);
    % keep only the data that occurs within the manually-approved awake regions
    [ACh_finalStimData_OLD,ACh_finalStimFileIDs_OLD,ACh_finalStimDuration_OLD,ACh_finalStimEventTimes_OLD] = RemoveInvalidData_IOS(ACh_stimData,ACh_stimFileIDs,ACh_stimDurations,ACh_stimEventTimes,ManualDecisions);
    [NE_finalStimData_OLD,NE_finalStimFileIDs_OLD,NE_finalStimDuration_OLD,NE_finalStimEventTimes_OLD] = RemoveInvalidData_IOS(NE_stimData,NE_stimFileIDs,NE_stimDurations,NE_stimEventTimes,ManualDecisions);
    % keep only the data that doesn't have sleep
    SleepFileIDs = [];
    SleepBinTimes = [];
    [ACh_finalStimData,ACh_finalStimFileIDs,~,~] = RemoveSleepData(animalID,ACh_finalStimData_OLD,ACh_finalStimFileIDs_OLD,ACh_finalStimDuration_OLD,ACh_finalStimEventTimes_OLD,SleepFileIDs,SleepBinTimes);
    [NE_finalStimData,NE_finalStimFileIDs,~,~] = RemoveSleepData(animalID,NE_finalStimData_OLD,NE_finalStimFileIDs_OLD,NE_finalStimDuration_OLD,NE_finalStimEventTimes_OLD,SleepFileIDs,SleepBinTimes);
    % filter [CBV] and mean-subtract 2 seconds prior to stimulus (left hem)
    for gg = 1:size(ACh_finalStimData,1)
        ACh_ProcStimData_temp = filtfilt(sos,g,ACh_finalStimData(gg,:));
        ACh_ProcStimData(gg,:) = ACh_ProcStimData_temp;% - mean(ACh_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    % filter [CBV] and mean-subtract 2 seconds prior to stimulus (right hem)
    for gg = 1:size(NE_finalStimData,1)
        NE_ProcStimData_temp = filtfilt(sos,g,NE_finalStimData(gg,:));
        NE_ProcStimData(gg,:) = NE_ProcStimData_temp;% - mean(NE_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    % take mean [CBV] 1-2 seconds after stimulation (left hem)
    for nn = 1:size(ACh_ProcStimData,1)
        ACh_stimCBVMean{nn,1} = mean(ACh_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        ACh_stimCBV{nn,1} = ACh_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    % take mean [CBV] 1-2 seconds after stimulation (right hem)
    for nn = 1:size(NE_ProcStimData,1)
        NE_stimCBVMean{nn,1} = mean(NE_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        NE_stimCBV{nn,1} = NE_ProcStimData(nn,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.MeanACh = cell2mat(ACh_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.MeanNE = cell2mat(NE_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.IndACh = ACh_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.IndNE = NE_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.ACh_FileIDs = ACh_finalStimFileIDs;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV.NE_FileIDs = NE_finalStimFileIDs;
% end
%% analyze [CBV] during periods of NREM sleep
% pull data from SleepData.mat structure
if firstHrs == "true"
    [AChremData,nremFileIDs,~] = NRemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NEremData,~,~] = NRemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
elseif firstHrs == "false"
    [AChremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NEremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
end
% filter and take mean [CBV] during NREM epochs
for nn = 1:length(AChremData)
    AChremCBVMean(nn,1) = mean(filtfilt(sos,g,AChremData{nn,1}(1:end)));
    NEremCBVMean(nn,1) = mean(filtfilt(sos,g,NEremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanCBV.NREM.CBV.MeanACh = AChremCBVMean;
AnalysisResults.(animalID).MeanCBV.NREM.CBV.MeanNE = NEremCBVMean;
AnalysisResults.(animalID).MeanCBV.NREM.CBV.IndACh = AChremData;
AnalysisResults.(animalID).MeanCBV.NREM.CBV.IndNE = NEremData;
AnalysisResults.(animalID).MeanCBV.NREM.CBV.FileIDs = nremFileIDs;
%% analyze [CBV] during periods of REM sleep

if isfield(SleepData.(modelType),'REM')
    % pull data from SleepData.mat structure
    [ACh_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [NE_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter and take mean [CBV] during REM epochs
    for nn = 1:length(ACh_remData)
        ACh_remCBVMean(nn,1) = mean(filtfilt(sos,g,ACh_remData{nn,1}(1:end)));
        NE_remCBVMean(nn,1) = mean(filtfilt(sos,g,NE_remData{nn,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.REM.CBV.MeanACh = ACh_remCBVMean;
    AnalysisResults.(animalID).MeanCBV.REM.CBV.MeanNE = NE_remCBVMean;
    AnalysisResults.(animalID).MeanCBV.REM.CBV.IndACh = ACh_remData;
    AnalysisResults.(animalID).MeanCBV.REM.CBV.IndNE = NE_remData;
    AnalysisResults.(animalID).MeanCBV.REM.CBV.FileIDs = remFileIDs; 
end
%}
%% save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end

end
