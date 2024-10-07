function [AnalysisResults] = AnalyzeMeanPupilDiameter_FP_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the hemodynamic signal  during different arousal states (IOS)
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
samplingRate = RestData.Pupil.zDiameter.CBVCamSamplingRate;
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
%% analyze  during periods of rest
% pull data from RestData.mat structure
[restLogical] = FilterEvents_IOS(RestData.Pupil.zDiameter,RestCriteria);
[puffLogical] = FilterEvents_IOS(RestData.Pupil.zDiameter,RestPuffCriteria);
combRestLogical = logical(restLogical.*puffLogical);
restFileIDs = RestData.Pupil.zDiameter.fileIDs(combRestLogical,:);
restEventTimes = RestData.Pupil.zDiameter.eventTimes(combRestLogical,:);
restDurations = RestData.Pupil.zDiameter.durations(combRestLogical,:);
Pupil_RestingData = RestData.Pupil.zDiameter.data(combRestLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[Pupil_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS(Pupil_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
% filter 
for gg = 1:length(Pupil_finalRestData)
    Pupil_ProcRestData{gg,1} = filtfilt(sos,g,Pupil_finalRestData{gg,1}); %#ok<*AGROW>
end
% take mean  during resting epochs
for nn = 1:length(Pupil_ProcRestData)
    Pupil_restPupilMean(nn,1) = mean(Pupil_ProcRestData{nn,1}(1:end));
end
% save results
AnalysisResults.(animalID).MeanPupil.Rest.Pupil.MeanzDiameter = Pupil_restPupilMean;
AnalysisResults.(animalID).MeanPupil.Rest.Pupil.IndzDiameter = Pupil_ProcRestData;
AnalysisResults.(animalID).MeanPupil.Rest.Pupil.FileIDs = finalRestFileIDs;
%% analyze  during periods of moderate whisking (2-5 seconds)
% pull data from EventData.mat structure
[whiskLogical] = FilterEvents_IOS(EventData.Pupil.zDiameter.whisk,WhiskCriteria);
[puffLogical] = FilterEvents_IOS(EventData.Pupil.zDiameter.whisk,WhiskPuffCriteria);
combWhiskLogical = logical(whiskLogical.*puffLogical);
whiskFileIDs = EventData.Pupil.zDiameter.whisk.fileIDs(combWhiskLogical,:);
whiskEventTimes = EventData.Pupil.zDiameter.whisk.eventTime(combWhiskLogical,:);
whiskDurations = EventData.Pupil.zDiameter.whisk.duration(combWhiskLogical,:);
Pupil_whiskData = EventData.Pupil.zDiameter.whisk.data(combWhiskLogical,:);
% keep only the data that occurs within the manually-approved awake regions
[Pupil_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS(Pupil_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
% filter  and mean-subtract 2 seconds prior to whisk
for gg = 1:size(Pupil_finalWhiskData,1)
    Pupil_ProcWhiskData_temp = filtfilt(sos,g,Pupil_finalWhiskData(gg,:));
    Pupil_ProcWhiskData(gg,:) = Pupil_ProcWhiskData_temp;% - mean(Pupil_ProcWhiskData_temp(1:params.Offset*samplingRate));
end
% take mean  during whisking epochs from onset through 5 seconds
for nn = 1:size(Pupil_ProcWhiskData,1)
    Pupil_whiskPupilMean(nn,1) = mean(Pupil_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
    Pupil_whiskPupil{nn,1} = Pupil_ProcWhiskData(nn,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
end
% save results
AnalysisResults.(animalID).MeanPupil.Whisk.Pupil.MeanzDiameter = (Pupil_whiskPupilMean);
AnalysisResults.(animalID).MeanPupil.Whisk.Pupil.IndzDiameter = Pupil_whiskPupil;
AnalysisResults.(animalID).MeanPupil.Whisk.Pupil.FileIDs = finalWhiskFileIDs;
%% analyze  during periods of stimulation
% pull data from EventData.mat structure
    if firstHrs == "true"
        Pupil_stimFilter = FilterEvents_IOS(EventData.Pupil.zDiameter.stim,StimCriteriaA);
        [Pupil_stimPupilData] = EventData.Pupil.zDiameter.stim.data(Pupil_stimFilter,:);
        [Pupil_stimFileIDs] = EventData.Pupil.zDiameter.stim.fileIDs(Pupil_stimFilter,:);
        [Pupil_stimEventTimes] = EventData.Pupil.zDiameter.stim.eventTime(Pupil_stimFilter,:);
        Pupil_stimDurations = zeros(length(Pupil_stimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [Pupil_finalStimData,Pupil_finalStimFileIDs,~,~] = RemoveInvalidData_IOS(Pupil_stimPupilData,Pupil_stimFileIDs,Pupil_stimDurations,Pupil_stimEventTimes,ManualDecisions);
        % filter  and mean-subtract 2 seconds prior to stimulus (left hem)
        for gg = 1:size(Pupil_finalStimData,1)
            Pupil_ProcStimData_temp = filtfilt(sos,g,Pupil_finalStimData(gg,:));
            Pupil_ProcStimData(gg,:) = Pupil_ProcStimData_temp;% - mean(Pupil_ProcStimData_temp(1:params.Offset*samplingRate));
        end
       
        % take mean  1-2 seconds after stimulation (left hem)
        for nn = 1:size(Pupil_ProcStimData,1)
            Pupil_stimPupilMean{nn,1} = mean(Pupil_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate),2);
            Pupil_stimPupil{nn,1} = Pupil_ProcStimData(nn,(params.Offset + 0)*samplingRate:params.minTime.Stim*samplingRate);
        end
       
        % save results
        AnalysisResults.(animalID).MeanPupil.Stim.Pupil.MeanzDiameter = cell2mat(Pupil_stimPupilMean);
        AnalysisResults.(animalID).MeanPupil.Stim.Pupil.IndzDiameter = Pupil_stimPupil;
        AnalysisResults.(animalID).MeanPupil.Stim.Pupil.Pupil_FileIDs = Pupil_finalStimFileIDs;
    end
%% analyze  during periods of NREM sleep
% pull data from SleepData.mat structure
    if firstHrs == "true"
    [zDiameterremData,nremFileIDs,~] = NRemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.zDiameter,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif firstHrs == "false"
    [zDiameterremData,nremFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.zDiameter,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
% filter and take mean  during NREM epochs
for nn = 1:length(zDiameterremData)
    zDiameterremPupilMean(nn,1) = mean(filtfilt(sos,g,zDiameterremData{nn,1}(1:end)));
end
% save results
AnalysisResults.(animalID).MeanPupil.NREM.Pupil.MeanzDiameter = zDiameterremPupilMean;
AnalysisResults.(animalID).MeanPupil.NREM.Pupil.IndzDiameter = zDiameterremData;
AnalysisResults.(animalID).MeanPupil.NREM.Pupil.FileIDs = nremFileIDs;
%% analyze  during periods of REM sleep
% pull data from SleepData.mat structure
if firstHrs == "false"
    [Pupil_remData,remFileIDs,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.zDiameter,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    % filter and take mean  during REM epochs
    for nn = 1:length(Pupil_remData)
        Pupil_remPupilMean(nn,1) = mean(filtfilt(sos,g,Pupil_remData{nn,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).MeanPupil.REM.Pupil.MeanzDiameter = Pupil_remPupilMean;
    AnalysisResults.(animalID).MeanPupil.REM.Pupil.IndzDiameter = Pupil_remData;
    AnalysisResults.(animalID).MeanPupil.REM.Pupil.FileIDs = remFileIDs;
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
