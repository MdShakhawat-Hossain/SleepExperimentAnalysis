function [AnalysisResults] = AnalyzeCorrCoeffs(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze Pearson's correlation coefficient between bilateral
%   hemodynamic [TRITC] and neural signals (IOS) and between bilateral GCaMP7s and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'TRITC','GCaMP7s','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.Whisk = 5;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
cd(dataLocation)
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
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
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% lowpass filter
samplingRate = RestData.TRITC.LH.TRITCCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for whisking
WhiskCriteria.Fieldname = {'duration','puffDistance'};
WhiskCriteria.Comparison = {'gt','gt'};
WhiskCriteria.Value = {5,5};
WhiskPuffCriteria.Fieldname = {'puffDistance'};
WhiskPuffCriteria.Comparison = {'gt'};
WhiskPuffCriteria.Value = {5};
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal-based correlation analysis
for a = 1:length(dataTypes)
    dataType = dataTypes{1,a};
    %% analyze Pearson's correlation coefficient during periods of rest
    % pull data from RestData.mat structure
    if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
        [restLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.(dataType).LH,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.(dataType).LH.fileIDs(combRestLogical,:);
        restEventTimes = RestData.(dataType).LH.eventTimes(combRestLogical,:);
        restDurations = RestData.(dataType).LH.durations(combRestLogical,:);
        LH_unstimRestingData = RestData.(dataType).LH.data(combRestLogical,:);
        RH_unstimRestingData = RestData.(dataType).RH.data(combRestLogical,:);
    else
        [restLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.cortical_LH.(dataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
        restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
        LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
        RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
    end
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    clear LH_ProcRestData RH_ProcRestData
    % filter, detrend, and truncate data to minimum length to match events
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,LH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant'); %#ok<*AGROW>
        RH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,RH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
    end
    % analyze correlation coefficient of resting epochs
    for n = 1:length(LH_ProcRestData)
        rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
        rest_R(n,1) = rest_CC(2,1);
    end
    meanRest_R = mean(rest_R);
    stdRest_R = std(rest_R,0,1);
    % save results
    AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).R = rest_R;
    AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).meanR = meanRest_R;
    AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).stdR = stdRest_R;
    %% analyze Pearson's correlation coefficient during periods of moderate whisking (2-5 seconds)
    % pull data from EventData.mat structure
    if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
        [whiskLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.(dataType).LH.whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFileIDs = EventData.(dataType).LH.whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.(dataType).LH.whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.(dataType).LH.whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.(dataType).LH.whisk.data(combWhiskLogical,:);
        RH_whiskData = EventData.(dataType).RH.whisk.data(combWhiskLogical,:);
   else
        [whiskLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.cortical_LH.(dataType).whisk,WhiskPuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        whiskFileIDs = EventData.cortical_LH.(dataType).whisk.fileIDs(combWhiskLogical,:);
        whiskEventTimes = EventData.cortical_LH.(dataType).whisk.eventTime(combWhiskLogical,:);
        whiskDurations = EventData.cortical_LH.(dataType).whisk.duration(combWhiskLogical,:);
        LH_whiskData = EventData.cortical_LH.(dataType).whisk.NormData(combWhiskLogical,:);
        RH_whiskData = EventData.cortical_RH.(dataType).whisk.NormData(combWhiskLogical,:);
    end
    % keep only the data that occurs within the manually-approved awake regions
    [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    clear LH_ProcWhiskData RH_ProcWhiskData
    % filter, detrend, and take data from whisk onset through 5 seconds
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        RH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
    end
    % analyze correlation coefficient of whisking epochs
    for n = 1:size(LH_ProcWhiskData,1)
        whisk_CC = corrcoef(LH_ProcWhiskData(n,:),RH_ProcWhiskData(n,:));
        whisk_R(n,1) = whisk_CC(2,1);
    end
    meanWhisk_R = mean(whisk_R);
    stdWhisk_R = std(whisk_R,0,1);
    % save results
    AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).R = whisk_R;
    AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).meanR = meanWhisk_R;
    AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).stdR = stdWhisk_R;
    %% analyze Pearson's correlation coefficient during periods of alert
    zz = 1;
    clear LH_AwakeData RH_AwakeData LH_ProcAwakeData RH_ProcAwakeData
    LH_AwakeData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 500   % 36 bins (180 total) or 3 minutes of sleep
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
                    LH_AwakeData{zz,1} = (ProcData.data.(dataType).LH - RestingBaselines.manualSelection.(dataType).LH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).LH.(strDay).std;
                    RH_AwakeData{zz,1} = (ProcData.data.(dataType).RH - RestingBaselines.manualSelection.(dataType).RH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).RH.(strDay).std;
                else
                    LH_AwakeData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                    RH_AwakeData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                end
                zz = zz + 1;
            end
        end
    end
    if isempty(LH_AwakeData) == false
        % filter and detrend data
        for gg = 1:length(LH_AwakeData)
            LH_ProcAwakeData{gg,1} = detrend(filtfilt(sos,g,LH_AwakeData{gg,1}),'constant');
            RH_ProcAwakeData{gg,1} = detrend(filtfilt(sos,g,RH_AwakeData{gg,1}),'constant');
        end
        % analyze correlation coefficient of alert epochs
        for n = 1:length(LH_ProcAwakeData)
            awake_CC = corrcoef(LH_ProcAwakeData{n,1},RH_ProcAwakeData{n,1});
            awake_R(n,1) = awake_CC(2,1);
        end
        meanAwake_R = mean(awake_R);
        stdAwake_R = std(awake_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).R = awake_R;
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).meanR = meanAwake_R;
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).stdR = stdAwake_R;
    else
        % save results
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).R = [];
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).meanR = [];
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).stdR = [];
    end
    %% analyze Pearson's correlation coefficient during periods of asleep
    zz = 1;
    clear LH_SleepData RH_SleepData LH_ProcSleepData RH_ProcSleepData
    LH_SleepData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 124   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
                    LH_SleepData{zz,1} = (ProcData.data.(dataType).LH - RestingBaselines.manualSelection.(dataType).LH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).LH.(strDay).std;
                    RH_SleepData{zz,1} = (ProcData.data.(dataType).RH - RestingBaselines.manualSelection.(dataType).RH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).RH.(strDay).std;
                else
                    LH_SleepData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                    RH_SleepData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
                end
                zz = zz + 1;
            end
        end
    end
    if isempty(LH_SleepData) == false
        % filter and detrend data
        for gg = 1:length(LH_SleepData)
            LH_ProcSleepData{gg,1} = detrend(filtfilt(sos,g,LH_SleepData{gg,1}),'constant');
            RH_ProcSleepData{gg,1} = detrend(filtfilt(sos,g,RH_SleepData{gg,1}),'constant');
        end
        % analyze correlation coefficient of asleep epochs
        for n = 1:length(LH_ProcSleepData)
            sleep_CC = corrcoef(LH_ProcSleepData{n,1},RH_ProcSleepData{n,1});
            sleep_R(n,1) = sleep_CC(2,1);
        end
        meanSleep_R = mean(sleep_R);
        stdSleep_R = std(sleep_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).R = sleep_R;
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).meanR = meanSleep_R;
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).stdR = stdSleep_R;
    else
        % save results
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).R = [];
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).meanR = [];
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).stdR = [];
    end
    %% analyze Pearson's correlation coefficient during periods of all data
    zz = 1;
    clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData
    LH_AllData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,allDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(allDataFileDate);
        load(procDataFileID,'-mat')
        puffs = ProcData.data.stimulations.LPadSol;
        % don't include trials with stimulation
        if isempty(puffs) == true
            if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
                LH_AllData{zz,1} = (ProcData.data.(dataType).LH - RestingBaselines.manualSelection.(dataType).LH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).LH.(strDay).std;
                RH_AllData{zz,1} = (ProcData.data.(dataType).RH - RestingBaselines.manualSelection.(dataType).RH.(strDay).mean)./RestingBaselines.manualSelection.(dataType).RH.(strDay).std;
            else
                LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay).mean;
                RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay).mean;
            end
            zz = zz + 1;
        end
    end
    % filter and detrend data
    for gg = 1:length(LH_AllData)
        LH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,LH_AllData{gg,1}),'constant');
        RH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,RH_AllData{gg,1}),'constant');
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:length(LH_ProcAllData)
        all_CC = corrcoef(LH_ProcAllData{n,1},RH_ProcAllData{n,1});
        all_R(n,1) = all_CC(2,1);
    end
    meanAll_R = mean(all_R);
    stdAll_R = std(all_R,0,1);
    % save results
    AnalysisResults.(animalID).CorrCoeff.All.(dataType).R = all_R;
    AnalysisResults.(animalID).CorrCoeff.All.(dataType).meanR = meanAll_R;
    AnalysisResults.(animalID).CorrCoeff.All.(dataType).stdR = stdAll_R;
    %% analyze Pearson's correlation coefficient during periods of NREM
    % pull data from SleepData.mat structure
    if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).LH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(dataType).RH,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
   else
        [LH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_LH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [RH_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for j = 1:length(LH_nremData)
        LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
    end
    % analyze correlation coefficient of NREM epochs
    for n = 1:length(LH_nremData)
        nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
        nrem_R(n,1) = nrem_CC(2,1);
    end
    meanNREM_R = mean(nrem_R);
    stdNREM_R = std(nrem_R,0,1);
    % save results
    AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).R = nrem_R;
    AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).meanR = meanNREM_R;
    AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).stdR = stdNREM_R;
    %% analyze Pearson's correlation coefficient during periods of REM
    % pull data from SleepData.mat structure
    if strcmp(dataType,'TRITC') == true || strcmp(dataType,'GCaMP7s') == true
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).LH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(dataType).RH,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
     else
        [LH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_LH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [RH_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    end
    % filter, detrend, and truncate data to data to minimum length to match events
    for m = 1:length(LH_remData)
        LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
    end
    % analyze correlation coefficient of REM epochs
    for n = 1:length(LH_remData)
        rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
        rem_R(n,1) = rem_CC(2,1);
    end
    meanREM_R = mean(rem_R);
    stdREM_R = std(rem_R,0,1);
    % save results
    AnalysisResults.(animalID).CorrCoeff.REM.(dataType).R = rem_R;
    AnalysisResults.(animalID).CorrCoeff.REM.(dataType).meanR = meanREM_R;
    AnalysisResults.(animalID).CorrCoeff.REM.(dataType).stdR = stdREM_R;
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
