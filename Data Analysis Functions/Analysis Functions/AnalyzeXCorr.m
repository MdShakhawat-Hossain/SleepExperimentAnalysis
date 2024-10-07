function [AnalysisResults] = AnalyzeXCorr(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the cross-correlation between neural activity and hemodynamics [TRITC] (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'LH','RH'};
modelType = 'Forest';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
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
% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStructC.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% lowpass filter
samplingRate = RestData.TRITC.LH.TRITCCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};RestData.TRITC
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal-based cross-correlation analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    neuralDataType = ['cortical_' dataType];
    % pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
    trialDuration_sec = RestData.TRITC.LH.trialDuration_sec;
    sleepBinWidth = 5;   % sec
    oneSecSpecFs = 30;   % Hz
    %% cross-correlation analysis for resting data
    % pull data from RestData.mat structure
    [restLogical] = FilterEvents_IOS(RestData.TRITC.(dataType),RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.TRITC.(dataType),RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.TRITC.(dataType).fileIDs(combRestLogical,:);
    restDurations = RestData.TRITC.(dataType).durations(combRestLogical,:);
    restEventTimes = RestData.TRITC.(dataType).eventTimes(combRestLogical,:);
    restingTRITCData = RestData.TRITC.(dataType).data(combRestLogical,:);
    restingMUAData = RestData.(neuralDataType).muaPower.NormData(combRestLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [restFinalRestTRITCData,restFinalFileIDs,restFinalDurations,restFinalEventTimes] = RemoveInvalidData_IOS(restingTRITCData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [restFinalRestMUAData,~,~,~] = RemoveInvalidData_IOS(restingMUAData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    cc = 1;
    for bb = 1:length(restFinalFileIDs)
        restFileID = restFinalFileIDs{bb,1};
        % check whether the event occurs in the appropriate time frame
        restStartTime = ceil(restFinalEventTimes(bb,1)*10)/10; % *10/10 used to round to first decimal place in a floor/ceil fashion.
        restDuration = floor(restFinalDurations(bb,1)*10)/10;
        if restStartTime >= 0.5 && (restStartTime + restDuration) <= (trialDuration_sec - 0.5)
            % remove the number of samples due to rounding up to start and rounding down to end. This is done to keep the TRITC/MUA vectores aligned positionally with the upcoming
            % spectral analysis which is at 10 Hz
            leadSamples = round((restStartTime - restFinalEventTimes(bb,1))*samplingRate);
            lagSamples = round((restFinalDurations(bb,1) - restDuration)*samplingRate);
            % load in TRITC from rest period
            restTRITC = restFinalRestTRITCData{bb,1};
            restMUA = restFinalRestMUAData{bb,1};
            % remove leading/lag samples due to rounding to nearest 0.1 up/0.1 down
            restSnipTRITC = restTRITC(1 + leadSamples:end - lagSamples);
            restSnipMUA = restMUA(1 + leadSamples:end - lagSamples);
            restFiltTRITC = filtfilt(sos,g,detrend(restSnipTRITC,'constant'));
            restFiltMUA = filtfilt(sos,g,detrend(restSnipMUA,'constant'));
            % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
            % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
            if length(restFiltTRITC) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltTRITC);
                restPadTRITC = (ones(1,restChunkSampleDiff))*restFiltTRITC(end);
                restPadMUA = (ones(1,restChunkSampleDiff))*restFiltMUA(end);
                restShortTRITC = horzcat(restFiltTRITC,restPadTRITC);
                restShortMUA = horzcat(restFiltMUA,restPadMUA);
            else
                restShortTRITC = restFiltTRITC(1:params.minTime.Rest*samplingRate);
                restShortMUA = restFiltMUA(1:params.minTime.Rest*samplingRate);
            end
            restProcData.TRITC{cc,1} = restShortTRITC;
            restProcData.MUA{cc,1} = restShortMUA;
            % extract LFP from spectrograms associated with the whisking indecies
            specDataFileID = [animalID '_' restFileID '_SpecDataC.mat'];
            clear S_data
            for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                    rest_S = AllSpecData.(neuralDataType).normS{g,1};
                    rest_T = round(AllSpecData.(neuralDataType).T{g,1},1);
                    rest_F = AllSpecData.(neuralDataType).F{g,1};
                end
            end
            restStartTimeIndex = find(rest_T == restStartTime);
            restStartTimeIndex = restStartTimeIndex(1);
            restDurationIndex = find(rest_T == round((restStartTime + restDuration),1));
            restDurationIndex = restDurationIndex(1);
            restS_Vals = rest_S(:,restStartTimeIndex:restDurationIndex);
            % only take the first min rest time in seconds
            shortRestS_Vals = restS_Vals(:,1:params.minTime.Rest*oneSecSpecFs);
            % mean subtract with detrend and lowpass filter each column
            restProcData.S{cc,1} = detrend(shortRestS_Vals','constant')';
            cc = cc + 1;
        end
        % set parameters for cross-correlation analysis
        restTRITCvLFPzhold = [];
        restLagTime = 5;   % seconds
        restFrequency = oneSecSpecFs;   % Hz
        restMaxLag = restLagTime*restFrequency;
        restTRITCvLFPxcVals = ones(length(rest_F),2*restMaxLag + 1);
        % run cross-correlation analysis - average through time
        for dd = 1:length(restProcData.TRITC)
            for ee = 1:size(restProcData.S{dd, 1}, 1)
                restTRITCarray = restProcData.TRITC{dd,1};
                restMUAarray = restProcData.MUA{dd,1};
                restNeuralArray = restProcData.S{dd,1}(ee,:);
                [restTRITCvLFPxcVals(ee,:),restLFP_lags] = xcorr(restTRITCarray,restNeuralArray,restMaxLag,'coeff');
            end
            [restTRITCvMUAxcVals(dd,:),restMUA_lags] = xcorr(restTRITCarray,restMUAarray,restMaxLag,'coeff'); %#ok<*AGROW>
            restTRITCvLFPzhold = cat(3,restTRITCvLFPzhold,restTRITCvLFPxcVals);
        end
        restMeanTRITCvLFPxcVals = mean(restTRITCvLFPzhold,3);
        restMeanTRITCvMUAxcVals = mean(restTRITCvMUAxcVals,1);
        restStdTRITCvMUAxcVals = std(restTRITCvMUAxcVals,0,1);
    end
    % save results
    AnalysisResults.(animalID).XCorr.Rest.(dataType).LFP_lags = restLFP_lags;
    AnalysisResults.(animalID).XCorr.Rest.(dataType).MUA_lags = restMUA_lags;
    AnalysisResults.(animalID).XCorr.Rest.(dataType).F = rest_F;
    AnalysisResults.(animalID).XCorr.Rest.(dataType).TRITCvLFPxcVals = restMeanTRITCvLFPxcVals;
    AnalysisResults.(animalID).XCorr.Rest.(dataType).TRITCvMUAxcVals = restMeanTRITCvMUAxcVals;
    AnalysisResults.(animalID).XCorr.Rest.(dataType).TRITCvMUAxcVals_std = restStdTRITCvMUAxcVals;
    %% cross-correlation analysis for NREM
    NREM_sleepTime = params.minTime.NREM;   % seconds
    [NREM_finalTRITC,NREM_allSleepFileIDs,NREM_finalBinTimes] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.TRITC.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    [NREM_finalMUA,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(neuralDataType).muaPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    NREM_uniqueSleepFileIDs = unique(NREM_allSleepFileIDs);
    jj = 1;
    for ff = 1:length(NREM_uniqueSleepFileIDs)
        % pull out the bin times (there may be multiple events) in each unique NREM sleep file
        NREM_uniqueSleepFileID = char(NREM_uniqueSleepFileIDs(ff));
        hh = 1;
        clear NREM_binTimes
        for gg = 1:length(NREM_allSleepFileIDs)
            NREM_sleepFileID = char(NREM_allSleepFileIDs(gg));
            if strcmp(NREM_uniqueSleepFileID,NREM_sleepFileID)
                NREM_binTimes{hh,1} = NREM_finalBinTimes{gg,1};
                hh = hh + 1;
            end
        end
        % pull out the Spectrogram data that matches the unique NREM sleep file
        NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecDataC.mat'];
        load(NREM_specDataFileID,'-mat')
        NREM_S_Data = SpecData.(neuralDataType).normS;
        for ii = 1:length(NREM_binTimes)
            NREM_Bins = NREM_binTimes{ii,1};
            NREM_startTime = NREM_Bins(1) - sleepBinWidth;
            NREM_endTime = NREM_Bins(end);
            if NREM_startTime > 5 && NREM_endTime < trialDuration_sec
                NREM_startTimeIndex = find(rest_T == NREM_startTime);
                NREM_startTimeIndex = NREM_startTimeIndex(1);
                NREM_durationIndex = find(rest_T == NREM_endTime);
                NREM_durationIndex = NREM_durationIndex(1);
                NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                editIndex{jj,1} = {'none'};
            elseif NREM_startTime == 5 && length(NREM_Bins) >= (params.minTime.NREM/sleepBinWidth + 1)
                NREM_startTime = NREM_Bins(2) - sleepBinWidth;
                NREM_endTime = NREM_Bins(end);
                NREM_startTimeIndex = find(rest_T == NREM_startTime);
                NREM_startTimeIndex = NREM_startTimeIndex(1);
                NREM_durationIndex = find(rest_T == NREM_endTime);
                NREM_durationIndex = NREM_durationIndex(1);
                NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                editIndex{jj,1} = {'leading'};
            elseif NREM_endTime == 900 && length(NREM_Bins) >= (params.minTime.NREM/sleepBinWidth + 1)
                NREM_startTime = NREM_Bins(1) - sleepBinWidth;
                NREM_endTime = NREM_Bins(end - 1);
                NREM_startTimeIndex = find(rest_T == NREM_startTime);
                NREM_startTimeIndex = NREM_startTimeIndex(1);
                NREM_durationIndex = find(rest_T == NREM_endTime);
                NREM_durationIndex = NREM_durationIndex(1);
                NREM_sleepNeuralVals{jj,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                editIndex{jj,1} = {'lagging'};
            else
                NREM_sleepNeuralVals{jj,1} = [];
                editIndex{jj,1} = {'delete'};
            end
            jj = jj + 1;
        end
    end
    % detrend spectrogram neural values
    for kk = 1:length(NREM_sleepNeuralVals)
        NREM_indSleepNeuralVals = NREM_sleepNeuralVals{kk,1};
        if isempty(NREM_indSleepNeuralVals) == false
            NREM_indSleepNeuralVals = NREM_indSleepNeuralVals(:,1:NREM_sleepTime*oneSecSpecFs)';
            NREM_dtSleepNeuralVals{kk,1} = detrend(NREM_indSleepNeuralVals,'constant')';
        else
            NREM_dtSleepNeuralVals{kk,1} = [];
        end
    end
    % adjust[TRITC] and MUA events to match the edits made to the length of each spectrogram
    mm = 1;
    for ll = 1:length(NREM_dtSleepNeuralVals)
        if isempty(NREM_dtSleepNeuralVals) == false
            NREM_finalSleepNeuralVals{mm,1} = NREM_dtSleepNeuralVals{ll,1};
            if strcmp(editIndex{ll,1},'none') == true
                NREM_TRITCVals = NREM_finalTRITC{ll,1}(1:NREM_sleepTime*samplingRate);
                NREM_MUAVals = NREM_finalMUA{ll,1}(1:NREM_sleepTime*samplingRate);
                NREM_finalTRITCVals{mm,1} = filtfilt(sos,g,detrend(NREM_TRITCVals,'constant'));
                NREM_finalMUAVals{mm,1} = filtfilt(sos,g,detrend(NREM_MUAVals,'constant'));
                mm = mm + 1;
            elseif strcmp(editIndex{ll,1},'leading') == true
                NREM_TRITCVals = NREM_finalTRITC{ll,1}((samplingRate*sleepBinWidth) + 1:(NREM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                NREM_MUAVals = NREM_finalMUA{ll,1}((samplingRate*sleepBinWidth) + 1:(NREM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                NREM_finalTRITCVals{mm,1} = filtfilt(sos,g,detrend(NREM_TRITCVals,'constant'));
                NREM_finalMUAVals{mm,1} = filtfilt(sos,g,detrend(NREM_MUAVals,'constant'));
                mm = mm + 1;
            elseif strcmp(editIndex{ll,1},'lagging') == true
                NREM_TRITCVals = NREM_finalTRITC{ll,1}(1:NREM_sleepTime*samplingRate);
                NREM_MUAVals = NREM_finalMUA{ll,1}(1:NREM_sleepTime*samplingRate);
                NREM_finalTRITCVals{mm,1} = filtfilt(sos,g,detrend(NREM_TRITCVals,'constant'));
                NREM_finalMUAVals{mm,1} = filtfilt(sos,g,detrend(NREM_MUAVals,'constant'));
                mm = mm + 1;
            elseif strcmp(editIndex{ll,1},'delete') == true
                % remove TRITC/MUA from final file
            end
        end
    end
    % run cross-correlation analysis - average through time
    NREM_F = SpecData.(neuralDataType).F;
    NREM_TRITCvLFPzHold = [];
    NREM_lagTime = 5;   % Seconds
    NREM_frequency = oneSecSpecFs;   % Hz
    NREM_maxLag = NREM_lagTime*NREM_frequency;
    NREM_TRITCvLFPxcVals = ones(size(NREM_indSleepNeuralVals,2),2*NREM_maxLag + 1);
    for nn = 1:length(NREM_finalSleepNeuralVals)
        for oo = 1:size(NREM_finalSleepNeuralVals{nn,1},1)
            NREM_TRITC_array = NREM_finalTRITCVals{nn,1};
            NREM_MUA_array = NREM_finalMUAVals{nn,1};
            NREM_Neural_array = NREM_finalSleepNeuralVals{nn,1}(oo,:);
            [NREM_TRITCvLFPxcVals(oo,:),NREM_LFP_lags] = xcorr(NREM_TRITC_array,NREM_Neural_array,NREM_maxLag,'coeff');
        end
        [NREM_TRITCvMUAxcVals(nn,:),NREM_MUA_lags] = xcorr(NREM_TRITC_array,NREM_MUA_array,NREM_maxLag,'coeff');
        NREM_TRITCvLFPzHold = cat(3,NREM_TRITCvLFPzHold,NREM_TRITCvLFPxcVals);
    end
    NREM_meanTRITCvLFPxcVals = mean(NREM_TRITCvLFPzHold,3);
    NREM_meanTRITCvMUAxcVals = mean(NREM_TRITCvMUAxcVals,1);
    NREM_stdTRITCvMUAxcVals = std(NREM_TRITCvMUAxcVals,0,1);
    % save results
    AnalysisResults.(animalID).XCorr.NREM.(dataType).LFP_lags = NREM_LFP_lags;
    AnalysisResults.(animalID).XCorr.NREM.(dataType).MUA_lags = NREM_MUA_lags;
    AnalysisResults.(animalID).XCorr.NREM.(dataType).F = NREM_F;
    AnalysisResults.(animalID).XCorr.NREM.(dataType).TRITCvLFPxcVals = NREM_meanTRITCvLFPxcVals;
    AnalysisResults.(animalID).XCorr.NREM.(dataType).TRITCvMUAxcVals = NREM_meanTRITCvMUAxcVals;
    AnalysisResults.(animalID).XCorr.NREM.(dataType).TRITCvMUAxcVals_std = NREM_stdTRITCvMUAxcVals;
    %% cross-correlation analysis for REM
    %{
    REM_sleepTime = params.minTime.REM;   % seconds
    [REM_finalTRITC,REM_allSleepFileIDs,REM_finalBinTimes] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.TRITC.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    [REM_finalMUA,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(neuralDataType).muaPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
    REM_uniqueSleepFileIDs = unique(REM_allSleepFileIDs);
    uu = 1;
    clear editIndex
    for qq = 1:length(REM_uniqueSleepFileIDs)
        % pull out the bin times (there may be multiple events) in each unique REM sleep file
        REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(qq));
        ss = 1;
        clear REM_binTimes
        for rr = 1:length(REM_allSleepFileIDs)
            REM_sleepFileID = char(REM_allSleepFileIDs(rr));
            if strcmp(REM_uniqueSleepFileID,REM_sleepFileID)
                REM_binTimes{ss,1} = REM_finalBinTimes{rr,1};
                ss = ss + 1;
            end
        end
        % pull out the Spectrogram data that matches the unique REM sleep file
        REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecDataC.mat'];
        load(REM_specDataFileID,'-mat')
        REM_S_Data = SpecData.(neuralDataType).normS;
        for tt = 1:length(REM_binTimes)
            REM_Bins = REM_binTimes{tt,1};
            REM_startTime = REM_Bins(1) - sleepBinWidth;
            REM_endTime = REM_Bins(end);
            if REM_startTime > 5 && REM_endTime < trialDuration_sec
                REM_startTimeIndex = find(rest_T == REM_startTime);
                REM_startTimeIndex = REM_startTimeIndex(1);
                REM_durationIndex = find(rest_T == REM_endTime);
                REM_durationIndex = REM_durationIndex(1);
                REM_sleepNeuralVals{uu,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                editIndex{uu,1} = {'none'};
            elseif REM_startTime == 5 && length(REM_Bins) >= (params.minTime.REM/sleepBinWidth + 1)
                REM_startTime = REM_Bins(2) - sleepBinWidth;
                REM_endTime = REM_Bins(end);
                REM_startTimeIndex = find(rest_T == REM_startTime);
                REM_startTimeIndex = REM_startTimeIndex(1);
                REM_durationIndex = find(rest_T == REM_endTime);
                REM_durationIndex = REM_durationIndex(1);
                REM_sleepNeuralVals{uu,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                editIndex{uu,1} = {'leading'};
            elseif REM_endTime == 900 && length(REM_Bins) >= (params.minTime.REM/sleepBinWidth + 1)
                REM_startTime = REM_Bins(1) - sleepBinWidth;
                REM_endTime = REM_Bins(end - 1);
                REM_startTimeIndex = find(rest_T == REM_startTime);
                REM_startTimeIndex = REM_startTimeIndex(1);
                REM_durationIndex = find(rest_T == REM_endTime);
                REM_durationIndex = REM_durationIndex(1);
                REM_sleepNeuralVals{uu,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                editIndex{uu,1} = {'lagging'};
            else
                REM_sleepNeuralVals{uu,1} = [];
                editIndex{uu,1} = {'delete'};
            end
            uu = uu + 1;
        end
    end
    % detrend spectrogram neural values
    for vv = 1:length(REM_sleepNeuralVals)
        REM_indSleepNeuralVals = REM_sleepNeuralVals{vv,1};
        if isempty(REM_indSleepNeuralVals) == false
            REM_indSleepNeuralVals = REM_indSleepNeuralVals(:,1:REM_sleepTime*oneSecSpecFs)';
            REM_dtSleepNeuralVals{vv,1} = detrend(REM_indSleepNeuralVals,'constant')';
        else
            REM_dtSleepNeuralVals{vv,1} = [];
        end
    end
    % adjust [TRITC] and MUA events to match the edits made to the length of each spectrogram
    xx = 1;
    for ww = 1:length(REM_dtSleepNeuralVals)
        if isempty(REM_dtSleepNeuralVals) == false
            REM_finalSleepNeuralVals{xx,1} = REM_dtSleepNeuralVals{ww,1};
            if strcmp(editIndex{ww,1},'none') == true
                REM_TRITCVals = REM_finalTRITC{ww,1}(1:REM_sleepTime*samplingRate);
                REM_MUAVals = REM_finalMUA{ww,1}(1:REM_sleepTime*samplingRate);
                REM_finalTRITCVals{xx,1} = filtfilt(sos,g,detrend(REM_TRITCVals,'constant'));
                REM_finalMUAVals{xx,1} = filtfilt(sos,g,detrend(REM_MUAVals,'constant'));
                xx = xx + 1;
            elseif strcmp(editIndex{ww,1},'leading') == true
                REM_TRITCVals = REM_finalTRITC{ww,1}((samplingRate*sleepBinWidth) + 1:(REM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                REM_MUAVals = REM_finalMUA{ww,1}((samplingRate*sleepBinWidth) + 1:(REM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                REM_finalTRITCVals{xx,1} = filtfilt(sos,g,detrend(REM_TRITCVals,'constant'));
                REM_finalMUAVals{xx,1} = filtfilt(sos,g,detrend(REM_MUAVals,'constant'));
                xx = xx + 1;
            elseif strcmp(editIndex{ww,1},'lagging') == true
                REM_TRITCVals = REM_finalTRITC{ww,1}(1:REM_sleepTime*samplingRate);
                REM_MUAVals = REM_finalMUA{ww,1}(1:REM_sleepTime*samplingRate);
                REM_finalTRITCVals{xx,1} = filtfilt(sos,g,detrend(REM_TRITCVals,'constant'));
                REM_finalMUAVals{xx,1} = filtfilt(sos,g,detrend(REM_MUAVals,'constant'));
                xx = xx + 1;
            elseif strcmp(editIndex{ww,1},'delete') == true
                % remove TRITC/MUA from final file
            end
        end
    end
    % run cross-correlation analysis - average through time
    REM_F = SpecData.(neuralDataType).F;
    REM_TRITCvLFPzHold = [];
    REM_lagTime = 5;   % Seconds
    REM_frequency = oneSecSpecFs;   % Hz
    REM_maxLag = REM_lagTime*REM_frequency;
    REM_TRITCvLFPxcVals = ones(size(REM_indSleepNeuralVals,2),2*REM_maxLag + 1);
    for yy = 1:length(REM_finalSleepNeuralVals)
        for zz = 1:size(REM_finalSleepNeuralVals{yy,1},1)
            REM_TRITC_array = REM_finalTRITCVals{yy,1};
            REM_MUA_array = REM_finalMUAVals{yy,1};
            REM_Neural_array = REM_finalSleepNeuralVals{yy,1}(zz,:);
            [REM_TRITCvLFPxcVals(zz,:),REM_LFP_lags] = xcorr(REM_TRITC_array,REM_Neural_array,REM_maxLag,'coeff');
        end
        [REM_TRITCvMUAxcVals(yy,:),REM_MUA_lags] = xcorr(REM_TRITC_array,REM_MUA_array,REM_maxLag,'coeff');
        REM_TRITCvLFPzHold = cat(3,REM_TRITCvLFPzHold,REM_TRITCvLFPxcVals);
    end
    REM_meanTRITCvLFPxcVals = mean(REM_TRITCvLFPzHold,3);
    REM_meanTRITCvMUAxcVals = mean(REM_TRITCvMUAxcVals,1);
    REM_stdTRITCvMUAxcVals = std(REM_TRITCvMUAxcVals,0,1);
    % save results
    AnalysisResults.(animalID).XCorr.REM.(dataType).LFP_lags = REM_LFP_lags;
    AnalysisResults.(animalID).XCorr.REM.(dataType).MUA_lags = REM_MUA_lags;
    AnalysisResults.(animalID).XCorr.REM.(dataType).F = REM_F;
    AnalysisResults.(animalID).XCorr.REM.(dataType).TRITCvLFPxcVals = REM_meanTRITCvLFPxcVals;
    AnalysisResults.(animalID).XCorr.REM.(dataType).TRITCvMUAxcVals = REM_meanTRITCvMUAxcVals;
    AnalysisResults.(animalID).XCorr.REM.(dataType).TRITCvMUAxcVals_std = REM_stdTRITCvMUAxcVals;
end
    %}
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
