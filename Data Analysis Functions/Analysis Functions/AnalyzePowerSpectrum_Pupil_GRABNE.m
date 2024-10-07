function [AnalysisResults] = AnalyzePowerSpectrum_Pupil_GRABNE(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP'};%{'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_gammaBandPower'};
modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs
delim = '\';
dataLocation = [rootFolder delim animalID delim 'CombinedImaging'];
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
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat strut
SleepDataFileStruct = dir('*_SleepData.mat');
SleepDataFile = {SleepDataFileStruct.name}';
SleepDataFileID = char(SleepDataFile);
load(SleepDataFileID,'-mat')
% scoring results
% find and load manual baseline event information
scoringResultsFileStruct = dir('*Manual_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% filter characteristics & resting criteria
samplingRate = RestData.Rhodamine.Z_Ach.RhodamineCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for behavior-based power spectrum analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze power spectra during periods of rest
    % pull data from RestData.mat structure   
    % get the rest data structure
    if strcmp(dataType,'Ach_Rhodamine') == true
        [restLogical] = FilterEvents_IOS(RestData.Rhodamine.Z_Ach,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.Rhodamine.Z_Ach,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Rhodamine.Z_Ach.fileIDs(combRestLogical,:);
        restEventTimes = RestData.Rhodamine.Z_Ach.eventTimes(combRestLogical,:);
        restDurations = RestData.Rhodamine.Z_Ach.durations(combRestLogical,:);
        restData = RestData.Rhodamine.Z_Ach.data(combRestLogical,:);
    elseif strcmp(dataType,'NE_Rhodamine') == true
        [restLogical] = FilterEvents_IOS(RestData.Rhodamine.Z_NE,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.Rhodamine.Z_NE,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Rhodamine.Z_NE.fileIDs(combRestLogical,:);
        restEventTimes = RestData.Rhodamine.Z_NE.eventTimes(combRestLogical,:);
        restDurations = RestData.Rhodamine.Z_NE.durations(combRestLogical,:);
        restData = RestData.Rhodamine.Z_NE.data(combRestLogical,:);
    elseif strcmp(dataType,'Ach_GFP') == true
        [restLogical] = FilterEvents_IOS(RestData.GFP.Z_Ach,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.GFP.Z_Ach,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.GFP.Z_Ach.fileIDs(combRestLogical,:);
        restEventTimes = RestData.GFP.Z_Ach.eventTimes(combRestLogical,:);
        restDurations = RestData.GFP.Z_Ach.durations(combRestLogical,:);
        restData = RestData.GFP.Z_Ach.data(combRestLogical,:);
    elseif strcmp(dataType,'NE_GFP') == true
        [restLogical] = FilterEvents_IOS(RestData.GFP.Z_NE,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.GFP.Z_NE,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.GFP.Z_NE.fileIDs(combRestLogical,:);
        restEventTimes = RestData.GFP.Z_NE.eventTimes(combRestLogical,:);
        restDurations = RestData.GFP.Z_NE.durations(combRestLogical,:);
        restData = RestData.GFP.Z_NE.data(combRestLogical,:);
    % elseif strcmp(dataType,'RH_thetaBandPower') == true
    %     [restLogical] = FilterEvents_IOS(RestData.cortical_RH.thetaBandPower,RestCriteria);
    %     [puffLogical] = FilterEvents_IOS(RestData.cortical_RH.thetaBandPower,RestPuffCriteria);
    %     combRestLogical = logical(restLogical.*puffLogical);
    %     restFileIDs = RestData.cortical_RH.thetaBandPower.fileIDs(combRestLogical,:);
    %     restEventTimes = RestData.cortical_RH.thetaBandPower.eventTimes(combRestLogical,:);
    %     restDurations = RestData.cortical_RH.thetaBandPower.durations(combRestLogical,:);
    %     restData = RestData.cortical_RH.thetaBandPower.data(combRestLogical,:);
    % elseif strcmp(dataType,'RH_gammaBandPower') == true
    %     [restLogical] = FilterEvents_IOS(RestData.cortical_RH.gammaBandPower,RestCriteria);
    %     [puffLogical] = FilterEvents_IOS(RestData.cortical_RH.gammaBandPower,RestPuffCriteria);
    %     combRestLogical = logical(restLogical.*puffLogical);
    %     restFileIDs = RestData.cortical_RH.gammaBandPower.fileIDs(combRestLogical,:);
    %     restEventTimes = RestData.cortical_RH.gammaBandPower.eventTimes(combRestLogical,:);
    %     restDurations = RestData.cortical_RH.gammaBandPower.durations(combRestLogical,:);
    %     restData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
    % else
    %     [restLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestCriteria);
    %     [puffLogical] = FilterEvents_IOS(RestData.Pupil.(dataType),RestPuffCriteria);
    %     combRestLogical = logical(restLogical.*puffLogical);
    %     restFileIDs = RestData.Pupil.(dataType).fileIDs(combRestLogical,:);
    %     restEventTimes = RestData.Pupil.(dataType).eventTimes(combRestLogical,:);
    %     restDurations = RestData.Pupil.(dataType).durations(combRestLogical,:);
    %     restData = RestData.Pupil.(dataType).data(combRestLogical,:);
    end

    % keep only the data that occurs within the manually-approved awake regions
    [finalRestData,~,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    clear procRestData
    if isempty(finalRestData) == false
        % detrend and truncate data to minimum length to match events
        for bb = 1:length(finalRestData)
            if length(finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{bb,1});
                restPad = (ones(1,restChunkSampleDiff))*finalRestData{bb,1}(end);
                procRestData{bb,1} = horzcat(finalRestData{bb,1},restPad); %#ok<*AGROW>
                procRestData{bb,1} = detrend(procRestData{bb,1},'constant');
            else
                procRestData{bb,1} = detrend(finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for cc = 1:length(procRestData)
            if sum(isnan(procRestData{cc,1})) == 0
                restDataMat(:,zz) = procRestData{cc,1};
                zz = zz + 1;
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [rest_S,rest_f,rest_sErr] = mtspectrumc(restDataMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).S = rest_S;
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).f = rest_f;
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).sErr = rest_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.Rest.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of alert
    zz = 1;
    clear awakeData procAwakeData
    awakeData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,fileDate,awakeDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(awakeDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) > 436   % 36 bins (180 total) or 3 minutes of Asleep
            load(procDataFileID,'-mat')
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if strcmp(dataType,'Ach_Rhodamine') == true
                        awakeData{zz,1} = ProcData.data.Rhodamine.Z_Ach;
                        zz = zz + 1;
                    elseif strcmp(dataType,'NE_Rhodamine') == true
                        awakeData{zz,1} = ProcData.data.Rhodamine.Z_NE;
                        zz = zz + 1;
                    elseif strcmp(dataType,'Ach_GFP') == true
                        awakeData{zz,1} = ProcData.data.GFP.Z_Ach;
                        zz = zz + 1;
                    elseif strcmp(dataType,'NE_GFP') == true
                        awakeData{zz,1} = ProcData.data.GFP.Z_NE;
                        zz = zz + 1;
                    % elseif strcmp(dataType,'RH_thetaBandPower') == true
                    %     awakeData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean;
                    %     zz = zz + 1;
                    % elseif strcmp(dataType,'RH_gammaBandPower') == true
                    %     awakeData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                    %     zz = zz + 1;
                    % else
                    %     if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                    %         awakeData{zz,1} = ProcData.data.Pupil.(dataType);
                    %         zz = zz + 1;
                    %     end
                    end
                end
            end
        end
    end
    if isempty(awakeData) == false
        % detrend data
        for bb = 1:length(awakeData)
            procAwakeData{bb,1} = detrend(awakeData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        awakeDataMat = zeros(length(procAwakeData{1,1}),length(procAwakeData));
        for cc = 1:length(procAwakeData)
            awakeDataMat(:,cc) = procAwakeData{cc,1}(1:length(procAwakeData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [1 1];   % Tapers [n, 2n - 1] [10,19];
        [awake_S,awake_f,awake_sErr] = mtspectrumc(awakeDataMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).S = awake_S;
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).f = awake_f;
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).sErr = awake_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.Awake.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of Asleep
    zz = 1;
    clear asleepData procAsleepData
    asleepData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,fileDate,asleepDataFileID] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        scoringLabels = [];
        for cc = 1:length(ScoringResults.fileIDs)
            if strcmp(asleepDataFileID,ScoringResults.fileIDs{cc,1}) == true
                scoringLabels = ScoringResults.labels{cc,1};
            end
        end
        % check labels to match arousal state
        if sum(strcmp(scoringLabels,'Not Sleep')) < 312   % 36 bins (180 total) or 3 minutes of awake
            load(procDataFileID,'-mat')
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end
                % don't include trials with stimulation
                if isempty(puffs) == true
                    if strcmp(dataType,'Ach_Rhodamine') == true
                        asleepData{zz,1} = ProcData.data.Rhodamine.Z_Ach;
                        zz = zz + 1;
                    elseif strcmp(dataType,'NE_Rhodamine') == true
                        asleepData{zz,1} = ProcData.data.Rhodamine.Z_NE;
                        zz = zz + 1;
                    elseif strcmp(dataType,'NE_GFP') == true
                        asleepData{zz,1} = ProcData.data.GFP.Z_NE;
                        zz = zz + 1;
                    elseif strcmp(dataType,'Ach_GFP') == true
                        asleepData{zz,1} = ProcData.data.GFP.Z_Ach;
                        zz = zz + 1;
                    % elseif strcmp(dataType,'RH_tehtaBandPower') == true
                    %     asleepData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean;
                    %     zz = zz + 1;
                    % elseif strcmp(dataType,'RH_gammaBandPower') == true
                    %     asleepData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                    %     zz = zz + 1;
                    % else
                    %     if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                    %         asleepData{zz,1} = ProcData.data.Pupil.(dataType);
                    %         zz = zz + 1;
                    %     end
                    end
                end
            end
        end
    end
    if isempty(asleepData) == false
        % detrend data
        for bb = 1:length(asleepData)
            procAsleepData{bb,1} = detrend(asleepData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        asleepDataMat = zeros(length(procAsleepData{1,1}),length(procAsleepData));
        for cc = 1:length(procAsleepData)
            asleepDataMat(:,cc) = procAsleepData{cc,1}(1:length(procAsleepData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        [asleep_S,asleep_f,asleep_sErr] = mtspectrumc(asleepDataMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).S = asleep_S;
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).f = asleep_f;
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).sErr = asleep_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.Asleep.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of all data
    zz = 1;
    clear allData procAllData
    allData = [];
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        [~,fileDate,~] = GetFileInfo_IOS(procDataFileID);
        strDay = ConvertDate_IOS(fileDate);
        load(procDataFileID,'-mat')
        if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
            try
                puffs = ProcData.data.stimulations.LPadSol;
            catch
                puffs = ProcData.data.solenoids.LPadSol;
            end
            % don't include trials with stimulation
            if isempty(puffs) == true
                if strcmp(dataType,'Ach_Rhodamine') == true
                    allData{zz,1} = ProcData.data.Rhodamine.Z_Ach;
                    zz = zz + 1;
                elseif strcmp(dataType,'NE_Rhodamine') == true
                    allData{zz,1} = ProcData.data.Rhodamine.Z_NE;
                    zz = zz + 1;
                elseif strcmp(dataType,'NE_GFP') == true
                    allData{zz,1} = ProcData.data.GFP.Z_NE;
                    zz = zz + 1;
                elseif strcmp(dataType,'Ach_GFP') == true
                    allData{zz,1} = ProcData.data.GFP.Z_Ach;
                    zz = zz + 1;
                % elseif strcmp(dataType,'RH_thetaBandPower') == true
                %     allData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean;
                %     zz = zz + 1;
                % elseif strcmp(dataType,'RH_gammaBandPower') == true
                %     allData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                %     zz = zz + 1;
                % else
                %     if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                %         allData{zz,1} = ProcData.data.Pupil.(dataType);
                %         zz = zz + 1;
                %     end
                end
            end
        end
    end
    if isempty(allData) == false
        % detrend data
        for bb = 1:length(allData)
            procAllData{bb,1} = detrend(allData{bb,1},'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        allDataMat = zeros(length(procAllData{1,1}),length(procAllData));
        for cc = 1:length(procAllData)
            allDataMat(:,cc) = procAllData{cc,1}(1:length(procAllData{1,1}));
        end
        % calculate the power spectra of the desired signals
        params.tapers = [10,19];   % Tapers [n, 2n - 1]
        [all_S,all_f,all_sErr] = mtspectrumc(allDataMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).S = all_S;
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).f = all_f;
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).sErr = all_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.All.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of NREM
      % pull data from SleepData.mat structure

    if strcmp(dataType,'Ach_Rhodamine') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.Z_Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'NE_Rhodamine') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.Z_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'Ach_GFP') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'NE_GFP') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Z_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % elseif strcmp(dataType,'RH_thetaBandPower') == true
    %     [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.thetaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % elseif strcmp(dataType,'RH_gammaBandPower') == true
    %     [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.gammaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % else
    %     [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    % 
    end
    
    if isempty(nremData) == false
        % detrend and truncate data to minimum length to match events
        for dd = 1:length(nremData)
            nremData{dd,1} = detrend(nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for ee = 1:length(nremData)
            if sum(isnan(nremData{ee,1})) == 0
                nremMat(:,zz) = nremData{ee,1};
                zz = zz + 1;
            end
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        [nrem_S,nrem_f,nrem_sErr] = mtspectrumc(nremMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).S = nrem_S;
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).f = nrem_f;
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).sErr = nrem_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.NREM.(dataType).sErr = [];
    end
    %% analyze power spectra during periods of REM
    % pull data from SleepData.mat structure
    if isempty(SleepData.(modelType).REM.data.Pupil) == false
        if strcmp(dataType,'Ach_Rhodamine') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.Z_Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'NE_Rhodamine') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.Z_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'Ach_GFP') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Z_Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'NE_GFP') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Z_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % elseif strcmp(dataType,'RH_thetaBandPower') == true
        %     [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.thetaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % elseif strcmp(dataType,'RH_gammaBandPower') == true
        %     [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.gammaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % else
        %     [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        end
    else
            remData = [];
    end
    if isempty(remData) == false
        % detrend and truncate data to minimum length to match events
        for dd = 1:length(remData)
            remData{dd,1} = detrend(remData{dd,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for ee = 1:length(remData)
            if sum(isnan(remData{ee,1})) == 0
                remMat(:,zz) = remData{ee,1};
                zz = zz + 1;
            end
        end
        % calculate the power spectra of the desired signals
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        [rem_S,rem_f,rem_sErr] = mtspectrumc(remMat,params);
        % save results
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).S = rem_S;
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).f = rem_f;
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).sErr = rem_sErr;
    else
        % save results
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).S = [];
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).f = [];
        AnalysisResults.(animalID).PowerSpectrum.REM.(dataType).sErr = [];
    end
end
%% save data
cd([rootFolder delim])
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end
