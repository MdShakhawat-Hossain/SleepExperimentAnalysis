function [AnalysisResults] = AnalyzePowerSpectrum_DataDivision(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral power of hemodynamic [HbT] and neural signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'ACh_CBV','NE_CBV','ACh_GFP','NE_GFP'};%{'zDiameter','ACh_CBV','NE_CBV','ACh_GFP','NE_GFP','RH_thetaBandPower','RH_gammaBandPower'};
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
samplingRate = RestData.CBV.P_ACh.CBVCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {10};%5
% go through each valid data type for behavior-based power spectrum analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    %% analyze power spectra during periods of rest
    % pull data from RestData.mat structure   
    % get the rest data structure
    if strcmp(dataType,'ACh_CBV') == true
        [restLogical] = FilterEvents_IOS(RestData.CBV.P_ACh,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.CBV.P_ACh,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.CBV.P_ACh.fileIDs(combRestLogical,:);
        restEventTimes = RestData.CBV.P_ACh.eventTimes(combRestLogical,:);
        restDurations = RestData.CBV.P_ACh.durations(combRestLogical,:);
        restData = RestData.CBV.P_ACh.data(combRestLogical,:);
    elseif strcmp(dataType,'NE_CBV') == true
        [restLogical] = FilterEvents_IOS(RestData.CBV.P_NE,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.CBV.P_NE,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.CBV.P_NE.fileIDs(combRestLogical,:);
        restEventTimes = RestData.CBV.P_NE.eventTimes(combRestLogical,:);
        restDurations = RestData.CBV.P_NE.durations(combRestLogical,:);
        restData = RestData.CBV.P_NE.data(combRestLogical,:);
    elseif strcmp(dataType,'ACh_GFP') == true
        [restLogical] = FilterEvents_IOS(RestData.GFP.P_ACh,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.GFP.P_ACh,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.GFP.P_ACh.fileIDs(combRestLogical,:);
        restEventTimes = RestData.GFP.P_ACh.eventTimes(combRestLogical,:);
        restDurations = RestData.GFP.P_ACh.durations(combRestLogical,:);
        restData = RestData.GFP.P_ACh.data(combRestLogical,:);
    elseif strcmp(dataType,'NE_GFP') == true
        [restLogical] = FilterEvents_IOS(RestData.GFP.P_NE,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.GFP.P_NE,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.GFP.P_NE.fileIDs(combRestLogical,:);
        restEventTimes = RestData.GFP.P_NE.eventTimes(combRestLogical,:);
        restDurations = RestData.GFP.P_NE.durations(combRestLogical,:);
        restData = RestData.GFP.P_NE.data(combRestLogical,:);
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
      for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,awakeDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(awakeDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
           % cut the data into three 10 minutes chunk
            ScoreSize = 10*12;% 10 minutes X 12 5 sec periods per minutes
            
            load(procDataFileID,'-mat')
            DataSize = 10*60*ProcData.notes.dsFs; % 15 minutes X 60 secs X 30 Hz;
                % only run on files with good pupil measurement
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                        % don't include trials with stimulation
                        if isempty(puffs) == true
                            if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                               for SSL =  1:1:5
                                ScoreLabels = scoringLabels(((SSL-1)*ScoreSize)+1:(SSL*ScoreSize));
                                  % check labels to match arousal state
                                    if sum(strcmp(ScoreLabels,'Not Sleep')) > 108 % 90% of the time awake                              
                                        % pull data based on data type
                                        if strcmp(dataType,'ACh_CBV') == true
                                            awakeData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'NE_CBV') == true
                                            awakeData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'ACh_GFP') == true
                                            awakeData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'NE_GFP') == true
                                            awakeData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        end
                                        zz = zz + 1;
                                    end
                               end
                            end
                        end
                end
      end
      % calculate the power spectra
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
        params.tapers = [5,9];   % Tapers [n, 2n - 1] [5,9];
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
      for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
           % cut the data into three 10 minutes chunk
            ScoreSize = 10*12;% 10 minutes X 12 5 sec periods per minutes
            
            load(procDataFileID,'-mat')
            DataSize = 10*60*ProcData.notes.dsFs; % 10 minutes X 60 secs X 30 Hz;
                % only run on files with good pupil measurement
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                        % don't include trials with stimulation
                        if isempty(puffs) == true
                            if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                               for SSL =  1:1:5
                                ScoreLabels = scoringLabels(((SSL-1)*ScoreSize)+1:(SSL*ScoreSize));
                                  % check labels to match arousal state
                                    if sum(strcmp(ScoreLabels,'Not Sleep')) < 36 % 70% of the time asleep                              
                                        % pull data based on data type
                                        if strcmp(dataType,'ACh_CBV') == true
                                            asleepData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'NE_CBV') == true
                                            asleepData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'ACh_GFP') == true
                                            asleepData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(dataType,'NE_GFP') == true
                                            asleepData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        end
                                        zz = zz + 1;
                                    end
                               end
                            end
                        end
                end
       end

    % calculate the power spectra
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
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
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

    for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,~] = GetFileInfo_JNeurosci2022(procDataFileID);
            load(procDataFileID,'-mat')

            DataSize = 10*60*ProcData.notes.dsFs; % 15 minutes X 60 secs X 30 Hz;
            % only run on files with good pupil measurement
            if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                try
                    puffs = ProcData.data.stimulations.LPadSol;
                catch
                    puffs = ProcData.data.solenoids.LPadSol;
                end                  
                    % don't include trials with stimulation
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                            % cut the data into three 10 minutes chunk
                             for SSL =  1:1:5
                                % pull data based on data type
                                if strcmp(dataType,'ACh_CBV') == true
                                    allData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(dataType,'NE_CBV') == true
                                    allData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(dataType,'ACh_GFP') == true
                                    allData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(dataType,'NE_GFP') == true
                                    allData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                end
                                zz = zz + 1;
                            end
                       end
                   end
           end 
    end

    % calculate the power spectra
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
        params.tapers = [5,9];   % Tapers [n, 2n - 1]
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

    if strcmp(dataType,'ACh_CBV') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'NE_CBV') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'ACh_GFP') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
    elseif strcmp(dataType,'NE_GFP') == true
        [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
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
        if strcmp(dataType,'ACh_CBV') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'NE_CBV') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'ACh_GFP') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        elseif strcmp(dataType,'NE_GFP') == true
            [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
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
