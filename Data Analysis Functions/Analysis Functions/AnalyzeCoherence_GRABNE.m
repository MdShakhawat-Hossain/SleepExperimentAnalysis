function [Results_Coherence] = AnalyzeCoherence_GRABNE(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral coherence between neural-hemodynamic [HbT] signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'Ach_Rhodamine','zDiameter','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_gammaBandPower'};
hemDataTypes = {'Ach_Rhodamine','zDiameter','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_gammaBandPower'};
modelType = 'Ensemble';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 30;
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
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load AsleepData.mat struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
%
scoringResultsFileStruct = dir('*Ensemble_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% lowpass filter
samplingRate = RestData.Rhodamine.Ach.RhodamineCamSamplingRate;
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal-based coherence analysis
for zzz = 1:length(hemDataTypes)
    hemDataType = hemDataTypes{1,zzz};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        %% analyze neural-hemo coherence during periods of rest
        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_IOS(RestData.Rhodamine.Ach,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.Rhodamine.Ach,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Rhodamine.Ach.fileIDs(combRestLogical,:);
        restEventTimes = RestData.Rhodamine.Ach.eventTimes(combRestLogical,:);
        restDurations = RestData.Rhodamine.Ach.durations(combRestLogical,:);
        
        if strcmp(hemDataType,'Ach_Rhodamine') == true
            hemDatatype_restData = RestData.Rhodamine.Ach.data(combRestLogical,1);
        elseif strcmp(hemDataType,'NE_Rhodamine') == true
            hemDatatype_restData = RestData.Rhodamine.NE.data(combRestLogical,:);
        elseif strcmp(hemDataType,'Ach_GFP') == true
            hemDatatype_restData = RestData.GFP.Ach.data(combRestLogical,:);
        elseif strcmp(hemDataType,'NE_GFP') == true
            hemDatatype_restData = RestData.GFP.NE.data(combRestLogical,:);
        elseif strcmp(hemDataType,'RH_thetaBandPower') == true
            hemDatatype_restData = RestData.cortical_RH.thetaBandPower.data(combRestLogical,:);
        elseif strcmp(hemDataType,'RH_gammaBandPower') == true
            hemDatatype_restData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
        elseif strcmp(hemDataType,'zDiameter') == true
            hemDatatype_restData = RestData.Pupil.zDiameter.data(combRestLogical,:);
        end

        if strcmp(datatype,'Ach_Rhodamine') == true
            Datatype_restData = RestData.Rhodamine.Ach.data(combRestLogical,1);
        elseif strcmp(datatype,'NE_Rhodamine') == true
            Datatype_restData = RestData.Rhodamine.NE.data(combRestLogical,:);
        elseif strcmp(datatype,'Ach_GFP') == true
            Datatype_restData = RestData.GFP.Ach.data(combRestLogical,:);
        elseif strcmp(datatype,'NE_GFP') == true
            Datatype_restData = RestData.GFP.NE.data(combRestLogical,:);
        elseif strcmp(datatype,'RH_thetaBandPower') == true
            Datatype_restData = RestData.cortical_RH.thetaBandPower.data(combRestLogical,:);
        elseif strcmp(datatype,'RH_gammaBandPower') == true
            Datatype_restData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
        elseif strcmp(datatype,'zDiameter') == true
            Datatype_restData = RestData.Pupil.zDiameter.data(combRestLogical,:);
        end

        % keep only the data that occurs within the manually-approved awake regions
        [hemDatatype_finalRestData,~,~,~] = RemoveInvalidData_IOS(hemDatatype_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Datatype_FinalRestData,~,~,~] = RemoveInvalidData_IOS(Datatype_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear hemDatatype_procRestData Datatype_procRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(hemDatatype_finalRestData)
            if length(hemDatatype_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(hemDatatype_finalRestData{bb,1});
                hemDatatype_restPad = (ones(1,restChunkSampleDiff))*hemDatatype_finalRestData{bb,1}(end);
                Datatype_restPad = (ones(1,restChunkSampleDiff))*Datatype_FinalRestData{bb,1}(end);
                hemDatatype_procRestData{bb,1} = horzcat(hemDatatype_finalRestData{bb,1},hemDatatype_restPad); %#ok<*AGROW>
                Datatype_procRestData{bb,1} = horzcat(Datatype_FinalRestData{bb,1},Datatype_restPad);
                hemDatatype_procRestData{bb,1} = detrend(hemDatatype_procRestData{bb,1},'constant');
                Datatype_procRestData{bb,1} = detrend(Datatype_procRestData{bb,1},'constant');
            else
                hemDatatype_procRestData{bb,1} = detrend(hemDatatype_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                Datatype_procRestData{bb,1} = detrend(Datatype_FinalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for cc = 1:length(hemDatatype_procRestData)
            if sum(isnan(Datatype_procRestData{cc,1})) == 0
                hemDatatype_restDataMat(:,zz) = hemDatatype_procRestData{cc,1};
                Datatype_restDataMat(:,zz) = Datatype_procRestData{cc,1};
                zz = zz + 1;
            end
        end
        % parameters for coherencyc - information available in function
        params.tapers = [1,1]; % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1]; % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(hemDatatype_restDataMat,Datatype_restDataMat,params);
        % save results
        Results_Coherence.(animalID).Rest.(dataType).(hemDataType).C = C_RestData;
        Results_Coherence.(animalID).Rest.(dataType).(hemDataType).f = f_RestData;
        Results_Coherence.(animalID).Rest.(dataType).(hemDataType).confC = confC_RestData;
        Results_Coherence.(animalID).Rest.(dataType).(hemDataType).cErr = cErr_RestData;
        %% analyze neural-hemo coherence during periods of alert
        zz = 1;
        clear hemDatatype_awakeData Datatype_awakeData hemDatatype_procAwakeData Datatype_procAwakeData
        hemDatatype_awakeData = []; Datatype_awakeData = [];
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
            if sum(strcmp(scoringLabels,'Not Sleep')) > 144   % 36 bins (180 total) or 3 minutes of asleep
                load(procDataFileID,'-mat')
                % don't include trials with stimulation
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemDataType,'Ach_Rhodamine') == true
                                hemDatatype_awakeData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(hemDataType,'NE_Rhodamine') == true
                                hemDatatype_awakeData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(hemDataType,'Ach_GFP') == true
                            hemDatatype_awakeData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(hemDataType,'NE_GFP') == true
                            hemDatatype_awakeData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(hemDataType,'RH_thetaBandPower') == true
                                hemDatatype_awakeData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                hemDatatype_awakeData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                            end
                            Datatype_awakeData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(hemDatatype_awakeData) == false
            for bb = 1:length(hemDatatype_awakeData)
                hemDatatype_procAwakeData{bb,1} = detrend(hemDatatype_awakeData{bb,1},'constant');
                Datatype_procAwakeData{bb,1} = detrend(Datatype_awakeData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            hemDatatype_awakeDataMat = zeros(length(hemDatatype_procAwakeData{1,1}),length(hemDatatype_procAwakeData));
            Datatype_awakeDataMat = zeros(length(Datatype_procAwakeData{1,1}),length(Datatype_procAwakeData));
            for cc = 1:length(hemDatatype_procAwakeData)
                hemDatatype_awakeDataMat(:,cc) = hemDatatype_procAwakeData{cc,1};
                Datatype_awakeDataMat(:,cc) = Datatype_procAwakeData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc(hemDatatype_awakeDataMat,Datatype_awakeDataMat,params);
            % save results
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).C = C_AwakeData;
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).f = f_AwakeData;
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).confC = confC_AwakeData;
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).cErr = cErr_AwakeData;
        else
            % save results
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).C = [];
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).f = [];
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).confC = [];
            Results_Coherence.(animalID).Awake.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of aasleep
        zz = 1;
        clear hemDatatype_asleepData Datatype_asleepData hemDatatype_procAsleepData Datatype_procAsleepData
        hemDatatype_asleepData = []; Datatype_asleepData = [];
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
            if sum(strcmp(scoringLabels,'Not Sleep')) < 36   % 36 bins (180 total) or 3 minutes of awake
                load(procDataFileID,'-mat')
                if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                    % don't include trials with stimulation
                    if isempty(puffs) == true
                        if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                            if strcmp(hemDataType,'Ach_Rhodamine') == true
                                hemDatatype_asleepData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(hemDataType,'NE_Rhodamine') == true
                                hemDatatype_asleepData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(hemDataType,'Ach_GFP') == true
                                hemDatatype_asleepData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(hemDataType,'NE_GFP') == true
                                hemDatatype_asleepData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(hemDataType,'RH_thetaBandPower') == true
                                hemDatatype_asleepData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_LH.thetaBandPower.(strDay).mean;
                            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                                hemDatatype_asleepData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                            end
                            Datatype_asleepData{zz,1} = ProcData.data.Pupil.(dataType);
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(hemDatatype_asleepData) == false
            for bb = 1:length(hemDatatype_asleepData)
                hemDatatype_procAsleepData{bb,1} = detrend(hemDatatype_asleepData{bb,1},'constant');
                Datatype_procAsleepData{bb,1} = detrend(Datatype_asleepData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            hemDatatype_asleepDataMat = zeros(length(hemDatatype_procAsleepData{1,1}),length(hemDatatype_procAsleepData));
            Datatype_asleepDataMat = zeros(length(Datatype_procAsleepData{1,1}),length(Datatype_procAsleepData));
            for cc = 1:length(hemDatatype_procAsleepData)
                hemDatatype_asleepDataMat(:,cc) = hemDatatype_procAsleepData{cc,1};
                Datatype_asleepDataMat(:,cc) = Datatype_procAsleepData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(hemDatatype_asleepDataMat,Datatype_asleepDataMat,params);
            % save results
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).C = C_AsleepData;
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).f = f_AsleepData;
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).confC = confC_AsleepData;
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).cErr = cErr_AsleepData;
        else
            % save results
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).C = [];
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).f = [];
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).confC = [];
            Results_Coherence.(animalID).Asleep.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of all data
        zz = 1;
        clear hemDatatype_allData Datatype_allData hemDatatype_allData Datatype_allData
        hemDatatype_allData = []; Datatype_allData = [];
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
                    if sum(isnan(ProcData.data.Pupil.(dataType))) == 0
                        if strcmp(hemDataType,'Ach_Rhodamine') == true
                            hemDatatype_allData{zz,1} = ProcData.data.Rhodamine.Ach;
                        elseif strcmp(hemDataType,'NE_Rhodamine') == true
                            hemDatatype_allData{zz,1} = ProcData.data.Rhodamine.NE;
                        elseif strcmp(hemDataType,'Ach_GFP') == true
                            hemDatatype_allData{zz,1} = ProcData.data.GFP.NE;
                        elseif strcmp(hemDataType,'NE_GFP') == true
                            hemDatatype_allData{zz,1} = ProcData.data.GFP.NE;
                        elseif strcmp(hemDataType,'RH_thetaBandPower') == true
                            hemDatatype_allData{zz,1} = (ProcData.data.cortical_RH.thetaBandPower - RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.thetaBandPower.(strDay).mean;
                        elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                            hemDatatype_allData{zz,1} = (ProcData.data.cortical_RH.gammaBandPower - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean)./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay).mean;
                        end
                        Datatype_allData{zz,1} = ProcData.data.Pupil.(dataType);
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(hemDatatype_allData) == false
            for bb = 1:length(hemDatatype_allData)
                hemDatatype_procAllData{bb,1} = detrend(hemDatatype_allData{bb,1},'constant');
                Datatype_procAllData{bb,1} = detrend(Datatype_allData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            hemDatatype_allDataMat = zeros(length(hemDatatype_procAllData{1,1}),length(hemDatatype_procAllData));
            Datatype_allDataMat = zeros(length(Datatype_procAllData{1,1}),length(Datatype_procAllData));
            for cc = 1:length(hemDatatype_procAllData)
                hemDatatype_allDataMat(:,cc) = hemDatatype_procAllData{cc,1};
                Datatype_allDataMat(:,cc) = Datatype_procAllData{cc,1};
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            [C_AllData,~,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(hemDatatype_allDataMat,Datatype_allDataMat,params);
            % save results
            Results_Coherence.(animalID).All.(dataType).(hemDataType).C = C_AllData;
            Results_Coherence.(animalID).All.(dataType).(hemDataType).f = f_AllData;
            Results_Coherence.(animalID).All.(dataType).(hemDataType).confC = confC_AllData;
            Results_Coherence.(animalID).All.(dataType).(hemDataType).cErr = cErr_AllData;
        else
            % save results
            Results_Coherence.(animalID).All.(dataType).(hemDataType).C = [];
            Results_Coherence.(animalID).All.(dataType).(hemDataType).f = [];
            Results_Coherence.(animalID).All.(dataType).(hemDataType).confC = [];
            Results_Coherence.(animalID).All.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of NREM
        % pull data from AsleepData.mat structure
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            if strcmp(hemDataType,'Ach_Rhodamine') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(hemDataType,'NE_Rhodamine') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(hemDataType,'Ach_GFP') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(hemDataType,'NE_GFP') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(hemDataType,'RH_thetaBandPower') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.thetaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                [hemDatatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.gammaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            end            
            [Datatype_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);            
        else
            hemDatatype_nremData = [];
            Datatype_nremData = [];
        end
        if isempty(hemDatatype_nremData) == false
            clear hemDatatype_procNremData Datatype_procNremData
            % filter, detrend, and truncate data to minimum length to match events
            for ee = 1:length(hemDatatype_nremData)
                hemDatatype_procNremData{ee,1} = detrend(hemDatatype_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
                Datatype_procNremData{ee,1} = detrend(Datatype_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            zz = 1;
            for ff = 1:length(hemDatatype_procNremData)
                if sum(isnan(Datatype_procNremData{ff,1})) == 0
                    hemDatatype_nremMat(:,zz) = hemDatatype_procNremData{ff,1};
                    Datatype_nremMat(:,zz) = Datatype_procNremData{ff,1};
                    zz = zz + 1;
                end
            end
            % calculate the coherence between desired signals
            params.tapers = [3,5]; % Tapers [n, 2n - 1]
            [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(hemDatatype_nremMat,Datatype_nremMat,params);
            % save results
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).C = C_nrem;
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).f = f_nrem;
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).confC = confC_nrem;
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).cErr = cErr_nrem;
        else
            % save results
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).C = [];
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).f = [];
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).confC = [];
            Results_Coherence.(animalID).NREM.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of REM
        % pull data from AsleepData.mat structure
        if isempty(SleepData.(modelType).REM.data.Pupil) == false
             if strcmp(hemDataType,'Ach_Rhodamine') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(hemDataType,'NE_Rhodamine') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(hemDataType,'Ach_GFP') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(hemDataType,'NE_GFP') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(hemDataType,'RH_thetaBandPower') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.thetaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(hemDataType,'RH_gammaBandPower') == true
                [hemDatatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.gammaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            end            
            [Datatype_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        else
            hemDatatype_remData = [];
            Datatype_remData = [];
        end
        if isempty(hemDatatype_remData) == false
            clear hemDatatype_procRemData Datatype_procRemData
            % filter, detrend, and truncate data to minimum length to match events
            for ee = 1:length(hemDatatype_remData)
                hemDatatype_procRemData{ee,1} = detrend(hemDatatype_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
                Datatype_procRemData{ee,1} = detrend(Datatype_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            zz = 1;
            for ff = 1:length(hemDatatype_procRemData)
                if sum(isnan(Datatype_procRemData{ff,1})) == 0
                    hemDatatype_remMat(:,zz) = hemDatatype_procRemData{ff,1};
                    Datatype_remMat(:,zz) = Datatype_procRemData{ff,1};
                    zz = zz + 1;
                end
            end
            % calculate the coherence between desired signals
            params.tapers = [5,9]; % Tapers [n, 2n - 1]
            [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(hemDatatype_remMat,Datatype_remMat,params);
            % save results
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).C = C_rem;
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).f = f_rem;
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).confC = confC_rem;
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).cErr = cErr_rem;
        else
            % save results
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).C = [];
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).f = [];
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).confC = [];
            Results_Coherence.(animalID).REM.(dataType).(hemDataType).cErr = [];
        end
    end
end
% save data
AnalysisResults.Coherence = Results_Coherence;
cd([rootFolder delim])
save('AnalysisResults.mat','AnalysisResults')

end
