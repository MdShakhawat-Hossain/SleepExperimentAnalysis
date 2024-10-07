function [Results_NeuralHemoCoher] = AnalyzeNeuralHemoCoherence(animalID,rootFolder,Results_NeuralHemoCoher)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral coherence between neural-hemodynamic [CBV] signals (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
hemDataTypes = {'LH','RH'};
modelType = 'Forest';
params.minTime.Rest = 10;
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
% find and load RestingBaselines.mat struct
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% find and load SleepData.mat struct
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% find and load Forest_ScoringResults.mat struct
forestScoringResultsFileID = [animalID '_Forest_ScoringResults.mat'];
load(forestScoringResultsFileID,'-mat')
% lowpass filter
samplingRate = RestData.HbT.LH.HbTCamSamplingRate;
% [z,p,k] = butter(4,1/(samplingRate/2),'low');
% [sos,g] = zp2sos(z,p,k);
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
        [restLogical] = FilterEvents_IOS(RestData.HbT.(hemDataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.HbT.(hemDataType),RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.HbT.(hemDataType).fileIDs(combRestLogical,:);
        restEventTimes = RestData.HbT.(hemDataType).eventTimes(combRestLogical,:);
        restDurations = RestData.HbT.(hemDataType).durations(combRestLogical,:);
        CBV_unstimRestingData = RestData.HbT.(hemDataType).data(combRestLogical,:);
        Gamma_unstimRestingData = RestData.(['cortical_' hemDataType(4:5)]).(dataType).NormData(combRestLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [CBV_finalRestData,~,~,~] = RemoveInvalidData_IOS(CBV_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Gamma_finalRestData,~,~,~] = RemoveInvalidData_IOS(Gamma_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear CBV_ProcRestData Gamma_ProcRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(CBV_finalRestData)
            if length(CBV_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(CBV_finalRestData{bb,1});
                CBV_restPad = (ones(1,restChunkSampleDiff))*CBV_finalRestData{bb,1}(end);
                Gamma_restPad = (ones(1,restChunkSampleDiff))*Gamma_finalRestData{bb,1}(end);
                CBV_ProcRestData{bb,1} = horzcat(CBV_finalRestData{bb,1},CBV_restPad); %#ok<*AGROW>
                Gamma_ProcRestData{bb,1} = horzcat(Gamma_finalRestData{bb,1},Gamma_restPad);
                % CBV_ProcRestData{bb,1} = filtfilt(sos,g,detrend(CBV_ProcRestData{bb,1},'constant'));
                % Gamma_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Gamma_ProcRestData{bb,1},'constant'));
                CBV_ProcRestData{bb,1} = detrend(CBV_ProcRestData{bb,1},'constant');
                Gamma_ProcRestData{bb,1} = detrend(Gamma_ProcRestData{bb,1},'constant');
            else
                % CBV_ProcRestData{bb,1} = filtfilt(sos,g,detrend(CBV_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                % Gamma_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Gamma_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                CBV_ProcRestData{bb,1} = detrend(CBV_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                Gamma_ProcRestData{bb,1} = detrend(Gamma_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        CBV_restData = zeros(length(CBV_ProcRestData{1,1}),length(CBV_ProcRestData));
        Gamma_restData = zeros(length(Gamma_ProcRestData{1,1}),length(Gamma_ProcRestData));
        for cc = 1:length(CBV_ProcRestData)
            CBV_restData(:,cc) = CBV_ProcRestData{cc,1};
            Gamma_restData(:,cc) = Gamma_ProcRestData{cc,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(CBV_restData,Gamma_restData,params);
        % save results
        Results_NeuralHemoCoher.(animalID).Rest.(dataType).(hemDataType).C = C_RestData;
        Results_NeuralHemoCoher.(animalID).Rest.(dataType).(hemDataType).f = f_RestData;
        Results_NeuralHemoCoher.(animalID).Rest.(dataType).(hemDataType).confC = confC_RestData;
        Results_NeuralHemoCoher.(animalID).Rest.(dataType).(hemDataType).cErr = cErr_RestData;
        %% analyze neural-hemo coherence during periods of alert
        zz = 1;
        clear CBV_AwakeData Gamma_AwakeData CBV_ProcAwakeData Gamma_ProcAwakeData
        CBV_AwakeData = [];
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
            if sum(strcmp(scoringLabels,'Not Sleep')) > 139   % 36 bins (180 total) or 3 minutes of sleep
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        CBV_AwakeData{zz,1} = ProcData.data.HbT.(hemDataType);
                        Gamma_AwakeData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(CBV_AwakeData) == false
            for bb = 1:length(CBV_AwakeData)
                % CBV_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(CBV_AwakeData{bb,1},'constant'));
                % Gamma_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(Gamma_AwakeData{bb,1},'constant'));
                CBV_ProcAwakeData{bb,1} = detrend(CBV_AwakeData{bb,1},'constant');
                Gamma_ProcAwakeData{bb,1} = detrend(Gamma_AwakeData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            CBV_awakeData = zeros(length(CBV_ProcAwakeData{1,1}),length(CBV_ProcAwakeData));
            Gamma_awakeData = zeros(length(Gamma_ProcAwakeData{1,1}),length(Gamma_ProcAwakeData));
            for cc = 1:length(CBV_ProcAwakeData)
                CBV_awakeData(:,cc) = CBV_ProcAwakeData{cc,1};
                Gamma_awakeData(:,cc) = Gamma_ProcAwakeData{cc,1};
            end
            % calculate the coherence between desired signals
            [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc(CBV_awakeData,Gamma_awakeData,params);
            % save results
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).C = C_AwakeData;
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).f = f_AwakeData;
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).confC = confC_AwakeData;
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).cErr = cErr_AwakeData;
        else
            % save results
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).C = [];
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).f = [];
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).confC = [];
            Results_NeuralHemoCoher.(animalID).Awake.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of asleep
        zz = 1;
        clear CBV_SleepData Gamma_SleepData CBV_ProcSleepData Gamma_ProcSleepData
        CBV_SleepData = [];
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
            if sum(strcmp(scoringLabels,'Not Sleep')) < 139   % 36 bins (180 total) or 3 minutes of awake
                load(procDataFileID,'-mat')
                puffs = ProcData.data.stimulations.LPadSol;
                % don't include trials with stimulation
                if isempty(puffs) == true
                    motionArtifact = ProcData.notes.motionArtifact;
                    if motionArtifact == false
                        CBV_SleepData{zz,1} = ProcData.data.HbT.(hemDataType);
                        Gamma_SleepData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
        end
        % filter and detrend data
        if isempty(CBV_SleepData) == false
            for bb = 1:length(CBV_SleepData)
                % CBV_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(CBV_SleepData{bb,1},'constant'));
                % Gamma_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(Gamma_SleepData{bb,1},'constant'));
                CBV_ProcSleepData{bb,1} = detrend(CBV_SleepData{bb,1},'constant');
                Gamma_ProcSleepData{bb,1} = detrend(Gamma_SleepData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            CBV_sleepData = zeros(length(CBV_ProcSleepData{1,1}),length(CBV_ProcSleepData));
            Gamma_sleepData = zeros(length(Gamma_ProcSleepData{1,1}),length(Gamma_ProcSleepData));
            for cc = 1:length(CBV_ProcSleepData)
                CBV_sleepData(:,cc) = CBV_ProcSleepData{cc,1};
                Gamma_sleepData(:,cc) = Gamma_ProcSleepData{cc,1};
            end
            % calculate the coherence between desired signals
            [C_SleepData,~,~,~,~,f_SleepData,confC_SleepData,~,cErr_SleepData] = coherencyc(CBV_sleepData,Gamma_sleepData,params);
            % save results
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).C = C_SleepData;
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).f = f_SleepData;
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).confC = confC_SleepData;
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).cErr = cErr_SleepData;
        else
            % save results
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).C = [];
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).f = [];
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).confC = [];
            Results_NeuralHemoCoher.(animalID).Sleep.(dataType).(hemDataType).cErr = [];
        end
        %% analyze neural-hemo coherence during periods of all data
        zz = 1;
        clear CBV_AllUnstimData Gamma_AllUnstimData CBV_ProcAllUnstimData Gamma_ProcAllUnstimData
        CBV_AllUnstimData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allUnstimDataFileDate,~] = GetFileInfo_IOS(procDataFileID);
            strDay = ConvertDate_IOS(allUnstimDataFileDate);
            load(procDataFileID,'-mat')
            puffs = ProcData.data.stimulations.LPadSol;
            % don't include trials with stimulation
            if isempty(puffs) == true
                motionArtifact = ProcData.notes.motionArtifact;
                if motionArtifact == false
                    CBV_AllUnstimData{zz,1} = ProcData.data.HbT.(hemDataType);
                    Gamma_AllUnstimData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                    zz = zz + 1;
                end
            end
        end
        % filter and detrend data
        if isempty(CBV_AllUnstimData) == false
            for bb = 1:length(CBV_AllUnstimData)
                % CBV_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(CBV_AllUnstimData{bb,1},'constant'));
                % Gamma_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(Gamma_AllUnstimData{bb,1},'constant'));
                CBV_ProcAllUnstimData{bb,1} = detrend(CBV_AllUnstimData{bb,1},'constant');
                Gamma_ProcAllUnstimData{bb,1} = detrend(Gamma_AllUnstimData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            CBV_allUnstimData = zeros(length(CBV_ProcAllUnstimData{1,1}),length(CBV_ProcAllUnstimData));
            Gamma_allUnstimData = zeros(length(Gamma_ProcAllUnstimData{1,1}),length(Gamma_ProcAllUnstimData));
            for cc = 1:length(CBV_ProcAllUnstimData)
                CBV_allUnstimData(:,cc) = CBV_ProcAllUnstimData{cc,1};
                Gamma_allUnstimData(:,cc) = Gamma_ProcAllUnstimData{cc,1};
            end
            % calculate the coherence between desired signals
            [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc(CBV_allUnstimData,Gamma_allUnstimData,params);
            % save results
            Results_NeuralHemoCoher.(animalID).All.(dataType).(hemDataType).C = C_AllUnstimData;
            Results_NeuralHemoCoher.(animalID).All.(dataType).(hemDataType).f = f_AllUnstimData;
            Results_NeuralHemoCoher.(animalID).All.(dataType).(hemDataType).confC = confC_AllUnstimData;
            Results_NeuralHemoCoher.(animalID).All.(dataType).(hemDataType).cErr = cErr_AllUnstimData;
        end
        %% analyze neural-hemo coherence during periods of NREM
        % pull data from SleepData.mat structure
        [HbTremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.(hemDataType(4:5)),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        [Gamma_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.(['cortical_' hemDataType(4:5)]).(dataType),SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
        % filter, detrend, and truncate data to minimum length to match events
        for ee = 1:length(HbTremData)
            % HbTremData{ee,1} = filtfilt(sos,g,detrend(HbTremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            % Gamma_nremData{ee,1} = filtfilt(sos,g,detrend(Gamma_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            HbTremData{ee,1} = detrend(HbTremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            Gamma_nremData{ee,1} = detrend(Gamma_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        HbTrem = zeros(length(HbTremData{1,1}),length(HbTremData));
        Gamma_nrem = zeros(length(Gamma_nremData{1,1}),length(Gamma_nremData));
        for ff = 1:length(HbTremData)
            HbTrem(:,ff) = HbTremData{ff,1};
            Gamma_nrem(:,ff) = Gamma_nremData{ff,1};
        end
        % calculate the coherence between desired signals
        [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(HbTrem,Gamma_nrem,params);
        % save results
        Results_NeuralHemoCoher.(animalID).NREM.(dataType).(hemDataType).C = C_nrem;
        Results_NeuralHemoCoher.(animalID).NREM.(dataType).(hemDataType).f = f_nrem;
        Results_NeuralHemoCoher.(animalID).NREM.(dataType).(hemDataType).confC = confC_nrem;
        Results_NeuralHemoCoher.(animalID).NREM.(dataType).(hemDataType).cErr = cErr_nrem;
        %% analyze neural-hemo coherence during periods of REM
        % pull data from SleepData.mat structure
        [CBV_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.(hemDataType(4:5)),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        [Gamma_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.(['cortical_' hemDataType(4:5)]).(dataType),SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
        % filter, detrend, and truncate data to minimum length to match events
        for gg = 1:length(CBV_remData)
            % CBV_remData{gg,1} = filtfilt(sos,g,detrend(CBV_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            % Gamma_remData{gg,1} = filtfilt(sos,g,detrend(Gamma_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            CBV_remData{gg,1} = detrend(CBV_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
            Gamma_remData{gg,1} = detrend(Gamma_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant');
        end
        % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        CBV_rem = zeros(length(CBV_remData{1,1}),length(CBV_remData));
        Gamma_rem = zeros(length(Gamma_remData{1,1}),length(Gamma_remData));
        for hh = 1:length(CBV_remData)
            CBV_rem(:,hh) = CBV_remData{hh,1};
            Gamma_rem(:,hh) = Gamma_remData{hh,1};
        end
        % calculate the coherence between desired signals
        [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(CBV_rem,Gamma_rem,params);
        % save results
        Results_NeuralHemoCoher.(animalID).REM.(dataType).(hemDataType).C = C_rem;
        Results_NeuralHemoCoher.(animalID).REM.(dataType).(hemDataType).f = f_rem;
        Results_NeuralHemoCoher.(animalID).REM.(dataType).(hemDataType).confC = confC_rem;
        Results_NeuralHemoCoher.(animalID).REM.(dataType).(hemDataType).cErr = cErr_rem;
    end
    
end
% save data
cd(rootFolder)
save('Results_NeuralHemoCoher.mat','Results_NeuralHemoCoher')

end
