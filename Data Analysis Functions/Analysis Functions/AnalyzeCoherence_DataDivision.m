function [AnalysisResults] = AnalyzeCoherence_DataDivision(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the spectral coherence between GRAB sensors - CBV
%________________________________________________________________________________________________________________________

%% function parameters
dataTypes = {'NE_GFP','ACh_GFP'};
pupilDataTypes = {'ACh_CBV','NE_CBV','ACh_GFP','NE_GFP'};
modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 60;
%% only run analysis for valid animal IDs

    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end

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
scoringResultsFileStruct = dir('*Manual_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% lowpass filter
try
    samplingRate = RestData.GFP.P_ACh.CBVCamSamplingRate;
catch
    samplingRate = RestData.GFP.P_ACh.RhodamineCamSamplingRate;
end
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {10};%5
% go through eACh valid data type for arousal-based coherence analysis
for zzz = 1:length(pupilDataTypes)
    pupilDataType = pupilDataTypes{1,zzz};
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        %% analyze GRAB sensors - CBV coherence during periods of rest

        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_IOS(RestData.GFP.P_NE,RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.GFP.P_NE,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.GFP.P_NE.fileIDs(combRestLogical,:);
        restEventTimes = RestData.GFP.P_NE.eventTimes(combRestLogical,:);
        restDurations = RestData.GFP.P_NE.durations(combRestLogical,:);

        if strcmp(dataType,'ACh_CBV') == true
             Pupil_restData = RestData.CBV.P_ACh.data(combRestLogical,1);
        elseif strcmp(dataType,'NE_CBV') == true
            Pupil_restData = RestData.CBV.P_NE.data(combRestLogical,:);
        elseif strcmp(dataType,'ACh_GFP') == true
            Pupil_restData = RestData.GFP.P_ACh.data(combRestLogical,:);
        elseif strcmp(dataType,'NE_GFP') == true
            Pupil_restData = RestData.GFP.P_NE.data(combRestLogical,:);
        end
        

        if strcmp(pupilDataType,'ACh_CBV') == true
             restData = RestData.CBV.P_ACh.data(combRestLogical,1);
        elseif strcmp(pupilDataType,'NE_CBV') == true
            restData = RestData.CBV.P_NE.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'ACh_GFP') == true
            restData = RestData.GFP.P_ACh.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'NE_GFP') == true
            restData = RestData.GFP.P_NE.data(combRestLogical,:);
        end

        % keep only the data that occurs within the manually-approved awake regions
        [finalRestData,~,~,~] = RemoveInvalidData_IOS(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [Pupil_FinalRestData,~,~,~] = RemoveInvalidData_IOS(Pupil_restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        clear procRestData Pupil_procRestData
        % filter, detrend, and truncate data to minimum length to match events
        for bb = 1:length(finalRestData)
            if length(finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{bb,1});
                restPad = (ones(1,restChunkSampleDiff))*finalRestData{bb,1}(end);
                Pupil_restPad = (ones(1,restChunkSampleDiff))*Pupil_FinalRestData{bb,1}(end);
                procRestData{bb,1} = horzcat(finalRestData{bb,1},restPad); %#ok<*AGROW>
                Pupil_procRestData{bb,1} = horzcat(Pupil_FinalRestData{bb,1},Pupil_restPad);
                procRestData{bb,1} = detrend(procRestData{bb,1},'constant');
                Pupil_procRestData{bb,1} = detrend(Pupil_procRestData{bb,1},'constant');
            else
                procRestData{bb,1} = detrend(finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
                Pupil_procRestData{bb,1} = detrend(Pupil_FinalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant');
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        zz = 1;
        for cc = 1:length(procRestData)
            if sum(isnan(Pupil_procRestData{cc,1})) == 0
                restDataMat(:,zz) = procRestData{cc,1};
                Pupil_restDataMat(:,zz) = Pupil_procRestData{cc,1};
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
            if size(restDataMat,1) > size(Pupil_restDataMat,1)
                restDataMat = restDataMat(1:size(Pupil_restDataMat,1),:);
            elseif size(restDataMat,1) < size(Pupil_restDataMat,1)
                Pupil_restDataMat = Pupil_restDataMat(1:size(restDataMat,1),:);
            end
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc(restDataMat,Pupil_restDataMat,params);
        % save results
        AnalysisResults.(animalID).Coherence.Rest.(dataType).(pupilDataType).C = C_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).(pupilDataType).f = f_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).(pupilDataType).confC = confC_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).(pupilDataType).cErr = cErr_RestData;
        %% analyze GRAB sensors - CBV coherence during periods of alert
        zz = 1;
        clear awakeData awakePupilData procAwakeData Pupil_procAwakeData
        awakeData = []; awakePupilData = [];
        
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,awakeDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(awakeDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
           % cut the data into three 5 minutes chunk
            ScoreSize = 10*12;% 10 minutes X 12 5 sec periods per minutes
           
            load(procDataFileID,'-mat')
             DataSize = 10*60*ProcData.notes.dsFs; % 10 minutes X 60 secs X 30 Hz;
                % only run on files with good pupil measurement
                % if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
                    try
                        puffs = ProcData.data.stimulations.LPadSol;
                    catch
                        puffs = ProcData.data.solenoids.LPadSol;
                    end
                        % don't include trials with stimulation
                        if isempty(puffs) == true
                            if sum(isnan(ProcData.data.Pupil.zDiameter)) == 0
                               for SSL =  1:1:5 % 10 for 5 seconds bin
                                ScoreLabels = scoringLabels(((SSL-1)*ScoreSize)+1:(SSL*ScoreSize));
                                  % check labels to match arousal state
                                    if sum(strcmp(ScoreLabels,'Not Sleep')) > 96 % 80% of the time awake                                
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
            
                                        % pull data based on pupil data type
                                        if strcmp(pupilDataType,'ACh_CBV') == true
                                            awakePupilData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'NE_CBV') == true
                                            awakePupilData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'ACh_GFP') == true
                                            awakePupilData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'NE_GFP') == true
                                            awakePupilData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        end
                                        zz = zz + 1;
                                    end
                               end
                            end
                        end
                % end
        end

        % filter and detrend data
        if isempty(awakeData) == false
            for bb = 1:length(awakeData)
                procAwakeData{bb,1} = detrend(awakeData{bb,1},'constant');
                Pupil_procAwakeData{bb,1} = detrend(awakePupilData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            awakeDataMat = zeros(length(procAwakeData{1,1}),length(procAwakeData));
            Pupil_awakeDataMat = zeros(length(Pupil_procAwakeData{1,1}),length(Pupil_procAwakeData));
            for cc = 1:length(procAwakeData)
                awakeDataMat(:,cc) = procAwakeData{cc,1}(1:length(procAwakeData{1,1}));
                Pupil_awakeDataMat(:,cc) = Pupil_procAwakeData{cc,1}(1:length(Pupil_procAwakeData{1,1}));
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            if size(awakeDataMat,1) > size(Pupil_awakeDataMat,1)
                awakeDataMat = awakeDataMat(1:size(Pupil_awakeDataMat,1),:);
            elseif size(awakeDataMat,1) < size(Pupil_awakeDataMat,1)
                Pupil_awakeDataMat = Pupil_awakeDataMat(1:size(awakeDataMat,1),:);
            end
            [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc(awakeDataMat,Pupil_awakeDataMat,params);
            % save results
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).C = C_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).f = f_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).confC = confC_AwakeData;
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).cErr = cErr_AwakeData;
        else
            % save results
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).C = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).f = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).confC = [];
            AnalysisResults.(animalID).Coherence.Awake.(dataType).(pupilDataType).cErr = [];
        end
        %% analyze GRAB sensors - CBVo coherence during periods of aasleep
        zz = 1;
        clear asleepData asleepPupilData procAsleepData Pupil_procAsleepData
        asleepData = []; asleepPupilData = [];

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
                % if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
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
            
                                        % pull data based on pupil data type
                                        if strcmp(pupilDataType,'ACh_CBV') == true
                                            asleepPupilData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'NE_CBV') == true
                                            asleepPupilData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'ACh_GFP') == true
                                            asleepPupilData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        elseif strcmp(pupilDataType,'NE_GFP') == true
                                            asleepPupilData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                        end
                                        zz = zz + 1;
                                    end
                               end
                            end
                        end
                % end
        end

        % filter and detrend data
        if isempty(asleepData) == false
            for bb = 1:length(asleepData)
                procAsleepData{bb,1} = detrend(asleepData{bb,1},'constant');
                Pupil_procAsleepData{bb,1} = detrend(asleepPupilData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            asleepDataMat = zeros(length(procAsleepData{1,1}),length(procAsleepData));
            Pupil_asleepDataMat = zeros(length(Pupil_procAsleepData{1,1}),length(Pupil_procAsleepData));
            for cc = 1:length(procAsleepData)
                asleepDataMat(:,cc) = procAsleepData{cc,1}(1:length(procAsleepData{1,1}));
                Pupil_asleepDataMat(:,cc) = Pupil_procAsleepData{cc,1}(1:length(Pupil_procAsleepData{1,1}));
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            if size(asleepDataMat,1) > size(Pupil_asleepDataMat,1)
                asleepDataMat = asleepDataMat(1:size(Pupil_asleepDataMat,1),:);
            elseif size(asleepDataMat,1) < size(Pupil_asleepDataMat,1)
                Pupil_asleepDataMat = Pupil_asleepDataMat(1:size(asleepDataMat,1),:);
            end
            [C_AsleepData,~,~,~,~,f_AsleepData,confC_AsleepData,~,cErr_AsleepData] = coherencyc(asleepDataMat,Pupil_asleepDataMat,params);
            % save results
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).C = C_AsleepData;
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).f = f_AsleepData;
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).confC = confC_AsleepData;
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).cErr = cErr_AsleepData;
        else
            % save results
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).C = [];
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).f = [];
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).confC = [];
            AnalysisResults.(animalID).Coherence.Asleep.(dataType).(pupilDataType).cErr = [];
        end
        %% analyze GRAB sensors - CBV coherence during periods of all data
        zz = 1;
        clear allData allPupilData procAllData Pupil_procAllData
        allData = []; allPupilData = [];


        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,~] = GetFileInfo_JNeurosci2022(procDataFileID);
            load(procDataFileID,'-mat')

            DataSize = 10*60*ProcData.notes.dsFs; % 10 minutes X 60 secs X 30 Hz;
            % only run on files with good pupil measurement
            % if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
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
    
                                % pull data based on data type
                                if strcmp(pupilDataType,'ACh_CBV') == true
                                    allPupilData{zz,1} = ProcData.data.CBV.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(pupilDataType,'NE_CBV') == true
                                    allPupilData{zz,1} = ProcData.data.CBV.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(pupilDataType,'ACh_GFP') == true
                                    allPupilData{zz,1} = ProcData.data.GFP.P_ACh(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                elseif strcmp(pupilDataType,'NE_GFP') == true
                                    allPupilData{zz,1} = ProcData.data.GFP.P_NE(((SSL-1)*DataSize)+1:(SSL*DataSize));
                                end
                                zz = zz + 1;
                            end
                       end
                   end
           % end 
        end

        % filter and detrend data
        if isempty(allData) == false
            for bb = 1:length(allData)
                procAllData{bb,1} = detrend(allData{bb,1},'constant');
                Pupil_procAllData{bb,1} = detrend(allPupilData{bb,1},'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            allDataMat = zeros(length(procAllData{1,1}),length(procAllData));
            Pupil_allDataMat = zeros(length(Pupil_procAllData{1,1}),length(Pupil_procAllData));
            for cc = 1:length(procAllData)
                allDataMat(:,cc) = procAllData{cc,1}(1:length(procAllData{1,1}));
                Pupil_allDataMat(:,cc) = Pupil_procAllData{cc,1}(1:length(Pupil_procAllData{1,1}));
            end
            % calculate the coherence between desired signals
            params.tapers = [10,19]; % Tapers [n, 2n - 1]
            if size(allDataMat,1) > size(Pupil_allDataMat,1)
                allDataMat = allDataMat(1:size(Pupil_allDataMat,1),:);
            elseif size(allDataMat,1) < size(Pupil_allDataMat,1)
                Pupil_allDataMat = Pupil_allDataMat(1:size(allDataMat,1),:);
            end
            
            [C_AllData,~,~,~,~,f_AllData,confC_AllData,~,cErr_AllData] = coherencyc(allDataMat,Pupil_allDataMat,params);
            % save results
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).C = C_AllData;
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).f = f_AllData;
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).confC = confC_AllData;
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).cErr = cErr_AllData;
        else
            % save results
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).C = [];
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).f = [];
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).confC = [];
            AnalysisResults.(animalID).Coherence.All.(dataType).(pupilDataType).cErr = [];
        end
        %% analyze GRAB sensors - CBV coherence during periods of NREM
        % pull data from AsleepData.mat structure
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            
            if strcmp(pupilDataType,'ACh_CBV') == true
                [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'NE_CBV') == true
                [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'ACh_GFP') == true
                [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'NE_GFP') == true
                [nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            end           

            if strcmp(dataType,'ACh_CBV') == true
                [Pupil_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'NE_CBV') == true
                [Pupil_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.CBV.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'ACh_GFP') == true
                [Pupil_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_ACh,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'NE_GFP') == true
                [Pupil_nremData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.P_NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            end        
        else
            nremData = [];
            Pupil_nremData = [];
        end
        if isempty(nremData) == false
            clear procNremData Pupil_procNremData
            % filter, detrend, and truncate data to minimum length to match events
            for ee = 1:length(nremData)
                procNremData{ee,1} = detrend(nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
                Pupil_procNremData{ee,1} = detrend(Pupil_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            end
            % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            zz = 1;
            for ff = 1:length(procNremData)
                if sum(isnan(Pupil_procNremData{ff,1})) == 0
                    nremMat(:,zz) = procNremData{ff,1};
                    Pupil_nremMat(:,zz) = Pupil_procNremData{ff,1};
                    zz = zz + 1;
                end
            end
            % calculate the coherence between desired signals
            params.tapers = [3,5]; % Tapers [n, 2n - 1]
            if size(nremMat,1) > size(Pupil_nremMat,1)
                nremMat = nremMat(1:size(Pupil_nremMat,1),:);
            elseif size(nremMat,1) < size(Pupil_nremMat,1)
                Pupil_nremMat = Pupil_nremMat(1:size(nremMat,1),:);
            end
            [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc(nremMat,Pupil_nremMat,params);
            % save results
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).C = C_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).f = f_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).confC = confC_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).cErr = cErr_nrem;
        else
            % save results
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).C = [];
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).f = [];
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).confC = [];
            AnalysisResults.(animalID).Coherence.NREM.(dataType).(pupilDataType).cErr = [];
        end
        %% analyze GRAB sensors - CBV coherence during periods of REM
        % pull data from AsleepData.mat structure
        
         if firstHrs == "false"
            if isfield(SleepData.(modelType),'REM')
                if strcmp(pupilDataType,'ACh_CBV') == true
                    [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(pupilDataType,'NE_CBV') == true
                    [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(pupilDataType,'ACh_GFP') == true
                    [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(pupilDataType,'NE_GFP') == true
                    [remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                end  
                if strcmp(dataType,'ACh_CBV') == true
                    [Pupil_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(dataType,'NE_CBV') == true
                    [Pupil_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.CBV.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(dataType,'ACh_GFP') == true
                    [Pupil_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_ACh,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                elseif strcmp(dataType,'NE_GFP') == true
                    [Pupil_remData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.P_NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
                end
            else
                remData = [];
                Pupil_remData = [];
            end
            if isempty(remData) == false
                clear procRemData Pupil_procRemData
                % filter, detrend, and truncate data to minimum length to match events
                for ee = 1:length(remData)
                    procRemData{ee,1} = detrend(remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
                    Pupil_procRemData{ee,1} = detrend(Pupil_remData{ee,1}(1:(params.minTime.REM*samplingRate)),'constant');
                end
                % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
                zz = 1;
                for ff = 1:length(procRemData)
                    if sum(isnan(Pupil_procRemData{ff,1})) == 0
                        remMat(:,zz) = procRemData{ff,1};
                        Pupil_remMat(:,zz) = Pupil_procRemData{ff,1};
                        zz = zz + 1;
                    end
                end
                % calculate the coherence between desired signals
                params.tapers = [5,9]; % Tapers [n, 2n - 1]
                if size(remMat,1) > size(Pupil_remMat,1)
                    remMat = remMat(1:size(Pupil_remMat,1),:);
                elseif size(remMat,1) < size(Pupil_remMat,1)
                    Pupil_remMat = Pupil_remMat(1:size(remMat,1),:);
                end
                [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc(remMat,Pupil_remMat,params);
                % save results
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).C = C_rem;
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).f = f_rem;
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).confC = confC_rem;
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).cErr = cErr_rem;
            else
                % save results
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).C = [];
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).f = [];
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).confC = [];
                AnalysisResults.(animalID).Coherence.REM.(dataType).(pupilDataType).cErr = [];
            end
         end
         %}
    end
end
 % save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end
end
