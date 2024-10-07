function [AnalysisResults] = AnalyzeCrossCorrelation_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Analyze the cross-correlation between neural activity/hemodynamics and pupil diameter
%________________________________________________________________________________________________________________________

% function parameters & data types
dataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};
pupilDataTypes = {'zDiameter','Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','RH_thetaBandPower','RH_alphaBandPower','RH_gammaBandPower'};
modelType = 'Manual';
params.minTime.Rest = 10;
params.minTime.NREM = 30;
params.minTime.REM = 30;
% go to animal's data location    
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end
    cd(dataLocation)
% identify and load RestData.mat struct
restDataFileStruct = dir('*_RestData.mat');
restDataFile = {restDataFileStruct.name}';
restDataFileID = char(restDataFile);
load(restDataFileID,'-mat')
% identify and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID,'-mat')
% identify and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID,'-mat')
% identify and load SleepData.mat strut
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
% identify and load ScoringResults.mat
scoringResultsFileStruct = dir('*Manual_ScoringResults.mat');
scoringResultsFile = {scoringResultsFileStruct.name}';
scoringResultsFileID = char(scoringResultsFile);
load(scoringResultsFileID,'-mat')
% character list of all ProcData file IDs
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% filter characteristics & resting criteria
samplingRate = RestData.Rhodamine.Ach.RhodamineCamSamplingRate;
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% criteria for resting
RestCriteria.Fieldname = {'durations'};
RestCriteria.Comparison = {'gt'};
RestCriteria.Value = {params.minTime.Rest};
RestPuffCriteria.Fieldname = {'puffDistances'};
RestPuffCriteria.Comparison = {'gt'};
RestPuffCriteria.Value = {5};
% go through each valid data type for arousal-dependent cross-correlation analysis
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    for bb = 1:length(pupilDataTypes)
        pupilDataType = pupilDataTypes{1,bb};
        %% cross-correlation analysis for resting data
        % pull data from RestData.mat structure
        [restLogical] = FilterEvents_JNeurosci2022(RestData.Pupil.zDiameter,RestCriteria);
        [puffLogical] = FilterEvents_JNeurosci2022(RestData.Pupil.zDiameter,RestPuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.Pupil.zDiameter.fileIDs(combRestLogical,:);
        restDurations = RestData.Pupil.zDiameter.durations(combRestLogical,:);
        restEventTimes = RestData.Pupil.zDiameter.eventTimes(combRestLogical,:);

        if strcmp(pupilDataType,'Ach_Rhodamine') == true
            pupilRestData = RestData.Rhodamine.Ach.data(combRestLogical,1);
        elseif strcmp(pupilDataType,'NE_Rhodamine') == true
            pupilRestData = RestData.Rhodamine.NE.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'Ach_GFP') == true
            pupilRestData = RestData.GFP.Ach.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'NE_GFP') == true
            pupilRestData = RestData.GFP.NE.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'zDiameter') == true
            pupilRestData = RestData.Pupil.zDiameter.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
            pupilRestData = RestData.cortical_RH.thetaBandPower.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
            pupilRestData = RestData.cortical_RH.alphaBandPower.data(combRestLogical,:);
        elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
            pupilRestData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
        end

        if strcmp(dataType,'Ach_Rhodamine') == true
            restData = RestData.Rhodamine.Ach.data(combRestLogical,1);
        elseif strcmp(dataType,'NE_Rhodamine') == true
            restData = RestData.Rhodamine.NE.data(combRestLogical,:);
        elseif strcmp(dataType,'Ach_GFP') == true
            restData = RestData.GFP.Ach.data(combRestLogical,:);
        elseif strcmp(dataType,'NE_GFP') == true
            restData = RestData.GFP.NE.data(combRestLogical,:);
        elseif strcmp(dataType,'zDiameter') == true
            restData = RestData.Pupil.zDiameter.data(combRestLogical,:);
        elseif strcmp(dataType,'RH_thetaBandPower') == true
            restData = RestData.cortical_RH.thetaBandPower.data(combRestLogical,:);
        elseif strcmp(dataType,'RH_alphaBandPower') == true
            restData = RestData.cortical_RH.alphaBandPower.data(combRestLogical,:);
        elseif strcmp(dataType,'RH_gammaBandPower') == true
            restData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
        end

        % keep only the data that occurs within the manually-approved alert regions
        [finalRestData,~,~,~] = RemoveInvalidData_JNeurosci2022(restData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [finalPupilRestData,~,~,~] = RemoveInvalidData_JNeurosci2022(pupilRestData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % process, filter + detrend each array
        catRestData = [];
        cc = 1;
        for dd = 1:length(finalRestData)
            if sum(isnan(finalPupilRestData{dd,1})) == 0
                if length(finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                    restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{bb,1});
                    restPad = (ones(1,restChunkSampleDiff))*finalRestData{bb,1}(end);
                    restPupilPad = (ones(1,restChunkSampleDiff))*finalPupilRestData{bb,1}(end);
                    procRestData = horzcat(finalRestData{bb,1},restPad);
                    procPupilRestData = horzcat(finalPupilRestData{bb,1},restPupilPad);
                    procRestData = filtfilt(sos,g,detrend(procRestData{bb,1},'constant'));
                    procPupilRestData = filtfilt(sos,g,detrend(procPupilRestData{bb,1},'constant'));
                else
                    procRestData = filtfilt(sos,g,detrend(finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                    procPupilRestData = filtfilt(sos,g,detrend(finalPupilRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                end
                catRestData.data{cc,1} = procRestData;
                catRestData.pupil{cc,1} = procPupilRestData;
                cc = cc + 1;
            end
        end
        % run cross correlation between data types
        restXcVals = [];
        if isempty(catRestData) == false
            restLagTime = 5; % seconds
            restFrequency = samplingRate; % Hz
            restMaxLag = restLagTime*restFrequency;
            % run cross-correlation analysis - average through time
            for ee = 1:length(catRestData.data)
                restArray = catRestData.data{ee,1};
                restPupilarray = catRestData.pupil{ee,1};
                [restXcVals(ee,:),restPupilLags] = xcorr(restArray,restPupilarray,restMaxLag,'coeff');
            end
            restMeanXcVals = mean(restXcVals,1);
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Rest.(dataType).(pupilDataType).lags = restPupilLags;
            AnalysisResults.(animalID).CrossCorrelation.Rest.(dataType).(pupilDataType).xcVals = restMeanXcVals;
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Rest.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.Rest.(dataType).(pupilDataType).xcVals = [];
        end
        %% cross-correlation analysis for alert data
        zz = 1;
        alertData = []; alertPupilData = []; alertProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,alertDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(alertDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) > 500 % 80% of the time awake % 144 % 12 minutes of awake
                load(procDataFileID,'-mat')
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
                            % pull data based on data type
                            if strcmp(dataType,'Ach_Rhodamine') == true
                                alertData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(dataType,'NE_Rhodamine') == true
                                alertData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(dataType,'Ach_GFP') == true
                                alertData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(dataType,'NE_GFP') == true
                                alertData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(dataType,'zDiameter') == true
                                alertData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(dataType,'RH_thetaBandPower') == true
                                alertData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(dataType,'RH_alphaBandPower') == true
                                alertData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(dataType,'RH_gammaBandPower') == true
                                alertData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end

                            % pull data based on pupil data type
                            if strcmp(pupilDataType,'Ach_Rhodamine') == true
                                alertPupilData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(pupilDataType,'NE_Rhodamine') == true
                                alertPupilData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(pupilDataType,'Ach_GFP') == true
                                alertPupilData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(pupilDataType,'NE_GFP') == true
                                alertPupilData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(pupilDataType,'zDiameter') == true
                                alertPupilData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
                                alertPupilData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
                                alertPupilData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
                                alertPupilData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % process, filter + detrend each array
        alertXcVals = [];
        if isempty(alertData) == false
            for ee = 1:length(alertData)
                alertProcData.data{ee,1} = filtfilt(sos,g,detrend(alertData{ee,1},'constant'));
                alertProcData.pupil{ee,1} = filtfilt(sos,g,detrend(alertPupilData{ee,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            alertLagTime = 30; % seconds
            alertFrequency = samplingRate; % Hz
            alertMaxLag = alertLagTime*alertFrequency;
            % run cross-correlation analysis - average through time
            for ff = 1:length(alertProcData.data)
                alertArray = alertProcData.data{ff,1};
                alertPupilarray = alertProcData.pupil{ff,1};
                if length(alertArray) ~=  length(alertPupilarray)
                    minlength = min(length(alertArray),length(alertPupilarray));
                    alertArray = alertArray(1:minlength);
                    alertPupilarray = alertPupilarray(1:minlength);
                end
                [alertXcVals(ff,:),alertPupilLags] = xcorr(alertArray,alertPupilarray,alertMaxLag,'coeff');
            end
            alertMeanXcVals = mean(alertXcVals,1);
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Alert.(dataType).(pupilDataType).lags = alertPupilLags;
            AnalysisResults.(animalID).CrossCorrelation.Alert.(dataType).(pupilDataType).xcVals = alertMeanXcVals;
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Alert.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.Alert.(dataType).(pupilDataType).xcVals = [];
        end
        %% cross-correlation analysis for asleep data
        zz = 1;
        asleepData = []; asleepPupilData = []; asleepProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,asleepDataFileID] = GetFileInfo_JNeurosci2022(procDataFileID);
            scoringLabels = [];
            for dd = 1:length(ScoringResults.fileIDs)
                if strcmp(asleepDataFileID,ScoringResults.fileIDs{dd,1}) == true
                    scoringLabels = ScoringResults.labels{dd,1};
                end
            end
            % check labels to match arousal state
            if sum(strcmp(scoringLabels,'Not Sleep')) < 124 % 36 % 80% of the time sleep 12 minutes of asleep
                load(procDataFileID,'-mat')
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
                            % pull data based on data type
                            if strcmp(dataType,'Ach_Rhodamine') == true
                                asleepData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(dataType,'NE_Rhodamine') == true
                                asleepData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(dataType,'Ach_GFP') == true
                                asleepData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(dataType,'NE_GFP') == true
                                asleepData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(dataType,'zDiameter') == true
                                asleepData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(dataType,'RH_thetaBandPower') == true
                                asleepData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(dataType,'RH_alphaBandPower') == true
                                asleepData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(dataType,'RH_gammaBandPower') == true
                                asleepData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end

                            % pull data based on pupil data type
                            if strcmp(pupilDataType,'Ach_Rhodamine') == true
                                asleepPupilData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(pupilDataType,'NE_Rhodamine') == true
                                asleepPupilData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(pupilDataType,'Ach_GFP') == true
                                asleepPupilData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(pupilDataType,'NE_GFP') == true
                                asleepPupilData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(pupilDataType,'zDiameter') == true
                                asleepPupilData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
                                asleepPupilData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
                                asleepPupilData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
                                asleepPupilData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end
                            zz = zz + 1;
                        end
                    end
                end
            end
        end
        % process, filter + detrend each array
        asleepXcVals = [];
        if isempty(asleepData) == false
            for ee = 1:length(asleepData)
                asleepProcData.data{ee,1} = filtfilt(sos,g,detrend(asleepData{ee,1},'constant'));
                asleepProcData.pupil{ee,1} = filtfilt(sos,g,detrend(asleepPupilData{ee,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            asleepLagTime = 30; % seconds
            asleepFrequency = samplingRate; % Hz
            asleepMaxLag = asleepLagTime*asleepFrequency;
            % run cross-correlation analysis - average through time
            for ff = 1:length(asleepProcData.data)
                asleepArray = asleepProcData.data{ff,1};
                asleepPupilarray = asleepProcData.pupil{ff,1};

                if length(asleepArray) ~=  length(asleepPupilarray)
                    minlength = min(length(asleepArray),length(asleepPupilarray));
                    asleepArray = asleepArray(1:minlength);
                    asleepPupilarray = asleepPupilarray(1:minlength);
                end
                [asleepXcVals(ff,:),asleepPupilLags] = xcorr(asleepArray,asleepPupilarray,asleepMaxLag,'coeff');
            end
            asleepMeanXcVals = mean(asleepXcVals,1);
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Asleep.(dataType).(pupilDataType).lags = asleepPupilLags;
            AnalysisResults.(animalID).CrossCorrelation.Asleep.(dataType).(pupilDataType).xcVals = asleepMeanXcVals;
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.Asleep.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.Asleep.(dataType).(pupilDataType).xcVals = [];
        end
        %% cross-correlation analysis for all data
        zz = 1;
        allData = []; allPupilData = []; allProcData = [];
        for cc = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(cc,:);
            [~,~,~] = GetFileInfo_JNeurosci2022(procDataFileID);
            load(procDataFileID,'-mat')
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
                            % pull data based on data type
                            if strcmp(dataType,'Ach_Rhodamine') == true
                                allData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(dataType,'NE_Rhodamine') == true
                                allData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(dataType,'Ach_GFP') == true
                                allData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(dataType,'NE_GFP') == true
                                allData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(dataType,'zDiameter') == true
                                allData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(dataType,'RH_thetaBandPower') == true
                                allData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(dataType,'RH_alphaBandPower') == true
                                allData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(dataType,'RH_gammaBandPower') == true
                                allData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end

                            % pull data based on pupil data type
                            if strcmp(pupilDataType,'Ach_Rhodamine') == true
                                allPupilData{zz,1} = ProcData.data.Rhodamine.Ach;
                            elseif strcmp(pupilDataType,'NE_Rhodamine') == true
                                allPupilData{zz,1} = ProcData.data.Rhodamine.NE;
                            elseif strcmp(pupilDataType,'Ach_GFP') == true
                                allPupilData{zz,1} = ProcData.data.GFP.Ach;
                            elseif strcmp(pupilDataType,'NE_GFP') == true
                                allPupilData{zz,1} = ProcData.data.GFP.NE;
                            elseif strcmp(pupilDataType,'zDiameter') == true
                                allPupilData{zz,1} = ProcData.data.Pupil.zDiameter;
                            elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
                                allPupilData{zz,1} = ProcData.data.cortical_RH.thetaBandPower;
                            elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
                                allPupilData{zz,1} = ProcData.data.cortical_RH.alphaBandPower;
                            elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
                                allPupilData{zz,1} = ProcData.data.cortical_RH.gammaBandPower;
                            end
                        zz = zz + 1;
                    end
                end
            end
        end
        % process, filter + detrend each array
        allXcVals = [];
        if isempty(allData) == false
            for dd = 1:length(allData)
                allProcData.data{dd,1} = filtfilt(sos,g,detrend(allData{dd,1},'constant'));
                allProcData.pupil{dd,1} = filtfilt(sos,g,detrend(allPupilData{dd,1},'constant'));
            end
            % set parameters for cross-correlation analysis
            allLagTime = 30; % seconds
            allFrequency = samplingRate; % Hz
            allMaxLag = allLagTime*allFrequency;
            % run cross-correlation analysis - average through time
            for ee = 1:length(allProcData.data)
                allHbTarray = allProcData.data{ee,1};
                allPupilarray = allProcData.pupil{ee,1};
                if length(allHbTarray) ~=  length(allPupilarray)
                    minlength = min(length(allHbTarray),length(allPupilarray));
                    allHbTarray = allHbTarray(1:minlength);
                    allPupilarray = allPupilarray(1:minlength);
                end
                [allXcVals(ee,:),allPupilLags] = xcorr(allHbTarray,allPupilarray,allMaxLag,'coeff');
            end
            allMeanXcVals = mean(allXcVals,1);
            % save results
            AnalysisResults.(animalID).CrossCorrelation.All.(dataType).(pupilDataType).lags = allPupilLags;
            AnalysisResults.(animalID).CrossCorrelation.All.(dataType).(pupilDataType).xcVals = allMeanXcVals;
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.All.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.All.(dataType).(pupilDataType).xcVals = [];
        end
        %% cross-correlation analysis for NREM data
        if isempty(SleepData.(modelType).NREM.data.Pupil) == false
            NREM_sleepTime = params.minTime.NREM; % seconds
            if strcmp(dataType,'Ach_Rhodamine') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'NE_Rhodamine') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'Ach_GFP') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'NE_GFP') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'zDiameter') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.zDiameter,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'RH_thetaBandPower') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.thetaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(dataType,'RH_alphaBandPower') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.alphaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);           
            elseif strcmp(dataType,'RH_gammaBandPower') == true
                [NREM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.gammaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            end

            if strcmp(pupilDataType,'Ach_Rhodamine') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'NE_Rhodamine') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Rhodamine.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'Ach_GFP') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.Ach,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'NE_GFP') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.GFP.NE,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'zDiameter') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.Pupil.zDiameter,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.thetaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.alphaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);           
            elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
                [NREM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).NREM.data.cortical_RH.gammaBandPower,SleepData.(modelType).NREM.FileIDs,SleepData.(modelType).NREM.BinTimes);
            end            
            
            NREM_finalVals = []; NREM_finalPupilVals = [];
            if isempty(NREM_finalData) == false
                % adjust events to match the edits made to the length of each spectrogram
                dd = 1;
                for cc = 1:length(NREM_finalData)
                    if sum(isnan(NREM_finalPupilData{cc,1})) == 0
                        NREM_vals = NREM_finalData{cc,1}(1:NREM_sleepTime*samplingRate);
                        NREM_pupilVals = NREM_finalPupilData{cc,1}(1:NREM_sleepTime*samplingRate);
                        NREM_finalVals{dd,1} = filtfilt(sos,g,detrend(NREM_vals,'constant'));
                        NREM_finalPupilVals{dd,1} = filtfilt(sos,g,detrend(NREM_pupilVals,'constant'));
                        dd = dd + 1;
                    end
                end
                % process, filter + detrend each array
                NREM_xcVals = [];
                if isempty(NREM_finalVals) == false
                    % run cross-correlation analysis - average through time
                    NREM_lagTime = 15; % seconds
                    NREM_frequency = samplingRate; % Hz
                    NREM_maxLag = NREM_lagTime*NREM_frequency;
                    for ee = 1:length(NREM_finalVals)
                        NREM_array = NREM_finalVals{ee,1};
                        NREM_pupilArray = NREM_finalPupilVals{ee,1};
                        if length(NREM_array) ~=  length(NREM_pupilArray)
                            minlength = min(length(NREM_array),length(NREM_pupilArray));
                            NREM_array = NREM_array(1:minlength);
                            NREM_pupilArray = NREM_pupilArray(1:minlength);
                        end
                        [NREM_xcVals(ee,:),NREM_PupilLags] = xcorr(NREM_array,NREM_pupilArray,NREM_maxLag,'coeff');
                    end
                    NREM_meanXcVals = mean(NREM_xcVals,1);
                    % save results
                    AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).lags = NREM_PupilLags;
                    AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).xcVals = NREM_meanXcVals;
                else
                    % save results
                    AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).lags = [];
                    AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).xcVals = [];
                end
            end
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.NREM.(dataType).(pupilDataType).xcVals = [];
        end
        %% cross-correlation analysis for REM
      if firstHrs == "false"   
        if isempty(SleepData.(modelType).REM.data.Pupil) == false
            REM_sleepTime = params.minTime.REM; % seconds
            if strcmp(dataType,'Ach_Rhodamine') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'NE_Rhodamine') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'Ach_GFP') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'NE_GFP') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'zDiameter') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.zDiameter,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'RH_thetaBandPower') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.thetaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(dataType,'RH_alphaBandPower') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.alphaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);           
            elseif strcmp(dataType,'RH_gammaBandPower') == true
                [REM_finalData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.gammaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            end

            if strcmp(pupilDataType,'Ach_Rhodamine') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'NE_Rhodamine') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Rhodamine.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'Ach_GFP') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.Ach,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'NE_GFP') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.GFP.NE,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'zDiameter') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.Pupil.zDiameter,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'RH_thetaBandPower') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.thetaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            elseif strcmp(pupilDataType,'RH_alphaBandPower') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.alphaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);           
            elseif strcmp(pupilDataType,'RH_gammaBandPower') == true
                [REM_finalPupilData,~,~] = RemoveStimSleepData_IOS(animalID,SleepData.(modelType).REM.data.cortical_RH.gammaBandPower,SleepData.(modelType).REM.FileIDs,SleepData.(modelType).REM.BinTimes);
            end  
            REM_finalVals = []; REM_finalPupilVals = [];
            if isempty(REM_finalData) == false
                % adjust events to match the edits made to the length of each spectrogram
                dd = 1;
                for cc = 1:length(REM_finalData)
                    if sum(isnan(REM_finalPupilData{cc,1})) == 0
                        REM_vals = REM_finalData{cc,1}(1:REM_sleepTime*samplingRate);
                        REM_pupilVals = REM_finalPupilData{cc,1}(1:REM_sleepTime*samplingRate);
                        REM_finalVals{dd,1} = filtfilt(sos,g,detrend(REM_vals,'constant'));
                        REM_finalPupilVals{dd,1} = filtfilt(sos,g,detrend(REM_pupilVals,'constant'));
                        dd = dd + 1;
                    end
                end
                % process, filter + detrend each array
                REM_xcVals = [];
                if isempty(REM_finalVals) == false
                    % run cross-correlation analysis - average through time
                    REM_lagTime = 15; % seconds
                    REM_frequency = samplingRate; % Hz
                    REM_maxLag = REM_lagTime*REM_frequency;
                    for ee = 1:length(REM_finalVals)
                        REM_array = REM_finalVals{ee,1};
                        REM_pupilArray = REM_finalPupilVals{ee,1};
                        if length(REM_array) ~=  length(REM_pupilArray)
                            minlength = min(length(REM_array),length(REM_pupilArray));
                            REM_array = REM_array(1:minlength);
                            REM_pupilArray = REM_pupilArray(1:minlength);
                        end
                        [REM_xcVals(ee,:),REM_PupilLags] = xcorr(REM_array,REM_pupilArray,REM_maxLag,'coeff');
                    end
                    REM_meanXcVals = mean(REM_xcVals,1);
                    % save results
                    AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).lags = REM_PupilLags;
                    AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).xcVals = REM_meanXcVals;
                else
                    % save results
                    AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).lags = [];
                    AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).xcVals = [];
                end
            end
        else
            % save results
            AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).lags = [];
            AnalysisResults.(animalID).CrossCorrelation.REM.(dataType).(pupilDataType).xcVals = [];
        end
      end
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
