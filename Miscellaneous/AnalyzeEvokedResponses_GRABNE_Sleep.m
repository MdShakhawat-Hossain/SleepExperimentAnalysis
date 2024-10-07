function [AnalysisResults] = AnalyzeEvokedResponses_GRABNE_Sleep(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
% extInd = strfind(group,delim);
% setName = group(extInd + 1:end);
% if strcmp(setName,'IOS Set A') == true
%     dataLocation = [rootFolder delim group delim animalID delim 'Bilateral Imaging'];
% elseif strcmp(setName,'IOS Set B') == true
%     dataLocation = [rootFolder delim group delim animalID delim 'Combined Imaging'];
% end
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end
% dataTypes = {'LH','RH'};
%% only run analysis for valid animal IDs
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID,'-mat')
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
% find and load AllSpecStruct.mat struct
allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
allSpecStructFile = {allSpecStructFileStruct.name}';
allSpecStructFileID = char(allSpecStructFile);
load(allSpecStructFileID,'-mat')
% criteria for whisking
WhiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaA.Comparison = {'gt','lt','gt'};
WhiskCriteriaA.Value = {0.5,2,5};
WhiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
WhiskCriteriaB.Comparison = {'gt','lt','gt'};
WhiskCriteriaB.Value = {2,5,5};
WhiskCriteriaC.Fieldname = {'duration','puffDistance'};
WhiskCriteriaC.Comparison = {'gt','gt'};
WhiskCriteriaC.Value = {5,5};
WhiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};

% criteria for movement
movementCriteriaA.Fieldname = {'duration','duration','puffDistance'};
movementCriteriaA.Comparison = {'gt','lt','gt'};
movementCriteriaA.Value = {2.5,5,5};
movementCriteriaB.Fieldname = {'duration','duration','puffDistance'};
movementCriteriaB.Comparison = {'gt','lt','gt'};
movementCriteriaB.Value = {5,10,5};
movementCriteriaC.Fieldname = {'duration','puffDistance'};
movementCriteriaC.Comparison = {'gt','gt'};
movementCriteriaC.Value = {10,5};
movementCriteriaNames = {'ShortMovement','IntermediateMovement','LongMovement'};

% criteria for stimulation
StimCriteriaA.Value = {'LPadSol'};
StimCriteriaA.Fieldname = {'solenoidName'};
StimCriteriaA.Comparison = {'equal'};
StimCriteriaB.Value = {'RPadSol'};
StimCriteriaB.Fieldname = {'solenoidName'};
StimCriteriaB.Comparison = {'equal'};
StimCriteriaC.Value = {'AudSol'};
StimCriteriaC.Fieldname = {'solenoidName'};
StimCriteriaC.Comparison = {'equal'};
stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
%% SleepData
sleepDataFileStruct = dir('*_SleepData.mat');
sleepDataFile = {sleepDataFileStruct.name}';
sleepDataFileID = char(sleepDataFile);
load(sleepDataFileID,'-mat')
if isfield (SleepData.Manual,'REM') == true
SleepFileIDs = [SleepData.Manual.NREM.FileIDs;SleepData.Manual.REM.FileIDs];
SleepBinTimes = [SleepData.Manual.NREM.BinTimes;SleepData.Manual.REM.BinTimes];
else
SleepFileIDs = [SleepData.Manual.NREM.FileIDs];
SleepBinTimes = [SleepData.Manual.NREM.BinTimes];    
end

%% analyze whisking-evoked responses
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.CBV.P_NE.whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.CBV.P_NE.whisk.trialDuration_sec;
    timeVector = (0:(EventData.CBV.P_NE.whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV.P_NE.whisk.epoch.offset;
    offset = EventData.CBV.P_NE.whisk.epoch.offset;

    if firstHrs == "false"
        for bb = 1:length(WhiskCriteriaNames)
            whiskCriteriaName = WhiskCriteriaNames{1,bb};
            if strcmp(whiskCriteriaName,'ShortWhisks') == true
                WhiskCriteria = WhiskCriteriaA;
            elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
                WhiskCriteria = WhiskCriteriaB;
            elseif strcmp(whiskCriteriaName,'LongWhisks') == true
                WhiskCriteria = WhiskCriteriaC;
            end
            % pull data from EventData.mat structure
            [whiskLogical] = FilterEvents_FP(EventData.CBV.P_NE.whisk,WhiskCriteria);
            combWhiskLogical = logical(whiskLogical);
            [AChWhiskCBVData] = EventData.CBV.P_ACh.whisk.NormData(combWhiskLogical,:);
            [NEWhiskCBVData] = EventData.CBV.P_NE.whisk.NormData(combWhiskLogical,:);
            [AChWhiskGFPData] = EventData.GFP.P_ACh.whisk.NormData(combWhiskLogical,:);
            [NEWhiskGFPData] = EventData.GFP.P_NE.whisk.NormData(combWhiskLogical,:);
    
            [allWhiskCorticalMUAData] = EventData.cortical_LH.corticalPower.whisk.NormData(combWhiskLogical,:);
    %         [allWhiskHippocampalMUAData] = EventData.hippocampus.corticalPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalGamData] = EventData.cortical_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
    %         [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskFileIDs] = EventData.CBV.P_NE.whisk.fileIDs(combWhiskLogical,:);
            [allWhiskEventTimes] = EventData.CBV.P_NE.whisk.eventTime(combWhiskLogical,:);
            allWhiskDurations = EventData.CBV.P_NE.whisk.duration(combWhiskLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [AChfinalWhiskCBVData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes] = RemoveInvalidData_IOS(AChWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [NEfinalWhiskCBVData_OLD,~,~,~] = RemoveInvalidData_IOS(NEWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);         
            [AChfinalWhiskGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(AChWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [NEfinalWhiskGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(NEWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalMUAData_OLD,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalGamData_OLD,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);

            % remove sleep data from stim
            [NEfinalWhiskCBVData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveSleepData(animalID,NEfinalWhiskCBVData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);
            [AChfinalWhiskCBVData,~,~,~] = RemoveSleepData(animalID,AChfinalWhiskCBVData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);
            [NEfinalWhiskGFPData,~,~,~] = RemoveSleepData(animalID,NEfinalWhiskGFPData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);        
            [AChfinalWhiskGFPData,~,~,~] = RemoveSleepData(animalID,AChfinalWhiskGFPData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalWhiskCorticalMUAData,~,~,~] = RemoveSleepData(animalID,finalWhiskCorticalMUAData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalWhiskCorticalGamData,~,~,~] = RemoveSleepData(animalID,finalWhiskCorticalGamData_OLD,finalWhiskFileIDs,finalWhiskDurations,finalWhiskFileEventTimes,SleepFileIDs,SleepBinTimes);

            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procWhiskGFPData procWhiskCBVData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
            dd = 1;
            for cc = 1:size(NEfinalWhiskCBVData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 0;%2;
                whiskEndTime = whiskStartTime + 15;
                finalWhiskFileID = finalWhiskFileIDs{cc,1};
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    NEwhiskCBVarray = NEfinalWhiskCBVData(cc,:);
                    AChwhiskCBVarray = AChfinalWhiskCBVData(cc,:);
    
                    NEwhiskGFParray = NEfinalWhiskGFPData(cc,:);
                    AChwhiskGFParray = AChfinalWhiskGFPData(cc,:);
    
                    whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
    %                 whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                    whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
    %                 whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                    
                    AChfiltWhiskCBVarray = sgolayfilt(AChwhiskCBVarray,3,9) - mean(AChwhiskCBVarray(1:(offset*samplingRate)));
                    NEfiltWhiskCBVarray = sgolayfilt(NEwhiskCBVarray,3,9) - mean(NEwhiskCBVarray(1:(offset*samplingRate)));
                    AChfiltWhiskGFParray = sgolayfilt(AChwhiskGFParray,3,9) - mean(AChwhiskGFParray(1:(offset*samplingRate)));
                    NEfiltWhiskGFParray = sgolayfilt(NEwhiskGFParray,3,9) - mean(NEwhiskGFParray(1:(offset*samplingRate)));
    
                    AChprocWhiskCBVData(dd,:) = AChfiltWhiskCBVarray;
                    NEprocWhiskCBVData(dd,:) = NEfiltWhiskCBVarray;
                    AChprocWhiskGFPData(dd,:) = AChfiltWhiskGFParray;
                    NEprocWhiskGFPData(dd,:) = NEfiltWhiskGFParray;
    
                    procWhiskCorticalMUAData(dd,:) = whiskCorticalMUAarray - mean(whiskCorticalMUAarray(1:(offset*samplingRate)));
                    procWhiskCorticalGamData(dd,:) = whiskCorticalGamArray - mean(whiskCorticalGamArray(1:(offset*samplingRate)));
                    finalWhiskStartTimes(dd,1) = whiskStartTime;
                    finalWhiskEndTimes(dd,1) = whiskEndTime;
                    finalWhiskFiles{dd,1} = finalWhiskFileID;
                    dd = dd + 1;
                end
            end
            AChmeanWhiskCBVData = mean(AChprocWhiskCBVData,1);%*100;
            AChstdWhiskCBVData = std(AChprocWhiskCBVData,0,1);%*100;

            NEmeanWhiskCBVData = mean(NEprocWhiskCBVData,1);%*100;
            NEstdWhiskCBVData = std(NEprocWhiskCBVData,0,1);%*100;
    
            AChmeanWhiskGFPData = mean(AChprocWhiskGFPData,1);%*100;
            AChstdWhiskGFPData = std(AChprocWhiskGFPData,0,1);%*100;

            NEmeanWhiskGFPData = mean(NEprocWhiskGFPData,1);%*100;
            NEstdWhiskGFPData = std(NEprocWhiskGFPData,0,1);%*100;
    
            meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
            stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
            meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
            stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
            % extract LFP from spectrograms associated with the whisking indecies
            whiskCorticalZhold = [];
            for ee = 1:length(finalWhiskFiles)
                % load normalized one-second bin data from each file
                whiskFileID = finalWhiskFiles{ee,1};
                whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
                whiskSpecField = 'cortical_LH';
                for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                    if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                        whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                        F = AllSpecData.(whiskSpecField).F{ff,1};
                        T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                    end
                end
                whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
                whiskStartTimeIndex = whiskStartTimeIndex(1);
                whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
                whiskDurationIndex = whiskDurationIndex(end);
                whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                whiskCorticalZhold = cat(3,whiskCorticalZhold,whiskCorticalS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanWhiskCortS = mean(whiskCorticalZhold,3);
            baseWhiskCortS_Vals = mean(meanWhiskCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixWhiskCortS_Vals = baseWhiskCortS_Vals.*ones(size(meanWhiskCortS));
            msStimWhiskS_Vals = (meanWhiskCortS - baseMatrixWhiskCortS_Vals);

            T2 = -20:(1/specSamplingRate):20;
            % save results
            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).CBV.CBV = AChmeanWhiskCBVData;
            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).CBV.CBVStD = AChstdWhiskCBVData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).CBV.CBV = NEmeanWhiskCBVData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).CBV.CBVStD = NEstdWhiskCBVData;
    
            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).GFP.GFP = AChmeanWhiskGFPData;
            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).GFP.GFPStD = AChstdWhiskGFPData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).GFP.GFP = NEmeanWhiskGFPData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).GFP.GFPStD = NEstdWhiskGFPData;

            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).CBV.CBVRaw = AChprocWhiskCBVData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).CBV.CBVRaw = NEprocWhiskCBVData;
    
            AnalysisResults.(animalID).Whisk.P_ACh.(whiskCriteriaName).GFP.GFPRaw = AChprocWhiskGFPData;
            AnalysisResults.(animalID).Whisk.P_NE.(whiskCriteriaName).GFP.GFPRaw = NEprocWhiskGFPData;

    
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.corticalS = msStimWhiskS_Vals;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.T = T2;
            AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.F = F;
        end
    end
    %% analyze movement evoked responses
    %{
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.CBV.P_NE.movement.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.CBV.P_NE.movement.trialDuration_sec;
    timeVector = (0:(EventData.CBV.P_NE.movement.epoch.duration*samplingRate))/samplingRate - EventData.CBV.P_NE.movement.epoch.offset;
    offset = EventData.CBV.P_NE.movement.epoch.offset;

    if firstHrs == "false"
        for bb = 1:length(movementCriteriaNames)
            movementCriteriaName = movementCriteriaNames{1,bb};
            if strcmp(movementCriteriaName,'ShortMovement') == true
                MovementCriteria = movementCriteriaA;
            % elseif strcmp(movementCriteriaName,'IntermediateMovement') == true
            %     MovementCriteria = movementCriteriaB;
            % elseif strcmp(movementCriteriaName,'LongMovement') == true
            %     MovementCriteria = movementCriteriaC;
            end
            % pull data from EventData.mat structure
            [movementLogical] = FilterEvents_FP(EventData.CBV.P_NE.movement,MovementCriteria);
            combMovementLogical = logical(movementLogical);
            [AChMovementCBVData] = EventData.CBV.P_ACh.movement.NormData(combMovementLogical,:);
            [NEMovementCBVData] = EventData.CBV.P_NE.movement.NormData(combMovementLogical,:);
            [AChMovementGFPData] = EventData.GFP.P_ACh.movement.NormData(combMovementLogical,:);
            [NEMovementGFPData] = EventData.GFP.P_NE.movement.NormData(combMovementLogical,:);
    
            [allMovementCorticalMUAData] = EventData.cortical_LH.corticalPower.movement.NormData(combMovementLogical,:);
            [allMovementCorticalGamData] = EventData.cortical_LH.gammaBandPower.movement.NormData(combMovementLogical,:);
            [allMovementFileIDs] = EventData.CBV.P_NE.movement.fileIDs(combMovementLogical,:);
            [allMovementEventTimes] = EventData.CBV.P_NE.movement.eventTime(combMovementLogical,:);
            allMovementDurations = EventData.CBV.P_NE.movement.duration(combMovementLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [AChfinalMovementCBVData_OLD,finalMovementFileIDs,~,finalMovementFileEventTimes] = RemoveInvalidData_IOS(AChMovementCBVData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [NEfinalMovementCBVData_OLD,~,~,~] = RemoveInvalidData_IOS(NEMovementCBVData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
                   
            [AChfinalMovementGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(AChMovementGFPData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [NEfinalMovementGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(NEMovementGFPData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
    
            [finalMovementCorticalMUAData_OLD,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalMUAData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [finalMovementCorticalGamData_OLD,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalGamData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            
            % remove sleep data from stim
            [NEfinalMovementCBVData,finalMovementFileIDs,~,finalMovementFileEventTimes] = RemoveSleepData(animalID,NEfinalMovementCBVData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);
            [AChfinalMovementCBVData,~,~,~] = RemoveSleepData(animalID,AChfinalMovementCBVData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);
            [NEfinalMovementGFPData,~,~,~] = RemoveSleepData(animalID,NEfinalMovementGFPData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);        
            [AChfinalMovementGFPData,~,~,~] = RemoveSleepData(animalID,AChfinalMovementGFPData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalMovementCorticalMUAData,~,~,~] = RemoveSleepData(animalID,finalMovementCorticalMUAData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalMovementCorticalGamData,~,~,~] = RemoveSleepData(animalID,finalMovementCorticalGamData_OLD,finalMovementFileIDs,finalMovementDurations,finalMovementFileEventTimes,SleepFileIDs,SleepBinTimes);

            % lowpass filter each movementing event and mean-subtract by the first 2 seconds
            clear procMovementGFPData procMovementCBVData procMovementCorticalMUAData procMovementHippocampalMUAData procMovementCorticalGamData procMovementHippocampalGamData finalMovementStartTimes finalMovementEndTimes finalMovementFiles
            dd = 1;
            for cc = 1:size(NEfinalMovementCBVData,1)
                movementStartTime = round(finalMovementFileEventTimes(cc,1),1) - 0;%2;
                movementEndTime = movementStartTime + 15;
                finalMovementFileID = finalMovementFileIDs{cc,1};
                if movementStartTime >= 0.5 && movementEndTime <= (trialDuration_sec - 0.5)
                    NEmovementCBVarray = NEfinalMovementCBVData(cc,:);
                    AChmovementCBVarray = AChfinalMovementCBVData(cc,:);
    
                    NEmovementGFParray = NEfinalMovementGFPData(cc,:);
                    AChmovementGFParray = AChfinalMovementGFPData(cc,:);
    
                    movementCorticalMUAarray = finalMovementCorticalMUAData(cc,:);
                    movementCorticalGamArray = finalMovementCorticalGamData(cc,:);
                    
                    AChfiltMovementCBVarray = sgolayfilt(AChmovementCBVarray,3,9) - mean(AChmovementCBVarray(1:(offset*samplingRate)));
                    NEfiltMovementCBVarray = sgolayfilt(NEmovementCBVarray,3,9) - mean(NEmovementCBVarray(1:(offset*samplingRate)));
                    AChfiltMovementGFParray = sgolayfilt(AChmovementGFParray,3,9) - mean(AChmovementGFParray(1:(offset*samplingRate)));
                    NEfiltMovementGFParray = sgolayfilt(NEmovementGFParray,3,9) - mean(NEmovementGFParray(1:(offset*samplingRate)));
    
                    AChprocMovementCBVData(dd,:) = AChfiltMovementCBVarray;
                    NEprocMovementCBVData(dd,:) = NEfiltMovementCBVarray;
                    AChprocMovementGFPData(dd,:) = AChfiltMovementGFParray;
                    NEprocMovementGFPData(dd,:) = NEfiltMovementGFParray;
    
                    procMovementCorticalMUAData(dd,:) = movementCorticalMUAarray - mean(movementCorticalMUAarray(1:(offset*samplingRate)));
                    procMovementCorticalGamData(dd,:) = movementCorticalGamArray - mean(movementCorticalGamArray(1:(offset*samplingRate)));
                    finalMovementStartTimes(dd,1) = movementStartTime;
                    finalMovementEndTimes(dd,1) = movementEndTime;
                    finalMovementFiles{dd,1} = finalMovementFileID;
                    dd = dd + 1;
                end
            end
            AChmeanMovementCBVData = mean(AChprocMovementCBVData,1);%*100;
            AChstdMovementCBVData = std(AChprocMovementCBVData,0,1);%*100;
            NEmeanMovementCBVData = mean(NEprocMovementCBVData,1);%*100;
            NEstdMovementCBVData = std(NEprocMovementCBVData,0,1);%*100;
    
            AChmeanMovementGFPData = mean(AChprocMovementGFPData,1);%*100;
            AChstdMovementGFPData = std(AChprocMovementGFPData,0,1);%*100;
            NEmeanMovementGFPData = mean(NEprocMovementGFPData,1);%*100;
            NEstdMovementGFPData = std(NEprocMovementGFPData,0,1);%*100;
    
            meanMovementCorticalMUAData = mean(procMovementCorticalMUAData,1)*100;
            stdMovementCorticalMUAData = std(procMovementCorticalMUAData,0,1)*100;
            meanMovementCorticalGamData = mean(procMovementCorticalGamData,1)*100;
            stdMovementCorticalGamData = std(procMovementCorticalGamData,0,1)*100;
            % extract LFP from spectrograms associated with the movementing indecies
            movementCorticalZhold = [];
            for ee = 1:length(finalMovementFiles)
                % load normalized one-second bin data from each file
                movementFileID = finalMovementFiles{ee,1};
                movementSpecDataFileID = [animalID '_' movementFileID '_SpecDataB.mat'];
                movementSpecField = 'cortical_LH';
                for ff = 1:length(AllSpecData.(movementSpecField).fileIDs)
                    if strcmp(AllSpecData.(movementSpecField).fileIDs{ff,1},movementSpecDataFileID) == true
                        movementCorticalS_Data = AllSpecData.(movementSpecField).normS{ff,1};
                        F = AllSpecData.(movementSpecField).F{ff,1};
                        T = round(AllSpecData.(movementSpecField).T{ff,1},1);
                    end
                end
                movementStartTimeIndex = find(T == round(finalMovementStartTimes(ee,1),1));
                movementStartTimeIndex = movementStartTimeIndex(1);
                movementDurationIndex = find(T == round(finalMovementEndTimes(ee,1),1));
                movementDurationIndex = movementDurationIndex(end);
                movementCorticalS_Vals = movementCorticalS_Data(:,movementStartTimeIndex:movementDurationIndex);
                movementCorticalZhold = cat(3,movementCorticalZhold,movementCorticalS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanMovementCortS = mean(movementCorticalZhold,3);
            baseMovementCortS_Vals = mean(meanMovementCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixMovementCortS_Vals = baseMovementCortS_Vals.*ones(size(meanMovementCortS));
            msStimMovementS_Vals = (meanMovementCortS - baseMatrixMovementCortS_Vals);
            T2 = -20:(1/specSamplingRate):20;
            % save results
            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).CBV.CBV = AChmeanMovementCBVData;
            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).CBV.CBVStD = AChstdMovementCBVData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).CBV.CBV = NEmeanMovementCBVData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).CBV.CBVStD = NEstdMovementCBVData;
    
            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).GFP.GFP = AChmeanMovementGFPData;
            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).GFP.GFPStD = AChstdMovementGFPData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).GFP.GFP = NEmeanMovementGFPData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).GFP.GFPStD = NEstdMovementGFPData;


            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).CBV.CBVRaw = AChprocMovementCBVData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).CBV.CBVRaw = NEprocMovementCBVData;
    
            AnalysisResults.(animalID).Movement.P_ACh.(movementCriteriaName).GFP.GFPRaw = AChprocMovementGFPData;
            AnalysisResults.(animalID).Movement.P_NE.(movementCriteriaName).GFP.GFPRaw = NEprocMovementGFPData;
    
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).MUA.corticalData = meanMovementCorticalMUAData;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).MUA.corticalStD = stdMovementCorticalMUAData;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).Gam.corticalData = meanMovementCorticalGamData;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).Gam.corticalStD = stdMovementCorticalGamData;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).LFP.corticalS = msStimMovementS_Vals;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).LFP.T = T2;
            AnalysisResults.(animalID).Movement.cortical.(movementCriteriaName).LFP.F = F;
        end
    end
    %}
    %% analyze stimulus-evoked responses
    if firstHrs == "true"
        for gg = 1:length(stimCriteriaNames)
            stimCriteriaName = stimCriteriaNames{1,gg};
            if strcmp(stimCriteriaName,'stimCriteriaA') == true
                StimCriteria = StimCriteriaA;
                solenoid = 'LPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
                StimCriteria = StimCriteriaB;
                solenoid = 'RPadSol';
            elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
                StimCriteria = StimCriteriaC;
                solenoid = 'AudSol';
            end
            % pull data from EventData.mat structure
            allStimFilter = FilterEvents_IOS(EventData.CBV.P_NE.stim,StimCriteria);
            [AChallStimCBVData] = EventData.CBV.P_ACh.stim.NormData(allStimFilter,:);
            [NEallStimCBVData] = EventData.CBV.P_NE.stim.NormData(allStimFilter,:);
            [AChallStimGFPData] = EventData.GFP.P_ACh.stim.NormData(allStimFilter,:);
            [NEallStimGFPData] = EventData.GFP.P_NE.stim.NormData(allStimFilter,:); 
    
            [allStimCortMUAData] = EventData.cortical_LH.corticalPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData] = EventData.cortical_LH.gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.CBV.P_NE.stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.CBV.P_NE.stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalStimCBVData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes] = RemoveInvalidData_IOS(NEallStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [AChfinalStimCBVData_OLD,~,~,~] = RemoveInvalidData_IOS(AChallStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [NEfinalStimGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(NEallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);        
            [AChfinalStimGFPData_OLD,~,~,~] = RemoveInvalidData_IOS(AChallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortMUAData_OLD,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData_OLD,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
             
            % remove sleep data from stim
            [NEfinalStimCBVData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveSleepData(animalID,NEfinalStimCBVData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);
            [AChfinalStimCBVData,~,~,~] = RemoveSleepData(animalID,AChfinalStimCBVData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);
            [NEfinalStimGFPData,~,~,~] = RemoveSleepData(animalID,NEfinalStimGFPData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);        
            [AChfinalStimGFPData,~,~,~] = RemoveSleepData(animalID,AChfinalStimGFPData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalStimCortMUAData,~,~,~] = RemoveSleepData(animalID,finalStimCortMUAData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);
            [finalStimCortGamData,~,~,~] = RemoveSleepData(animalID,finalStimCortGamData_OLD,finalStimFileIDs,finalStimDurations,finalStimFileEventTimes,SleepFileIDs,SleepBinTimes);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procStimGFPData procStimCBVData procStimCBVData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
            ii = 1;
            for hh = 1:size(NEfinalStimCBVData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 0;
                stimEndTime = stimStartTime + 15;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    AChstimCBVarray = AChfinalStimCBVData(hh,:);
                    NEstimCBVarray = NEfinalStimCBVData(hh,:);
                    AChstimGFParray = AChfinalStimGFPData(hh,:);
                    NEstimGFParray = NEfinalStimGFPData(hh,:);
    
                    stimCortMUAarray = finalStimCortMUAData(hh,:);
                    stimCortGamArray = finalStimCortGamData(hh,:);
    
                    AChfiltStimCBVarray = sgolayfilt(AChstimCBVarray,3,9);
                    NEfiltStimCBVarray = sgolayfilt(NEstimCBVarray,3,9);
                    AChfiltStimGFParray = sgolayfilt(AChstimGFParray,3,9);
                    NEfiltStimGFParray = sgolayfilt(NEstimGFParray,3,9);
    
                    AChprocStimCBVData(hh,:) = AChfiltStimCBVarray - mean(AChfiltStimCBVarray(1:(offset*samplingRate)));
                    NEprocStimCBVData(hh,:) = NEfiltStimCBVarray - mean(NEfiltStimCBVarray(1:(offset*samplingRate)));
                    AChprocStimGFPData(hh,:) = AChfiltStimGFParray - mean(AChfiltStimGFParray(1:(offset*samplingRate)));
                    NEprocStimGFPData(hh,:) = NEfiltStimGFParray - mean(NEfiltStimGFParray(1:(offset*samplingRate)));
    
                    procStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
                    procStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
                    finalStimStartTimes(ii,1) = stimStartTime;
                    finalStimEndTimes(ii,1) = stimEndTime;
                    finalStimFiles{ii,1} = finalStimFileID;
                    ii = ii + 1;
                end
            end
            AChmeanStimCBVData = mean(AChprocStimCBVData,1);%*100;
            AChstdStimCBVData = std(AChprocStimCBVData,0,1);%*100;
            NEmeanStimCBVData = mean(NEprocStimCBVData,1);%*100;
            NEstdStimCBVData = std(NEprocStimCBVData,0,1);%*100;
            AChmeanStimGFPData = mean(AChprocStimGFPData,1);%*100;
            AChstdStimGFPData = std(AChprocStimGFPData,0,1);%*100;
            NEmeanStimGFPData = mean(NEprocStimGFPData,1);%*100;
            NEstdStimGFPData = std(NEprocStimGFPData,0,1);%*100;
    
            meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
            stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
            meanStimCortGamData = mean(procStimCortGamData,1)*100;
            stdStimCortGamData = std(procStimCortGamData,0,1)*100;
            % extract LFP from spectrograms associated with the stimuli indecies
            stimCortZhold = [];
%             stimHipZhold = [];
            for jj = 1:length(finalStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_LH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                        F = AllSpecData.(stimSpecField).F{kk,1};
                        T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanStimCortS = mean(stimCortZhold,3);
            baseStimCortS_Vals = mean(meanStimCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixStimCortS_Vals = baseStimCortS_Vals.*ones(size(meanStimCortS));
            msStimCortS_Vals = (meanStimCortS - baseMatrixStimCortS_Vals);
            % save results
            T2 = -20:(1/specSamplingRate):20;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).count = size(procStimCortMUAData,1);
    
            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).CBV.CBV = AChmeanStimCBVData;
            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).CBV.CBVStD = AChstdStimCBVData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).CBV.CBV = NEmeanStimCBVData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).CBV.CBVStD = NEstdStimCBVData;
            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).GFP.GFP= AChmeanStimGFPData;
            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).GFP.GFPStD = AChstdStimGFPData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).GFP.GFP= NEmeanStimGFPData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).GFP.GFPStD = NEstdStimGFPData;

            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).CBV.CBVRaw = AChprocStimCBVData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).CBV.CBVRaw = NEprocStimCBVData;
            AnalysisResults.(animalID).Stim.P_ACh.(solenoid).GFP.GFPRaw= AChprocStimGFPData;
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).GFP.GFPRaw= NEprocStimGFPData;

    
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalStD = stdStimCortMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalData = meanStimCortGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalStD = stdStimCortGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.corticalS = msStimCortS_Vals;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.F = F;
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

    %% 

%     %% [1-S2o] CBV auditory stim
%     figure;
% ax15 = subplot(1,2,1);
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.CBV.CBV,'-','color',colors('indian red'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.CBV.CBV   + AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.CBV.CBVStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.CBV.CBV   - AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.CBV.CBVStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% title('Auditory stim Blood Volume')
% ylabel('Z \DeltaF/F (ACh)')
% ax15.YLim = [-3 6];
% 
% yyaxis right
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.CBV.CBV  ,'-','color',colors('army green'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.CBV.CBV   + AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.CBV.CBVStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.CBV.CBV   - AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.CBV.CBVStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% ylabel('Z \DeltaF/F (NE)')
% ax15.YAxis(1).Color = colors('indian red');
% ax15.YAxis(2).Color = colors('army green');
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax15.TickLength = [0.03,0.03];
% ax15.YLim = [-3 6];
% 
% %% [1-S2r] GFP auditory stim
% ax18 = subplot(1,2,2);
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.GFP.GFP  ,'-','color',colors('indian red'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.GFP.GFP   + AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.GFP.GFPStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.GFP.GFP   - AnalysisResults_firstHrs.NEACh002.Stim.P_ACh.AudSol.GFP.GFPStD  ,'-','color',colors('indian red'),'LineWidth',0.10)
% title('Auditory stim GFP')
% ylabel('Z \DeltaF/F GRAB ACh')
% ax18.YLim = [-3 6];
% 
% yyaxis right
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.GFP.GFP  ,'-','color',colors('army green'),'LineWidth',2)
% hold on
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.GFP.GFP   + AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.GFP.GFPStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% plot(AnalysisResults_firstHrs.NEACh002.Stim.cortical.AudSol.timeVector  ,AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.GFP.GFP   - AnalysisResults_firstHrs.NEACh002.Stim.P_NE.AudSol.GFP.GFPStD  ,'-','color',colors('army green'),'LineWidth',0.10)
% title('Auditory stim GFP')
% ylabel('Z \DeltaF/F GRAB NE')
% ax18.YAxis(1).Color = colors('indian red');
% ax18.YAxis(2).Color = colors('army green');
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax18.TickLength = [0.03,0.03];
% ax18.YLim = [-3 6];

end

