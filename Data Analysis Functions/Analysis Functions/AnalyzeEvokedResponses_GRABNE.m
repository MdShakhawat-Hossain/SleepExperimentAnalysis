function [AnalysisResults] = AnalyzeEvokedResponses_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
%________________________________________________________________________________________________________________________

%% function parameters
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
    end
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
            % elseif strcmp(whiskCriteriaName,'LongWhisks') == true
            %     WhiskCriteria = WhiskCriteriaC;
            end
            % pull data from EventData.mat structure
            [whiskLogical] = FilterEvents_FP(EventData.CBV.P_NE.whisk,WhiskCriteria);
            combWhiskLogical = logical(whiskLogical);
            [AChWhiskCBVData] = EventData.CBV.P_ACh.whisk.NormData(combWhiskLogical,:);
            [NEWhiskCBVData] = EventData.CBV.P_NE.whisk.NormData(combWhiskLogical,:);
            [AChWhiskGFPData] = EventData.GFP.P_ACh.whisk.NormData(combWhiskLogical,:);
            [NEWhiskGFPData] = EventData.GFP.P_NE.whisk.NormData(combWhiskLogical,:);

            [allWhiskPupilData] = EventData.Pupil.Diameter.whisk.NormData(combWhiskLogical,:);
    
            [allWhiskCorticalMUAData_LH] = EventData.cortical_LH.corticalPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalGamData_LH] = EventData.cortical_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
            % [allWhiskCorticalMUAData_RH] = EventData.cortical_RH.corticalPower.whisk.NormData(combWhiskLogical,:);
            % [allWhiskCorticalGamData_RH] = EventData.cortical_RH.gammaBandPower.whisk.NormData(combWhiskLogical,:);

            [allWhiskFileIDs] = EventData.CBV.P_NE.whisk.fileIDs(combWhiskLogical,:);
            [allWhiskEventTimes] = EventData.CBV.P_NE.whisk.eventTime(combWhiskLogical,:);
            allWhiskDurations = EventData.CBV.P_NE.whisk.duration(combWhiskLogical,:);
            % keep only the data that occurs within the manually-approved awake regions
            [AChfinalWhiskCBVData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS(AChWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [NEfinalWhiskCBVData,~,~,~] = RemoveInvalidData_IOS(NEWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
                   
            [AChfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(AChWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [NEfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(NEWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            
            [finalWhiskPupilData,~,~,~] = RemoveInvalidData_IOS(allWhiskPupilData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
  
            [finalWhiskCorticalMUAData_LH,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData_LH,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalGamData_LH,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData_LH,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % [finalWhiskCorticalMUAData_RH,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData_LH,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % [finalWhiskCorticalGamData_RH,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData_LH,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles procWhiskPupilData NEprocWhiskGFPData AChprocWhiskGFPData NEprocWhiskCBVData AChprocWhiskCBVData procWhiskCorticalMUAData_LH procWhiskCorticalGamData_LH procWhiskCorticalMUAData_RH procWhiskCorticalGamData_RH
            dd = 1;
            for cc = 1:size(NEfinalWhiskCBVData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 3;
                whiskEndTime = whiskStartTime + 15;
                finalWhiskFileID = finalWhiskFileIDs{cc,1};
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    NEwhiskCBVarray = NEfinalWhiskCBVData(cc,:);
                    AChwhiskCBVarray = AChfinalWhiskCBVData(cc,:);
    
                    NEwhiskGFParray = NEfinalWhiskGFPData(cc,:);
                    AChwhiskGFParray = AChfinalWhiskGFPData(cc,:);

                    whiskPupilarray = finalWhiskPupilData(cc,:);
    
                    whiskCorticalMUAarray_LH = finalWhiskCorticalMUAData_LH(cc,:);
                    whiskCorticalGamArray_LH = finalWhiskCorticalGamData_LH(cc,:);

                    % whiskCorticalMUAarray_RH = finalWhiskCorticalMUAData_RH(cc,:);
                    % whiskCorticalGamArray_RH = finalWhiskCorticalGamData_RH(cc,:);
                    
                    AChfiltWhiskCBVarray = sgolayfilt(AChwhiskCBVarray,3,9) - mean(AChwhiskCBVarray(1:(offset*samplingRate)));
                    NEfiltWhiskCBVarray = sgolayfilt(NEwhiskCBVarray,3,9) - mean(NEwhiskCBVarray(1:(offset*samplingRate)));
                    AChfiltWhiskGFParray = sgolayfilt(AChwhiskGFParray,3,9) - mean(AChwhiskGFParray(1:(offset*samplingRate)));
                    NEfiltWhiskGFParray = sgolayfilt(NEwhiskGFParray,3,9) - mean(NEwhiskGFParray(1:(offset*samplingRate)));             
                    filtWhiskPupilarray = sgolayfilt(whiskPupilarray,3,9) - mean(whiskPupilarray(1:(offset*samplingRate)));

                    AChprocWhiskCBVData(dd,:) = AChfiltWhiskCBVarray;
                    NEprocWhiskCBVData(dd,:) = NEfiltWhiskCBVarray;
                    AChprocWhiskGFPData(dd,:) = AChfiltWhiskGFParray;
                    NEprocWhiskGFPData(dd,:) = NEfiltWhiskGFParray;
                    procWhiskPupilData(dd,:) = filtWhiskPupilarray;

                    procWhiskCorticalMUAData_LH(dd,:) = whiskCorticalMUAarray_LH - mean(whiskCorticalMUAarray_LH(1:(offset*samplingRate)));
                    procWhiskCorticalGamData_LH(dd,:) = whiskCorticalGamArray_LH - mean(whiskCorticalGamArray_LH(1:(offset*samplingRate)));
                    % procWhiskCorticalMUAData_RH(dd,:) = whiskCorticalMUAarray_RH - mean(whiskCorticalMUAarray_RH(1:(offset*samplingRate)));
                    % procWhiskCorticalGamData_RH(dd,:) = whiskCorticalGamArray_RH - mean(whiskCorticalGamArray_RH(1:(offset*samplingRate)));

                    finalWhiskStartTimes(dd,1) = whiskStartTime;
                    finalWhiskEndTimes(dd,1) = whiskEndTime;
                    finalWhiskFiles{dd,1} = finalWhiskFileID;
                    dd = dd + 1;
                end
            end
            AChmeanWhiskCBVData = mean(AChprocWhiskCBVData,1);
            AChstdWhiskCBVData = std(AChprocWhiskCBVData,0,1);

            NEmeanWhiskCBVData = mean(NEprocWhiskCBVData,1);
            NEstdWhiskCBVData = std(NEprocWhiskCBVData,0,1);
    
            AChmeanWhiskGFPData = mean(AChprocWhiskGFPData,1);
            AChstdWhiskGFPData = std(AChprocWhiskGFPData,0,1);

            NEmeanWhiskGFPData = mean(NEprocWhiskGFPData,1);
            NEstdWhiskGFPData = std(NEprocWhiskGFPData,0,1);

            meanWhiskPupilData = mean(procWhiskPupilData,1);
            stdWhiskPupilData = std(procWhiskPupilData,0,1);
    
            meanWhiskCorticalMUAData_LH = mean(procWhiskCorticalMUAData_LH,1)*100;
            stdWhiskCorticalMUAData_LH = std(procWhiskCorticalMUAData_LH,0,1)*100;
            meanWhiskCorticalGamData_LH = mean(procWhiskCorticalGamData_LH,1)*100;
            stdWhiskCorticalGamData_LH = std(procWhiskCorticalGamData_LH,0,1)*100;

            % meanWhiskCorticalMUAData_RH = mean(procWhiskCorticalMUAData_RH,1)*100;
            % stdWhiskCorticalMUAData_RH = std(procWhiskCorticalMUAData_RH,0,1)*100;
            % meanWhiskCorticalGamData_RH = mean(procWhiskCorticalGamData_RH,1)*100;
            % stdWhiskCorticalGamData_RH = std(procWhiskCorticalGamData_RH,0,1)*100;            
            %% extract ECoG spectrograms associated with the whisking indecies
            % left cortex
            whiskCorticalZhold_LH = [];
            for ee = 1:length(finalWhiskFiles)
                % load normalized one-second bin data from each file
                whiskFileID = finalWhiskFiles{ee,1};
                whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
                whiskSpecField = 'cortical_LH';
                for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                    if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                        whiskCorticalS_Data_LH = AllSpecData.(whiskSpecField).normS{ff,1};
                        F = AllSpecData.(whiskSpecField).F{ff,1};
                        T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                    end
                end
                whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
                whiskStartTimeIndex = whiskStartTimeIndex(1);
                whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
                whiskDurationIndex = whiskDurationIndex(end);
                whiskCorticalS_Vals_LH = whiskCorticalS_Data_LH(:,whiskStartTimeIndex:whiskDurationIndex);
                whiskCorticalZhold_LH = cat(3,whiskCorticalZhold_LH,whiskCorticalS_Vals_LH);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanWhiskCortS_LH = mean(whiskCorticalZhold_LH,3);
            baseWhiskCortS_Vals_LH = mean(meanWhiskCortS_LH(:,1:1.5*specSamplingRate),2);
            baseMatrixWhiskCortS_Vals_LH = baseWhiskCortS_Vals_LH.*ones(size(meanWhiskCortS_LH));
            msStimWhiskS_Vals_LH = (meanWhiskCortS_LH - baseMatrixWhiskCortS_Vals_LH);

            % right cortex
            % whiskCorticalZhold_RH = [];
            % for ee = 1:length(finalWhiskFiles)
            %     % load normalized one-second bin data from each file
            %     whiskFileID = finalWhiskFiles{ee,1};
            %     whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            %     whiskSpecField = 'cortical_RH';
            %     for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
            %         if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
            %             whiskCorticalS_Data_RH = AllSpecData.(whiskSpecField).normS{ff,1};
            %             F = AllSpecData.(whiskSpecField).F{ff,1};
            %             T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
            %         end
            %     end
            %     whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            %     whiskStartTimeIndex = whiskStartTimeIndex(1);
            %     whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            %     whiskDurationIndex = whiskDurationIndex(end);
            %     whiskCorticalS_Vals_RH = whiskCorticalS_Data_RH(:,whiskStartTimeIndex:whiskDurationIndex);
            %     whiskCorticalZhold_RH = cat(3,whiskCorticalZhold_RH,whiskCorticalS_Vals_RH);
            % end
            % % cortical mean-subtract by first 2 seconds prior to stimulus
            % meanWhiskCortS_RH = mean(whiskCorticalZhold_RH,3);
            % baseWhiskCortS_Vals_RH = mean(meanWhiskCortS_RH(:,1:1.5*specSamplingRate),2);
            % baseMatrixWhiskCortS_Vals_RH = baseWhiskCortS_Vals_RH.*ones(size(meanWhiskCortS_RH));
            % msStimWhiskS_Vals_RH = (meanWhiskCortS_RH - baseMatrixWhiskCortS_Vals_RH);

            T2 = -5:(1/specSamplingRate):15;
            %% save results
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

            AnalysisResults.(animalID).Whisk.Pupil.(whiskCriteriaName).Diameter.Diameter = meanWhiskPupilData;
            AnalysisResults.(animalID).Whisk.Pupil.(whiskCriteriaName).Diameter.DiameterStD = stdWhiskPupilData;
            AnalysisResults.(animalID).Whisk.Pupil.(whiskCriteriaName).Diameter.DiameterRaw = procWhiskPupilData;

            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData_LH;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData_LH;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData_LH;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData_LH;
            
            % AnalysisResults.(animalID).Whisk.cortical_RH.(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData_RH;
            % AnalysisResults.(animalID).Whisk.cortical_RH.(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData_RH;
            % AnalysisResults.(animalID).Whisk.cortical_RH.(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData_RH;
            % AnalysisResults.(animalID).Whisk.cortical_RH.(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData_RH;

            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).LFP.corticalS_LH = msStimWhiskS_Vals_LH;
            % AnalysisResults.(animalID).Whisk.cortical_RH.(whiskCriteriaName).LFP.corticalS_RH = msStimWhiskS_Vals_RH;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).LFP.T = T2;
            AnalysisResults.(animalID).Whisk.cortical_LH.(whiskCriteriaName).LFP.F = F;
        end
    end
    %% analyze movement evoked responses
    
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
                        
            [MovementPupilData] = EventData.Pupil.Diameter.movement.NormData(combMovementLogical,:);
    
            [allMovementCorticalMUAData_LH] = EventData.cortical_LH.corticalPower.movement.NormData(combMovementLogical,:);
            [allMovementCorticalGamData_LH] = EventData.cortical_LH.gammaBandPower.movement.NormData(combMovementLogical,:);

            % [allMovementCorticalMUAData_RH] = EventData.cortical_RH.corticalPower.movement.NormData(combMovementLogical,:);
            % [allMovementCorticalGamData_RH] = EventData.cortical_RH.gammaBandPower.movement.NormData(combMovementLogical,:);

            [allMovementFileIDs] = EventData.CBV.P_NE.movement.fileIDs(combMovementLogical,:);
            [allMovementEventTimes] = EventData.CBV.P_NE.movement.eventTime(combMovementLogical,:);
            allMovementDurations = EventData.CBV.P_NE.movement.duration(combMovementLogical,:);

            % keep only the data that occurs within the manually-approved awake regions
            [AChfinalMovementCBVData,finalMovementFileIDs,~,finalMovementFileEventTimes] = RemoveInvalidData_IOS(AChMovementCBVData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [NEfinalMovementCBVData,~,~,~] = RemoveInvalidData_IOS(NEMovementCBVData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
                   
            [AChfinalMovementGFPData,~,~,~] = RemoveInvalidData_IOS(AChMovementGFPData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [NEfinalMovementGFPData,~,~,~] = RemoveInvalidData_IOS(NEMovementGFPData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
    
            [finalMovementPupilData,~,~,~] = RemoveInvalidData_IOS(MovementPupilData,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);

            [finalMovementCorticalMUAData_LH,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalMUAData_LH,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            [finalMovementCorticalGamData_LH,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalGamData_LH,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);

            % [finalMovementCorticalMUAData_RH,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalMUAData_RH,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
            % [finalMovementCorticalGamData_RH,~,~,~] = RemoveInvalidData_IOS(allMovementCorticalGamData_RH,allMovementFileIDs,allMovementDurations,allMovementEventTimes,ManualDecisions);
                        
            % lowpass filter each movementing event and mean-subtract by the first 2 seconds
            clear finalMovementStartTimes finalMovementEndTimes finalMovementFiles procMovementPupilData NEprocMovementGFPData AChprocMovementGFPData NEprocMovementCBVData AChprocMovementCBVData procMovementCorticalMUAData_LH procMovementCorticalGamData_LH procMovementCorticalMUAData_RH procMovementCorticalGamData_RH
            dd = 1;
            for cc = 1:size(NEfinalMovementCBVData,1)
                movementStartTime = round(finalMovementFileEventTimes(cc,1),1) - 3;
                movementEndTime = movementStartTime + 15;
                finalMovementFileID = finalMovementFileIDs{cc,1};
                if movementStartTime >= 0.5 && movementEndTime <= (trialDuration_sec - 0.5)
                    NEmovementCBVarray = NEfinalMovementCBVData(cc,:);
                    AChmovementCBVarray = AChfinalMovementCBVData(cc,:);
    
                    NEmovementGFParray = NEfinalMovementGFPData(cc,:);
                    AChmovementGFParray = AChfinalMovementGFPData(cc,:);
                                        
                    movementPupilarray = finalMovementPupilData(cc,:);
    
                    movementCorticalMUAarray_LH = finalMovementCorticalMUAData_LH(cc,:);
                    movementCorticalGamArray_LH = finalMovementCorticalGamData_LH(cc,:);

                    % movementCorticalMUAarray_RH = finalMovementCorticalMUAData_RH(cc,:);
                    % movementCorticalGamArray_RH = finalMovementCorticalGamData_RH(cc,:);
                    
                    AChfiltMovementCBVarray = sgolayfilt(AChmovementCBVarray,3,9) - mean(AChmovementCBVarray(1:(offset*samplingRate)));
                    NEfiltMovementCBVarray = sgolayfilt(NEmovementCBVarray,3,9) - mean(NEmovementCBVarray(1:(offset*samplingRate)));
                    AChfiltMovementGFParray = sgolayfilt(AChmovementGFParray,3,9) - mean(AChmovementGFParray(1:(offset*samplingRate)));
                    NEfiltMovementGFParray = sgolayfilt(NEmovementGFParray,3,9) - mean(NEmovementGFParray(1:(offset*samplingRate)));
                    filtMovementPupilarray = sgolayfilt(movementPupilarray,3,9) - mean(movementPupilarray(1:(offset*samplingRate)));
        
                    AChprocMovementCBVData(dd,:) = AChfiltMovementCBVarray;
                    NEprocMovementCBVData(dd,:) = NEfiltMovementCBVarray;
                    AChprocMovementGFPData(dd,:) = AChfiltMovementGFParray;
                    NEprocMovementGFPData(dd,:) = NEfiltMovementGFParray;
                    procMovementPupilData(dd,:) = filtMovementPupilarray;

                    procMovementCorticalMUAData_LH(dd,:) = movementCorticalMUAarray_LH - mean(movementCorticalMUAarray_LH(1:(offset*samplingRate)));
                    procMovementCorticalGamData_LH(dd,:) = movementCorticalGamArray_LH - mean(movementCorticalGamArray_LH(1:(offset*samplingRate)));
                    % procMovementCorticalMUAData_RH(dd,:) = movementCorticalMUAarray_RH - mean(movementCorticalMUAarray_RH(1:(offset*samplingRate)));
                    % procMovementCorticalGamData_RH(dd,:) = movementCorticalGamArray_RH - mean(movementCorticalGamArray_RH(1:(offset*samplingRate)));

                    finalMovementStartTimes(dd,1) = movementStartTime;
                    finalMovementEndTimes(dd,1) = movementEndTime;
                    finalMovementFiles{dd,1} = finalMovementFileID;
                    dd = dd + 1;
                end
            end
            AChmeanMovementCBVData = mean(AChprocMovementCBVData,1);
            AChstdMovementCBVData = std(AChprocMovementCBVData,0,1);
            NEmeanMovementCBVData = mean(NEprocMovementCBVData,1);
            NEstdMovementCBVData = std(NEprocMovementCBVData,0,1);
    
            AChmeanMovementGFPData = mean(AChprocMovementGFPData,1);
            AChstdMovementGFPData = std(AChprocMovementGFPData,0,1);
            NEmeanMovementGFPData = mean(NEprocMovementGFPData,1);
            NEstdMovementGFPData = std(NEprocMovementGFPData,0,1);
    
            meanMovementPupilData = mean(procMovementPupilData,1);
            stdMovementPupilData = std(procMovementPupilData,0,1);

            meanMovementCorticalMUAData_LH = mean(procMovementCorticalMUAData_LH,1)*100;
            stdMovementCorticalMUAData_LH = std(procMovementCorticalMUAData_LH,0,1)*100;
            meanMovementCorticalGamData_LH = mean(procMovementCorticalGamData_LH,1)*100;
            stdMovementCorticalGamData_LH = std(procMovementCorticalGamData_LH,0,1)*100;

            % meanMovementCorticalMUAData_RH = mean(procMovementCorticalMUAData_RH,1)*100;
            % stdMovementCorticalMUAData_RH = std(procMovementCorticalMUAData_RH,0,1)*100;
            % meanMovementCorticalGamData_RH = mean(procMovementCorticalGamData_RH,1)*100;
            % stdMovementCorticalGamData_RH = std(procMovementCorticalGamData_RH,0,1)*100;
            %% extract ECoG spectrograms associated with the movementing indecies
            % LH            
            movementCorticalZhold_LH = [];
            for ee = 1:length(finalMovementFiles)
                % load normalized one-second bin data from each file
                movementFileID = finalMovementFiles{ee,1};
                movementSpecDataFileID = [animalID '_' movementFileID '_SpecDataB.mat'];
                movementSpecField = 'cortical_LH';
                for ff = 1:length(AllSpecData.(movementSpecField).fileIDs)
                    if strcmp(AllSpecData.(movementSpecField).fileIDs{ff,1},movementSpecDataFileID) == true
                        movementCorticalS_Data_LH = AllSpecData.(movementSpecField).normS{ff,1};
                        F = AllSpecData.(movementSpecField).F{ff,1};
                        T = round(AllSpecData.(movementSpecField).T{ff,1},1);
                    end
                end
                movementStartTimeIndex = find(T == round(finalMovementStartTimes(ee,1),1));
                movementStartTimeIndex = movementStartTimeIndex(1);
                movementDurationIndex = find(T == round(finalMovementEndTimes(ee,1),1));
                movementDurationIndex = movementDurationIndex(end);
                movementCorticalS_Vals_LH = movementCorticalS_Data_LH(:,movementStartTimeIndex:movementDurationIndex);
                movementCorticalZhold_LH = cat(3,movementCorticalZhold_LH,movementCorticalS_Vals_LH);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanMovementCortS_LH = mean(movementCorticalZhold_LH,3);
            baseMovementCortS_Vals_LH = mean(meanMovementCortS_LH(:,1:1.5*specSamplingRate),2);
            baseMatrixMovementCortS_Vals_LH = baseMovementCortS_Vals_LH.*ones(size(meanMovementCortS_LH));
            msStimMovementS_Vals_LH = (meanMovementCortS_LH - baseMatrixMovementCortS_Vals_LH);
                        
            % RH            
            % movementCorticalZhold_RH = [];
            % for ee = 1:length(finalMovementFiles)
            %     % load normalized one-second bin data from each file
            %     movementFileID = finalMovementFiles{ee,1};
            %     movementSpecDataFileID = [animalID '_' movementFileID '_SpecDataB.mat'];
            %     movementSpecField = 'cortical_RH';
            %     for ff = 1:length(AllSpecData.(movementSpecField).fileIDs)
            %         if strcmp(AllSpecData.(movementSpecField).fileIDs{ff,1},movementSpecDataFileID) == true
            %             movementCorticalS_Data_RH = AllSpecData.(movementSpecField).normS{ff,1};
            %             F = AllSpecData.(movementSpecField).F{ff,1};
            %             T = round(AllSpecData.(movementSpecField).T{ff,1},1);
            %         end
            %     end
            %     movementStartTimeIndex = find(T == round(finalMovementStartTimes(ee,1),1));
            %     movementStartTimeIndex = movementStartTimeIndex(1);
            %     movementDurationIndex = find(T == round(finalMovementEndTimes(ee,1),1));
            %     movementDurationIndex = movementDurationIndex(end);
            %     movementCorticalS_Vals_RH = movementCorticalS_Data_RH(:,movementStartTimeIndex:movementDurationIndex);
            %     movementCorticalZhold_RH = cat(3,movementCorticalZhold_RH,movementCorticalS_Vals_RH);
            % end
            % % cortical mean-subtract by first 2 seconds prior to stimulus
            % meanMovementCortS_RH = mean(movementCorticalZhold_RH,3);
            % baseMovementCortS_Vals_RH = mean(meanMovementCortS_RH(:,1:1.5*specSamplingRate),2);
            % baseMatrixMovementCortS_Vals_RH = baseMovementCortS_Vals_RH.*ones(size(meanMovementCortS_RH));
            % msStimMovementS_Vals_RH = (meanMovementCortS_RH - baseMatrixMovementCortS_Vals_RH);
     
            T2 = -5:(1/specSamplingRate):15;
            %% save results
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

            AnalysisResults.(animalID).Movement.Pupil.(movementCriteriaName).Diameter.Diamter = meanMovementPupilData;
            AnalysisResults.(animalID).Movement.Pupil.(movementCriteriaName).Diameter.DiameterStD = stdMovementPupilData;
            AnalysisResults.(animalID).Movement.Pupil.(movementCriteriaName).Diameter.DiameterRaw = procMovementPupilData;
    
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).MUA.corticalData = meanMovementCorticalMUAData_LH;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).MUA.corticalStD = stdMovementCorticalMUAData_LH;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).Gam.corticalData = meanMovementCorticalGamData_LH;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).Gam.corticalStD = stdMovementCorticalGamData_LH;

            % AnalysisResults.(animalID).Movement.cortical_RH.(movementCriteriaName).MUA.corticalData = meanMovementCorticalMUAData_RH;
            % AnalysisResults.(animalID).Movement.cortical_RH.(movementCriteriaName).MUA.corticalStD = stdMovementCorticalMUAData_RH;
            % AnalysisResults.(animalID).Movement.cortical_RH.(movementCriteriaName).Gam.corticalData = meanMovementCorticalGamData_RH;
            % AnalysisResults.(animalID).Movement.cortical_RH.(movementCriteriaName).Gam.corticalStD = stdMovementCorticalGamData_RH;

            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).LFP.corticalS_LH = msStimMovementS_Vals_LH;
            % AnalysisResults.(animalID).Movement.cortical_RH.(movementCriteriaName).LFP.corticalS_RH = msStimMovementS_Vals_RH;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).LFP.T = T2;
            AnalysisResults.(animalID).Movement.cortical_LH.(movementCriteriaName).LFP.F = F;
           
        end
    end
     %}
    %% analyze stimulus-evoked responses
  %{  
    if firstHrs == "true" || firstHrs == "false"
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
            [allStimPupilData] = EventData.Pupil.Diameter.stim.NormData(allStimFilter,:); 
    
            [allStimCortMUAData_LH] = EventData.cortical_LH.corticalPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData_LH] = EventData.cortical_LH.gammaBandPower.stim.NormData(allStimFilter,:);
            % [allStimCortMUAData_RH] = EventData.cortical_RH.corticalPower.stim.NormData(allStimFilter,:);
            % [allStimCortGamData_RH] = EventData.cortical_RH.gammaBandPower.stim.NormData(allStimFilter,:);

            [allStimFileIDs] = EventData.CBV.P_NE.stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.CBV.P_NE.stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalStimCBVData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS(NEallStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [AChfinalStimCBVData,~,~,~] = RemoveInvalidData_IOS(AChallStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [NEfinalStimGFPData,~,~,~] = RemoveInvalidData_IOS(NEallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);        
            [AChfinalStimGFPData,~,~,~] = RemoveInvalidData_IOS(AChallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimPupilData,~,~,~] = RemoveInvalidData_IOS(allStimPupilData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            
            [finalStimCortMUAData_LH,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData_LH,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData_LH,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData_LH,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % [finalStimCortMUAData_RH,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData_RH,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % [finalStimCortGamData_RH,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData_RH,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear finalStimStartTimes finalStimEndTimes finalStimFiles procStimPupilData NEprocStimGFPData AChprocStimGFPData NEprocStimCBVData AChprocStimCBVData procStimCorticalMUAData_LH procStimCorticalGamData_LH procStimCorticalMUAData_RH procStimCorticalGamData_RH
            ii = 1;
            for hh = 1:size(NEfinalStimCBVData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 3;
                stimEndTime = stimStartTime + 15;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    AChstimCBVarray = AChfinalStimCBVData(hh,:);
                    NEstimCBVarray = NEfinalStimCBVData(hh,:);
                    AChstimGFParray = AChfinalStimGFPData(hh,:);
                    NEstimGFParray = NEfinalStimGFPData(hh,:);
                    stimPupilarray = finalStimPupilData(hh,:);
  
                    stimCortMUAarray_LH = finalStimCortMUAData_LH(hh,:);
                    stimCortGamArray_LH = finalStimCortGamData_LH(hh,:);
                    % stimCortMUAarray_RH = finalStimCortMUAData_RH(hh,:);
                    % stimCortGamArray_RH = finalStimCortGamData_RH(hh,:);
    
                    AChfiltStimCBVarray = sgolayfilt(AChstimCBVarray,3,9);
                    NEfiltStimCBVarray = sgolayfilt(NEstimCBVarray,3,9);
                    AChfiltStimGFParray = sgolayfilt(AChstimGFParray,3,9);
                    NEfiltStimGFParray = sgolayfilt(NEstimGFParray,3,9);
                    filtStimPupilarray = sgolayfilt(stimPupilarray,3,9);
    
                    AChprocStimCBVData(hh,:) = AChfiltStimCBVarray - mean(AChfiltStimCBVarray(1:(offset*samplingRate)));
                    NEprocStimCBVData(hh,:) = NEfiltStimCBVarray - mean(NEfiltStimCBVarray(1:(offset*samplingRate)));
                    AChprocStimGFPData(hh,:) = AChfiltStimGFParray - mean(AChfiltStimGFParray(1:(offset*samplingRate)));
                    NEprocStimGFPData(hh,:) = NEfiltStimGFParray - mean(NEfiltStimGFParray(1:(offset*samplingRate)));
                    procStimPupilData(hh,:) = filtStimPupilarray - mean(filtStimPupilarray(1:(offset*samplingRate)));
 
                    procStimCortMUAData_LH(hh,:) = stimCortMUAarray_LH - mean(stimCortMUAarray_LH(1:(offset*samplingRate)));
                    procStimCortGamData_LH(hh,:) = stimCortGamArray_LH - mean(stimCortGamArray_LH(1:(offset*samplingRate)));
                    % procStimCortMUAData_RH(hh,:) = stimCortMUAarray_RH - mean(stimCortMUAarray_RH(1:(offset*samplingRate)));
                    % procStimCortGamData_RH(hh,:) = stimCortGamArray_RH - mean(stimCortGamArray_RH(1:(offset*samplingRate)));

                    finalStimStartTimes(ii,1) = stimStartTime;
                    finalStimEndTimes(ii,1) = stimEndTime;
                    finalStimFiles{ii,1} = finalStimFileID;
                    ii = ii + 1;
                end
            end
            AChmeanStimCBVData = mean(AChprocStimCBVData,1);
            AChstdStimCBVData = std(AChprocStimCBVData,0,1);
            NEmeanStimCBVData = mean(NEprocStimCBVData,1);
            NEstdStimCBVData = std(NEprocStimCBVData,0,1);
            AChmeanStimGFPData = mean(AChprocStimGFPData,1);
            AChstdStimGFPData = std(AChprocStimGFPData,0,1);
            NEmeanStimGFPData = mean(NEprocStimGFPData,1);
            NEstdStimGFPData = std(NEprocStimGFPData,0,1);
            meanStimPupilData = mean(procStimPupilData,1);
            stdStimPupilData = std(procStimPupilData,0,1);
    
            meanStimCortMUAData_LH = mean(procStimCortMUAData_LH,1)*100;
            stdStimCortMUAData_LH = std(procStimCortMUAData_LH,0,1)*100;
            meanStimCortGamData_LH = mean(procStimCortGamData_LH,1)*100;
            stdStimCortGamData_LH = std(procStimCortGamData_LH,0,1)*100;
            %% extract ECoG spectrograms associated with the stimuli indecies
            % left
            stimCortZhold_LH = [];
            for jj = 1:length(finalStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_LH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data_LH = AllSpecData.(stimSpecField).normS{kk,1};
                        F = AllSpecData.(stimSpecField).F{kk,1};
                        T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals_LH = stimCorticalS_Data_LH(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold_LH = cat(3,stimCortZhold_LH,stimCortS_Vals_LH);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanStimCortS_LH = mean(stimCortZhold_LH,3);
            baseStimCortS_Vals_LH = mean(meanStimCortS_LH(:,1:1.5*specSamplingRate),2);
            baseMatrixStimCortS_Vals_LH = baseStimCortS_Vals_LH.*ones(size(meanStimCortS_LH));
            msStimCortS_Vals_LH = (meanStimCortS_LH - baseMatrixStimCortS_Vals_LH);

            % right
            %             stimCortZhold_RH = [];
            % for jj = 1:length(finalStimFiles)
            %     % load normalized one-second bin data from each file
            %     stimFileID = finalStimFiles{jj,1};
            %     stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
            %     stimSpecField = 'cortical_RH';
            %     for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
            %         if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
            %             stimCorticalS_Data_RH = AllSpecData.(stimSpecField).normS{kk,1};
            %             F = AllSpecData.(stimSpecField).F{kk,1};
            %             T = round(AllSpecData.(stimSpecField).T{kk,1},1);
            %         end
            %     end
            %     stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
            %     stimStartTimeIndex = stimStartTimeIndex(1);
            %     stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
            %     stimDurationIndex = stimDurationIndex(end);
            %     stimCortS_Vals_RH = stimCorticalS_Data_RH(:,stimStartTimeIndex:stimDurationIndex);
            %     stimCortZhold_RH = cat(3,stimCortZhold_RH,stimCortS_Vals_RH);
            % end
            % % cortical mean-subtract by first 2 seconds prior to stimulus
            % meanStimCortS_RH = mean(stimCortZhold_RH,3);
            % baseStimCortS_Vals_RH = mean(meanStimCortS_RH(:,1:1.5*specSamplingRate),2);
            % baseMatrixStimCortS_Vals_RH = baseStimCortS_Vals_RH.*ones(size(meanStimCortS_RH));
            % msStimCortS_Vals_RH = (meanStimCortS_RH - baseMatrixStimCortS_Vals_RH);
                        
            T2 = -5:(1/specSamplingRate):15;
            %% save results
            AnalysisResults.(animalID).Stim.P_NE.(solenoid).count = size(procStimCortMUAData_LH,1);
    
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

            AnalysisResults.(animalID).Stim.Pupil.(solenoid).Diameter.Diameter= meanStimPupilData;
            AnalysisResults.(animalID).Stim.Pupil.(solenoid).Diameter.DiameterStD = stdStimPupilData;                
            AnalysisResults.(animalID).Stim.Pupil.(solenoid).Diameter.DiameterRaw= procStimPupilData;

            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).MUA.corticalData = meanStimCortMUAData_LH;
            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).MUA.corticalStD = stdStimCortMUAData_LH;
            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).Gam.corticalData = meanStimCortGamData_LH;
            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).Gam.corticalStD = stdStimCortGamData_LH;
            % 
            % AnalysisResults.(animalID).Stim.cortical_RH.(solenoid).MUA.corticalData = meanStimCortMUAData_RH;
            % AnalysisResults.(animalID).Stim.cortical_RH.(solenoid).MUA.corticalStD = stdStimCortMUAData_RH;
            % AnalysisResults.(animalID).Stim.cortical_RH.(solenoid).Gam.corticalData = meanStimCortGamData_RH;
            % AnalysisResults.(animalID).Stim.cortical_RH.(solenoid).Gam.corticalStD = stdStimCortGamData_RH;

            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).LFP.corticalS_LH = msStimCortS_Vals_LH;
            % AnalysisResults.(animalID).Stim.cortical_RH.(solenoid).LFP.corticalS_RH = msStimCortS_Vals_RH;

            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).Stim.cortical_LH.(solenoid).LFP.F = F;
        end
    end
    %}
    %}
    % save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end
end

