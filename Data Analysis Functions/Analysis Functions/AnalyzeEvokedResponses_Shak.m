function [AnalysisResults] = AnalyzeEvokedResponses_Shak(animalID,rootFolder,AnalysisResults)
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
dataLocation = [rootFolder '\' animalID '\CombinedImaging'];
dataTypes = {'LH','RH'};
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
for aa = 1:length(dataTypes)
    dataType = dataTypes{1,aa};
    neuralDataType = ['cortical_' dataType];
%     if strcmp(setName,'IOS Set B') == true
%         dataType = 'adjBarrels';
%     end
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.Rhodamine.(dataType).whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.Rhodamine.(dataType).whisk.trialDuration_sec;
    timeVector = (0:(EventData.Rhodamine.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.Rhodamine.(dataType).whisk.epoch.offset;
    offset = EventData.Rhodamine.(dataType).whisk.epoch.offset;
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
        [whiskLogical] = FilterEvents_IOS(EventData.Rhodamine.(dataType).whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [allWhiskRhodamineData] = EventData.Rhodamine.(dataType).whisk.NormData(combWhiskLogical,:);
        [allWhiskGCaMP7sData] = EventData.GCaMP7s.(dataType).whisk.NormData(combWhiskLogical,:);

        [allWhiskCorticalMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskCorticalGamData] = EventData.(neuralDataType).gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskFileIDs] = EventData.Rhodamine.(dataType).whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.Rhodamine.(dataType).whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.Rhodamine.(dataType).whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [finalWhiskRhodamineData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS(allWhiskRhodamineData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskGCaMP7sData,~,~,~] = RemoveInvalidData_IOS(allWhiskGCaMP7sData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);

        [finalWhiskCorticalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskHippocampalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskCorticalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskHippocampalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskGCaMP7sData procWhiskRhodamineData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
        dd = 1;
        for cc = 1:size(finalWhiskRhodamineData,1)
            whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 5;%2;
            whiskEndTime = whiskStartTime + 15;
            finalWhiskFileID = finalWhiskFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                whiskRhodaminearray = finalWhiskRhodamineData(cc,:);
                whiskGCaMP7sarray = finalWhiskGCaMP7sData(cc,:);

                whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
                whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
                whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                
                filtWhiskRhodaminearray = sgolayfilt(whiskRhodaminearray,3,17);
                filtWhiskGCaMP7sarray = sgolayfilt(whiskGCaMP7sarray,3,17);

                procWhiskRhodamineData(dd,:) = filtWhiskRhodaminearray;%- mean(filtWhiskRhodaminearray(1:(offset*samplingRate))); %#ok<*AGROW>
                procWhiskGCaMP7sData(dd,:) = filtWhiskGCaMP7sarray;% - mean(filtWhiskGCaMP7sarray(1:(offset*samplingRate))); %#ok<*AGROW>

                procWhiskCorticalMUAData(dd,:) = whiskCorticalMUAarray - mean(whiskCorticalMUAarray(1:(offset*samplingRate)));
                procWhiskHippocampalMUAData(dd,:) = whiskHippocampalMUAarray - mean(whiskHippocampalMUAarray(1:(offset*samplingRate)));
                procWhiskCorticalGamData(dd,:) = whiskCorticalGamArray - mean(whiskCorticalGamArray(1:(offset*samplingRate)));
                procWhiskHippocampalGamData(dd,:) = whiskHippocampalGamArray - mean(whiskHippocampalGamArray(1:(offset*samplingRate)));
                finalWhiskStartTimes(dd,1) = whiskStartTime;
                finalWhiskEndTimes(dd,1) = whiskEndTime;
                finalWhiskFiles{dd,1} = finalWhiskFileID;
                dd = dd + 1;
            end
        end
        meanWhiskRhodamineData = mean(procWhiskRhodamineData,1);%*100;
        stdWhiskRhodamineData = std(procWhiskRhodamineData,0,1);%*100;
        meanWhiskGCaMP7sData = mean(procWhiskGCaMP7sData,1);%*100;
        stdWhiskGCaMP7sData = std(procWhiskGCaMP7sData,0,1);%*100;

        meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
        stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
        meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
        stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
        meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
        stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
        meanWhiskHippocampalGamData = mean(procWhiskHippocampalGamData,1)*100;
        stdWhiskHippocampalGamData = std(procWhiskHippocampalGamData,0,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
        whiskHippocampalZhold = [];
        for ee = 1:length(finalWhiskFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = neuralDataType;
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                    whiskHippocampalS_Data = AllSpecData.hippocampus.normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskCorticalZhold = cat(3,whiskCorticalZhold,whiskCorticalS_Vals);
            whiskHippocampalZhold = cat(3,whiskHippocampalZhold,whiskHippocampalS_Vals);
        end
        % cortical mean-subtract by first 2 seconds prior to stimulus
        meanWhiskCortS = mean(whiskCorticalZhold,3);
        baseWhiskCortS_Vals = mean(meanWhiskCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixWhiskCortS_Vals = baseWhiskCortS_Vals.*ones(size(meanWhiskCortS));
        msStimWhiskS_Vals = (meanWhiskCortS - baseMatrixWhiskCortS_Vals);
        % hippocampal mean-subtract by first 2 seconds prior to stimulus
        meanWhiskHipS = mean(whiskHippocampalZhold,3);
        baseWhiskHipS_Vals = mean(meanWhiskHipS(:,1:1.5*specSamplingRate),2);
        baseMatrixWhiskHipS_Vals = baseWhiskHipS_Vals.*ones(size(meanWhiskHipS));
        msWhiskHipS_Vals = (meanWhiskHipS - baseMatrixWhiskHipS_Vals);
        T2 = -5:(1/specSamplingRate):10;
        % save results
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Rhodamine.Rhodamine = meanWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Rhodamine.RhodamineStD = stdWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).GCaMP7s.GCaMP7s = meanWhiskGCaMP7sData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).GCaMP7s.GCaMP7sStD = stdWhiskGCaMP7sData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalData = meanWhiskHippocampalGamData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalStD = stdWhiskHippocampalGamData;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.corticalS = msStimWhiskS_Vals;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.hippocampalS = msWhiskHipS_Vals;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.T = T2;
        AnalysisResults.(animalID).Whisk.(dataType).(whiskCriteriaName).LFP.F = F;
    end
    %% analyze stimulus-evoked responses
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
        allStimFilter = FilterEvents_IOS(EventData.Rhodamine.(dataType).stim,StimCriteria);
        [allStimRhodamineData] = EventData.Rhodamine.(dataType).stim.NormData(allStimFilter,:);
        [allStimGCaMP7sData] = EventData.GCaMP7s.(dataType).stim.NormData(allStimFilter,:); % EventData.GCaMP7s.(dataType).stim.NormData(allStimFilter,:); 

        [allStimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(allStimFilter,:);
        [allStimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(allStimFilter,:);
        [allStimCortGamData] = EventData.(neuralDataType).gammaBandPower.stim.NormData(allStimFilter,:);
        [allStimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(allStimFilter,:);
        [allStimFileIDs] = EventData.Rhodamine.(dataType).stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.Rhodamine.(dataType).stim.eventTime(allStimFilter,:);
        allStimDurations = zeros(length(allStimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [finalStimRhodamineData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS(allStimRhodamineData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimGCaMP7sData,~,~,~] = RemoveInvalidData_IOS(allStimGCaMP7sData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimHipMUAData,~,~,~] = RemoveInvalidData_IOS(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimHipGamData,~,~,~] = RemoveInvalidData_IOS(allStimHipGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        % lowpass filter each stim event and mean-subtract by the first 2 seconds
        clear procStimGCaMP7sData procStimRhodamineData procStimRhodamineData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
        ii = 1;
        for hh = 1:size(finalStimRhodamineData,1)
            stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 5;%2;
            stimEndTime = stimStartTime + 15;
            finalStimFileID = finalStimFileIDs{hh,1};
            if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                stimRhodaminearray = finalStimRhodamineData(hh,:);
                stimGCaMP7sarray = finalStimGCaMP7sData(hh,:);

                stimCortMUAarray = finalStimCortMUAData(hh,:);
                stimHipMUAarray = finalStimHipMUAData(hh,:);
                stimCortGamArray = finalStimCortGamData(hh,:);
                stimHipGamArray = finalStimHipGamData(hh,:);

                filtStimRhodaminearray = sgolayfilt(stimRhodaminearray,3,17);
                filtStimGCaMP7sarray = sgolayfilt(stimGCaMP7sarray,3,17);

                procStimRhodamineData(hh,:) = filtStimRhodaminearray- mean(filtStimRhodaminearray(1:(offset*samplingRate)));
                procStimGCaMP7sData(hh,:) = filtStimGCaMP7sarray- mean(filtStimGCaMP7sarray(1:(offset*samplingRate)));

                procStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
                procStimHipMUAData(hh,:) = stimHipMUAarray - mean(stimHipMUAarray(1:(offset*samplingRate)));
                procStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
                procStimHipGamData(hh,:) = stimHipGamArray - mean(stimHipGamArray(1:(offset*samplingRate)));
                finalStimStartTimes(ii,1) = stimStartTime;
                finalStimEndTimes(ii,1) = stimEndTime;
                finalStimFiles{ii,1} = finalStimFileID;
                ii = ii + 1;
            end
        end
        meanStimRhodamineData = mean(procStimRhodamineData,1);%*100;
        stdStimRhodamineData = std(procStimRhodamineData,0,1);%*100;
        meanStimGCaMP7sData = mean(procStimGCaMP7sData,1);%*100;
        stdStimGCaMP7sData = std(procStimGCaMP7sData,0,1);%*100;

        meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
        stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
        meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
        stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
        meanStimCortGamData = mean(procStimCortGamData,1)*100;
        stdStimCortGamData = std(procStimCortGamData,0,1)*100;
        meanStimHipGamData = mean(procStimHipGamData,1)*100;
        stdStimHipGamData = std(procStimHipGamData,0,1)*100;
        % extract LFP from spectrograms associated with the stimuli indecies
        stimCortZhold = [];
        stimHipZhold = [];
        for jj = 1:length(finalStimFiles)
            % load normalized one-second bin data from each file
            stimFileID = finalStimFiles{jj,1};
            stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
            stimSpecField = neuralDataType;
            for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                    stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                    stimHippocampalS_Data = AllSpecData.hippocampus.normS{kk,1};
                end
            end
            stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
            stimStartTimeIndex = stimStartTimeIndex(1);
            stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
            stimDurationIndex = stimDurationIndex(end);
            stimCortS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
            stimHipS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
            stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
            stimHipZhold = cat(3,stimHipZhold,stimHipS_Vals);
        end
        % cortical mean-subtract by first 2 seconds prior to stimulus
        meanStimCortS = mean(stimCortZhold,3);
        baseStimCortS_Vals = mean(meanStimCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixStimCortS_Vals = baseStimCortS_Vals.*ones(size(meanStimCortS));
        msStimCortS_Vals = (meanStimCortS - baseMatrixStimCortS_Vals);
        % hippocampal mean-subtract by first 2 seconds prior to stimulus
        meanStimHipS = mean(stimHipZhold,3);
        baseStimHipS_Vals = mean(meanStimHipS(:,1:1.5*specSamplingRate),2);
        baseMatrixStimHipS_Vals = baseStimHipS_Vals.*ones(size(meanStimHipS));
        msStimHipS_Vals = (meanStimHipS - baseMatrixStimHipS_Vals);
        % save results
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).count = size(procStimHipMUAData,1);
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Rhodamine.Rhodamine = meanStimRhodamineData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Rhodamine.RhodamineStD = stdStimRhodamineData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).GCaMP7s.GCaMP7s= meanStimGCaMP7sData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).GCaMP7s.GCaMP7sStD = stdStimGCaMP7sData;

        AnalysisResults.(animalID).Stim.(dataType).(solenoid).MUA.corticalData = meanStimCortMUAData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).MUA.corticalStD = stdStimCortMUAData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).MUA.hippocampalData = meanStimHipMUAData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Gam.corticalData = meanStimCortGamData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Gam.corticalStD = stdStimCortGamData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Gam.hippocampalData = meanStimHipGamData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).Gam.hippocampalStD = stdStimHipGamData;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).timeVector = timeVector;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).LFP.corticalS = msStimCortS_Vals;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).LFP.hippocampalS = msStimHipS_Vals;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).LFP.T = T2;
        AnalysisResults.(animalID).Stim.(dataType).(solenoid).LFP.F = F;
    end
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end

