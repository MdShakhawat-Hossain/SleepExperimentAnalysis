function [AnalysisResults] = AnalyzeEvokedResponses_EEG(animalID,rootFolder,AnalysisResults)
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
    samplingRate = EventData.Rhodamine.RH.whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.Rhodamine.RH.whisk.trialDuration_sec;
    timeVector = (0:(EventData.Rhodamine.RH.whisk.epoch.duration*samplingRate))/samplingRate - EventData.Rhodamine.RH.whisk.epoch.offset;
    offset = EventData.Rhodamine.RH.whisk.epoch.offset;
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
        [whiskLogical] = FilterEvents_IOS(EventData.Rhodamine.RH.whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [RHWhiskRhodamineData] = EventData.Rhodamine.RH.whisk.NormData(combWhiskLogical,:);
        [RHWhiskGFPData] = EventData.GFP.RH.whisk.NormData(combWhiskLogical,:);

        [allWhiskEEGMUAData] = EventData.EEG_LH.EEGPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskEEGGamData] = EventData.EEG_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);

        [allWhiskFileIDs] = EventData.Rhodamine.RH.whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.Rhodamine.RH.whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.Rhodamine.RH.whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [RHfinalWhiskRhodamineData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS(RHWhiskRhodamineData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
               
        [RHfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(RHWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);

        [finalWhiskEEGMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskEEGMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskEEGGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskEEGGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskGFPData procWhiskRhodamineData procWhiskEEGMUAData  procWhiskEEGGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
        dd = 1;
        for cc = 1:size(RHfinalWhiskRhodamineData,1)
            whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 5;%2;
            whiskEndTime = whiskStartTime + 15;
            finalWhiskFileID = finalWhiskFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                RHwhiskRhodaminearray = RHfinalWhiskRhodamineData(cc,:);

                RHwhiskGFParray = RHfinalWhiskGFPData(cc,:);

                whiskEEGMUAarray = finalWhiskEEGMUAData(cc,:);
                whiskEEGGamArray = finalWhiskEEGGamData(cc,:);
                
                RHfiltWhiskRhodaminearray = sgolayfilt(RHwhiskRhodaminearray,3,17);
                RHfiltWhiskGFParray = sgolayfilt(RHwhiskGFParray,3,17);

                RHprocWhiskRhodamineData(dd,:) = RHfiltWhiskRhodaminearray - mean(RHfiltWhiskRhodaminearray(1:(offset*samplingRate)));
                RHprocWhiskGFPData(dd,:) = RHfiltWhiskGFParray - mean(RHfiltWhiskGFParray(1:(offset*samplingRate)));

                procWhiskEEGMUAData(dd,:) = whiskEEGMUAarray - mean(whiskEEGMUAarray(1:(offset*samplingRate)));
                procWhiskEEGGamData(dd,:) = whiskEEGGamArray - mean(whiskEEGGamArray(1:(offset*samplingRate)));
                finalWhiskStartTimes(dd,1) = whiskStartTime;
                finalWhiskEndTimes(dd,1) = whiskEndTime;
                finalWhiskFiles{dd,1} = finalWhiskFileID;
                dd = dd + 1;
            end
        end
        RHmeanWhiskRhodamineData = mean(RHprocWhiskRhodamineData,1);%*100;
        RHstdWhiskRhodamineData = std(RHprocWhiskRhodamineData,0,1);%*100;

        RHmeanWhiskGFPData = mean(RHprocWhiskGFPData,1);%*100;
        RHstdWhiskGFPData = std(RHprocWhiskGFPData,0,1);%*100;

        meanWhiskEEGMUAData = mean(procWhiskEEGMUAData,1)*100;
        stdWhiskEEGMUAData = std(procWhiskEEGMUAData,0,1)*100;
        meanWhiskEEGGamData = mean(procWhiskEEGGamData,1)*100;
        stdWhiskEEGGamData = std(procWhiskEEGGamData,0,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskEEGZhold = [];
        for ee = 1:length(finalWhiskFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = 'EEG_LH';
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskEEGS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskEEGS_Vals = whiskEEGS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskEEGZhold = cat(3,whiskEEGZhold,whiskEEGS_Vals);
        end
        % EEG mean-subtract by first 2 seconds prior to stimulus
        meanWhiskCortS = mean(whiskEEGZhold,3);
        baseWhiskCortS_Vals = mean(meanWhiskCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixWhiskCortS_Vals = baseWhiskCortS_Vals.*ones(size(meanWhiskCortS));
        msStimWhiskS_Vals = (meanWhiskCortS - baseMatrixWhiskCortS_Vals);
        T2 = -20:(1/specSamplingRate):20;
        % save results
        AnalysisResults.(animalID).Whisk.RH.(whiskCriteriaName).Rhodamine.Rhodamine = RHmeanWhiskRhodamineData;
        AnalysisResults.(animalID).Whisk.RH.(whiskCriteriaName).Rhodamine.RhodamineStD = RHstdWhiskRhodamineData;

        AnalysisResults.(animalID).Whisk.RH.(whiskCriteriaName).GFP.GFP = RHmeanWhiskGFPData;
        AnalysisResults.(animalID).Whisk.RH.(whiskCriteriaName).GFP.GFPStD = RHstdWhiskGFPData;

        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).MUA.EEGData = meanWhiskEEGMUAData;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).MUA.EEGStD = stdWhiskEEGMUAData;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).Gam.EEGData = meanWhiskEEGGamData;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).Gam.EEGStD = stdWhiskEEGGamData;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).LFP.EEGS = msStimWhiskS_Vals;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).LFP.T = T2;
        AnalysisResults.(animalID).Whisk.EEG.(whiskCriteriaName).LFP.F = F;
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
        allStimFilter = FilterEvents_IOS(EventData.Rhodamine.RH.stim,StimCriteria);
        [RHallStimRhodamineData] = EventData.Rhodamine.RH.stim.NormData(allStimFilter,:);
        [RHallStimGFPData] = EventData.GFP.RH.stim.NormData(allStimFilter,:); 

        [allStimCortMUAData] = EventData.EEG_LH.EEGPower.stim.NormData(allStimFilter,:);
        [allStimCortGamData] = EventData.EEG_LH.gammaBandPower.stim.NormData(allStimFilter,:);
        [allStimFileIDs] = EventData.Rhodamine.RH.stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.Rhodamine.RH.stim.eventTime(allStimFilter,:);
        allStimDurations = zeros(length(allStimEventTimes),1);
        % keep only the data that occurs within the manually-approved awake regions
        [RHfinalStimRhodamineData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS(RHallStimRhodamineData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [RHfinalStimGFPData,~,~,~] = RemoveInvalidData_IOS(RHallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);        
        [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        % lowpass filter each stim event and mean-subtract by the first 2 seconds
        clear procStimGFPData procStimRhodamineData procStimRhodamineData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
        ii = 1;
        for hh = 1:size(RHfinalStimRhodamineData,1)
            stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 5;%2;
            stimEndTime = stimStartTime + 15;
            finalStimFileID = finalStimFileIDs{hh,1};
            if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                RHstimRhodaminearray = RHfinalStimRhodamineData(hh,:);
                RHstimGFParray = RHfinalStimGFPData(hh,:);

                stimCortMUAarray = finalStimCortMUAData(hh,:);
                stimCortGamArray = finalStimCortGamData(hh,:);

                RHfiltStimRhodaminearray = sgolayfilt(RHstimRhodaminearray,3,17);
                RHfiltStimGFParray = sgolayfilt(RHstimGFParray,3,17);

                RHprocStimRhodamineData(hh,:) = RHfiltStimRhodaminearray- mean(RHfiltStimRhodaminearray(1:(offset*samplingRate)));
                RHprocStimGFPData(hh,:) = RHfiltStimGFParray- mean(RHfiltStimGFParray(1:(offset*samplingRate)));

                procStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
                procStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
                finalStimStartTimes(ii,1) = stimStartTime;
                finalStimEndTimes(ii,1) = stimEndTime;
                finalStimFiles{ii,1} = finalStimFileID;
                ii = ii + 1;
            end
        end
        RHmeanStimRhodamineData = mean(RHprocStimRhodamineData,1);%*100;
        RHstdStimRhodamineData = std(RHprocStimRhodamineData,0,1);%*100;
        RHmeanStimGFPData = mean(RHprocStimGFPData,1);%*100;
        RHstdStimGFPData = std(RHprocStimGFPData,0,1);%*100;

        meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
        stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
        meanStimCortGamData = mean(procStimCortGamData,1)*100;
        stdStimCortGamData = std(procStimCortGamData,0,1)*100;
        % extract LFP from spectrograms associated with the stimuli indecies
        stimCortZhold = [];
        for jj = 1:length(finalStimFiles)
            % load normalized one-second bin data from each file
            stimFileID = finalStimFiles{jj,1};
            stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
            stimSpecField = 'EEG_LH';
            for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                    stimEEGS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                end
            end
            stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
            stimStartTimeIndex = stimStartTimeIndex(1);
            stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
            stimDurationIndex = stimDurationIndex(end);
            stimCortS_Vals = stimEEGS_Data(:,stimStartTimeIndex:stimDurationIndex);
            stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
        end
        % EEG mean-subtract by first 2 seconds prior to stimulus
        meanStimCortS = mean(stimCortZhold,3);
        baseStimCortS_Vals = mean(meanStimCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixStimCortS_Vals = baseStimCortS_Vals.*ones(size(meanStimCortS));
        msStimCortS_Vals = (meanStimCortS - baseMatrixStimCortS_Vals);
        % save results
        AnalysisResults.(animalID).Stim.RH.(solenoid).count = size(procStimCortMUAData,1);

        AnalysisResults.(animalID).Stim.RH.(solenoid).Rhodamine.Rhodamine = RHmeanStimRhodamineData;
        AnalysisResults.(animalID).Stim.RH.(solenoid).Rhodamine.RhodamineStD = RHstdStimRhodamineData;
        AnalysisResults.(animalID).Stim.RH.(solenoid).GFP.GFP= RHmeanStimGFPData;
        AnalysisResults.(animalID).Stim.RH.(solenoid).GFP.GFPStD = RHstdStimGFPData;

        AnalysisResults.(animalID).Stim.EEG.(solenoid).MUA.EEGData = meanStimCortMUAData;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).MUA.EEGStD = stdStimCortMUAData;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).Gam.EEGData = meanStimCortGamData;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).Gam.EEGStD = stdStimCortGamData;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).timeVector = timeVector;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).LFP.EEGS = msStimCortS_Vals;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).LFP.T = T2;
        AnalysisResults.(animalID).Stim.EEG.(solenoid).LFP.F = F;
    end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults','-v7.3')

end

