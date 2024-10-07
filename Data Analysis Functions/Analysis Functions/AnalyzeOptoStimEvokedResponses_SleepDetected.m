function [AnalysisResults] = AnalyzeOptoStimEvokedResponses_SleepDetected(animalID,rootFolder,AnalysisResults,firstHrs,saveFigs)
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

% criteria for optostimulation
StimCriteriaO.Value = {'OptoStim'};
StimCriteriaO.Fieldname = {'solenoidName'};
StimCriteriaO.Comparison = {'equal'};
OptostimCriteriaNames = {'stimCriteriaO'};
%%
samplingRate = EventData.CBV.P_NE.whisk.samplingRate;
specSamplingRate = 10;
trialDuration_sec = EventData.CBV.P_NE.whisk.trialDuration_sec;
timeVector = (0:(EventData.CBV.P_NE.whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV.P_NE.whisk.epoch.offset;
offset = EventData.CBV.P_NE.whisk.epoch.offset;
%% analyze opto stimulus-evoked responses
    if firstHrs == "false"
        for gg = 1:length(OptostimCriteriaNames)
            stimCriteriaName = OptostimCriteriaNames{1,gg};
            if strcmp(stimCriteriaName,'stimCriteriaO') == true
                OptoStimCriteria = StimCriteriaO;
                solenoid = 'OptoStim';
            end
            % pull data from EventData.mat structure
            allOptoStimFilter = FilterEvents_IOS(EventData.CBV.P_NE.stim,OptoStimCriteria);
            [AChallOptoStimCBVData] = EventData.CBV.P_ACh.stim.NormData(allOptoStimFilter,:);
            [NEallOptoStimCBVData] = EventData.CBV.P_NE.stim.NormData(allOptoStimFilter,:);
            [AChallOptoStimGFPData] = EventData.GFP.P_ACh.stim.NormData(allOptoStimFilter,:);
            [NEallOptoStimGFPData] = EventData.GFP.P_NE.stim.NormData(allOptoStimFilter,:); 
    
            [allOptoStimCortMUAData_LH] = EventData.cortical_LH.corticalPower.stim.NormData(allOptoStimFilter,:);
            [allOptoStimCortGamData_LH] = EventData.cortical_LH.gammaBandPower.stim.NormData(allOptoStimFilter,:);

            [allOptoStimCortMUAData_RH] = EventData.cortical_RH.corticalPower.stim.NormData(allOptoStimFilter,:);
            [allOptoStimCortGamData_RH] = EventData.cortical_RH.gammaBandPower.stim.NormData(allOptoStimFilter,:);
    
            [allOptoStimPupilData] = EventData.Pupil.Diameter.stim.NormData(allOptoStimFilter,:);

            [allOptoStimFileIDs] = EventData.CBV.P_NE.stim.fileIDs(allOptoStimFilter,:);
            [allOptoStimEventTimes] = EventData.CBV.P_NE.stim.eventTime(allOptoStimFilter,:);
            allOptoStimDurations = zeros(length(allOptoStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalOptoStimCBVData,finalOptoStimFileIDs,~,finalOptoStimFileEventTimes] = RemoveInvalidData_IOS(NEallOptoStimCBVData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [AChfinalOptoStimCBVData,~,~,~] = RemoveInvalidData_IOS(AChallOptoStimCBVData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [NEfinalOptoStimGFPData,~,~,~] = RemoveInvalidData_IOS(NEallOptoStimGFPData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);        
            [AChfinalOptoStimGFPData,~,~,~] = RemoveInvalidData_IOS(AChallOptoStimGFPData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [finalOptoStimCortMUAData_LH,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortMUAData_LH,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [finalOptoStimCortGamData_LH,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortGamData_LH,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            
            [finalOptoStimCortMUAData_RH,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortMUAData_RH,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            [finalOptoStimCortGamData_RH,~,~,~] = RemoveInvalidData_IOS(allOptoStimCortGamData_RH,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);

            [finalOptoStimPupilData,~,~,~] = RemoveInvalidData_IOS(allOptoStimPupilData,allOptoStimFileIDs,allOptoStimDurations,allOptoStimEventTimes,ManualDecisions);
            %% detect if the optostim cause any arousal or if the optostim occured during arousal
            [NoSleepFileTag,CauseArousalFileTag,CauseSleepFileTag] = RemoveSleepStimData(finalOptoStimFileIDs,finalOptoStimFileEventTimes,animalID);
            % separate stim data based on arousal state
            % OptoStim During Arousal
            NE_NoSleep_OptoStimCBVData = NEfinalOptoStimCBVData(NoSleepFileTag,:);
            NE_NoSleep_OptoStimGFPData = NEfinalOptoStimGFPData(NoSleepFileTag,:);
            ACh_NoSleep_OptoStimCBVData = AChfinalOptoStimCBVData(NoSleepFileTag,:);
            ACh_NoSleep_OptoStimGFPData = AChfinalOptoStimGFPData(NoSleepFileTag,:);

            NoSleep_OptoStimPupilData = finalOptoStimPupilData(NoSleepFileTag,:);

            NoSleep_OptoStimCortMUAData_LH = finalOptoStimCortMUAData_LH(NoSleepFileTag,:);
            NoSleep_OptoStimCortGamData_LH = finalOptoStimCortGamData_LH(NoSleepFileTag,:);
            NoSleep_OptoStimCortMUAData_RH = finalOptoStimCortMUAData_RH(NoSleepFileTag,:);
            NoSleep_OptoStimCorGamData_RH = finalOptoStimCortGamData_RH(NoSleepFileTag,:);

            NoSleep_OptoStimFileIDs = finalOptoStimFileIDs(NoSleepFileTag);
            NoSleep_OptoStimFileEventTimes = finalOptoStimFileEventTimes(NoSleepFileTag,1);

            % OptoStim Cause Arousal
            NE_CauseArousal_OptoStimCBVData = NEfinalOptoStimCBVData(CauseArousalFileTag,:);
            NE_CauseArousal_OptoStimGFPData = NEfinalOptoStimGFPData(CauseArousalFileTag,:);
            ACh_CauseArousal_OptoStimCBVData = AChfinalOptoStimCBVData(CauseArousalFileTag,:);
            ACh_CauseArousal_OptoStimGFPData = AChfinalOptoStimGFPData(CauseArousalFileTag,:);

            CauseArousal_OptoStimPupilData = finalOptoStimPupilData(CauseArousalFileTag,:);

            CauseArousal_OptoStimCortMUAData_LH = finalOptoStimCortMUAData_LH(CauseArousalFileTag,:);
            CauseArousal_OptoStimCortGamData_LH = finalOptoStimCortGamData_LH(CauseArousalFileTag,:);
            CauseArousal_OptoStimCortMUAData_RH = finalOptoStimCortMUAData_RH(CauseArousalFileTag,:);
            CauseArousal_OptoStimCorGamData_RH = finalOptoStimCortGamData_RH(CauseArousalFileTag,:);

            CauseArousal_OptoStimFileIDs = finalOptoStimFileIDs(CauseArousalFileTag);
            CauseArousal_OptoStimFileEventTimes = finalOptoStimFileEventTimes(CauseArousalFileTag,1);

            % OptoStim Does not Cause Arousal
            NE_CauseSleep_OptoStimCBVData = NEfinalOptoStimCBVData(CauseSleepFileTag,:);
            NE_CauseSleep_OptoStimGFPData = NEfinalOptoStimGFPData(CauseSleepFileTag,:);
            ACh_CauseSleep_OptoStimCBVData = AChfinalOptoStimCBVData(CauseSleepFileTag,:);
            ACh_CauseSleep_OptoStimGFPData = AChfinalOptoStimGFPData(CauseSleepFileTag,:);

            CauseSleep_OptoStimPupilData = finalOptoStimPupilData(CauseSleepFileTag,:);

            CauseSleep_OptoStimCortMUAData_LH = finalOptoStimCortMUAData_LH(CauseSleepFileTag,:);
            CauseSleep_OptoStimCortGamData_LH = finalOptoStimCortGamData_LH(CauseSleepFileTag,:);
            CauseSleep_OptoStimCortMUAData_RH = finalOptoStimCortMUAData_RH(CauseSleepFileTag,:);
            CauseSleep_OptoStimCorGamData_RH = finalOptoStimCortGamData_RH(CauseSleepFileTag,:);

            CauseSleep_OptoStimFileIDs = finalOptoStimFileIDs(CauseSleepFileTag);
            CauseSleep_OptoStimFileEventTimes = finalOptoStimFileEventTimes(CauseSleepFileTag,1);
            %% lowpass filter each stim event and mean-subtract by the first 2 seconds
            [AllOptoStimData] = OptoStimDataOrganize(NEfinalOptoStimCBVData, NEfinalOptoStimGFPData, AChfinalOptoStimCBVData, AChfinalOptoStimGFPData, finalOptoStimPupilData, finalOptoStimCortMUAData_LH, finalOptoStimCortGamData_LH, finalOptoStimCortMUAData_RH, finalOptoStimCortGamData_RH, finalOptoStimFileIDs, finalOptoStimFileEventTimes,solenoid,trialDuration_sec,offset,samplingRate,animalID,AllSpecData,specSamplingRate,timeVector);
            [NoSleep_OptoStimData] = OptoStimDataOrganize(NE_NoSleep_OptoStimCBVData, NE_NoSleep_OptoStimGFPData, ACh_NoSleep_OptoStimCBVData, ACh_NoSleep_OptoStimGFPData, NoSleep_OptoStimPupilData, NoSleep_OptoStimCortMUAData_LH, NoSleep_OptoStimCortGamData_LH, NoSleep_OptoStimCortMUAData_RH, NoSleep_OptoStimCorGamData_RH, NoSleep_OptoStimFileIDs, NoSleep_OptoStimFileEventTimes,solenoid,trialDuration_sec,offset,samplingRate,animalID,AllSpecData,specSamplingRate,timeVector);
            [CauseArousal_OptoStimData] = OptoStimDataOrganize(NE_CauseArousal_OptoStimCBVData, NE_CauseArousal_OptoStimGFPData, ACh_CauseArousal_OptoStimCBVData, ACh_CauseArousal_OptoStimGFPData, CauseArousal_OptoStimPupilData, CauseArousal_OptoStimCortMUAData_LH, CauseArousal_OptoStimCortGamData_LH, CauseArousal_OptoStimCortMUAData_RH, CauseArousal_OptoStimCorGamData_RH, CauseArousal_OptoStimFileIDs, CauseArousal_OptoStimFileEventTimes,solenoid,trialDuration_sec,offset,samplingRate,animalID,AllSpecData,specSamplingRate,timeVector);
            [CauseSleep_OptoStimData] = OptoStimDataOrganize(NE_CauseSleep_OptoStimCBVData, NE_CauseSleep_OptoStimGFPData, ACh_CauseSleep_OptoStimCBVData, ACh_CauseSleep_OptoStimGFPData, CauseSleep_OptoStimPupilData, CauseSleep_OptoStimCortMUAData_LH, CauseSleep_OptoStimCortGamData_LH, CauseSleep_OptoStimCortMUAData_RH, CauseSleep_OptoStimCorGamData_RH, CauseSleep_OptoStimFileIDs, CauseSleep_OptoStimFileEventTimes,solenoid,trialDuration_sec,offset,samplingRate,animalID,AllSpecData,specSamplingRate,timeVector);
            %% save results
            AnalysisResults.(animalID).OptoStim.AllOptoStim = AllOptoStimData;
            AnalysisResults.(animalID).OptoStim.NoSleep_OptoStim = NoSleep_OptoStimData;
            AnalysisResults.(animalID).OptoStim.CauseArousal_OptoStim = CauseArousal_OptoStimData;
            AnalysisResults.(animalID).OptoStim.CauseSleep_OptoStim = CauseSleep_OptoStimData;
        end
    end
 %% save data
    cd(rootFolder)
    if firstHrs == "false"
        save('AnalysisResults.mat','AnalysisResults','-v7.3')
    elseif firstHrs == "true"
        AnalysisResults_firstHrs = AnalysisResults;
        save('AnalysisResults_firstHrs','AnalysisResults_firstHrs','-v7.3')
    end
end

