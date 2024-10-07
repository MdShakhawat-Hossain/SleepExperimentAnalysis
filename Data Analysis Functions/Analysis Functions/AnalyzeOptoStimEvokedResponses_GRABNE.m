function [AnalysisResults] = AnalyzeOptoStimEvokedResponses_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs,saveFigs)
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
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procOptoStimGFPData procOptoStimCBVData procOptoStimCBVData procOptoStimCortMUAData_LH procOptoStimCortMUAData_RH procOptoStimPupilData procOptoStimCortGamData_LH procOptoStimCortGamData_RH finalOptoStimStartTimes finalOptoStimEndTimes finalOptoStimFiles
            ii = 1;
            for hh = 1:size(NEfinalOptoStimCBVData,1)
                stimStartTime = round(finalOptoStimFileEventTimes(hh,1),1) - 2;
                stimEndTime = stimStartTime + 15;
                finalOptoStimFileID = finalOptoStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    AChstimCBVarray = AChfinalOptoStimCBVData(hh,:);
                    NEstimCBVarray = NEfinalOptoStimCBVData(hh,:);
                    AChstimGFParray = AChfinalOptoStimGFPData(hh,:);
                    NEstimGFParray = NEfinalOptoStimGFPData(hh,:);
    
                    stimCortMUAarray_LH = finalOptoStimCortMUAData_LH(hh,:);
                    stimCortGamArray_LH = finalOptoStimCortGamData_LH(hh,:);
                    stimCortMUAarray_RH = finalOptoStimCortMUAData_RH(hh,:);
                    stimCortGamArray_RH = finalOptoStimCortGamData_RH(hh,:);

                    stimPupilarray = finalOptoStimPupilData(hh,:);
    
                    AChfiltOptoStimCBVarray = sgolayfilt(AChstimCBVarray,3,9);
                    NEfiltOptoStimCBVarray = sgolayfilt(NEstimCBVarray,3,9);
                    AChfiltOptoStimGFParray = sgolayfilt(AChstimGFParray,3,9);
                    NEfiltOptoStimGFParray = sgolayfilt(NEstimGFParray,3,9);
                    PupilfiltOptoStimarray = sgolayfilt(stimPupilarray,3,9);
    
                    AChprocOptoStimCBVData(hh,:) = AChfiltOptoStimCBVarray - mean(AChfiltOptoStimCBVarray(1:(offset*samplingRate)));
                    NEprocOptoStimCBVData(hh,:) = NEfiltOptoStimCBVarray - mean(NEfiltOptoStimCBVarray(1:(offset*samplingRate)));
                    AChprocOptoStimGFPData(hh,:) = AChfiltOptoStimGFParray - mean(AChfiltOptoStimGFParray(1:(offset*samplingRate)));
                    NEprocOptoStimGFPData(hh,:) = NEfiltOptoStimGFParray - mean(NEfiltOptoStimGFParray(1:(offset*samplingRate)));
                    procOptoStimPupilData(hh,:) = PupilfiltOptoStimarray - mean(PupilfiltOptoStimarray(1:(offset*samplingRate)));
  
                    procOptoStimCortMUAData_LH(hh,:) = stimCortMUAarray_LH - mean(stimCortMUAarray_LH(1:(offset*samplingRate)));
                    procOptoStimCortGamData_LH(hh,:) = stimCortGamArray_LH - mean(stimCortGamArray_LH(1:(offset*samplingRate)));
                    
                    procOptoStimCortMUAData_RH(hh,:) = stimCortMUAarray_RH - mean(stimCortMUAarray_RH(1:(offset*samplingRate)));
                    procOptoStimCortGamData_RH(hh,:) = stimCortGamArray_RH - mean(stimCortGamArray_RH(1:(offset*samplingRate)));
                    
                    finalOptoStimStartTimes(ii,1) = stimStartTime;
                    finalOptoStimEndTimes(ii,1) = stimEndTime;
                    finalOptoStimFiles{ii,1} = finalOptoStimFileID;
                    ii = ii + 1;
                end
            end
            AChmeanOptoStimCBVData = mean(AChprocOptoStimCBVData,1);
            AChstdOptoStimCBVData = std(AChprocOptoStimCBVData,0,1);
            NEmeanOptoStimCBVData = mean(NEprocOptoStimCBVData,1);
            NEstdOptoStimCBVData = std(NEprocOptoStimCBVData,0,1);
            AChmeanOptoStimGFPData = mean(AChprocOptoStimGFPData,1);
            AChstdOptoStimGFPData = std(AChprocOptoStimGFPData,0,1);
            NEmeanOptoStimGFPData = mean(NEprocOptoStimGFPData,1);
            NEstdOptoStimGFPData = std(NEprocOptoStimGFPData,0,1);
    
            meanOptoStimCortMUAData_LH = mean(procOptoStimCortMUAData_LH,1)*100;
            stdOptoStimCortMUAData_LH = std(procOptoStimCortMUAData_LH,0,1)*100;
            meanOptoStimCortGamData_LH = mean(procOptoStimCortGamData_LH,1)*100;
            stdOptoStimCortGamData_LH = std(procOptoStimCortGamData_LH,0,1)*100;

            meanOptoStimCortMUAData_RH = mean(procOptoStimCortMUAData_RH,1)*100;
            stdOptoStimCortMUAData_RH = std(procOptoStimCortMUAData_RH,0,1)*100;
            meanOptoStimCortGamData_RH = mean(procOptoStimCortGamData_RH,1)*100;
            stdOptoStimCortGamData_RH = std(procOptoStimCortGamData_RH,0,1)*100;

            meanOptoStimPupilData = mean(procOptoStimPupilData,1);
            stdOptoStimPupilData = std(procOptoStimPupilData,0,1);
            %% extract ECoG spectrograms associated with the stimuli indecies
            % left hemispheres
            stimCortZhold_LH = [];
            for jj = 1:length(finalOptoStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalOptoStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_LH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data_LH = AllSpecData.(stimSpecField).normS{kk,1};
                        F = AllSpecData.(stimSpecField).F{kk,1};
                        T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                    end
                end
                stimStartTimeIndex = find(T == round(finalOptoStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalOptoStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals_LH = stimCorticalS_Data_LH(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold_LH = cat(3,stimCortZhold_LH,stimCortS_Vals_LH);
            end
            % cortical mean-subtract by first 3 seconds prior to stimulus
            meanOptoStimCortS_LH = mean(stimCortZhold_LH,3);
            baseOptoStimCortS_Vals_LH = mean(meanOptoStimCortS_LH(:,1:1.5*specSamplingRate),2);
            baseMatrixOptoStimCortS_Vals_LH = baseOptoStimCortS_Vals_LH.*ones(size(meanOptoStimCortS_LH));
            msOptoStimCortS_Vals_LH = (meanOptoStimCortS_LH - baseMatrixOptoStimCortS_Vals_LH);  
            
            % right hemispheres
            stimCortZhold_RH = [];
            for jj = 1:length(finalOptoStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalOptoStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = 'cortical_RH';
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data_RH = AllSpecData.(stimSpecField).normS{kk,1};
                        F = AllSpecData.(stimSpecField).F{kk,1};
                        T = round(AllSpecData.(stimSpecField).T{kk,1},1);
                    end
                end
                stimStartTimeIndex = find(T == round(finalOptoStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalOptoStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals_RH = stimCorticalS_Data_RH(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold_RH = cat(3,stimCortZhold_RH,stimCortS_Vals_RH);
            end
            % cortical mean-subtract by first 3 seconds prior to stimulus
            meanOptoStimCortS_RH = mean(stimCortZhold_RH,3);
            baseOptoStimCortS_Vals_RH = mean(meanOptoStimCortS_RH(:,1:1.5*specSamplingRate),2);
            baseMatrixOptoStimCortS_Vals_RH = baseOptoStimCortS_Vals_RH.*ones(size(meanOptoStimCortS_RH));
            msOptoStimCortS_Vals_RH = (meanOptoStimCortS_RH - baseMatrixOptoStimCortS_Vals_RH);   
            
            T2 = -5:(1/specSamplingRate):15;
            %% save results
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).count = size(procOptoStimCortMUAData_LH,1);
    
            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).CBV.CBV = AChmeanOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).CBV.CBVStD = AChstdOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).CBV.CBV = NEmeanOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).CBV.CBVStD = NEstdOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).GFP.GFP= AChmeanOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).GFP.GFPStD = AChstdOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).GFP.GFP= NEmeanOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).GFP.GFPStD = NEstdOptoStimGFPData;

            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).CBV.CBVRaw = AChprocOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).CBV.CBVRaw = NEprocOptoStimCBVData;
            AnalysisResults.(animalID).OptoStim.P_ACh.(solenoid).GFP.GFPRaw = AChprocOptoStimGFPData;
            AnalysisResults.(animalID).OptoStim.P_NE.(solenoid).GFP.GFPRaw = NEprocOptoStimGFPData;
    
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).MUA.corticalData = meanOptoStimCortMUAData_LH;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).MUA.corticalStD = stdOptoStimCortMUAData_LH;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).Gam.corticalData = meanOptoStimCortGamData_LH;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).Gam.corticalStD = stdOptoStimCortGamData_LH;

            AnalysisResults.(animalID).OptoStim.cortical_RH.(solenoid).MUA.corticalData = meanOptoStimCortMUAData_RH;
            AnalysisResults.(animalID).OptoStim.cortical_RH.(solenoid).MUA.corticalStD = stdOptoStimCortMUAData_RH;
            AnalysisResults.(animalID).OptoStim.cortical_RH.(solenoid).Gam.corticalData = meanOptoStimCortGamData_RH;
            AnalysisResults.(animalID).OptoStim.cortical_RH.(solenoid).Gam.corticalStD = stdOptoStimCortGamData_RH;

            AnalysisResults.(animalID).OptoStim.Pupil.(solenoid).Diameter.DiameterData = meanOptoStimPupilData;
            AnalysisResults.(animalID).OptoStim.Pupil.(solenoid).Diameter.DiameterStD = stdOptoStimPupilData;
            AnalysisResults.(animalID).OptoStim.Pupil.(solenoid).Diameter.DiameterRaw = procOptoStimPupilData;

            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).ECoG.corticalS_LH = msOptoStimCortS_Vals_LH;
            AnalysisResults.(animalID).OptoStim.cortical_RH.(solenoid).ECoG.corticalS_RH = msOptoStimCortS_Vals_RH;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).ECoG.T = T2;
            AnalysisResults.(animalID).OptoStim.cortical_LH.(solenoid).ECoG.F = F;
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

