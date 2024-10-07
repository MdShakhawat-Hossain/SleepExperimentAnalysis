function [AnalysisResults] = AnalyzeEvokedResponses_GRABNE_2Photon(animalID,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by MD Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the stimulus-evoked and whisker-evoked responses (2 photon)
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
    samplingRate = EventData.GFP.Z_NE.whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = EventData.GFP.Z_NE.whisk.trialDuration_sec;
    timeVector = (0:(EventData.GFP.Z_NE.whisk.epoch.duration*samplingRate))/samplingRate - EventData.GFP.Z_NE.whisk.epoch.offset;
    offset = EventData.GFP.Z_NE.whisk.epoch.offset;
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
        [whiskLogical] = FilterEvents_FP(EventData.GFP.Z_NE.whisk,WhiskCriteria);
        combWhiskLogical = logical(whiskLogical);
        [NEWhiskGFPData] = EventData.GFP.Z_NE.whisk.NormData(combWhiskLogical,:);

        [allWhiskCorticalMUAData] = EventData.cortical_LH.corticalPower.whisk.NormData(combWhiskLogical,:);
%         [allWhiskHippocampalMUAData] = EventData.hippocampus.corticalPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskCorticalGamData] = EventData.cortical_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
%         [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
        [allWhiskFileIDs] = EventData.GFP.Z_NE.whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.GFP.Z_NE.whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.GFP.Z_NE.whisk.duration(combWhiskLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
               
        [NEfinalWhiskGFPData,~,~,~] = RemoveInvalidData_IOS(NEWhiskGFPData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);

        [finalWhiskCorticalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
%         [finalWhiskHippocampalMUAData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        [finalWhiskCorticalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
%         [finalWhiskHippocampalGamData,~,~,~] = RemoveInvalidData_IOS(allWhiskHippocampalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskGFPData procWhiskRhodamineData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
        dd = 1;
        for cc = 1:size(NEfinalWhiskRhodamineData,1)
            whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 0;%2;
            whiskEndTime = whiskStartTime + 15;
            finalWhiskFileID = finalWhiskFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                NEwhiskGFParray = NEfinalWhiskGFPData(cc,:);

                whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
%                 whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
%                 whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                
                NEfiltWhiskGFParray = sgolayfilt(NEwhiskGFParray,3,9) - mean(NEwhiskGFParray(1:(offset*samplingRate)));
                NEprocWhiskGFPData(dd,:) = NEfiltWhiskGFParray;
                procWhiskCorticalMUAData(dd,:) = whiskCorticalMUAarray - mean(whiskCorticalMUAarray(1:(offset*samplingRate)));
%                 procWhiskHippocampalMUAData(dd,:) = whiskHippocampalMUAarray - mean(whiskHippocampalMUAarray(1:(offset*samplingRate)));
                procWhiskCorticalGamData(dd,:) = whiskCorticalGamArray - mean(whiskCorticalGamArray(1:(offset*samplingRate)));
%                 procWhiskHippocampalGamData(dd,:) = whiskHippocampalGamArray - mean(whiskHippocampalGamArray(1:(offset*samplingRate)));
                finalWhiskStartTimes(dd,1) = whiskStartTime;
                finalWhiskEndTimes(dd,1) = whiskEndTime;
                finalWhiskFiles{dd,1} = finalWhiskFileID;
                dd = dd + 1;
            end
        end

        NEmeanWhiskGFPData = mean(NEprocWhiskGFPData,1);%*100;
        NEstdWhiskGFPData = std(NEprocWhiskGFPData,0,1);%*100;

        meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
        stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
%         meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
%         stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
        meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
        stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
%         meanWhiskHippocampalGamData = mean(procWhiskHippocampalGamData,1)*100;
%         stdWhiskHippocampalGamData = std(procWhiskHippocampalGamData,0,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
%         whiskHippocampalZhold = [];
        for ee = 1:length(finalWhiskFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalWhiskFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = 'cortical_LH';
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
%                     whiskHippocampalS_Data = AllSpecData.hippocampus.normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
%             whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskCorticalZhold = cat(3,whiskCorticalZhold,whiskCorticalS_Vals);
%             whiskHippocampalZhold = cat(3,whiskHippocampalZhold,whiskHippocampalS_Vals);
        end
        % cortical mean-subtract by first 2 seconds prior to stimulus
        meanWhiskCortS = mean(whiskCorticalZhold,3);
        baseWhiskCortS_Vals = mean(meanWhiskCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixWhiskCortS_Vals = baseWhiskCortS_Vals.*ones(size(meanWhiskCortS));
        msStimWhiskS_Vals = (meanWhiskCortS - baseMatrixWhiskCortS_Vals);
        % hippocampal mean-subtract by first 2 seconds prior to stimulus
%         meanWhiskHipS = mean(whiskHippocampalZhold,3);
%         baseWhiskHipS_Vals = mean(meanWhiskHipS(:,1:1.5*specSamplingRate),2);
%         baseMatrixWhiskHipS_Vals = baseWhiskHipS_Vals.*ones(size(meanWhiskHipS));
%         msWhiskHipS_Vals = (meanWhiskHipS - baseMatrixWhiskHipS_Vals);
        T2 = -20:(1/specSamplingRate):20;
        % save results
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).GFP.GFP = NEmeanWhiskGFPData;
        AnalysisResults.(animalID).Whisk.Z_NE.(whiskCriteriaName).GFP.GFPStD = NEstdWhiskGFPData;

        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.hippocampalData = meanWhiskHippocampalGamData;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).Gam.hippocampalStD = stdWhiskHippocampalGamData;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.corticalS = msStimWhiskS_Vals;
%         AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.hippocampalS = msWhiskHipS_Vals;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.T = T2;
        AnalysisResults.(animalID).Whisk.cortical.(whiskCriteriaName).LFP.F = F;
    end
    %% analyze stimulus-evoked responses
    % if firstHrs == "true"
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
            allStimFilter = FilterEvents_IOS(EventData.GFP.Z_NE.stim,StimCriteria);
            [NEallStimGFPData] = EventData.GFP.Z_NE.stim.NormData(allStimFilter,:); 
    
            [allStimCortMUAData] = EventData.cortical_LH.corticalPower.stim.NormData(allStimFilter,:);
%             [allStimHipMUAData] = EventData.hippocampus.corticalPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData] = EventData.cortical_LH.gammaBandPower.stim.NormData(allStimFilter,:);
%             [allStimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.Rhodamine.Z_NE.stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.Rhodamine.Z_NE.stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % keep only the data that occurs within the manually-approved awake regions
            [NEfinalStimGFPData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS(NEallStimGFPData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);   
            [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
%             [finalStimHipMUAData,~,~,~] = RemoveInvalidData_IOS(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
%             [finalStimHipGamData,~,~,~] = RemoveInvalidData_IOS(allStimHipGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % lowpass filter each stim event and mean-subtract by the first 2 seconds
            clear procStimGFPData procStimRhodamineData procStimRhodamineData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
            ii = 1;
            for hh = 1:size(NEfinalStimGFPData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 0;%2;
                stimEndTime = stimStartTime + 15;%15;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    NEstimGFParray = NEfinalStimGFPData(hh,:);
    
                    stimCortMUAarray = finalStimCortMUAData(hh,:);
%                     stimHipMUAarray = finalStimHipMUAData(hh,:);
                    stimCortGamArray = finalStimCortGamData(hh,:);
%                     stimHipGamArray = finalStimHipGamData(hh,:);
    
                    NEfiltStimGFParray = sgolayfilt(NEstimGFParray,3,9);
    
                    NEprocStimGFPData(hh,:) = NEfiltStimGFParray - mean(NEfiltStimGFParray(1:(offset*samplingRate)));
    
                    procStimCortMUAData(hh,:) = stimCortMUAarray - mean(stimCortMUAarray(1:(offset*samplingRate)));
%                     procStimHipMUAData(hh,:) = stimHipMUAarray - mean(stimHipMUAarray(1:(offset*samplingRate)));
                    procStimCortGamData(hh,:) = stimCortGamArray - mean(stimCortGamArray(1:(offset*samplingRate)));
%                     procStimHipGamData(hh,:) = stimHipGamArray - mean(stimHipGamArray(1:(offset*samplingRate)));
                    finalStimStartTimes(ii,1) = stimStartTime;
                    finalStimEndTimes(ii,1) = stimEndTime;
                    finalStimFiles{ii,1} = finalStimFileID;
                    ii = ii + 1;
                end
            end

            NEmeanStimGFPData = mean(NEprocStimGFPData,1);%*100;
            NEstdStimGFPData = std(NEprocStimGFPData,0,1);%*100;
    
            meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
            stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
%             meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
%             stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
            meanStimCortGamData = mean(procStimCortGamData,1)*100;
            stdStimCortGamData = std(procStimCortGamData,0,1)*100;
%             meanStimHipGamData = mean(procStimHipGamData,1)*100;
%             stdStimHipGamData = std(procStimHipGamData,0,1)*100;
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
%                         stimHippocampalS_Data = AllSpecData.hippocampus.normS{kk,1};
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCortS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
%                 stimHipS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimCortZhold = cat(3,stimCortZhold,stimCortS_Vals);
%                 stimHipZhold = cat(3,stimHipZhold,stimHipS_Vals);
            end
            % cortical mean-subtract by first 2 seconds prior to stimulus
            meanStimCortS = mean(stimCortZhold,3);
            baseStimCortS_Vals = mean(meanStimCortS(:,1:1.5*specSamplingRate),2);
            baseMatrixStimCortS_Vals = baseStimCortS_Vals.*ones(size(meanStimCortS));
            msStimCortS_Vals = (meanStimCortS - baseMatrixStimCortS_Vals);
            % hippocampal mean-subtract by first 2 seconds prior to stimulus
%             meanStimHipS = mean(stimHipZhold,3);
%             baseStimHipS_Vals = mean(meanStimHipS(:,1:1.5*specSamplingRate),2);
%             baseMatrixStimHipS_Vals = baseStimHipS_Vals.*ones(size(meanStimHipS));
%             msStimHipS_Vals = (meanStimHipS - baseMatrixStimHipS_Vals);
            % save results
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).count = size(procStimCortMUAData,1);
    
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).GFP.GFP= NEmeanStimGFPData;
            AnalysisResults.(animalID).Stim.Z_NE.(solenoid).GFP.GFPStD = NEstdStimGFPData;
    
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.corticalStD = stdStimCortMUAData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.hippocampalData = meanStimHipMUAData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalData = meanStimCortGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.corticalStD = stdStimCortGamData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.hippocampalData = meanStimHipGamData;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).Gam.hippocampalStD = stdStimHipGamData;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.corticalS = msStimCortS_Vals;
%             AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.hippocampalS = msStimHipS_Vals;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).Stim.cortical.(solenoid).LFP.F = F;
        end
    % end
 % save data
 save('AnalysisResults.mat','AnalysisResults','-v7.3')

end

