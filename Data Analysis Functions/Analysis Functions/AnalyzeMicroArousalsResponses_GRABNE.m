function [AnalysisResults] = AnalyzeMicroArousalsResponses_GRABNE(animalID,rootFolder,AnalysisResults,firstHrs)
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
% find and load MicroArousalsData.data.mat struct
eventDataFileStruct = dir('*_MicroArousalsData.data.mat');
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
ArousalsCriteriaA.Fieldname = {'duration','duration','puffDistance'};
ArousalsCriteriaA.Comparison = {'gt','lt','gt'};
ArousalsCriteriaA.Value = {0.5,2,5};
ArousalsCriteriaB.Fieldname = {'duration','duration','puffDistance'};
ArousalsCriteriaB.Comparison = {'gt','lt','gt'};
ArousalsCriteriaB.Value = {2,5,5};
ArousalsCriteriaC.Fieldname = {'duration','puffDistance'};
ArousalsCriteriaC.Comparison = {'gt','gt'};
ArousalsCriteriaC.Value = {5,5};
ArousalsCriteriaNames = {'ShortArousalss','IntermediateArousalss','LongArousalss'};
% criteria for stimulation
% StimCriteriaA.Value = {'LPadSol'};
% StimCriteriaA.Fieldname = {'solenoidName'};
% StimCriteriaA.Comparison = {'equal'};
% StimCriteriaB.Value = {'RPadSol'};
% StimCriteriaB.Fieldname = {'solenoidName'};
% StimCriteriaB.Comparison = {'equal'};
% StimCriteriaC.Value = {'AudSol'};
% StimCriteriaC.Fieldname = {'solenoidName'};
% StimCriteriaC.Comparison = {'equal'};
% stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
%% analyze whisking-evoked responses
    % pull a few necessary numbers from the MicroArousalsData.data.mat struct such as trial duration and sampling rate
    samplingRate = MicroArousalsData.data.Rhodamine.Z_NE.whisk.samplingRate;
    specSamplingRate = 10;
    trialDuration_sec = MicroArousalsData.data.Rhodamine.Z_NE.whisk.trialDuration_sec;
    timeVector = (0:(MicroArousalsData.data.Rhodamine.Z_NE.whisk.epoch.duration*samplingRate))/samplingRate - MicroArousalsData.data.Rhodamine.Z_NE.whisk.epoch.offset;
    offset = MicroArousalsData.data.Rhodamine.Z_NE.whisk.epoch.offset;
    for bb = 1:length(ArousalsCriteriaNames)
        whiskCriteriaName = ArousalsCriteriaNames{1,bb};
        if strcmp(whiskCriteriaName,'ShortArousalss') == true
            ArousalsCriteria = ArousalsCriteriaA;
        elseif strcmp(whiskCriteriaName,'IntermediateArousalss') == true
            ArousalsCriteria = ArousalsCriteriaB;
        elseif strcmp(whiskCriteriaName,'LongArousalss') == true
            ArousalsCriteria = ArousalsCriteriaC;
        end
        % pull data from MicroArousalsData.data.mat structure
        [whiskLogical] = FilterEvents_IOS(MicroArousalsData.data.Rhodamine.Z_NE.whisk,ArousalsCriteria);
        combArousalsLogical = logical(whiskLogical);
        [AchArousalsRhodamineData] = MicroArousalsData.data.Rhodamine.Z_Ach.whisk.NormData(combArousalsLogical,:);
        [NEArousalsRhodamineData] = MicroArousalsData.data.Rhodamine.Z_NE.whisk.NormData(combArousalsLogical,:);
        [AchArousalsGFPData] = MicroArousalsData.data.GFP.Z_Ach.whisk.NormData(combArousalsLogical,:);
        [NEArousalsGFPData] = MicroArousalsData.data.GFP.Z_NE.whisk.NormData(combArousalsLogical,:);

        [allArousalsCorticalMUAData] = MicroArousalsData.data.cortical_LH.corticalPower.whisk.NormData(combArousalsLogical,:);
        [allArousalsCorticalGamData] = MicroArousalsData.data.cortical_LH.gammaBandPower.whisk.NormData(combArousalsLogical,:);
        [allArousalsFileIDs] = MicroArousalsData.data.Rhodamine.Z_NE.whisk.fileIDs(combArousalsLogical,:);
        [allArousalsEventTimes] = MicroArousalsData.data.Rhodamine.Z_NE.whisk.eventTime(combArousalsLogical,:);
        allArousalsDurations = MicroArousalsData.data.Rhodamine.Z_NE.whisk.duration(combArousalsLogical,:);
        % keep only the data that occurs within the manually-approved awake regions
        [AchfinalArousalsRhodamineData,finalArousalsFileIDs,~,finalArousalsFileEventTimes] = RemoveInvalidData_IOS(AchArousalsRhodamineData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);
        [NEfinalArousalsRhodamineData,~,~,~] = RemoveInvalidData_IOS(NEArousalsRhodamineData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);
               
        [AchfinalArousalsGFPData,~,~,~] = RemoveInvalidData_IOS(AchArousalsGFPData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);
        [NEfinalArousalsGFPData,~,~,~] = RemoveInvalidData_IOS(NEArousalsGFPData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);

        [finalArousalsCorticalMUAData,~,~,~] = RemoveInvalidData_IOS(allArousalsCorticalMUAData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);
        [finalArousalsCorticalGamData,~,~,~] = RemoveInvalidData_IOS(allArousalsCorticalGamData,allArousalsFileIDs,allArousalsDurations,allArousalsEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procArousalsGFPData procArousalsRhodamineData procArousalsCorticalMUAData procArousalsHippocampalMUAData procArousalsCorticalGamData procArousalsHippocampalGamData finalArousalsStartTimes finalArousalsEndTimes finalArousalsFiles
        dd = 1;
        for cc = 1:size(NEfinalArousalsRhodamineData,1)
            whiskStartTime = round(finalArousalsFileEventTimes(cc,1),1) - 0;%2;
            whiskEndTime = whiskStartTime + 15;
            finalArousalsFileID = finalArousalsFileIDs{cc,1};
            if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                NEwhiskRhodaminearray = NEfinalArousalsRhodamineData(cc,:);
                AchwhiskRhodaminearray = AchfinalArousalsRhodamineData(cc,:);

                NEwhiskGFParray = NEfinalArousalsGFPData(cc,:);
                AchwhiskGFParray = AchfinalArousalsGFPData(cc,:);

                whiskCorticalMUAarray = finalArousalsCorticalMUAData(cc,:);
                whiskCorticalGamArray = finalArousalsCorticalGamData(cc,:);
                
                AchfiltArousalsRhodaminearray = sgolayfilt(AchwhiskRhodaminearray,3,17) - mean(AchwhiskRhodaminearray(1:(offset*samplingRate)));
                NEfiltArousalsRhodaminearray = sgolayfilt(NEwhiskRhodaminearray,3,17) - mean(NEwhiskRhodaminearray(1:(offset*samplingRate)));
                AchfiltArousalsGFParray = sgolayfilt(AchwhiskGFParray,3,17) - mean(AchwhiskGFParray(1:(offset*samplingRate)));
                NEfiltArousalsGFParray = sgolayfilt(NEwhiskGFParray,3,17) - mean(NEwhiskGFParray(1:(offset*samplingRate)));

                AchprocArousalsRhodamineData(dd,:) = AchfiltArousalsRhodaminearray;
                NEprocArousalsRhodamineData(dd,:) = NEfiltArousalsRhodaminearray;
                AchprocArousalsGFPData(dd,:) = AchfiltArousalsGFParray;
                NEprocArousalsGFPData(dd,:) = NEfiltArousalsGFParray;

                procArousalsCorticalMUAData(dd,:) = whiskCorticalMUAarray - mean(whiskCorticalMUAarray(1:(offset*samplingRate)));
                procArousalsCorticalGamData(dd,:) = whiskCorticalGamArray - mean(whiskCorticalGamArray(1:(offset*samplingRate)));
                finalArousalsStartTimes(dd,1) = whiskStartTime;
                finalArousalsEndTimes(dd,1) = whiskEndTime;
                finalArousalsFiles{dd,1} = finalArousalsFileID;
                dd = dd + 1;
            end
        end
        AchmeanArousalsRhodamineData = mean(AchprocArousalsRhodamineData,1);%*100;
        AchstdArousalsRhodamineData = std(AchprocArousalsRhodamineData,0,1);%*100;
        NEmeanArousalsRhodamineData = mean(NEprocArousalsRhodamineData,1);%*100;
        NEstdArousalsRhodamineData = std(NEprocArousalsRhodamineData,0,1);%*100;

        AchmeanArousalsGFPData = mean(AchprocArousalsGFPData,1);%*100;
        AchstdArousalsGFPData = std(AchprocArousalsGFPData,0,1);%*100;
        NEmeanArousalsGFPData = mean(NEprocArousalsGFPData,1);%*100;
        NEstdArousalsGFPData = std(NEprocArousalsGFPData,0,1);%*100;

        meanArousalsCorticalMUAData = mean(procArousalsCorticalMUAData,1)*100;
        stdArousalsCorticalMUAData = std(procArousalsCorticalMUAData,0,1)*100;
        meanArousalsCorticalGamData = mean(procArousalsCorticalGamData,1)*100;
        stdArousalsCorticalGamData = std(procArousalsCorticalGamData,0,1)*100;
        % extract LFP from spectrograms associated with the whisking indecies
        whiskCorticalZhold = [];
        for ee = 1:length(finalArousalsFiles)
            % load normalized one-second bin data from each file
            whiskFileID = finalArousalsFiles{ee,1};
            whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
            whiskSpecField = 'cortical_LH';
            for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                    whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                    F = AllSpecData.(whiskSpecField).F{ff,1};
                    T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                end
            end
            whiskStartTimeIndex = find(T == round(finalArousalsStartTimes(ee,1),1));
            whiskStartTimeIndex = whiskStartTimeIndex(1);
            whiskDurationIndex = find(T == round(finalArousalsEndTimes(ee,1),1));
            whiskDurationIndex = whiskDurationIndex(end);
            whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
            whiskCorticalZhold = cat(3,whiskCorticalZhold,whiskCorticalS_Vals);
        end
        % cortical mean-subtract by first 2 seconds prior to stimulus
        meanArousalsCortS = mean(whiskCorticalZhold,3);
        baseArousalsCortS_Vals = mean(meanArousalsCortS(:,1:1.5*specSamplingRate),2);
        baseMatrixArousalsCortS_Vals = baseArousalsCortS_Vals.*ones(size(meanArousalsCortS));
        msStimArousalsS_Vals = (meanArousalsCortS - baseMatrixArousalsCortS_Vals);
        T2 = -20:(1/specSamplingRate):20;
        % save results
        AnalysisResults.(animalID).MicroArousals.Z_Ach.(whiskCriteriaName).Rhodamine.Rhodamine = AchmeanArousalsRhodamineData;
        AnalysisResults.(animalID).MicroArousals.Z_Ach.(whiskCriteriaName).Rhodamine.RhodamineStD = AchstdArousalsRhodamineData;
        AnalysisResults.(animalID).MicroArousals.Z_NE.(whiskCriteriaName).Rhodamine.Rhodamine = NEmeanArousalsRhodamineData;
        AnalysisResults.(animalID).MicroArousals.Z_NE.(whiskCriteriaName).Rhodamine.RhodamineStD = NEstdArousalsRhodamineData;

        AnalysisResults.(animalID).MicroArousals.Z_Ach.(whiskCriteriaName).GFP.GFP = AchmeanArousalsGFPData;
        AnalysisResults.(animalID).MicroArousals.Z_Ach.(whiskCriteriaName).GFP.GFPStD = AchstdArousalsGFPData;
        AnalysisResults.(animalID).MicroArousals.Z_NE.(whiskCriteriaName).GFP.GFP = NEmeanArousalsGFPData;
        AnalysisResults.(animalID).MicroArousals.Z_NE.(whiskCriteriaName).GFP.GFPStD = NEstdArousalsGFPData;

        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).MUA.corticalData = meanArousalsCorticalMUAData;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).MUA.corticalStD = stdArousalsCorticalMUAData;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).Gam.corticalData = meanArousalsCorticalGamData;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).Gam.corticalStD = stdArousalsCorticalGamData;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).timeVector = timeVector;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).LFP.corticalS = msStimArousalsS_Vals;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).LFP.T = T2;
        AnalysisResults.(animalID).MicroArousals.cortical.(whiskCriteriaName).LFP.F = F;
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

