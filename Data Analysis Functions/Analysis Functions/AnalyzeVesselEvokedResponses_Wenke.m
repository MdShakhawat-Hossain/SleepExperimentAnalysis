function [AnalysisResults_Wenke] = AnalyzeVesselEvokedResponses_Wenke(animalID,group,rootFolder,AnalysisResults_Wenke)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the whisking-evoked arteriole D/D responses (2PLSM)
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
dataLocation = [rootFolder '\' group '\' animalID '\2P Data\'];
cd(dataLocation)
% find and load EventData.mat struct
eventDataFileStruct = dir('*_EventData.mat');
eventDataFile = {eventDataFileStruct.name}';
eventDataFileID = char(eventDataFile);
load(eventDataFileID)
% find and load manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)
% find and load RestingBaselines.mat strut
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% criteria for whisking
whiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
whiskCriteriaA.Comparison = {'gt','lt','gt'};
whiskCriteriaA.Value = {0.5,2,5};
whiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
whiskCriteriaB.Comparison = {'gt','lt','gt'};
whiskCriteriaB.Value = {2,5,5};
whiskCriteriaC.Fieldname = {'duration','puffDistance'};
whiskCriteriaC.Comparison = {'gt','gt'};
whiskCriteriaC.Value = {5,5};
whiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
%% analyze whisking-evoked responses
for qq = 1:length(whiskCriteriaNames)
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = 5;
    offset = EventData.vesselDiameter.data.whisk.epoch.offset;
    timeVector = (0:(EventData.vesselDiameter.data.whisk.epoch.duration*samplingRate))/samplingRate - EventData.vesselDiameter.data.whisk.epoch.offset;
    whiskCriteriaName = whiskCriteriaNames{1,qq};
    if strcmp(whiskCriteriaName,'ShortWhisks') == true
        WhiskCriteria = whiskCriteriaA;
    elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
        WhiskCriteria = whiskCriteriaB;
    elseif strcmp(whiskCriteriaName,'LongWhisks') == true
        WhiskCriteria = whiskCriteriaC;
    end
    % pull data from EventData.mat structure
    [whiskLogical] = FilterEvents_2P(EventData.vesselDiameter.data.whisk,WhiskCriteria);
    whiskLogical = logical(whiskLogical);
    whiskingData = EventData.vesselDiameter.data.whisk.data(whiskLogical,:);
    whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(whiskLogical,:);
    whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(whiskLogical,:);
    whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(whiskLogical,:);
    whiskDurations = EventData.vesselDiameter.data.whisk.duration(whiskLogical,:);
    % keep only the data that occurs within the manually-approved awake regions
    [finalWhiskData,~,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P(whiskingData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    clear procWhiskData
    % filter and detrend data
    for aa = 1:size(finalWhiskData,1)
%         whiskStrDay = ConvertDate_2P(finalWhiskFileIDs{aa,1}(1:6));
        normWhiskData = (finalWhiskData(aa,:));% - RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay);
        filtWhiskData = sgolayfilt(normWhiskData,3,17);
        procWhiskData{aa,1} = filtWhiskData - mean(filtWhiskData(1:(offset*samplingRate))); %#ok<*AGROW>
    end
    % input data as time (1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedWhiskData = zeros(length(procWhiskData{1,1}),length(procWhiskData));
    for bb = 1:length(procWhiskData)
        reshapedWhiskData(:,bb) = procWhiskData{bb,1};
    end
    % remove veins from the artery list
    uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
    dd = 1;
    clear whiskArterioleIDs
    for cc = 1:length(uniqueWhiskVesselIDs)
        if strcmp(uniqueWhiskVesselIDs{cc,1}(1),'V') == false
            whiskArterioleIDs{dd,1} = uniqueWhiskVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    clear whiskArterioleEvoked
    for ee = 1:length(whiskArterioleIDs)
        whiskArterioleID = whiskArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalWhiskVesselIDs)
            if strcmp(whiskArterioleID,finalWhiskVesselIDs{ff,1}) == true
                whiskArterioleEvoked.(whiskArterioleID)(gg,:) = reshapedWhiskData(:,ff);
                gg = gg + 1;
            end
        end
    end
    % take mean/std of each arteriole's whisk-evoked response
    for aa = 1:length(whiskArterioleIDs)
        vID = whiskArterioleIDs{aa,1};
        meanWhiskEvokedDiam.(vID) = mean(whiskArterioleEvoked.(vID),1);
        stdWhiskEvokedDiam.(vID) = std(whiskArterioleEvoked.(vID),0,1);
        % save results
        AnalysisResults_Wenke.(group).(animalID).EvokedAvgs.(whiskCriteriaName).(vID).mean = meanWhiskEvokedDiam.(vID);
        AnalysisResults_Wenke.(group).(animalID).EvokedAvgs.(whiskCriteriaName).(vID).StD = stdWhiskEvokedDiam.(vID);
        AnalysisResults_Wenke.(group).(animalID).EvokedAvgs.(whiskCriteriaName).(vID).timeVector = timeVector;
    end
end
% save data
cd(rootFolder)
save('AnalysisResults_Wenke.mat','AnalysisResults_Wenke')
end
