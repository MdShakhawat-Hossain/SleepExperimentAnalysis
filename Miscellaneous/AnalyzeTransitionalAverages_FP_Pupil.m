function [AnalysisResults] = AnalyzeTransitionalAverages_FP_Pupil(animalID,saveFigs,rootFolder,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the transitions between different arousal-states (FP)
%________________________________________________________________________________________________________________________

%% only run analysis for valid animal IDs
% load model
% modelDirectory = [rootFolder '\' animalID '\CombinedImaging\Figures\Sleep Models\'];
% cd(modelDirectory)
% modelName = [animalID '_FP_RF_SleepScoringModel.mat'];
% load(modelName)
% go to data and load the model files
    if firstHrs == "false"
         dataLocation = [rootFolder '\' animalID '\CombinedImaging\'];
         transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    elseif firstHrs == "true"
        dataLocation = [rootFolder '\' animalID '\FirstHours\'];
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};%,'NREMtoREM','REMtoAWAKE'};

    end
cd(dataLocation)
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFile = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFile);

baselineFileStruct = dir('*_RestingBaselines.mat');
baselineFile = {baselineFileStruct.name}';
baselineFileID = char(baselineFile);
load(baselineFileID)
samplingRate = 30;
specSamplingRate = 10;
fileDates = fieldnames(RestingBaselines.manualSelection.cortical_LH.deltaBandPower);


    for a = 1:size(modelDataFileIDs,1)
        modelDataFileID = modelDataFileIDs(a,:);
        mdlIdx = strfind(modelDataFileID,'_');
        trainDataFileID = [modelDataFileID(1:mdlIdx(end)) 'TrainingData.mat'];
        if a == 1
            load(trainDataFileID)
            dataLength = size(trainingTable,1);
            joinedTable = trainingTable;
            joinedFileList = cell(size(trainingTable,1),1);
            joinedFileList(:) = {modelDataFileID};
        else
            load(trainDataFileID)
            fileIDCells = cell(size(trainingTable,1),1);
            fileIDCells(:) = {modelDataFileID};
            joinedTable = vertcat(joinedTable,trainingTable); %#ok<*AGROW>
            joinedFileList = vertcat(joinedFileList,fileIDCells);
        end
    end
    scoringTable = table2cell(joinedTable(:,end));
    labels = scoringTable;
    % go through each file and sleep score the data
% for a = 1:size(modelDataFileIDs,1)
%     modelDataFileID = modelDataFileIDs(a,:);
%     if a == 1
%         load(modelDataFileID)
%         dataLength = size(paramsTable,1);
%         joinedTable = paramsTable;
%         joinedFileList = cell(size(paramsTable,1),1);
%         joinedFileList(:) = {modelDataFileID};
%     else
%         load(modelDataFileID)
%         fileIDCells = cell(size(paramsTable,1),1);
%         fileIDCells(:) = {modelDataFileID};
%         joinedTable = vertcat(joinedTable,paramsTable); %#ok<*AGROW>
%         joinedFileList = vertcat(joinedFileList,fileIDCells);
%     end
% end
% scoringTable = joinedTable;


% [labels,~] = predict(RF_MDL,scoringTable);
% % apply a logical patch on the REM events
% REMindex = strcmp(labels,'REM Sleep');
% numFiles = length(labels)/dataLength;
% reshapedREMindex = reshape(REMindex,dataLength,numFiles);
% patchedREMindex = [];
% % patch missing REM indeces due to theta band falling off
% for b = 1:size(reshapedREMindex,2)
%     remArray = reshapedREMindex(:,b);
%     patchedREMarray = LinkBinaryEvents_FP(remArray',[5,0]);
%     patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
% end
% % change labels for each event
% for c = 1:length(labels)
%     if patchedREMindex(c,1) == 1
%         labels{c,1} = 'REM Sleep';
%     end
% end


% convert strings to numbers for easier comparisons
labelNumbers = zeros(length(labels),1);
for d = 1:length(labels)
    if strcmp(labels{d,1},'Not Sleep') == true
        labelNumbers(d,1) = 1;
    elseif strcmp(labels{d,1},'NREM Sleep') == true
        labelNumbers(d,1) = 2;
    elseif strcmp(labels{d,1},'REM Sleep') == true
        labelNumbers(d,1) = 3;
    end
end
% reshape
fileIDs = unique(joinedFileList);
fileLabels = reshape(labelNumbers,dataLength,size(modelDataFileIDs,1))';
stringArray = zeros(1,12);
for d = 1:length(transitions)
    transition = transitions{1,d};
    if strcmp(transition,'AWAKEtoNREM') == true
        for e = 1:12
            if e <= 6
                stringArray(1,e) = 1;
            else
                stringArray(1,e) = 2;
            end
        end
    elseif strcmp(transition,'NREMtoAWAKE') == true
        for e = 1:12
            if e <= 6
                stringArray(1,e) = 2;
            else
                stringArray(1,e) = 1;
            end
        end
    elseif strcmp(transition,'NREMtoREM') == true
        for e = 1:12
            if e <= 6
                stringArray(1,e) = 2;
            else
                stringArray(1,e) = 3;
            end
        end
    elseif strcmp(transition,'REMtoAWAKE') == true
        for e = 1:12
            if e <= 6
                stringArray(1,e) = 3;
            else
                stringArray(1,e) = 1;
            end
        end
    end
    % go through and pull out all indeces of matching strings
    idx = 1;
    for f = 1:length(fileIDs)
        fileID = fileIDs{f,1};
        labelArray = fileLabels(f,:);
        indeces = strfind(labelArray,stringArray);
        if isempty(indeces) == false
            for g1 = 1:length(indeces)
                data.(transition).files{idx,1} = fileID;
                data.(transition).startInd(idx,1) = indeces(1,g1);
                idx = idx + 1;
            end
        end
    end
end
% extract data
for h = 1:length(transitions)
    transition = transitions{1,h};
    iqx = 1;
    for i = 1:length(data.(transition).files)
        file = data.(transition).files{i,1};
        startBin = data.(transition).startInd(i,1);
        if startBin > 1 && startBin < (624- 12) % max bin - used before
            [animalID,fileDate,fileID] = GetFileInfo_FP(file);
            strDay = ConvertDate_FP(fileDate);
            procDataFileID = [animalID '_' fileID '_ProcData.mat'];
            load(procDataFileID)
            startTime = (startBin - 1)*5;   % sec
            endTime = startTime + (12*5);   % sec
            % whisking data
            [z1,p1,k1] = butter(4,10/(samplingRate/2),'low');
            [sos1,g1] = zp2sos(z1,p1,k1);
            filtWhiskAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle(startTime*samplingRate + 1:endTime*samplingRate));
            % motion  data
            filtforce = filtfilt(sos1,g1,ProcData.data.forceSensor(startTime*samplingRate + 1:endTime*samplingRate));
            binforce = (ProcData.data.binForceSensor(startTime*samplingRate + 1:endTime*samplingRate-1));
            binwhisk = (ProcData.data.binWhiskerAngle(startTime*samplingRate + 1:endTime*samplingRate-2));
            % EMG
            EMG = (ProcData.data.EMG.emg(startTime*samplingRate + 1:endTime*samplingRate)); %- RestingBaselines.manualSelection.EMG.emg.(strDay).mean))./RestingBaselines.manualSelection.EMG.emg.(strDay).mean;
            % Rhodamine data
            [z2,p2,k2] = butter(4,1/(samplingRate/2),'low');
            [sos2,g2] = zp2sos(z2,p2,k2);
            Ach_Rhodamine = ProcData.data.Rhodamine.Z_Ach;
            NE_Rhodamine = ProcData.data.Rhodamine.Z_NE;
 
            filtAch_Rhodamine = (filtfilt(sos2,g2,(Ach_Rhodamine(startTime*samplingRate + 1:endTime*samplingRate))));%- RestingBaselines.manualSelection.Rhodamine.Z_Ach.(strDay).mean);%./RestingBaselines.manualSelection.Rhodamine.Z_Ach.(strDay).std;
            filtNE_Rhodamine = (filtfilt(sos2,g2,(NE_Rhodamine(startTime*samplingRate + 1:endTime*samplingRate))));% - RestingBaselines.manualSelection.Rhodamine.Z_NE.(strDay).mean);%./RestingBaselines.manualSelection.Rhodamine.Z_NE.(strDay).std;
            % GFP data
            GRAB_ACh = ProcData.data.GFP.Z_Ach;
            GRAB_NE = ProcData.data.GFP.Z_NE;

            filtGRAB_ACh = (filtfilt(sos2,g2,(GRAB_ACh(startTime*samplingRate + 1:endTime*samplingRate))));% - RestingBaselines.manualSelection.GFP.Z_Ach.(strDay).mean);%./RestingBaselines.manualSelection.GFP.Z_Ach.(strDay).std;
            filtGRAB_NE = (filtfilt(sos2,g2,(GRAB_NE(startTime*samplingRate + 1:endTime*samplingRate))));% - RestingBaselines.manualSelection.GFP.Z_NE.(strDay).mean);%./RestingBaselines.manualSelection.GFP.Z_NE.(strDay).std;
            % Pupil data
            PupilDiameter = ProcData.data.Pupil.zDiameter;

            filtPupilDiameter = (filtfilt(sos2,g2,(PupilDiameter(startTime*samplingRate + 1:endTime*samplingRate))));% - RestingBaselines.manualSelection.GFP.Z_Ach.(strDay).mean);%./RestingBaselines.manualSelection.GFP.Z_Ach.(strDay).std;

            data.(transition).fileDate{iqx,1} = strDay;
            data.(transition).whisk(iqx,:) = filtWhiskAngle;
            data.(transition).binforce(iqx,:) = binforce;
            data.(transition).binwhisk(iqx,:) = binwhisk;
            data.(transition).force(iqx,:) = filtforce;
            data.(transition).EMG(iqx,:) = EMG;
            data.(transition).Ach_Rhodamine(iqx,:) = filtAch_Rhodamine;
            data.(transition).NE_Rhodamine(iqx,:) = filtNE_Rhodamine;
            data.(transition).GRAB_ACh(iqx,:) = filtGRAB_ACh;
            data.(transition).GRAB_NE(iqx,:) = filtGRAB_NE;
            data.(transition).PupilDiameter(iqx,:) = filtPupilDiameter;

            iqx = iqx + 1;
        end
    end
end
%% take averages of each behavior
movecond = {'move_0s','move_1s','move_2s','nomove_0s','nomove_1s','nomove_2s'};
for d = 1:length(transitions)
    transition = transitions{1,d};
    % find if there was movement in the transitions and separate those
    % movements into separate data
    % forcemovement = sum(data.(transition).binforce,2);
    % arrayIdx = 1:length(forcemovement);
    % forceidx_0s = find(forcemovement==0);  % no movement
    % forceidx_1s = find(forcemovement>30); % More than 1s movement
    % forceidx_2s = find(forcemovement>60); % More than 2s movement
    % whiskmovement = sum(data.(transition).binwhisk,2);
    % whiskidx_1s = find(whiskmovement>30); % More than 1s movement
    % whiskidx_2s = find(whiskmovement>60); % More than 2s movement
    % whiskidx_0s = find(whiskmovement==0); % no movement
    % data.(transition).nomove_0s = union(forceidx_0s,whiskidx_0s);
    % data.(transition).move_1s = union(forceidx_1s,whiskidx_1s);
    % data.(transition).move_2s = union(forceidx_2s,whiskidx_2s);
    % data.(transition).move_0s = setdiff(arrayIdx,data.(transition).nomove_0s);
    % data.(transition).nomove_1s = setdiff(arrayIdx,data.(transition).move_1s);
    % data.(transition).nomove_2s = setdiff(arrayIdx,data.(transition).move_2s); 
    % %% lets get to the separation
    % for mvn =  1:length(movecond)
    %     movementname = movecond{1,mvn};
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).Idx = data.(transition).(movementname);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).whisk = nanmean(data.(transition).whisk(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).force = nanmean(data.(transition).force(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).EMG = nanmean(data.(transition).EMG(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).T = data.(transition).T_short;
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).F = data.(transition).F;
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).indFileDate = data.(transition).fileDate(data.(transition).(movementname));
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).fileDates = fileDates;
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).LH_Cort = nanmean(data.(transition).LH_Cort(:,:,data.(transition).(movementname)),3);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).Ach_Rhodamine = nanmean(data.(transition).Ach_Rhodamine(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).NE_Rhodamine = nanmean(data.(transition).NE_Rhodamine(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).GRAB_ACh = nanmean(data.(transition).GRAB_ACh(data.(transition).(movementname),:),1);
    %     AnalysisResults.(animalID).Transitions.(transition).(movementname).GRAB_NE = nanmean(data.(transition).GRAB_NE(data.(transition).(movementname),:),1);
    % end
    %% this is the whole results
    AnalysisResults.(animalID).Transitions.(transition).whisk = mean(data.(transition).whisk,1);
    AnalysisResults.(animalID).Transitions.(transition).force = mean(data.(transition).force,1);
    AnalysisResults.(animalID).Transitions.(transition).whisk_std = std(data.(transition).whisk,1);
    AnalysisResults.(animalID).Transitions.(transition).force_std = std(data.(transition).force,1);
    
    AnalysisResults.(animalID).Transitions.(transition).EMG = mean(data.(transition).EMG,1);
    AnalysisResults.(animalID).Transitions.(transition).EMG_std = std(data.(transition).EMG,1);

    AnalysisResults.(animalID).Transitions.(transition).PupilDiameter = mean(data.(transition).PupilDiameter,1);
    AnalysisResults.(animalID).Transitions.(transition).PupilDiameter_std = std(data.(transition).PupilDiameter,1);
    AnalysisResults.(animalID).Transitions.(transition).PupilDiameterRaw = data.(transition).PupilDiameter;


    AnalysisResults.(animalID).Transitions.(transition).Ach_Rhodamine = mean(data.(transition).Ach_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).NE_Rhodamine = mean(data.(transition).NE_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).GRAB_ACh = mean(data.(transition).GRAB_ACh,1);
    AnalysisResults.(animalID).Transitions.(transition).GRAB_NE = mean(data.(transition).GRAB_NE,1);

    AnalysisResults.(animalID).Transitions.(transition).Ach_RhodamineRaw = data.(transition).Ach_Rhodamine;
    AnalysisResults.(animalID).Transitions.(transition).NE_RhodamineRaw = data.(transition).NE_Rhodamine;
    AnalysisResults.(animalID).Transitions.(transition).GRAB_AChRaw = data.(transition).GRAB_ACh;
    AnalysisResults.(animalID).Transitions.(transition).GRAB_NERaw = data.(transition).GRAB_NE;

    AnalysisResults.(animalID).Transitions.(transition).Ach_Rhodamine_std = std(data.(transition).Ach_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).NE_Rhodamine_std = std(data.(transition).NE_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).GRAB_ACh_std = std(data.(transition).GRAB_ACh,1);
    AnalysisResults.(animalID).Transitions.(transition).GRAB_NE_std = std(data.(transition).GRAB_NE,1);
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
