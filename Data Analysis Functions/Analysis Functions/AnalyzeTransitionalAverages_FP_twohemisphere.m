function [AnalysisResults] = AnalyzeTransitionalAverages_FP_twohemisphere(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the transitions between different arousal-states (FP)
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs  = {'T282','T285'};%{'T281','T282','T284','T285'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% only run analysis for valid animal IDs
% load model
modelDirectory = [rootFolder '\' animalID '\CombinedImaging\Figures\Sleep Models\'];
cd(modelDirectory)
modelName = [animalID '_FP_RF_SleepScoringModel.mat'];
load(modelName)
% go to data and load the model files
dataLocation = [rootFolder '/' animalID '/CombinedImaging/'];
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
% go through each file and sleep score the data
for a = 1:size(modelDataFileIDs,1)
    modelDataFileID = modelDataFileIDs(a,:);
    if a == 1
        load(modelDataFileID)
        dataLength = size(paramsTable,1);
        joinedTable = paramsTable;
        joinedFileList = cell(size(paramsTable,1),1);
        joinedFileList(:) = {modelDataFileID};
    else
        load(modelDataFileID)
        fileIDCells = cell(size(paramsTable,1),1);
        fileIDCells(:) = {modelDataFileID};
        joinedTable = vertcat(joinedTable,paramsTable); %#ok<*AGROW>
        joinedFileList = vertcat(joinedFileList,fileIDCells);
    end
end
scoringTable = joinedTable;
[labels,~] = predict(RF_MDL,scoringTable);
% apply a logical patch on the REM events
REMindex = strcmp(labels,'REM Sleep');
numFiles = length(labels)/dataLength;
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for b = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,b);
    patchedREMarray = LinkBinaryEvents_FP(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for c = 1:length(labels)
    if patchedREMindex(c,1) == 1
        labels{c,1} = 'REM Sleep';
    end
end
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
            specDataFileID = [animalID '_' fileID '_SpecDataB.mat'];
            load(specDataFileID)
            startTime = (startBin - 1)*5;   % sec
            endTime = startTime + (12*5);   % sec
            % whisking data
            [z1,p1,k1] = butter(4,10/(samplingRate/2),'low');
            [sos1,g1] = zp2sos(z1,p1,k1);
            filtWhiskAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle(startTime*samplingRate + 1:endTime*samplingRate));
            % EMG
            EMG = (ProcData.data.EMG.emg(startTime*samplingRate + 1:endTime*samplingRate)); %- RestingBaselines.manualSelection.EMG.emg.(strDay).mean))./RestingBaselines.manualSelection.EMG.emg.(strDay).mean;
            % spectrogram data
            cortical_LHnormS = SpecData.cortical_LH.normS;%S;%
            cortical_RHnormS = SpecData.cortical_RH.normS;%S;%
            hippocampusNormS = SpecData.hippocampus.normS;%S;%
            T = round(SpecData.cortical_LH.T,1);
            F = SpecData.cortical_LH.F;
            specStartIndex = find(T == startTime);
            specStartIndex = specStartIndex(1);
            specEndIndex = find(T == endTime);
            specEndIndex = specEndIndex(end);
            LH_CortSpec = cortical_LHnormS(:,specStartIndex + 1:specEndIndex);
            RH_CortSpec = cortical_RHnormS(:,specStartIndex + 1:specEndIndex);
            Hip_spec = hippocampusNormS(:,specStartIndex + 1:specEndIndex);
            T_short = T(1:size(LH_CortSpec,2));
            % Rhodamine data
            [z2,p2,k2] = butter(4,1/(samplingRate/2),'low');
            [sos2,g2] = zp2sos(z2,p2,k2);
            LH_Rhodamine = ProcData.data.Rhodamine.LH;
            RH_Rhodamine = ProcData.data.Rhodamine.RH;
 
            filtLH_Rhodamine = filtfilt(sos2,g2,(LH_Rhodamine(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.Rhodamine.LH.(strDay).mean)./RestingBaselines.manualSelection.Rhodamine.LH.(strDay).std);
            filtRH_Rhodamine = filtfilt(sos2,g2,(RH_Rhodamine(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.Rhodamine.RH.(strDay).mean)./RestingBaselines.manualSelection.Rhodamine.RH.(strDay).std);
            % GCaMP7s data
            LH_GCaMP7s = ProcData.data.GCaMP7s.LH;
            RH_GCaMP7s = ProcData.data.GCaMP7s.RH;

            filtLH_GCaMP7s = filtfilt(sos2,g2,((LH_GCaMP7s(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.GCaMP7s.LH.(strDay).mean)./RestingBaselines.manualSelection.GCaMP7s.LH.(strDay).std));
            filtRH_GCaMP7s = filtfilt(sos2,g2,((RH_GCaMP7s(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.GCaMP7s.RH.(strDay).mean)./RestingBaselines.manualSelection.GCaMP7s.RH.(strDay).std));

            data.(transition).fileDate{iqx,1} = strDay;
            data.(transition).whisk(iqx,:) = filtWhiskAngle;
            data.(transition).EMG(iqx,:) = EMG;
            data.(transition).LH_Cort(:,:,iqx) = LH_CortSpec(:,1:specSamplingRate*60);
            data.(transition).RH_Cort(:,:,iqx) = RH_CortSpec(:,1:specSamplingRate*60);
            data.(transition).Hip(:,:,iqx) = Hip_spec(:,1:specSamplingRate*60);
            data.(transition).T_short = T_short(1:specSamplingRate*60);
            data.(transition).F = F;
            data.(transition).LH_Rhodamine(iqx,:) = filtLH_Rhodamine;
            data.(transition).RH_Rhodamine(iqx,:) = filtRH_Rhodamine;
            data.(transition).LH_GCaMP7s(iqx,:) = filtLH_GCaMP7s;
            data.(transition).RH_GCaMP7s(iqx,:) = filtRH_GCaMP7s;
            iqx = iqx + 1;
        end
    end
end
% take averages of each behavior
for d = 1:length(transitions)
    transition = transitions{1,d};
    % save results
    AnalysisResults.(animalID).Transitions.(transition).whisk = mean(data.(transition).whisk,1);
    AnalysisResults.(animalID).Transitions.(transition).EMG = mean(data.(transition).EMG,1);
    AnalysisResults.(animalID).Transitions.(transition).Hip = mean(data.(transition).Hip,3);
    AnalysisResults.(animalID).Transitions.(transition).T = data.(transition).T_short;
    AnalysisResults.(animalID).Transitions.(transition).F = data.(transition).F;
    AnalysisResults.(animalID).Transitions.(transition).indFileDate = data.(transition).fileDate;
    AnalysisResults.(animalID).Transitions.(transition).fileDates = fileDates;
    allCort = cat(3,data.(transition).LH_Cort,data.(transition).RH_Cort);
    allRhodamine = cat(1,data.(transition).LH_Rhodamine,data.(transition).RH_Rhodamine);
    allGCaMP7s = cat(1,data.(transition).LH_GCaMP7s,data.(transition).RH_GCaMP7s);
    AnalysisResults.(animalID).Transitions.(transition).Cort = mean(allCort,3);
    AnalysisResults.(animalID).Transitions.(transition).LH_Cort = mean(data.(transition).LH_Cort,3);
    AnalysisResults.(animalID).Transitions.(transition).RH_Cort = mean(data.(transition).RH_Cort,3);
    AnalysisResults.(animalID).Transitions.(transition).Rhodamine = mean(allRhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).LH_Rhodamine = mean(data.(transition).LH_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).RH_Rhodamine = mean(data.(transition).RH_Rhodamine,1);
    AnalysisResults.(animalID).Transitions.(transition).GCaMP7s = mean(allGCaMP7s,1);
    AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s = mean(data.(transition).LH_GCaMP7s,1);
    AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s = mean(data.(transition).RH_GCaMP7s,1);
end
% save figures if desired
if strcmp(saveFigs,'y') == true
    for e = 1:length(transitions)
        transition = transitions{1,e};
         % Plot transitions
        transitionFig = figure;
        sgtitle([animalID ' average ' transition])
        % behavior
        subplot(4,1,1)
        p1 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).Rhodamine))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).Rhodamine,'color','r','LineWidth',2);
        ylabel('Zscored \DeltaRhodamine')
        xlim([0,samplingRate*60])
        yyaxis right
        p2 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).EMG))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).EMG,'color','k','LineWidth',2);
        ylabel('EMG (Volts^2)')
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        axis tight
        legend([p1,p2],'Rhodamine','EMG')
        % GCaMP7s
        subplot(4,1,2)
        p1 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).GCaMP7s))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).GCaMP7s,'color','g','LineWidth',2);
        ylabel('Zscored \DeltaGCaMP7s')
        xlim([0,samplingRate*60])
        yyaxis right
        p2 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).EMG))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).EMG,'color','k','LineWidth',2);
        ylabel('EMG (Volts^2)')
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        axis tight
        legend([p1,p2],'GCaMP7s','EMG')
        % Cort neural
        subplot(4,1,3)
        Semilog_ImageSC(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).Cort,'y')
        axis xy
%          caxis([-1,2])
        ylabel('Frequency (Hz)')
        set(gca,'Yticklabel','10^1')
        xlim([0,60])
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        yyaxis right
        ylabel('cortical LFP')
        % hip
        subplot(4,1,4)
        Semilog_ImageSC(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).Hip,'y')
%         caxis([-1,2])
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)')
        set(gca,'Yticklabel','10^1')
        xlim([0,60])
        set(gca,'TickLength',[0,0])
        set(gca,'box','off')
        yyaxis right
        ylabel('Hippocampal LFP')
        % save figure
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Transitions/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(transitionFig,[dirpath animalID '_' transition '_Transition']);
        saveas(transitionFig,[dirpath animalID '_' transition '_Transition.tiff'],'tiff');
        close(transitionFig)
        % Plot bilateral transitions
        transitionFig_Bilateral = figure;
        sgtitle([animalID ' average ' transition])
        % behavior
        subplot(5,1,1)
        p1 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).LH_Rhodamine))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).LH_Rhodamine,'r-.','LineWidth',2);
        ylabel('Zscored \DeltaRhodamine')
        xlim([0,samplingRate*60])
        hold on
        p2 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).RH_Rhodamine))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).RH_Rhodamine,'color','r','LineWidth',2);
        yyaxis right
        p3 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).EMG))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).EMG,'color','k','LineWidth',2);
        ylabel('EMG (Volts^2)')
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        axis tight
        legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG')

        % GCaMP7s
        subplot(5,1,2)
        p1 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s,'g--','LineWidth',2);
        ylabel('Zscored \DeltaGCaMP7s')
        xlim([0,samplingRate*60])
        hold on
        p2 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s,'color','g','LineWidth',2);
        yyaxis right
        p3 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).EMG))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).EMG,'color','k','LineWidth',2);
        ylabel('EMG (Volts^2)')
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        axis tight
        legend([p1,p2, p3],'LH_GCaMP7s','RH_GCaMP7s','EMG')
        
        % Cort neural LH
        subplot(5,1,3)
        Semilog_ImageSC(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).LH_Cort,'y')
        axis xy
%         caxis([-1,2])
        ylabel('Frequency (Hz)')
        set(gca,'Yticklabel','10^1')
        xlim([0,60])
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        yyaxis right
        ylabel('LH cortical LFP')

        % Cort neural RH
        subplot(5,1,4)
        Semilog_ImageSC(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).RH_Cort,'y')
        axis xy
%         caxis([-1,2])
        ylabel('Frequency (Hz)')
        set(gca,'Yticklabel','10^1')
        xlim([0,60])
        set(gca,'TickLength',[0,0])
        set(gca,'Xticklabel',[])
        set(gca,'box','off')
        yyaxis right
        ylabel('RH cortical LFP')
        
        % hip
        subplot(5,1,5)
        Semilog_ImageSC(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).Hip,'y')
%         caxis([-1,2])
        xlabel('Time (sec)')
        ylabel('Frequency (Hz)')
        set(gca,'Yticklabel','10^1')
        xlim([0,60])
        set(gca,'TickLength',[0,0])
        set(gca,'box','off')
        yyaxis right
        ylabel('Hippocampal LFP')
        % save figure
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Transitions/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(transitionFig_Bilateral,[dirpath animalID '_' transition '_BilateralTransition']);
        saveas(transitionFig_Bilateral,[dirpath animalID '_' transition '_BilateralTransition.tiff'],'tiff');
        close(transitionFig_Bilateral)
    end
end
% save data
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')

end
