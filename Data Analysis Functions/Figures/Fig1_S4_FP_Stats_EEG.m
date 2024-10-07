function [AnalysisResults] = Fig1_S4_FP_Stats_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorWhisk = [(0/256),(166/256),(81/256)];
colorRest = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorRest = [(255/256),(191/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
FP_animalIDs = {'SHF025'};
behavFields = {'Rest','Whisk','Stim','NREM','REM'};
%% Rhodamine comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.IndRH)
            data.Rhodamine.(animalID).(behavField).meanRH(cc,1) = mean(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.IndRH{cc,1});
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.Rhodamine,behavField) == false
            data.Rhodamine.(behavField).catmeanRH = [];
            data.Rhodamine.(behavField).animalID = {};
            data.Rhodamine.(behavField).behavior = {};
            data.Rhodamine.(behavField).hemisphere = {};
        end
        data.Rhodamine.(behavField).catmeanRH = cat(1,data.Rhodamine.(behavField).catmeanRH,mean(data.Rhodamine.(animalID).(behavField).meanRH));
        data.Rhodamine.(behavField).animalID = cat(1,data.Rhodamine.(behavField).animalID,animalID,animalID);
        data.Rhodamine.(behavField).behavior = cat(1,data.Rhodamine.(behavField).behavior,behavField,behavField);
        data.Rhodamine.(behavField).hemisphere = cat(1,data.Rhodamine.(behavField).hemisphere,'RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.Rhodamine.(behavField).meanmeanRH = mean(data.Rhodamine.(behavField).catmeanRH,1);
    data.Rhodamine.(behavField).stdmeanRH = std(data.Rhodamine.(behavField).catmeanRH,0,1);
end
%% GFP comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndRH)
            data.GFP.(animalID).(behavField).meanRH(cc,1) = mean(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndRH{cc,1});
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GFP,behavField) == false
            data.GFP.(behavField).catmeanRH = [];
            data.GFP.(behavField).animalID = {};
            data.GFP.(behavField).behavior = {};
            data.GFP.(behavField).hemisphere = {};
        end
        data.GFP.(behavField).catmeanRH = cat(1,data.GFP.(behavField).catmeanRH,mean(data.GFP.(animalID).(behavField).meanRH));
        data.GFP.(behavField).animalID = cat(1,data.GFP.(behavField).animalID,animalID,animalID);
        data.GFP.(behavField).behavior = cat(1,data.GFP.(behavField).behavior,behavField,behavField);
        data.GFP.(behavField).hemisphere = cat(1,data.GFP.(behavField).hemisphere,'RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.GFP.(behavField).meanmeanRH = mean(data.GFP.(behavField).catmeanRH,1);
    data.GFP.(behavField).stdmeanRH = std(data.GFP.(behavField).catmeanRH,0,1);
end
%% statistics - generalized linear mixed-effects model
% Rhodamine
% tableSize = cat(1,data.Rhodamine.Rest.animalID,data.Rhodamine.Stim.animalID,data.Rhodamine.Whisk.animalID,data.Rhodamine.NREM.animalID);
% Rhodamine_meanTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Mean'});
% Rhodamine_meanTable.Mouse = cat(1,data.Rhodamine.Rest.animalID,data.Rhodamine.Stim.animalID,data.Rhodamine.Whisk.animalID,data.Rhodamine.NREM.animalID,data.Rhodamine.REM.animalID);
% Rhodamine_meanTable.Hemisphere = cat(1,data.Rhodamine.Rest.hemisphere,data.Rhodamine.Stim.hemisphere,data.Rhodamine.Whisk.hemisphere,data.Rhodamine.NREM.hemisphere,data.Rhodamine.REM.hemisphere);
% Rhodamine_meanTable.Behavior = cat(1,data.Rhodamine.Rest.behavior,data.Rhodamine.Stim.behavior,data.Rhodamine.Whisk.behavior,data.Rhodamine.NREM.behavior,data.Rhodamine.REM.behavior);
% Rhodamine_meanTable.Mean = cat(1,data.Rhodamine.Rest.catmean,data.Rhodamine.Stim.catmean,data.Rhodamine.Whisk.catmean,data.Rhodamine.NREM.catmean,data.Rhodamine.REM.catmean);
% Rhodamine_meanFitFormula = 'Mean ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% Rhodamine_meanStats = fitglme(Rhodamine_meanTable,Rhodamine_meanFitFormula)

% GFP
% tableSize = cat(1,data.GFP.Rest.animalID,data.GFP.Stim.animalID,data.GFP.Whisk.animalID,data.GFP.NREM.animalID);
% GFP_meanTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Mean'});
% GFP_meanTable.Mouse = cat(1,data.GFP.Rest.animalID,data.GFP.Stim.animalID,data.GFP.Whisk.animalID,data.GFP.NREM.animalID,data.GFP.REM.animalID);
% GFP_meanTable.Hemisphere = cat(1,data.GFP.Rest.hemisphere,data.GFP.Stim.hemisphere,data.GFP.Whisk.hemisphere,data.GFP.NREM.hemisphere,data.GFP.REM.hemisphere);
% GFP_meanTable.Behavior = cat(1,data.GFP.Rest.behavior,data.GFP.Stim.behavior,data.GFP.Whisk.behavior,data.GFP.NREM.behavior,data.GFP.REM.behavior);
% GFP_meanTable.Mean = cat(1,data.GFP.Rest.catmean,data.GFP.Stim.catmean,data.GFP.Whisk.catmean,data.GFP.NREM.catmean,data.GFP.REM.catmean);
% GFP_meanFitFormula = 'Mean ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% GFP_meanStats = fitglme(GFP_meanTable,GFP_meanFitFormula)

%% Fig. 1-S4
summaryFigure = figure('Name','Fig1-S4 Stats');
sgtitle('Mean hemodynamic changes')
%% Mean Rhodamine RH
ax1 = subplot(1,2,1);
xInds = ones(1,length(FP_animalIDs));
s1=scatter(xInds*2,data.Rhodamine.Whisk.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(2,data.Rhodamine.Whisk.meanmeanRH,data.Rhodamine.Whisk.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

s2=scatter(xInds*3,data.Rhodamine.Stim.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(3,data.Rhodamine.Stim.meanmeanRH,data.Rhodamine.Stim.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

s3=scatter(xInds*1,data.Rhodamine.Rest.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e3 = errorbar(1,data.Rhodamine.Rest.meanmeanRH,data.Rhodamine.Rest.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

s4=scatter(xInds*4,data.Rhodamine.NREM.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.Rhodamine.NREM.meanmeanRH,data.Rhodamine.NREM.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

s5=scatter(xInds*5,data.Rhodamine.REM.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Rhodamine.REM.meanmeanRH,data.Rhodamine.REM.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;

title({'Mean Zscored \DeltaRhodamine RH'})
ylabel('Mean  \DeltaRhodamine RH')
legend([s3,s1,s2,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];

%% Mean GFP RH
ax3 = subplot(1,2,2);
xInds = ones(1,length(FP_animalIDs));
scatter(xInds*2,data.GFP.Whisk.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(2,data.GFP.Whisk.meanmeanRH,data.GFP.Whisk.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

scatter(xInds*3,data.GFP.Stim.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(3,data.GFP.Stim.meanmeanRH,data.GFP.Stim.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

scatter(xInds*1,data.GFP.Rest.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e3 = errorbar(1,data.GFP.Rest.meanmeanRH,data.GFP.Rest.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

scatter(xInds*4,data.GFP.NREM.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.GFP.NREM.meanmeanRH,data.GFP.NREM.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

scatter(xInds*5,data.GFP.REM.catmeanRH,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.GFP.REM.meanmeanRH,data.GFP.REM.stdmeanRH,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;

title({'Mean Zscored \DeltaGFP RH'})
ylabel('Mean Zscored \DeltaGFP RH')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1-S4-Stats']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S4-Stats'])
    close 
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S4_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     Mean-to-Mean Rhodamine statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for Mean-to-Mean Rhodamine during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(Rhodamine_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [Rhodamine]: ' num2str(round(data.Rhodamine.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.Rhodamine.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [Rhodamine]: ' num2str(round(data.Rhodamine.NREM.meanP2P,1)) ' +/- ' num2str(round(data.Rhodamine.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [Rhodamine]: ' num2str(round(data.Rhodamine.REM.meanP2P,1)) ' +/- ' num2str(round(data.Rhodamine.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     Mean GFP statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for Mean GFP during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GFP_meanStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Mean [GFP]: ' num2str(round(data.GFP.Whisk.meanmean,1)) ' +/- ' num2str(round(data.GFP.Whisk.stdmean,1))]); disp(' ')
%     disp(['NREM  Mean [GFP]: ' num2str(round(data.GFP.NREM.meanmean,1)) ' +/- ' num2str(round(data.GFP.NREM.stdmean,1))]); disp(' ')
%     disp(['REM  Mean [GFP]: ' num2str(round(data.GFP.REM.meanmean,1)) ' +/- ' num2str(round(data.GFP.REM.stdmean,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
% 
%     Mean-to-Mean GFP statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for Mean-to-Mean GFP during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GFP_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [GFP]: ' num2str(round(data.GFP.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.GFP.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [GFP]: ' num2str(round(data.GFP.NREM.meanP2P,1)) ' +/- ' num2str(round(data.GFP.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [GFP]: ' num2str(round(data.GFP.REM.meanP2P,1)) ' +/- ' num2str(round(data.GFP.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     Mean GFP statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for Mean GFP during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GFP_meanStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Mean [GFP]: ' num2str(round(data.GFP.Whisk.meanmean,1)) ' +/- ' num2str(round(data.GFP.Whisk.stdmean,1))]); disp(' ')
%     disp(['NREM  Mean [GFP]: ' num2str(round(data.GFP.NREM.meanmean,1)) ' +/- ' num2str(round(data.GFP.NREM.stdmean,1))]); disp(' ')
%     disp(['REM  Mean [GFP]: ' num2str(round(data.GFP.REM.meanmean,1)) ' +/- ' num2str(round(data.GFP.REM.stdmean,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     
%     diary off
end

end
