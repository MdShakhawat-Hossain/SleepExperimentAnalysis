function [AnalysisResults] = Fig1_S4_FP(rootFolder,saveFigs,delim,AnalysisResults)
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
FP_animalIDs = {'T282','T285'};%{'T281','T282','T284','T285'};
% TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
behavFields = {'Whisk','Stim','Rest','NREM','REM'};
%% TRITC comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:min(length(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH),length(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH))
            data.TRITC.(animalID).(behavField).meanLH(cc,1) = mean(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1});
            data.TRITC.(animalID).(behavField).meanRH(cc,1) = mean(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1});
%             data.TRITC.(animalID).(behavField).p2pLH(cc,1) = abs(mean(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1}));
%             data.TRITC.(animalID).(behavField).p2pRH(cc,1) = abs(mean(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1}));
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.TRITC,behavField) == false
            data.TRITC.(behavField).catmean = [];
%             data.TRITC.(behavField).catP2P = [];
            data.TRITC.(behavField).animalID = {};
            data.TRITC.(behavField).behavior = {};
            data.TRITC.(behavField).hemisphere = {};
        end
        data.TRITC.(behavField).catmean = cat(1,data.TRITC.(behavField).catmean,mean(data.TRITC.(animalID).(behavField).meanLH),mean(data.TRITC.(animalID).(behavField).meanRH));
%         data.TRITC.(behavField).catP2P = cat(1,data.TRITC.(behavField).catP2P,mean(data.TRITC.(animalID).(behavField).p2pLH),mean(data.TRITC.(animalID).(behavField).p2pRH));
        data.TRITC.(behavField).animalID = cat(1,data.TRITC.(behavField).animalID,animalID,animalID);
        data.TRITC.(behavField).behavior = cat(1,data.TRITC.(behavField).behavior,behavField,behavField);
        data.TRITC.(behavField).hemisphere = cat(1,data.TRITC.(behavField).hemisphere,'LH','RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.TRITC.(behavField).meanmean = mean(data.TRITC.(behavField).catmean,1);
    data.TRITC.(behavField).stdmean = std(data.TRITC.(behavField).catmean,0,1);
%     data.TRITC.(behavField).meanP2P = mean(data.TRITC.(behavField).catP2P,1);
%     data.TRITC.(behavField).stdP2P = std(data.TRITC.(behavField).catP2P,0,1);
end
%% GCaMP7s comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:min(length(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH),length(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH))
            data.GCaMP7s.(animalID).(behavField).meanLH(cc,1) = mean(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1});
            data.GCaMP7s.(animalID).(behavField).meanRH(cc,1) = mean(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1});
%             data.GCaMP7s.(animalID).(behavField).p2pLH(cc,1) = abs(mean(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1}));
%             data.GCaMP7s.(animalID).(behavField).p2pRH(cc,1) = abs(mean(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1}));
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GCaMP7s,behavField) == false
            data.GCaMP7s.(behavField).catmean = [];
%             data.GCaMP7s.(behavField).catP2P = [];
            data.GCaMP7s.(behavField).animalID = {};
            data.GCaMP7s.(behavField).behavior = {};
            data.GCaMP7s.(behavField).hemisphere = {};
        end
        data.GCaMP7s.(behavField).catmean = cat(1,data.GCaMP7s.(behavField).catmean,mean(data.GCaMP7s.(animalID).(behavField).meanLH),mean(data.GCaMP7s.(animalID).(behavField).meanRH));
%         data.GCaMP7s.(behavField).catP2P = cat(1,data.GCaMP7s.(behavField).catP2P,mean(data.GCaMP7s.(animalID).(behavField).p2pLH),mean(data.GCaMP7s.(animalID).(behavField).p2pRH));
        data.GCaMP7s.(behavField).animalID = cat(1,data.GCaMP7s.(behavField).animalID,animalID,animalID);
        data.GCaMP7s.(behavField).behavior = cat(1,data.GCaMP7s.(behavField).behavior,behavField,behavField);
        data.GCaMP7s.(behavField).hemisphere = cat(1,data.GCaMP7s.(behavField).hemisphere,'LH','RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.GCaMP7s.(behavField).meanmean = mean(data.GCaMP7s.(behavField).catmean,1);
    data.GCaMP7s.(behavField).stdmean = std(data.GCaMP7s.(behavField).catmean,0,1);
%     data.GCaMP7s.(behavField).meanP2P = mean(data.GCaMP7s.(behavField).catP2P,1);
%     data.GCaMP7s.(behavField).stdP2P = std(data.GCaMP7s.(behavField).catP2P,0,1);
end
%% statistics - generalized linear mixed-effects model
% TRITC
tableSize = cat(1,data.TRITC.Rest.animalID,data.TRITC.Stim.animalID,data.TRITC.Whisk.animalID,data.TRITC.NREM.animalID,data.TRITC.REM.animalID);
TRITC_meanTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Mean'});
TRITC_meanTable.Mouse = cat(1,data.TRITC.Rest.animalID,data.TRITC.Stim.animalID,data.TRITC.Whisk.animalID,data.TRITC.NREM.animalID,data.TRITC.REM.animalID);
TRITC_meanTable.Hemisphere = cat(1,data.TRITC.Rest.hemisphere,data.TRITC.Stim.hemisphere,data.TRITC.Whisk.hemisphere,data.TRITC.NREM.hemisphere,data.TRITC.REM.hemisphere);
TRITC_meanTable.Behavior = cat(1,data.TRITC.Rest.behavior,data.TRITC.Stim.behavior,data.TRITC.Whisk.behavior,data.TRITC.NREM.behavior,data.TRITC.REM.behavior);
TRITC_meanTable.Mean = cat(1,data.TRITC.Rest.catmean,data.TRITC.Stim.catmean,data.TRITC.Whisk.catmean,data.TRITC.NREM.catmean,data.TRITC.REM.catmean);
TRITC_meanFitFormula = 'Mean ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
TRITC_meanStats = fitglme(TRITC_meanTable,TRITC_meanFitFormula)

% GCaMP7s
tableSize = cat(1,data.GCaMP7s.Rest.animalID,data.GCaMP7s.Stim.animalID,data.GCaMP7s.Whisk.animalID,data.GCaMP7s.NREM.animalID,data.GCaMP7s.REM.animalID);
GCaMP7s_meanTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Mean'});
GCaMP7s_meanTable.Mouse = cat(1,data.GCaMP7s.Rest.animalID,data.GCaMP7s.Stim.animalID,data.GCaMP7s.Whisk.animalID,data.GCaMP7s.NREM.animalID,data.GCaMP7s.REM.animalID);
GCaMP7s_meanTable.Hemisphere = cat(1,data.GCaMP7s.Rest.hemisphere,data.GCaMP7s.Stim.hemisphere,data.GCaMP7s.Whisk.hemisphere,data.GCaMP7s.NREM.hemisphere,data.GCaMP7s.REM.hemisphere);
GCaMP7s_meanTable.Behavior = cat(1,data.GCaMP7s.Rest.behavior,data.GCaMP7s.Stim.behavior,data.GCaMP7s.Whisk.behavior,data.GCaMP7s.NREM.behavior,data.GCaMP7s.REM.behavior);
GCaMP7s_meanTable.Mean = cat(1,data.GCaMP7s.Rest.catmean,data.GCaMP7s.Stim.catmean,data.GCaMP7s.Whisk.catmean,data.GCaMP7s.NREM.catmean,data.GCaMP7s.REM.catmean);
GCaMP7s_meanFitFormula = 'Mean ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
GCaMP7s_meanStats = fitglme(GCaMP7s_meanTable,GCaMP7s_meanFitFormula)

%% Fig. 1-S4
summaryFigure = figure('Name','Fig1-S4');
sgtitle('Mean hemodynamic changes')
%% Mean TRITC
ax2 = subplot(1,2,1);
xInds = ones(1,length(FP_animalIDs)*2);
s1=scatter(xInds*3,data.TRITC.Whisk.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(3,data.TRITC.Whisk.meanmean,data.TRITC.Whisk.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

s2=scatter(xInds*2,data.TRITC.Stim.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.TRITC.Stim.meanmean,data.TRITC.Stim.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

s3=scatter(xInds*1,data.TRITC.Rest.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e3 = errorbar(1,data.TRITC.Rest.meanmean,data.TRITC.Rest.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

s4=scatter(xInds*4,data.TRITC.NREM.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.TRITC.NREM.meanmean,data.TRITC.NREM.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

s5=scatter(xInds*5,data.TRITC.REM.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.TRITC.REM.meanmean,data.TRITC.REM.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;

title({'Mean Zscored \DeltaTRITC '})
ylabel('Mean  \DeltaTRITC ')
legend([s3,s2,s1,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% Mean GCaMP7s
ax4 = subplot(1,2,2);
xInds = ones(1,length(FP_animalIDs)*2);
scatter(xInds*3,data.GCaMP7s.Whisk.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(3,data.GCaMP7s.Whisk.meanmean,data.GCaMP7s.Whisk.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

scatter(xInds*2,data.GCaMP7s.Stim.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.GCaMP7s.Stim.meanmean,data.GCaMP7s.Stim.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

scatter(xInds*1,data.GCaMP7s.Rest.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e3 = errorbar(1,data.GCaMP7s.Rest.meanmean,data.GCaMP7s.Rest.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

scatter(xInds*4,data.GCaMP7s.NREM.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.GCaMP7s.NREM.meanmean,data.GCaMP7s.NREM.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

scatter(xInds*5,data.GCaMP7s.REM.catmean,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.GCaMP7s.REM.meanmean,data.GCaMP7s.REM.stdmean,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;


title({'Mean Zscored \DeltaGCaMP7s '})
ylabel('Mean Zscored \DeltaGCaMP7s ')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1-S4']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S4'])
    close 
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S4_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     Mean-to-Mean TRITC statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for Mean-to-Mean TRITC during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(TRITC_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [TRITC]: ' num2str(round(data.TRITC.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [TRITC]: ' num2str(round(data.TRITC.NREM.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [TRITC]: ' num2str(round(data.TRITC.REM.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     Mean GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for Mean GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_meanStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdmean,1))]); disp(' ')
%     disp(['NREM  Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdmean,1))]); disp(' ')
%     disp(['REM  Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdmean,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
% 
%     Mean-to-Mean GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for Mean-to-Mean GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     Mean GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for Mean GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_meanStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdmean,1))]); disp(' ')
%     disp(['NREM  Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdmean,1))]); disp(' ')
%     disp(['REM  Mean [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanmean,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdmean,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     
%     diary off
end

end
