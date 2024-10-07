function [AnalysisResults] = Fig1_S4_FP_Sleep(rootFolder,saveFigs,delim,AnalysisResults)
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
colorRest = [(0/256),(166/256),(81/256)];
% colorRest = [(0/256),(256/256),(256/256)];
% colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorRest = [(255/256),(191/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
FP_animalIDs = {'T282','T284','T285'};
behavFields = {'Rest','NREM','REM'};
%% TRITC comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH)
            data.TRITC.(animalID).(behavField).maxLH(cc,1) = max(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1});
            data.TRITC.(animalID).(behavField).maxRH(cc,1) = max(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1});
            data.TRITC.(animalID).(behavField).p2pLH(cc,1) = abs(max(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndLH{cc,1}));
            data.TRITC.(animalID).(behavField).p2pRH(cc,1) = abs(max(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanTRITC.(behavField).TRITC.IndRH{cc,1}));
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.TRITC,behavField) == false
            data.TRITC.(behavField).catMax = [];
            data.TRITC.(behavField).catP2P = [];
            data.TRITC.(behavField).animalID = {};
            data.TRITC.(behavField).behavior = {};
            data.TRITC.(behavField).hemisphere = {};
        end
        data.TRITC.(behavField).catMax = cat(1,data.TRITC.(behavField).catMax,mean(data.TRITC.(animalID).(behavField).maxLH),mean(data.TRITC.(animalID).(behavField).maxRH));
        data.TRITC.(behavField).catP2P = cat(1,data.TRITC.(behavField).catP2P,mean(data.TRITC.(animalID).(behavField).p2pLH),mean(data.TRITC.(animalID).(behavField).p2pRH));
        data.TRITC.(behavField).animalID = cat(1,data.TRITC.(behavField).animalID,animalID,animalID);
        data.TRITC.(behavField).behavior = cat(1,data.TRITC.(behavField).behavior,behavField,behavField);
        data.TRITC.(behavField).hemisphere = cat(1,data.TRITC.(behavField).hemisphere,'LH','RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.TRITC.(behavField).meanMax = mean(data.TRITC.(behavField).catMax,1);
    data.TRITC.(behavField).stdMax = std(data.TRITC.(behavField).catMax,0,1);
    data.TRITC.(behavField).meanP2P = mean(data.TRITC.(behavField).catP2P,1);
    data.TRITC.(behavField).stdP2P = std(data.TRITC.(behavField).catP2P,0,1);
end
%% GCaMP7s comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH)
            data.GCaMP7s.(animalID).(behavField).maxLH(cc,1) = max(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1});
            data.GCaMP7s.(animalID).(behavField).maxRH(cc,1) = max(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1});
            data.GCaMP7s.(animalID).(behavField).p2pLH(cc,1) = abs(max(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndLH{cc,1}));
            data.GCaMP7s.(animalID).(behavField).p2pRH(cc,1) = abs(max(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanGCaMP7s.(behavField).GCaMP7s.IndRH{cc,1}));
        end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GCaMP7s,behavField) == false
            data.GCaMP7s.(behavField).catMax = [];
            data.GCaMP7s.(behavField).catP2P = [];
            data.GCaMP7s.(behavField).animalID = {};
            data.GCaMP7s.(behavField).behavior = {};
            data.GCaMP7s.(behavField).hemisphere = {};
        end
        data.GCaMP7s.(behavField).catMax = cat(1,data.GCaMP7s.(behavField).catMax,mean(data.GCaMP7s.(animalID).(behavField).maxLH),mean(data.GCaMP7s.(animalID).(behavField).maxRH));
        data.GCaMP7s.(behavField).catP2P = cat(1,data.GCaMP7s.(behavField).catP2P,mean(data.GCaMP7s.(animalID).(behavField).p2pLH),mean(data.GCaMP7s.(animalID).(behavField).p2pRH));
        data.GCaMP7s.(behavField).animalID = cat(1,data.GCaMP7s.(behavField).animalID,animalID,animalID);
        data.GCaMP7s.(behavField).behavior = cat(1,data.GCaMP7s.(behavField).behavior,behavField,behavField);
        data.GCaMP7s.(behavField).hemisphere = cat(1,data.GCaMP7s.(behavField).hemisphere,'LH','RH');
    end
end
% take mean/StD
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.GCaMP7s.(behavField).meanMax = mean(data.GCaMP7s.(behavField).catMax,1);
    data.GCaMP7s.(behavField).stdMax = std(data.GCaMP7s.(behavField).catMax,0,1);
    data.GCaMP7s.(behavField).meanP2P = mean(data.GCaMP7s.(behavField).catP2P,1);
    data.GCaMP7s.(behavField).stdP2P = std(data.GCaMP7s.(behavField).catP2P,0,1);
end
%% statistics - generalized linear mixed-effects model
% TRITC
tableSize = cat(1,data.TRITC.Rest.animalID,data.TRITC.NREM.animalID,data.TRITC.REM.animalID);
TRITC_maxTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Peak'});
TRITC_maxTable.Mouse = cat(1,data.TRITC.Rest.animalID,data.TRITC.NREM.animalID,data.TRITC.REM.animalID);
TRITC_maxTable.Hemisphere = cat(1,data.TRITC.Rest.hemisphere,data.TRITC.NREM.hemisphere,data.TRITC.REM.hemisphere);
TRITC_maxTable.Behavior = cat(1,data.TRITC.Rest.behavior,data.TRITC.NREM.behavior,data.TRITC.REM.behavior);
TRITC_maxTable.Peak = cat(1,data.TRITC.Rest.catMax,data.TRITC.NREM.catMax,data.TRITC.REM.catMax);
TRITC_maxFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
TRITC_maxStats = fitglme(TRITC_maxTable,TRITC_maxFitFormula)

% GCaMP7s
tableSize = cat(1,data.GCaMP7s.Rest.animalID,data.GCaMP7s.NREM.animalID,data.GCaMP7s.REM.animalID);
GCaMP7s_maxTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Peak'});
GCaMP7s_maxTable.Mouse = cat(1,data.GCaMP7s.Rest.animalID,data.GCaMP7s.NREM.animalID,data.GCaMP7s.REM.animalID);
GCaMP7s_maxTable.Hemisphere = cat(1,data.GCaMP7s.Rest.hemisphere,data.GCaMP7s.NREM.hemisphere,data.GCaMP7s.REM.hemisphere);
GCaMP7s_maxTable.Behavior = cat(1,data.GCaMP7s.Rest.behavior,data.GCaMP7s.NREM.behavior,data.GCaMP7s.REM.behavior);
GCaMP7s_maxTable.Peak = cat(1,data.GCaMP7s.Rest.catMax,data.GCaMP7s.NREM.catMax,data.GCaMP7s.REM.catMax);
GCaMP7s_maxFitFormula = 'Peak ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
GCaMP7s_maxStats = fitglme(GCaMP7s_maxTable,GCaMP7s_maxFitFormula)

%% Fig. 1-S4
summaryFigure = figure('Name','Fig1-S4');
sgtitle('Peak hemodynamic changes')
%% peak TRITC
ax1 = subplot(1,2,1);
xInds = ones(1,length(FP_animalIDs)*2);

s1=scatter(xInds*1,data.TRITC.Rest.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e1 = errorbar(1,data.TRITC.Rest.meanMax,data.TRITC.Rest.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

s2=scatter(xInds*2,data.TRITC.NREM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.TRITC.NREM.meanMax,data.TRITC.NREM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

s3=scatter(xInds*3,data.TRITC.REM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.TRITC.REM.meanMax,data.TRITC.REM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

title({'Mean Peak Zscored \DeltaTRITC '})
ylabel('Peak  \DeltaTRITC ')
legend([s1,s2,s3],'Rest','NREM','REM','Location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% peak GCaMP7s
ax2 = subplot(1,2,2);
xInds = ones(1,length(FP_animalIDs)*2);
scatter(xInds*3,data.GCaMP7s.Whisk.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk,'jitter','on','jitterAmount',0.25);
hold on

scatter(xInds*1,data.GCaMP7s.Rest.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest,'jitter','on','jitterAmount',0.25);
e1 = errorbar(1,data.GCaMP7s.Rest.meanMax,data.GCaMP7s.Rest.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

scatter(xInds*2,data.GCaMP7s.NREM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.GCaMP7s.NREM.meanMax,data.GCaMP7s.NREM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;

scatter(xInds*3,data.GCaMP7s.REM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.GCaMP7s.REM.meanMax,data.GCaMP7s.REM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;


title({'Mean Peak Zscored \DeltaGCaMP7s '})
ylabel('Peak Zscored \DeltaGCaMP7s ')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
axis tight
xlim([0,length(behavFields) + 1])
% ylim([0,140])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];

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
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S4_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     peak-to-peak TRITC statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for peak-to-peak TRITC during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(TRITC_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [TRITC]: ' num2str(round(data.TRITC.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [TRITC]: ' num2str(round(data.TRITC.NREM.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [TRITC]: ' num2str(round(data.TRITC.REM.meanP2P,1)) ' +/- ' num2str(round(data.TRITC.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     peak GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for peak GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_maxStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdMax,1))]); disp(' ')
%     disp(['NREM  Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdMax,1))]); disp(' ')
%     disp(['REM  Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdMax,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
% 
%     peak-to-peak GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for peak-to-peak GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanP2P,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdP2P,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     peak GCaMP7s statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4b] Generalized linear mixed-effects model statistics for peak GCaMP7s during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(GCaMP7s_maxStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.Whisk.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.Whisk.stdMax,1))]); disp(' ')
%     disp(['NREM  Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.NREM.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.NREM.stdMax,1))]); disp(' ')
%     disp(['REM  Peak [GCaMP7s]: ' num2str(round(data.GCaMP7s.REM.meanMax,1)) ' +/- ' num2str(round(data.GCaMP7s.REM.stdMax,1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     
%     diary off
end

end
