function [AnalysisResults] = Fig1_S4_FP_Stats_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FP_animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);
%%
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

    % behavFields = {'Rest','NREM','REM','Stim', 'Whisk'}; %{'Rest','NREM','Stim', 'Whisk'}; %
% elseif firstHrs == "false" 
    % behavFields = {'Rest','NREM','REM', 'Whisk'};
% end
%% CBV comparison between behaviors
% pre-allocate the date for eACh day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};

    if ~isfield(AnalysisResults.(animalID).MeanCBV,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(AnalysisResults.(animalID).MeanCBV,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
            data.CBV.(animalID).(behavField).meanACh = AnalysisResults.(animalID).MeanCBV.(behavField).CBV.MeanACh;
            data.CBV.(animalID).(behavField).meanNE = AnalysisResults.(animalID).MeanCBV.(behavField).CBV.MeanNE;
    end
end


% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};

    if ~isfield(AnalysisResults.(animalID).MeanCBV,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(AnalysisResults.(animalID).MeanCBV,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.CBV,behavField) == false
            data.CBV.(behavField).catmeanACh = [];
            data.CBV.(behavField).catmeanNE = [];
            data.CBV.(behavField).animalID = {};
            data.CBV.(behavField).behavior = {};
            data.CBV.(behavField).hemisphere = {};
        end
        data.CBV.(behavField).catmeanACh = cat(1,data.CBV.(behavField).catmeanACh,mean(data.CBV.(animalID).(behavField).meanACh));
        data.CBV.(behavField).catmeanNE = cat(1,data.CBV.(behavField).catmeanNE,mean(data.CBV.(animalID).(behavField).meanNE));        
        data.CBV.(behavField).animalID = cat(1,data.CBV.(behavField).animalID,animalID,animalID);
        data.CBV.(behavField).behavior = cat(1,data.CBV.(behavField).behavior,behavField,behavField);
        data.CBV.(behavField).hemisphere = cat(1,data.CBV.(behavField).hemisphere,'ACh','NE');
    end
end
%% take mean/StD
    if ~isfield(data.CBV,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(data.CBV,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.CBV.(behavField).meanmeanACh = mean(data.CBV.(behavField).catmeanACh,1);
    data.CBV.(behavField).stdmeanACh = std(data.CBV.(behavField).catmeanACh,0,1);
    data.CBV.(behavField).meanmeanNE = mean(data.CBV.(behavField).catmeanNE,1);
    data.CBV.(behavField).stdmeanNE = std(data.CBV.(behavField).catmeanNE,0,1);    
end
%% combining CBV data from bilateral hemispheres
    % for ee = 1:length(behavFields)
    %     behavField = behavFields{1,ee};
    %     data.CBV.(behavField).catmeanCBV = [];
    %     data.CBV.(behavField).catmeanCBV = cat(1,data.CBV.(behavField).catmeanACh,data.CBV.(behavField).catmeanNE);
    % end
% take mean/StD
% for ff = 1:length(behavFields)
%     behavField = behavFields{1,ff};
%     data.CBV.(behavField).meanmeanCBV = mean(data.CBV.(behavField).catmeanCBV,1);
%     data.CBV.(behavField).stdmeanCBV = std(data.CBV.(behavField).catmeanCBV,0,1);
% end
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        data.CBV.(behavField).catmeanCBV = [];
        data.CBV.(behavField).catmeanCBV = mean([data.CBV.(behavField).catmeanACh,data.CBV.(behavField).catmeanNE],2);
    end
    % mean and std
    for ff = 1:length(behavFields)
        behavField = behavFields{1,ff};
        data.CBV.(behavField).meanmeanCBV = mean(data.CBV.(behavField).catmeanCBV,1);
        data.CBV.(behavField).stdmeanCBV = std(data.CBV.(behavField).catmeanCBV,0,1);
    end
%% GFP comparison between behaviors
% pre-allocate the date for eACh day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};

    if ~isfield(AnalysisResults.(animalID).MeanCBV,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(AnalysisResults.(animalID).MeanCBV,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
            data.GFP.(animalID).(behavField).meanACh = AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanACh;
            data.GFP.(animalID).(behavField).meanNE = AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE;

%         for cc = 1:length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndACh)
%             data.GFP.(animalID).(behavField).meanACh(cc,1) = mean(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndACh{cc,1});
%         end
%         for cc = 1:length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndNE)
%             data.GFP.(animalID).(behavField).meanNE(cc,1) = mean(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndNE{cc,1});
%         end
    end
end
% put data into arrays and prep for stats
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};

    if ~isfield(AnalysisResults.(animalID).MeanCBV,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(AnalysisResults.(animalID).MeanCBV,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GFP,behavField) == false
            data.GFP.(behavField).catmeanACh = [];
            data.GFP.(behavField).catmeanNE = [];
            data.GFP.(behavField).animalID = {};
            data.GFP.(behavField).behavior = {};
            data.GFP.(behavField).hemisphere = {};
        end
        data.GFP.(behavField).catmeanACh = cat(1,data.GFP.(behavField).catmeanACh,mean(data.GFP.(animalID).(behavField).meanACh));
        data.GFP.(behavField).catmeanNE = cat(1,data.GFP.(behavField).catmeanNE,mean(data.GFP.(animalID).(behavField).meanNE));        
        data.GFP.(behavField).animalID = cat(1,data.GFP.(behavField).animalID,animalID,animalID);
        data.GFP.(behavField).behavior = cat(1,data.GFP.(behavField).behavior,behavField,behavField);
        data.GFP.(behavField).hemisphere = cat(1,data.GFP.(behavField).hemisphere,'ACh','NE');
    end
end
%% take mean/StD
    if ~isfield(data.GFP,'REM') % no REM data available
    behavFields = {'Rest','NREM','Stim', 'Whisk'};
    end

    if isfield(data.GFP,'REM') % REM data available
    behavFields = {'Rest','NREM','REM','Stim', 'Whisk'};
    end

for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.GFP.(behavField).meanmeanACh = mean(data.GFP.(behavField).catmeanACh,1);
    data.GFP.(behavField).stdmeanACh = std(data.GFP.(behavField).catmeanACh,0,1);
    data.GFP.(behavField).meanmeanNE = mean(data.GFP.(behavField).catmeanNE,1);
    data.GFP.(behavField).stdmeanNE = std(data.GFP.(behavField).catmeanNE,0,1);    
end
%% statistics - generalized linear mixed-effects model
% CBV
% tableSize = cat(1,data.CBV.Rest.animalID,data.CBV.Stim.animalID,data.CBV.Whisk.animalID,data.CBV.NREM.animalID);
% CBV_meanTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','string','string','double'},'VariableNames',{'Mouse','Hemisphere','Behavior','Mean'});
% CBV_meanTable.Mouse = cat(1,data.CBV.Rest.animalID,data.CBV.Stim.animalID,data.CBV.Whisk.animalID,data.CBV.NREM.animalID,data.CBV.REM.animalID);
% CBV_meanTable.Hemisphere = cat(1,data.CBV.Rest.hemisphere,data.CBV.Stim.hemisphere,data.CBV.Whisk.hemisphere,data.CBV.NREM.hemisphere,data.CBV.REM.hemisphere);
% CBV_meanTable.Behavior = cat(1,data.CBV.Rest.behavior,data.CBV.Stim.behavior,data.CBV.Whisk.behavior,data.CBV.NREM.behavior,data.CBV.REM.behavior);
% CBV_meanTable.Mean = cat(1,data.CBV.Rest.catmean,data.CBV.Stim.catmean,data.CBV.Whisk.catmean,data.CBV.NREM.catmean,data.CBV.REM.catmean);
% CBV_meanFitFormula = 'Mean ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
% CBV_meanStats = fitglme(CBV_meanTable,CBV_meanFitFormula)

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
summaryFigure = figure('Name','Stats');
sgtitle('Mean hemodynamic changes')
%% Mean CBV Bilateral
ax1 = subplot(2,2,1);
xInds = ones(1,length(FP_animalIDs));
s1=scatter(xInds*2,data.CBV.Whisk.catmeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
hold on
e1 = errorbar(2,data.CBV.Whisk.meanmeanCBV,data.CBV.Whisk.stdmeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'k';
e1.MarkerSize = 10;
e1.CapSize = 10;

% if firstHrs == "true"
    s2=scatter(xInds*3,data.CBV.Stim.catmeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    e2 = errorbar(3,data.CBV.Stim.meanmeanCBV,data.CBV.Stim.stdmeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e2.Color = 'k';
    e2.MarkerSize = 10;
    e2.CapSize = 10;
% end

s3=scatter(xInds*1,data.CBV.Rest.catmeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);
e3 = errorbar(1,data.CBV.Rest.meanmeanCBV,data.CBV.Rest.stdmeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'k';
e3.MarkerSize = 10;
e3.CapSize = 10;
hold on
s4=scatter(xInds*4,data.CBV.NREM.catmeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
e4 = errorbar(4,data.CBV.NREM.meanmeanCBV,data.CBV.NREM.stdmeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'k';
e4.MarkerSize = 10;
e4.CapSize = 10;

if isfield(data.GFP,'REM') % REM data available
    s5=scatter(xInds*5,data.CBV.REM.catmeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    e5 = errorbar(5,data.CBV.REM.meanmeanCBV,data.CBV.REM.stdmeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e5.Color = 'k';
    e5.MarkerSize = 10;
    e5.CapSize = 10;
end

title({'Mean \Delta F/F (%)CBV'})
ylabel('Mean  \DeltaCBV LH')

if isfield(data.GFP,'REM') % REM data available
    legend([s3,s1,s2,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','best')
elseif ~sum(data.GFP,'REM') % check if there is no REM data
    legend([s3,s1,s2,s4],'Rest','Whisk','Stim','NREM','Location','best')
end

set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
% axis tight
 xlim([0,6])
% ylim([-1.5,1.5])
set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
%{
%% Mean CBV ACh
ax1 = subplot(2,2,1);
xInds = ones(1,length(FP_animalIDs));
s1=scatter(xInds*2,data.CBV.Whisk.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
hold on
e1 = errorbar(2,data.CBV.Whisk.meanmeanACh,data.CBV.Whisk.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'k';
e1.MarkerSize = 10;
e1.CapSize = 10;

% if firstHrs == "true"
    s2=scatter(xInds*3,data.CBV.Stim.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    e2 = errorbar(3,data.CBV.Stim.meanmeanACh,data.CBV.Stim.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e2.Color = 'k';
    e2.MarkerSize = 10;
    e2.CapSize = 10;
% end

s3=scatter(xInds*1,data.CBV.Rest.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);
e3 = errorbar(1,data.CBV.Rest.meanmeanACh,data.CBV.Rest.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'k';
e3.MarkerSize = 10;
e3.CapSize = 10;
hold on
s4=scatter(xInds*4,data.CBV.NREM.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
e4 = errorbar(4,data.CBV.NREM.meanmeanACh,data.CBV.NREM.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'k';
e4.MarkerSize = 10;
e4.CapSize = 10;

% if firstHrs == "false"
    s5=scatter(xInds*5,data.CBV.REM.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    e5 = errorbar(5,data.CBV.REM.meanmeanACh,data.CBV.REM.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e5.Color = 'k';
    e5.MarkerSize = 10;
    e5.CapSize = 10;
% end

title({'Mean \Delta F/F (%)CBV LH'})
ylabel('Mean  \DeltaCBV LH')
% if firstHrs == "false"
legend([s3,s1,s2,s4,s5],'Rest','Whisk','Stim','NREM','REM','Location','best')
% elseif firstHrs == "true"
% legend([s3,s1,s2,s4],'Rest','Whisk','Stim','NREM','Location','best')
% end
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
% axis tight
 xlim([0,6])
% ylim([-1.5,1.5])
set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
%% mScarlet NE
ax2 = subplot(2,2,2);
xInds = ones(1,length(FP_animalIDs));
s1=scatter(xInds*2,data.CBV.Whisk.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
hold on
e1 = errorbar(2,data.CBV.Whisk.meanmeanNE,data.CBV.Whisk.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

% if firstHrs == "true"
    s2=scatter(xInds*3,data.CBV.Stim.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    e2 = errorbar(3,data.CBV.Stim.meanmeanNE,data.CBV.Stim.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e2.Color = 'black';
    e2.MarkerSize = 10;
    e2.CapSize = 10;
% end

s3=scatter(xInds*1,data.CBV.Rest.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);
e3 = errorbar(1,data.CBV.Rest.meanmeanNE,data.CBV.Rest.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
hold on

s4=scatter(xInds*4,data.CBV.NREM.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
e4 = errorbar(4,data.CBV.NREM.meanmeanNE,data.CBV.NREM.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

% if firstHrs == "false"
    s5=scatter(xInds*5,data.CBV.REM.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    e5 = errorbar(5,data.CBV.REM.meanmeanNE,data.CBV.REM.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e5.Color = 'black';
    e5.MarkerSize = 10;
    e5.CapSize = 10;
% end

title({'Mean \Delta F/F (%)CBV RH'})
ylabel('Mean  \DeltaCBV RH')
% legend([s3,s1,s4,s5],'Rest','Whisk','NREM','REM','Location','best')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
% axis tight
 xlim([0,6])
% ylim([-1.5,1.5])
set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
%}
%% Mean GFP ACh
ax3 = subplot(2,2,3);
xInds = ones(1,length(FP_animalIDs));
scatter(xInds*2,data.GFP.Whisk.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
hold on
ex1 = errorbar(2,data.GFP.Whisk.meanmeanACh,data.GFP.Whisk.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
% ex1.MarkerSize = 10;
% ex1.CapSize = 10;

% if firstHrs == "true"
scatter(xInds*3,data.GFP.Stim.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
e2 = errorbar(3,data.GFP.Stim.meanmeanACh,data.GFP.Stim.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
% e2.MarkerSize = 10;
% e2.CapSize = 10;
% end

scatter(xInds*1,data.GFP.Rest.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);
hold on
e3 = errorbar(1,data.GFP.Rest.meanmeanACh,data.GFP.Rest.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
% e3.MarkerSize = 10;
% e3.CapSize = 10;

hold on

scatter(xInds*4,data.GFP.NREM.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
e4 = errorbar(4,data.GFP.NREM.meanmeanACh,data.GFP.NREM.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
% e4.MarkerSize = 10;
% e4.CapSize = 10;

if isfield(data.GFP,'REM') % REM data available
    scatter(xInds*5,data.GFP.REM.catmeanACh,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    e5 = errorbar(5,data.GFP.REM.meanmeanACh,data.GFP.REM.stdmeanACh,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
end

title({'Mean \Delta F/F (%) ACh'})
ylabel('Mean \Delta F/F (%) ACh')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
 xlim([0,6])
% ylim([-1.5,1.5])
set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
%% Mean GFP NE
ax4 = subplot(2,2,4);
xInds = ones(1,length(FP_animalIDs));
scatter(xInds*2,data.GFP.Whisk.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
hold on
e1 = errorbar(2,data.GFP.Whisk.meanmeanNE,data.GFP.Whisk.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;

% if firstHrs == "true"
scatter(xInds*3,data.GFP.Stim.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
e2 = errorbar(3,data.GFP.Stim.meanmeanNE,data.GFP.Stim.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
% end

scatter(xInds*1,data.GFP.Rest.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);
hold on
e3 = errorbar(1,data.GFP.Rest.meanmeanNE,data.GFP.Rest.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;

hold on

scatter(xInds*4,data.GFP.NREM.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
e4 = errorbar(4,data.GFP.NREM.meanmeanNE,data.GFP.NREM.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
% e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;

if isfield(data.GFP,'REM') % REM data available
    scatter(xInds*5,data.GFP.REM.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    e5 = errorbar(5,data.GFP.REM.meanmeanNE,data.GFP.REM.stdmeanNE,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
    % e5.Color = 'black';
    e5.MarkerSize = 10;
    e5.CapSize = 10;
end

title({'Mean \Delta F/F (%) NE'})
ylabel('Mean \Delta F/F (%) NE')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
% axis tight
 xlim([0,6])
% ylim([-1.5,1.5])
set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath ManipulationType '_Stats_CBV_ACh_NE']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath ManipulationType '_Stats_CBV_ACh_NE'])
    close 
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S4_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     Mean-to-Mean CBV statistical diary
%     disp('======================================================================================================================')
%     disp('[1-S4a] Generalized linear mixed-effects model statistics for Mean-to-Mean CBV during Whisk, NREM, and REM')
%     disp('======================================================================================================================')
%     disp(CBV_p2pStats)
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     disp(['Whisk P2P [CBV]: ' num2str(round(data.CBV.Whisk.meanP2P,1)) ' +/- ' num2str(round(data.CBV.Whisk.stdP2P,1))]); disp(' ')
%     disp(['NREM P2P [CBV]: ' num2str(round(data.CBV.NREM.meanP2P,1)) ' +/- ' num2str(round(data.CBV.NREM.stdP2P,1))]); disp(' ')
%     disp(['REM P2P [CBV]: ' num2str(round(data.CBV.REM.meanP2P,1)) ' +/- ' num2str(round(data.CBV.REM.stdP2P,1))]); disp(' ')
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
