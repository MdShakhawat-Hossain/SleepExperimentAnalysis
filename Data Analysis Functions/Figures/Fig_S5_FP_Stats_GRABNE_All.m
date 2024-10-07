function [AnalysisResults] = Fig_S5_FP_Stats_GRABNE_All(rootFolder,saveFigs,delim,AnalysisResults,FP_animalIDs)
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
behavFields = {'Rest','NREM','REM'}; 
%% Rhodamine comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};    
            data.Rhodamine.(animalID).(behavField).meanAch = mean([AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch; AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE]);        
            data.Rhodamine.(animalID).(behavField).stdAch = std([AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch; AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE])/sqrt(length([AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch; AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE]));        
            % data.Rhodamine.(animalID).(behavField).stdNE = std(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE)/sqrt(length(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE));       
        
    end
end

for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.Rhodamine,behavField) == false
            data.Rhodamine.(behavField).catmeanAch = [];
            % data.Rhodamine.(behavField).catmeanNE = [];
            data.Rhodamine.(behavField).catstdAch = [];
            % data.Rhodamine.(behavField).catstdNE = [];
        end
        data.Rhodamine.(behavField).catmeanAch = cat(1,data.Rhodamine.(behavField).catmeanAch,(data.Rhodamine.(animalID).(behavField).meanAch));
        % data.Rhodamine.(behavField).catmeanNE = cat(1,data.Rhodamine.(behavField).catmeanNE,(data.Rhodamine.(animalID).(behavField).meanNE)); 
        data.Rhodamine.(behavField).catstdAch = cat(1,data.Rhodamine.(behavField).catstdAch,(data.Rhodamine.(animalID).(behavField).stdAch));
        % data.Rhodamine.(behavField).catstdNE = cat(1,data.Rhodamine.(behavField).catstdNE,(data.Rhodamine.(animalID).(behavField).stdNE));  
    end
end
%% GFP comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
            data.GFP.(animalID).(behavField).meanAch = mean(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch);
            data.GFP.(animalID).(behavField).meanNE = mean(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE);  
            data.GFP.(animalID).(behavField).stdAch = std(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch)/sqrt(length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch));
            data.GFP.(animalID).(behavField).stdNE = std(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE)/sqrt(length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE));  
    end
end
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GFP,behavField) == false
            data.GFP.(behavField).catmeanAch = [];
            data.GFP.(behavField).catmeanNE = [];
            data.GFP.(behavField).catstdAch = [];
            data.GFP.(behavField).catstdNE = [];
        end
        data.GFP.(behavField).catmeanAch = cat(1,data.GFP.(behavField).catmeanAch,(data.GFP.(animalID).(behavField).meanAch));
        data.GFP.(behavField).catmeanNE = cat(1,data.GFP.(behavField).catmeanNE,(data.GFP.(animalID).(behavField).meanNE));     
        data.GFP.(behavField).catstdAch = cat(1,data.GFP.(behavField).catstdAch,(data.GFP.(animalID).(behavField).stdAch));
        data.GFP.(behavField).catstdNE = cat(1,data.GFP.(behavField).catstdNE,(data.GFP.(animalID).(behavField).stdNE));  
    end
end
%% Fig. 1-S4
summaryFigure = figure('Name','Fig1-S5 Stats');
sgtitle('hemodynamic changes')
%% Mean CBV LH
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    
    ax1 = subplot(2,2,1);
    xMeans = ones(1,length(data.Rhodamine.Rest.catmeanAch));
    % s1=scatter(xMeans*2,data.Rhodamine.Whisk.catmeanAch,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
    
    
    % s2=scatter(ones(length(data.GFP.Rest.catmeanNE),1),data.Rhodamine.Stim.catmeanAch,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    
   
    plot(2.98,data.Rhodamine.Rest.catmeanAch(1),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
     hold on
    e3 = errorbar(2.98,data.Rhodamine.Rest.catmeanAch(1),data.Rhodamine.Rest.catstdAch(1),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;

    plot(3.02,data.Rhodamine.Rest.catmeanAch(2),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3 = errorbar(3.02,data.Rhodamine.Rest.catmeanAch(2),data.Rhodamine.Rest.catstdAch(2),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;
    
    plot(3.98,data.Rhodamine.NREM.catmeanAch(1),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(3.98,data.Rhodamine.NREM.catmeanAch(1),data.Rhodamine.NREM.catstdAch(1),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;

    plot(4.02,data.Rhodamine.NREM.catmeanAch(2),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(4.02,data.Rhodamine.NREM.catmeanAch(2),data.Rhodamine.NREM.catstdAch(2),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;
    
    plot(4.98,data.Rhodamine.REM.catmeanAch(1),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(4.98,data.Rhodamine.REM.catmeanAch(1),data.Rhodamine.REM.catstdAch(1),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;

    plot(5.02,data.Rhodamine.REM.catmeanAch(2),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(5.02,data.Rhodamine.REM.catmeanAch(2),data.Rhodamine.REM.catstdAch(2),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;
    
    title({'Zscored \DeltaCBV'})
    ylabel('Zscored \DeltaCBV')
    % legend([s3, s4, s5],'Rest','NREM','REM')
        set(gca,'xtick',[3, 4, 5])
    set(gca,'xticklabel',{'Rest','NREM','REM'});
    axis square
    axis tight
    xlim([2 6])
    ylim([-.5,1.5])
    set(gca,'box','off')
    ax1.TickLength = [0.03,0.03];
     %% CBV RH
    % ax1 = subplot(2,2,2);
    % % xMeans = ones(1,length(FP_animalIDs));
    % hold on
    % % s1=scatter(xMeans*2,data.Rhodamine.Whisk.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
    % 
    % 
    % % s2=scatter(ones(length(data.GFP.Rest.catmeanNE),1),data.Rhodamine.Stim.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    % 
    % 
    % plot(2.98,data.Rhodamine.Rest.catmeanNE(1),'o','MarkerEdgeColor','k','MarkerFaceColor',colorRest);
    % hold on
    % e3 = errorbar(2.98,data.Rhodamine.Rest.catmeanNE(1),data.Rhodamine.Rest.catstdNE(1),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    % e3.Color = colorRest;
    % 
    % plot(3.02,data.Rhodamine.Rest.catmeanNE(2),'o','MarkerEdgeColor','k','MarkerFaceColor',colorRest);
    % e3 = errorbar(3.02,data.Rhodamine.Rest.catmeanNE(2),data.Rhodamine.Rest.catstdNE(2),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    % e3.Color = colorRest;
    % 
    % plot(3.98,data.Rhodamine.NREM.catmeanNE(1),'o','MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
    % e4 = errorbar(3.98,data.Rhodamine.NREM.catmeanNE(1),data.Rhodamine.NREM.catstdNE(1),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    % e4.Color = colorNREM;
    % 
    % plot(4.02,data.Rhodamine.NREM.catmeanNE(2),'o','MarkerEdgeColor','k','MarkerFaceColor',colorNREM);
    % e4 = errorbar(4.02,data.Rhodamine.NREM.catmeanNE(2),data.Rhodamine.NREM.catstdNE(2),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    % e4.Color = colorNREM;
    % 
    % plot(4.98,data.Rhodamine.REM.catmeanNE(1),'o','MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    % e5 = errorbar(4.98,data.Rhodamine.REM.catmeanNE(1),data.Rhodamine.REM.catstdNE(1),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    % e5.Color = colorREM;
    % 
    % plot(5.02,data.Rhodamine.REM.catmeanNE(2),'o','MarkerEdgeColor','k','MarkerFaceColor',colorREM);
    % e5 = errorbar(5.02,data.Rhodamine.REM.catmeanNE(2),data.Rhodamine.REM.catstdNE(2),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    % e5.Color = colorREM;
    % 
    % 
    % title({'Zscored \DeltaCBV RH'})
    % ylabel('Zscored \DeltaCBV RH')
    % % legend([s3,s1,s4,s5],'Rest','Whisk','NREM','REM','Location','best')
    % set(gca,'xtick',[3, 4, 5])
    % set(gca,'xticklabel',{'Rest','NREM','REM'});
    % axis square
    % axis tight
    % xlim([2 6])
    % ylim([-.5,2])
    % set(gca,'box','off')
    % ax1.TickLength = [0.03,0.03];
    %% Mean ACh
    ax3 = subplot(2,2,3);
    % xMeans = ones(1,length(FP_animalIDs));
    hold on
    % scatter(xMeans*2,data.GFP.Whisk.catmeanAch,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
    
    
    plot(2.98,data.GFP.Rest.catmeanAch(1),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
     hold on
    e3 = errorbar(2.98,data.GFP.Rest.catmeanAch(1),data.GFP.Rest.catstdAch(1),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;

    plot(3.02,data.GFP.Rest.catmeanAch(2),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3 = errorbar(3.02,data.GFP.Rest.catmeanAch(2),data.GFP.Rest.catstdAch(2),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;
    
    plot(3.98,data.GFP.NREM.catmeanAch(1),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(3.98,data.GFP.NREM.catmeanAch(1),data.GFP.NREM.catstdAch(1),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;

    plot(4.02,data.GFP.NREM.catmeanAch(2),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(4.02,data.GFP.NREM.catmeanAch(2),data.GFP.NREM.catstdAch(2),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;
    
    plot(4.98,data.GFP.REM.catmeanAch(1),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(4.98,data.GFP.REM.catmeanAch(1),data.GFP.REM.catstdAch(1),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;

    plot(5.02,data.GFP.REM.catmeanAch(2),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(5.02,data.GFP.REM.catmeanAch(2),data.GFP.REM.catstdAch(2),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;
    
    %'jitter','on','jitterAmount',0.25
    title({'Zscored \DeltaACh'})
    ylabel('Zscored \DeltaACh')
    set(gca,'xtick',[3, 4, 5])
    set(gca,'xticklabel',{'Rest','NREM','REM'});
    axis square
    axis tight
    xlim([2 6])
    ylim([-3,0.5])
    set(gca,'box','off')
    ax3.TickLength = [0.03,0.03];
    %% Mean NE
    ax3 = subplot(2,2,4);
    % xMeans = ones(1,length(FP_animalIDs));
    hold on
    % scatter(xMeans*2,data.GFP.Whisk.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorWhisk);
    
    
    % scatter(ones(length(data.GFP.Rest.catmeanNE),1),data.GFP.Stim.catmeanNE,75,'MarkerEdgeColor','k','MarkerFaceColor',colorStim);
    
    
       plot(2.98,data.GFP.Rest.catmeanNE(1),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
     hold on
    e3 = errorbar(2.98,data.GFP.Rest.catmeanNE(1),data.GFP.Rest.catstdNE(1),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;

    plot(3.02,data.GFP.Rest.catmeanNE(2),'o','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3 = errorbar(3.02,data.GFP.Rest.catmeanNE(2),data.GFP.Rest.catstdNE(2),'d','MarkerEdgeColor',colorRest,'MarkerFaceColor',colorRest);
    e3.Color = colorRest;
    
    plot(3.98,data.GFP.NREM.catmeanNE(1),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(3.98,data.GFP.NREM.catmeanNE(1),data.GFP.NREM.catstdNE(1),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;

    plot(4.02,data.GFP.NREM.catmeanNE(2),'o','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4 = errorbar(4.02,data.GFP.NREM.catmeanNE(2),data.GFP.NREM.catstdNE(2),'d','MarkerEdgeColor',colorNREM,'MarkerFaceColor',colorNREM);
    e4.Color = colorNREM;
    
    plot(4.98,data.GFP.REM.catmeanNE(1),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(4.98,data.GFP.REM.catmeanNE(1),data.GFP.REM.catstdNE(1),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;

    plot(5.02,data.GFP.REM.catmeanNE(2),'o','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5 = errorbar(5.02,data.GFP.REM.catmeanNE(2),data.GFP.REM.catstdNE(2),'d','MarkerEdgeColor',colorREM,'MarkerFaceColor',colorREM);
    e5.Color = colorREM;
    
    
    title({'Zscored \DeltaNE'})
    ylabel('Zscored \DeltaNE')
        set(gca,'xtick',[3, 4, 5])
    set(gca,'xticklabel',{'Rest','NREM','REM'});
    axis square
    axis tight
    xlim([2 6])
    ylim([-3,0.5])
    set(gca,'box','off')
    ax3.TickLength = [0.03,0.03];
end
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1-S5-Stats_SEM_']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S5-Stats_SEM_'])
    close 
end

end
