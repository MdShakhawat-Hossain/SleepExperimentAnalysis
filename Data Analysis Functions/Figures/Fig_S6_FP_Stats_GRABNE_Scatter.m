function [AnalysisResults] = Fig_S6_FP_Stats_GRABNE_Scatter(rootFolder,saveFigs,delim,AnalysisResults,FP_animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorWhisk = [(0/256),(166/256),(81/256)];
% colorRest = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
colorRest = [(255/256),(191/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
behavFields = {'Rest','Whisk','NREM','REM'}; 
%% Rhodamine comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};    
            data.Rhodamine.(animalID).(behavField).meanAch = mean([AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch,AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch],2);        
            % data.Rhodamine.(animalID).(behavField).meanNE = mean(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanNE);
            % data.Rhodamine.(animalID).(behavField).stdAch = std(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch)/sqrt(length(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.MeanAch));        
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
            % data.Rhodamine.(behavField).catstdAch = [];
            % data.Rhodamine.(behavField).catstdNE = [];
        end
        % data.Rhodamine.(behavField).catmeanAch = cat(1,data.Rhodamine.(behavField).catmeanAch,(data.Rhodamine.(animalID).(behavField).meanAch));
        % data.Rhodamine.(behavField).catmeanNE = cat(1,data.Rhodamine.(behavField).catmeanNE,(data.Rhodamine.(animalID).(behavField).meanNE)); 
        % data.Rhodamine.(behavField).catstdAch = cat(1,data.Rhodamine.(behavField).catstdAch,(data.Rhodamine.(animalID).(behavField).stdAch));
        % data.Rhodamine.(behavField).catstdNE = cat(1,data.Rhodamine.(behavField).catstdNE,(data.Rhodamine.(animalID).(behavField).stdNE));  
        data.Rhodamine.(behavField).catmeanAch = cat(1,data.Rhodamine.(behavField).catmeanAch,(data.Rhodamine.(animalID).(behavField).meanAch));
        % data.Rhodamine.(behavField).catmeanNE = cat(1,data.Rhodamine.(behavField).catmeanNE,(data.Rhodamine.(animalID).(behavField).meanNE)); 
        % data.Rhodamine.(behavField).catstdAch = cat(1,data.Rhodamine.(behavField).catstdAch,(data.Rhodamine.(animalID).(behavField).stdAch));
        % data.Rhodamine.(behavField).catstdNE = cat(1,data.Rhodamine.(behavField).catstdNE,(data.Rhodamine.(animalID).(behavField).stdNE));  
    end
end
%% GFP comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
            data.GFP.(animalID).(behavField).meanAch = (AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch);
            data.GFP.(animalID).(behavField).meanNE = (AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE);  
            % data.GFP.(animalID).(behavField).stdAch = std(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch)/sqrt(length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanAch));
            % data.GFP.(animalID).(behavField).stdNE = std(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE)/sqrt(length(AnalysisResults.(animalID).MeanGFP.(behavField).GFP.MeanNE));  
    end
end
for dd = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.GFP,behavField) == false
            data.GFP.(behavField).catmeanAch = [];
            data.GFP.(behavField).catmeanNE = [];
            % data.GFP.(behavField).catstdAch = [];
            % data.GFP.(behavField).catstdNE = [];
        end
        data.GFP.(behavField).catmeanAch = cat(1,data.GFP.(behavField).catmeanAch,(data.GFP.(animalID).(behavField).meanAch));
        data.GFP.(behavField).catmeanNE = cat(1,data.GFP.(behavField).catmeanNE,(data.GFP.(animalID).(behavField).meanNE));     
        % data.GFP.(behavField).catstdAch = cat(1,data.GFP.(behavField).catstdAch,(data.GFP.(animalID).(behavField).stdAch));
        % data.GFP.(behavField).catstdNE = cat(1,data.GFP.(behavField).catstdNE,(data.GFP.(animalID).(behavField).stdNE));  
    end
end
%% Fig. 1-S4
summaryFigure = figure('Name','Fig-S6 Scatter');
sgtitle('Relation Betweeen CBV, ACh, and NE')
%% Mean CBV LH

    ax1 =  subplot(1,2,1);
       
    p1= scatter3(data.GFP.Rest.catmeanNE,data.Rhodamine.Rest.catmeanAch,data.GFP.Rest.catmeanAch,40,'MarkerEdgeColor','k','MarkerFaceColor','k');
    hold on
    p4= scatter3(data.GFP.Whisk.catmeanNE,data.Rhodamine.Whisk.catmeanAch,data.GFP.Whisk.catmeanAch,40,'MarkerEdgeColor',colorWhisk,'MarkerFaceColor',colorWhisk);
    p2 = scatter3(data.GFP.NREM.catmeanNE,data.Rhodamine.NREM.catmeanAch,data.GFP.NREM.catmeanAch,40,'MarkerEdgeColor','b','MarkerFaceColor','b');
    p3 = scatter3(data.GFP.REM.catmeanNE,data.Rhodamine.REM.catmeanAch,data.GFP.REM.catmeanAch,40,'MarkerEdgeColor','r','MarkerFaceColor','r');
    
    % title({'3D Scatter plot of Zscored \DeltaF/F'})
    ylabel('CBV (Zscored \DeltaF/F))')
    zlabel('ACh (Zscored \DeltaF/F)')
    xlabel('NE (Zscored \DeltaF/F)')
    axis square
    axis tight
    legend([p1,p2,p3,p4],'Rest','NREM','REM','Whisk','Location','best')
    xlim([-3 3])
    ylim([-4 4])
    zlim([-3 3])
    set(gca,'box','off')
    ax1.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1-S6-Scatter']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S6-Scatter'])
    close 
end

end
