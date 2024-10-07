% function [AnalysisResults] = Fig1_S3_Whisk_GRABNE_PharacologyCompare(rootFolder,saveFigs,delim,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
load('S:\NEACh\NEACh001\Tradozone\\AnalysisResults_firstHrs.mat');
AnalysisResults = AnalysisResults_firstHrs;
%% set-up and process data
animalIDs = {'NEACh001'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        data.(whiskDataType).Z_Ach.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.Z_Ach.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).Z_Ach.GFP(:,aa) = AnalysisResults.(animalID).Whisk.Z_Ach.(whiskDataType).GFP.GFP;
        data.(whiskDataType).Z_NE.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.Z_NE.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).Z_NE.GFP(:,aa) = AnalysisResults.(animalID).Whisk.Z_NE.(whiskDataType).GFP.GFP;

        data.(whiskDataType).cortical.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).cortical.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).Gam.corticalData;
        data.(whiskDataType).cortical.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).cortical.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.(whiskDataType).cortical.cortT(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.T;
        data.(whiskDataType).cortical.cortF(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.F;
        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).Z_AchRhodamine = data.(whiskDataType).Z_Ach.Rhodamine;
    data.(whiskDataType).Z_AchGFP = data.(whiskDataType).Z_Ach.GFP;
    data.(whiskDataType).Z_NERhodamine = data.(whiskDataType).Z_NE.Rhodamine;
    data.(whiskDataType).Z_NEGFP = data.(whiskDataType).Z_NE.GFP;

    data.(whiskDataType).cortMUA = data.(whiskDataType).cortical.cortMUA;
    data.(whiskDataType).cortGam = data.(whiskDataType).cortical.cortGam;
    data.(whiskDataType).cortS = data.(whiskDataType).cortical.cortS;
    data.(whiskDataType).cortS_Gam = data.(whiskDataType).cortical.cortS_Gam;
    data.(whiskDataType).cortT = data.(whiskDataType).cortical.cortT;
    data.(whiskDataType).cortF = data.(whiskDataType).cortical.cortF;
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).Z_AchmeanRhodamine = mean(data.(whiskDataType).Z_AchRhodamine,2);
    data.(whiskDataType).Z_AchstdRhodamine = std(data.(whiskDataType).Z_AchRhodamine,0,2);
    data.(whiskDataType).Z_NEmeanRhodamine = mean(data.(whiskDataType).Z_NERhodamine,2);
    data.(whiskDataType).Z_NEstdRhodamine = std(data.(whiskDataType).Z_NERhodamine,0,2);

    data.(whiskDataType).Z_AchmeanGFP = mean(data.(whiskDataType).Z_AchGFP,2);
    data.(whiskDataType).Z_AchstdGFP = std(data.(whiskDataType).Z_AchGFP,0,2);
    data.(whiskDataType).Z_NEmeanGFP = mean(data.(whiskDataType).Z_NEGFP,2);
    data.(whiskDataType).Z_NEstdGFP = std(data.(whiskDataType).Z_NEGFP,0,2);

    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortGam = mean(data.(whiskDataType).cortGam,2);
    data.(whiskDataType).stdCortGam = std(data.(whiskDataType).cortGam,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    data.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),3);
    data.(whiskDataType).std_CortS_Gam = std(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),0,3);
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end
tradozone_data = data;
clear data AnalysisResults
%%
load('S:\NEACh\AnalysisResults_firstHrs.mat');
AnalysisResults = AnalysisResults_firstHrs;

%% set-up and process data
animalIDs = {'NEACh001'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        data.(whiskDataType).Z_Ach.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.Z_Ach.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).Z_Ach.GFP(:,aa) = AnalysisResults.(animalID).Whisk.Z_Ach.(whiskDataType).GFP.GFP;
        data.(whiskDataType).Z_NE.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.Z_NE.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).Z_NE.GFP(:,aa) = AnalysisResults.(animalID).Whisk.Z_NE.(whiskDataType).GFP.GFP;

        data.(whiskDataType).cortical.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).cortical.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).Gam.corticalData;
        data.(whiskDataType).cortical.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).cortical.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.(whiskDataType).cortical.cortT(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.T;
        data.(whiskDataType).cortical.cortF(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.F;
        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).Z_AchRhodamine = data.(whiskDataType).Z_Ach.Rhodamine;
    data.(whiskDataType).Z_AchGFP = data.(whiskDataType).Z_Ach.GFP;
    data.(whiskDataType).Z_NERhodamine = data.(whiskDataType).Z_NE.Rhodamine;
    data.(whiskDataType).Z_NEGFP = data.(whiskDataType).Z_NE.GFP;

    data.(whiskDataType).cortMUA = data.(whiskDataType).cortical.cortMUA;
    data.(whiskDataType).cortGam = data.(whiskDataType).cortical.cortGam;
    data.(whiskDataType).cortS = data.(whiskDataType).cortical.cortS;
    data.(whiskDataType).cortS_Gam = data.(whiskDataType).cortical.cortS_Gam;
    data.(whiskDataType).cortT = data.(whiskDataType).cortical.cortT;
    data.(whiskDataType).cortF = data.(whiskDataType).cortical.cortF;
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).Z_AchmeanRhodamine = mean(data.(whiskDataType).Z_AchRhodamine,2);
    data.(whiskDataType).Z_AchstdRhodamine = std(data.(whiskDataType).Z_AchRhodamine,0,2);
    data.(whiskDataType).Z_NEmeanRhodamine = mean(data.(whiskDataType).Z_NERhodamine,2);
    data.(whiskDataType).Z_NEstdRhodamine = std(data.(whiskDataType).Z_NERhodamine,0,2);

    data.(whiskDataType).Z_AchmeanGFP = mean(data.(whiskDataType).Z_AchGFP,2);
    data.(whiskDataType).Z_AchstdGFP = std(data.(whiskDataType).Z_AchGFP,0,2);
    data.(whiskDataType).Z_NEmeanGFP = mean(data.(whiskDataType).Z_NEGFP,2);
    data.(whiskDataType).Z_NEstdGFP = std(data.(whiskDataType).Z_NEGFP,0,2);

    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortGam = mean(data.(whiskDataType).cortGam,2);
    data.(whiskDataType).stdCortGam = std(data.(whiskDataType).cortGam,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    data.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),3);
    data.(whiskDataType).std_CortS_Gam = std(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),0,3);
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end

%% plot the comparison of GCaMP and Rhodamine with GRABNE and Rhodamine
summaryFigureN = figure('Name','Fig1-S3 Stim');
sgtitle('whisking evoked responses in fiber photometry signals')


% Ach intermediate
ax2 = subplot(2,3,1);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(tradozone_data.IntermediateWhisks.meanTimeVector,tradozone_data.IntermediateWhisks.Z_AchmeanRhodamine,'-','color',colors('electric purple'),'LineWidth',2)

% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine + data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine - data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Moderate whisk Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax2.YLim = [-3 6];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(tradozone_data.IntermediateWhisks.meanTimeVector,tradozone_data.IntermediateWhisks.Z_AchmeanGFP,'-','color',colors('international klein blue'),'LineWidth',2)

% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP + data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP - data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')
ax2.YAxis(1).Color = colors('indian red');
ax2.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-3 6];
xlim([-5 10])
legend('CBV','CBV-T','ACh','ACh-T')




% NE intermediate whisk
ax5 = subplot(2,3,4);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(tradozone_data.IntermediateWhisks.meanTimeVector,tradozone_data.IntermediateWhisks.Z_NEmeanRhodamine,'-','color',colors('electric purple'),'LineWidth',2)

% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine + data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine - data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Moderate whisk NE fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax5.YLim = [-3 6];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(tradozone_data.IntermediateWhisks.meanTimeVector,tradozone_data.IntermediateWhisks.Z_NEmeanGFP,'-','color',colors('international klein blue'),'LineWidth',2)

% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP + data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP - data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
ax5.YLim = [-3 6];
xlim([-5 10])

legend('CBV','CBV-T','ACh','ACh-T')

% % save figure(s)
% if strcmp(saveFigs,'y') == true
%     if firstHrs == "false"
%         dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
%     elseif firstHrs == "true"
%         dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
%     end
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigureN,[dirpath 'Fig1-S3-Whisk-FiberSignals']);
%     set(summaryFigureN,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S3-Whisk-FiberSignals'])
%     close 
% end
% end
