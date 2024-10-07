function [AnalysisResults] = Fig1_S3_Whisk_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
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
        % hippocampal
%         data.(whiskDataType).Hip.hipMUA(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).MUA.hippocampalData;
%         data.(whiskDataType).Hip.hipGam(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).Gam.hippocampalData;
%         data.(whiskDataType).Hip.hipS(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.hippocampalS;
%         data.(whiskDataType).Hip.hipS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.hippocampalS(49:end,20:23);
%         data.(whiskDataType).Hip.hipT(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.T;
%         data.(whiskDataType).Hip.hipF(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.F;
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
%     data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
%     data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
%     data.(whiskDataType).meanHipGam = mean(data.(whiskDataType).Hip.hipGam,2);
%     data.(whiskDataType).stdHipGam = std(data.(whiskDataType).Hip.hipGam,0,2);
%     data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3).*100;
%     data.(whiskDataType).mean_HipS_Gam = mean(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),3);
%     data.(whiskDataType).std_HipS_Gam = std(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),0,3);
%     data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
%     data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end
%% Fig. 1-S3
%{
summaryFigure = figure('Name','Whisk');
sgtitle('Whisk evoked responses')
%% [1-S3a] brief whisks cortical MUA
ax1 = subplot(6,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA + data.ShortWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA - data.ShortWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title(' Brief whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S3b] moderate whisks cortical MUA
ax2 = subplot(6,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA + data.IntermediateWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA - data.IntermediateWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title('Moderate whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S3c] extended whisks cortical MUA
ax3 = subplot(6,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA + data.LongWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA - data.LongWhisks.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title('Extended whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S3d] brief whisks cortical LFP
ax4 = subplot(6,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title(' Brief whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S3e] moderate whisks cortical LFP
ax5 = subplot(6,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title(' Moderate whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S3f] extended whisks cortical LFP
ax6 = subplot(6,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('Extended whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S3g] brief whisks hippocampal MUA
% ax7 = subplot(6,3,7);
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA,'-','color',colors('rich black'),'LineWidth',2);
% hold on
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA + data.ShortWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA - data.ShortWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Brief whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax7.TickLength = [0.03,0.03];
%% [1-S3h] moderate whisks hippocampal MUA
% ax8 = subplot(6,3,8);
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA,'-','color',colors('rich black'),'LineWidth',2);
% hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA + data.IntermediateWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA - data.IntermediateWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Moderate whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];
%% [1-S3i] extended whisks hippocampal MUA
% ax9 = subplot(6,3,9);
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA,'-','color',colors('rich black'),'LineWidth',2);
% hold on
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA + data.LongWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA - data.LongWhisks.stdHipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Extended whisk hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-whisk time (s)')
% axis square
% set(gca,'box','off')
% ax9.TickLength = [0.03,0.03];
%% [1-S3j] brief whisks hippocampal LFP
% ax10 = subplot(6,3,10);
% imagesc(data.ShortWhisks.meanHipT,data.ShortWhisks.meanHipF,data.ShortWhisks.meanHipS)
% title('Brief whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c10 = colorbar;
% ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,100])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax10.TickLength = [0.03,0.03];
%% [1-S3k] moderate whisks hippocampal LFP
% ax11 = subplot(6,3,11);
% imagesc(data.IntermediateWhisks.meanHipT,data.IntermediateWhisks.meanHipF,data.IntermediateWhisks.meanHipS)
% title(' Moderate whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c11 = colorbar;
% ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,100])
% set(gca,'Ticklength',[0 0])
% axis square
% axis xy
% set(gca,'box','off')
% ax11.TickLength = [0.03,0.03];
%% [1-S3l] extended whisks hippocampal LFP
% ax12 = subplot(6,3,12);
% imagesc(data.LongWhisks.meanHipT,data.LongWhisks.meanHipF,data.LongWhisks.meanHipS)
% title('Extended whisk hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-whisk time (s)')
% c12 = colorbar;
% ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,100])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];
%% [1-S3m] brief whisks Rhodamine
ax13 = subplot(6,3,13);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine + data.ShortWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine - data.ShortWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Brief whisk Blood Volume (Z)')
ylabel('Zscored \DeltaF/F Ach')
ax13.YLim = [-1 1];

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine + data.ShortWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine - data.ShortWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax13.YAxis(1).Color = colors('indian red');
ax13.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
ax13.YLim = [-1 1];
%% [1-S3n] moderate whisks Rhodamine
ax14 = subplot(6,3,14);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine + data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine - data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Moderate whisk Blood Volume (Z)')
ylabel('Zscored \DeltaF/F Ach')
ax14.YLim = [-1 1];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine + data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine - data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax14.YAxis(1).Color = colors('indian red');
ax14.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
ax14.YLim = [-1 1];
%% [1-S3o] extended whisks Rhodamine
ax15 = subplot(6,3,15);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine + data.LongWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine - data.LongWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F Ach')
ax15.YLim = [-1 1];

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine + data.LongWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine - data.LongWhisks.Z_NEstdRhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax15.YAxis(1).Color = colors('indian red');
ax15.YAxis(2).Color = colors('army green');
title('Extended whisk Blood Volume (Z)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
ax15.YLim = [-1 1];
%% [1-S3p] brief whisks GFP
ax16 = subplot(6,3,16);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP + data.ShortWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP - data.ShortWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F Ach')
ax16.YLim = [-1 1];

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP + data.ShortWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP - data.ShortWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax16.YAxis(1).Color = colors('indian red');
ax16.YAxis(2).Color = colors('army green');
title('Brief whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
ax16.YLim = [-1 1];
%% [1-S3q] moderate whisks GFP
ax17 = subplot(6,3,17);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP + data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP - data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F Ach')
ax17.YLim = [-1 1];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP + data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP - data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax17.YAxis(1).Color = colors('indian red');
ax17.YAxis(2).Color = colors('army green');
title('Moderate whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
ax17.YLim = [-1 1];
%% [1-S3r] extended whisks GFP
ax18 = subplot(6,3,18);

plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP + data.LongWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP - data.LongWhisks.Z_AchstdGFP,'-','color',colors('indian red'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F Ach')
ax18.YLim = [-1 1];

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP + data.LongWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP - data.LongWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax18.YAxis(1).Color = colors('indian red');
ax18.YAxis(2).Color = colors('army green');
title('Extended whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
ax18.YLim = [-1 1];
%% axes positions
linkaxes([ax1,ax2,ax3],'xy')
linkaxes([ax4,ax5,ax6],'xy')
linkaxes([ax13,ax14,ax15],'xy')
linkaxes([ax16,ax17,ax18],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
% ax10Pos(3:4) = ax1Pos(3:4);
% ax11Pos(3:4) = ax2Pos(3:4);
% ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
% set(ax10,'position',ax10Pos);
% set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);

%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Whisk']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Whisk'])
    close 
end
%}
%{
%% plot the comparison of GCaMP and Rhodamine with GRABNE and Rhodamine
summaryFigureN = figure('Name','Whisk');
sgtitle('Whisking evoked responses in fiber photometry signals')

% Ach contra stim
ax1 = subplot(3,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine + data.ShortWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine - data.ShortWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Brief whisk Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax1.YLim = [-1 1];

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP + data.ShortWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP - data.ShortWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')

ax1.YAxis(1).Color = colors('indian red');
ax1.YAxis(2).Color = colors('army green');

xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.YLim = [-1 1];
xlim([-5 10])

% Ach ispi stim
ax2 = subplot(3,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine + data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine - data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Moderate whisk Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax2.YLim = [-1 1];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP + data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP - data.IntermediateWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')
ax2.YAxis(1).Color = colors('indian red');
ax2.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-1 1];
xlim([-5 10])

% Ach auditory stim
ax3 = subplot(3,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine + data.LongWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine - data.LongWhisks.Z_AchstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Long whisk Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax3.YLim = [-1 1];

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP + data.LongWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP - data.LongWhisks.Z_AchstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')
ax3.YAxis(1).Color = colors('indian red');
ax3.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YLim = [-1 1];
xlim([-5 10])

% NE contra stim
ax4 = subplot(3,3,4);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine + data.ShortWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine - data.ShortWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Brief whisk NE fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax4.YLim = [-1 1];

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP + data.ShortWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP - data.ShortWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax4.YAxis(1).Color = colors('indian red');
ax4.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
ax4.YLim = [-1 1];
xlim([-5 10])

% NE ispi stim
ax5 = subplot(3,3,5);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine + data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine - data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Moderate whisk NE fibers')
ylabel('\DeltaF/F Blood Volume (Z)')
ax5.YLim = [-1 1];

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP + data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP - data.IntermediateWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
ax5.YLim = [-1 1];
xlim([-5 10])

% NE auditory stim
ax6 = subplot(3,3,6);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine + data.LongWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine - data.LongWhisks.Z_NEstdRhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Long whisk NE fibers')
ylabel('\DeltaF/F Blood Volume (Z)')
ax6.YLim = [-1 1];

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP + data.LongWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP - data.LongWhisks.Z_NEstdGFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax6.YAxis(1).Color = colors('indian red');
ax6.YAxis(2).Color = colors('army green');
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
ax6.YLim = [-1 1];
xlim([-5 10])

% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'Whisk-FiberSignals']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Whisk-FiberSignals'])
    close 
end
%}

%% Only the fiber signals consolidated
summaryFigureN = figure('Name','Whisk_fiber_consolidated');
sgtitle('Whisking evoked responses in fiber photometry signals')

% Short Movement
ax1 = subplot(2,2,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine + data.ShortWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanRhodamine - data.ShortWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Brief Whisks')

plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine + data.ShortWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanRhodamine - data.ShortWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('\DeltaF/F Blood Volume (Z)')
ax1.YLim = [-2 6];

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP + data.ShortWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_AchmeanGFP - data.ShortWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
%
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP + data.ShortWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP - data.ShortWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

ylabel('\DeltaF/F GRAB(Z)')

ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-Whisks time (s)')
axis square
set(gca,'box','off')
ax1.YLim = [-2 6];
xlim([-5 15])


% Intermediate Movement
ax2 = subplot(2,2,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine + data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanRhodamine - data.IntermediateWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine + data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanRhodamine - data.IntermediateWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

title('Moderate Whisks')
ylabel('\DeltaF/F (Z)')
% ax2.YLim = [-2 6];

% yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP + data.IntermediateWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_AchmeanGFP - data.IntermediateWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)

plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP + data.IntermediateWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP - data.IntermediateWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('\DeltaF/F GRAB(Z)')
% ax2.YAxis(1).Color = 'k';%[0.8500 0.3250 0.0980];
% ax2.YAxis(2).Color = 'k';%[;%0 0.4470 0.7410];
xlabel('Peri-Whisks time (s)')
axis square
set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ax2.YLim = [-2 6];
ylim([-2 6])
xlim([-5 15])
legend('CBV-LH','CBV-RH', 'ACh', 'NE')

% Long Movement
ax3 = subplot(2,2,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine + data.LongWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanRhodamine - data.LongWhisks.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine + data.LongWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanRhodamine - data.LongWhisks.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

title('Long Whisks')
ylabel('\DeltaF/F Blood Volume (Z)')
ax3.YLim = [-2 6];

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP + data.LongWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_AchmeanGFP - data.LongWhisks.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)

plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP + data.LongWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP - data.LongWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

ylabel('\DeltaF/F GRAB(Z)')
ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax3.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-Whisks time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YLim = [-2 6];
xlim([-5 15])

% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'Whisk-FiberSignals-consolidated']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Whisk-FiberSignals-consolidated'])
    close 
end
end
