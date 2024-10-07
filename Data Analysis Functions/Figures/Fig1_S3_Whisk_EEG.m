function [AnalysisResults] = Fig1_S3_Whisk_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'SHF025'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        data.(whiskDataType).RH.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).RH.GFP(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).GFP.GFP;

        data.(whiskDataType).EEG.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).MUA.EEGData;
        data.(whiskDataType).EEG.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).Gam.EEGData;
        data.(whiskDataType).EEG.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).LFP.EEGS;
        data.(whiskDataType).EEG.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).LFP.EEGS(49:end,20:23);
        data.(whiskDataType).EEG.cortT(:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).LFP.T;
        data.(whiskDataType).EEG.cortF(:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).LFP.F;

        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.EEG.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).RHRhodamine = data.(whiskDataType).RH.Rhodamine;
    data.(whiskDataType).RHGFP = data.(whiskDataType).RH.GFP;

    data.(whiskDataType).cortMUA = data.(whiskDataType).EEG.cortMUA;
    data.(whiskDataType).cortGam = data.(whiskDataType).EEG.cortGam;
    data.(whiskDataType).cortS = data.(whiskDataType).EEG.cortS;
    data.(whiskDataType).cortS_Gam = data.(whiskDataType).EEG.cortS_Gam;
    data.(whiskDataType).cortT = data.(whiskDataType).EEG.cortT;
    data.(whiskDataType).cortF = data.(whiskDataType).EEG.cortF;
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).RHmeanRhodamine = mean(data.(whiskDataType).RHRhodamine,2);
    data.(whiskDataType).RHstdRhodamine = std(data.(whiskDataType).RHRhodamine,0,2);

    data.(whiskDataType).RHmeanGFP = mean(data.(whiskDataType).RHGFP,2);
    data.(whiskDataType).RHstdGFP = std(data.(whiskDataType).RHGFP,0,2);

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
%% Fig. 1-S3
summaryFigure = figure('Name','Fig1-S3 Whisk');
sgtitle('Whisk evoked responses')
%% [1-S3a] brief whisks EEG MUA
ax1 = subplot(4,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA + data.ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA - data.ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title(' Brief whisk EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S3b] moderate whisks EEG MUA
ax2 = subplot(4,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA + data.IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA - data.IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title('Moderate whisk EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S3c] extended whisks EEG MUA
ax3 = subplot(4,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA + data.LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA - data.LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title('Extended whisk EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S3d] brief whisks EEG LFP
ax4 = subplot(4,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title(' Brief whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S3e] moderate whisks EEG LFP
ax5 = subplot(4,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title(' Moderate whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S3f] extended whisks EEG LFP
ax6 = subplot(4,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('Extended whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

%% [1-S3m] brief whisks Rhodamine
ax13 = subplot(4,3,7);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine + data.ShortWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine - data.ShortWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Brief whisk Rhodamine')
ylabel('Zscored \DeltaF/F RH')


xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S3n] moderate whisks Rhodamine
ax14 = subplot(4,3,8);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine + data.IntermediateWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine - data.IntermediateWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Moderate whisk Rhodamine')
ylabel('Zscored \DeltaF/F RH')


xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S3o] extended whisks Rhodamine
ax15 = subplot(4,3,9);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine + data.LongWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine - data.LongWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
ylabel('Zscored \DeltaF/F RH')


title('Extended whisk Rhodamine')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S3p] brief whisks GFP
ax16 = subplot(4,3,10);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP + data.ShortWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP - data.ShortWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
ylabel('Zscored \DeltaF/F RH')

title('Brief whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S3q] moderate whisks GFP
ax17 = subplot(4,3,11);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP + data.IntermediateWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP - data.IntermediateWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
ylabel('Zscored \DeltaF/F RH')

title('Moderate whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S3r] extended whisks GFP
ax18 = subplot(4,3,12);

plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP + data.LongWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP - data.LongWhisks.RHstdGFP,'color',colors('indian red'),'LineWidth',0.25)
ylabel('Zscored \DeltaF/F RH')

title('Extended whisk GFP')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
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
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S3-Whisk']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S3-Whisk'])
    close 
end
%% plot the comparison of GCaMP and Rhodamine with GRABNE and Rhodamine
summaryFigureN = figure('Name','Fig1-S2 Stim');
sgtitle('Stimulus evoked responses in fiber photometry signals')

% RH contra stim
ax1 = subplot(1,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine + data.ShortWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanRhodamine - data.ShortWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Brief whisk RH fiber')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP + data.ShortWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.RHmeanGFP - data.ShortWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GCaMP7s')

ax1.YAxis(1).Color = colors('indian red');
ax1.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];

% RH ispi stim
ax2 = subplot(1,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine + data.IntermediateWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanRhodamine - data.IntermediateWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Moderate whisk RH fiber')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP + data.IntermediateWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.RHmeanGFP - data.IntermediateWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GCaMP7s')
ax2.YAxis(1).Color = colors('indian red');
ax2.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% RH auditory stim
ax3 = subplot(1,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine + data.LongWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanRhodamine - data.LongWhisks.RHstdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Long whisk RH fiber')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP + data.LongWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.RHmeanGFP - data.LongWhisks.RHstdGFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GCaMP7s')
ax3.YAxis(1).Color = colors('indian red');
ax3.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'Fig1-S3-Whisk-FiberSignals']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2-Whisk-FiberSignals'])
    close 
end
end
