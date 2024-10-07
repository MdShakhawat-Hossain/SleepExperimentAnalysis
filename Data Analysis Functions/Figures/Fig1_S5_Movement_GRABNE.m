function [AnalysisResults] = Fig1_S5_Movement_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%________________________________________________________________________________________________________________________

%% set-up and process data
movementDataTypes = {'ShortMovement'};%,'IntermediateMovement'};%,'LongMovement'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(movementDataTypes)
        movementDataType = movementDataTypes{1,cc};
        data.(movementDataType).Z_Ach.Rhodamine(:,aa) = AnalysisResults.(animalID).Movement.Z_Ach.(movementDataType).Rhodamine.Rhodamine;
        data.(movementDataType).Z_Ach.GFP(:,aa) = AnalysisResults.(animalID).Movement.Z_Ach.(movementDataType).GFP.GFP;
        data.(movementDataType).Z_NE.Rhodamine(:,aa) = AnalysisResults.(animalID).Movement.Z_NE.(movementDataType).Rhodamine.Rhodamine;
        data.(movementDataType).Z_NE.GFP(:,aa) = AnalysisResults.(animalID).Movement.Z_NE.(movementDataType).GFP.GFP;

        data.(movementDataType).cortical.cortMUA(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).MUA.corticalData;
        data.(movementDataType).cortical.cortGam(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).Gam.corticalData;
        data.(movementDataType).cortical.cortS(:,:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).LFP.corticalS;
        data.(movementDataType).cortical.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).LFP.corticalS(49:end,20:23);
        data.(movementDataType).cortical.cortT(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).LFP.T;
        data.(movementDataType).cortical.cortF(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).LFP.F;
        % time vector
        data.(movementDataType).timeVector(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).timeVector;
    end
end
% 
for ee = 1:length(movementDataTypes)
    movementDataType = movementDataTypes{1,ee};
    data.(movementDataType).Z_AchRhodamine = data.(movementDataType).Z_Ach.Rhodamine;
    data.(movementDataType).Z_AchGFP = data.(movementDataType).Z_Ach.GFP;
    data.(movementDataType).Z_NERhodamine = data.(movementDataType).Z_NE.Rhodamine;
    data.(movementDataType).Z_NEGFP = data.(movementDataType).Z_NE.GFP;

    data.(movementDataType).cortMUA = data.(movementDataType).cortical.cortMUA;
    data.(movementDataType).cortGam = data.(movementDataType).cortical.cortGam;
    data.(movementDataType).cortS = data.(movementDataType).cortical.cortS;
    data.(movementDataType).cortS_Gam = data.(movementDataType).cortical.cortS_Gam;
    data.(movementDataType).cortT = data.(movementDataType).cortical.cortT;
    data.(movementDataType).cortF = data.(movementDataType).cortical.cortF;
end
% mean
for ee = 1:length(movementDataTypes)
    movementDataType = movementDataTypes{1,ee};
    data.(movementDataType).Z_AchmeanRhodamine = mean(data.(movementDataType).Z_AchRhodamine,2);
    data.(movementDataType).Z_AchstdRhodamine = std(data.(movementDataType).Z_AchRhodamine,0,2);
    data.(movementDataType).Z_NEmeanRhodamine = mean(data.(movementDataType).Z_NERhodamine,2);
    data.(movementDataType).Z_NEstdRhodamine = std(data.(movementDataType).Z_NERhodamine,0,2);

    data.(movementDataType).Z_AchmeanGFP = mean(data.(movementDataType).Z_AchGFP,2);
    data.(movementDataType).Z_AchstdGFP = std(data.(movementDataType).Z_AchGFP,0,2);
    data.(movementDataType).Z_NEmeanGFP = mean(data.(movementDataType).Z_NEGFP,2);
    data.(movementDataType).Z_NEstdGFP = std(data.(movementDataType).Z_NEGFP,0,2);

    data.(movementDataType).meanCortMUA = mean(data.(movementDataType).cortMUA,2);
    data.(movementDataType).stdCortMUA = std(data.(movementDataType).cortMUA,0,2);
    data.(movementDataType).meanCortGam = mean(data.(movementDataType).cortGam,2);
    data.(movementDataType).stdCortGam = std(data.(movementDataType).cortGam,0,2);
    data.(movementDataType).meanCortS = mean(data.(movementDataType).cortS,3).*100;
    data.(movementDataType).mean_CortS_Gam = mean(mean(mean(data.(movementDataType).cortS_Gam.*100,2),1),3);
    data.(movementDataType).std_CortS_Gam = std(mean(mean(data.(movementDataType).cortS_Gam.*100,2),1),0,3);
    data.(movementDataType).meanCortT = mean(data.(movementDataType).cortT,2);
    data.(movementDataType).meanCortF = mean(data.(movementDataType).cortF,2);
    data.(movementDataType).meanTimeVector = mean(data.(movementDataType).timeVector(:,aa),2);
end
%% Fig. 1-S3
%{
summaryFigure = figure('Name','Movement');
sgtitle('Movement evoked responses')
%% [1-S3a] brief movements cortical MUA
ax1 = subplot(6,3,1);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.meanCortMUA + data.ShortMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.meanCortMUA - data.ShortMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title(' Brief movement cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S3b] moderate movements cortical MUA
% ax2 = subplot(6,3,2);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.meanCortMUA + data.IntermediateMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.meanCortMUA - data.IntermediateMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Moderate movement cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
%% [1-S3c] extended movements cortical MUA
% ax3 = subplot(6,3,3);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.meanCortMUA,'-','color',colors('rich black'),'LineWidth',2);
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.meanCortMUA + data.LongMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.meanCortMUA - data.LongMovement.stdCortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Extended movement cortical MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
%% [1-S3d] brief movements cortical LFP
ax4 = subplot(6,3,4);
imagesc(data.ShortMovement.meanCortT,data.ShortMovement.meanCortF,data.ShortMovement.meanCortS)
title(' Brief movement cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-movement time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
clim([-100,100])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S3e] moderate movements cortical LFP
% ax5 = subplot(6,3,5);
% imagesc(data.IntermediateMovement.meanCortT,data.IntermediateMovement.meanCortF,data.IntermediateMovement.meanCortS)
% title(' Moderate movement cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-movement time (s)')
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,100])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
%% [1-S3f] extended movements cortical LFP
% ax6 = subplot(6,3,6);
% imagesc(data.LongMovement.meanCortT,data.LongMovement.meanCortF,data.LongMovement.meanCortS)
% title('Extended movement cortical LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-movement time (s)')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% clim([-100,100])
% set(gca,'Ticklength',[0,0])
% axis square
% axis xy
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];

%% [1-S3m] brief movements Rhodamine
ax13 = subplot(6,3,13);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine + data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine - data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Brief movement Blood Volume (Z)')
ylabel('Zscored \DeltaF/F Ach')
ax13.YLim = [-2 6];

yyaxis right
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine + data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine - data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax13.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax13.YAxis(2).Color = [0.4660 0.6740 0.1880];
xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
ax13.YLim = [-2 6];
%% [1-S3n] moderate movements Rhodamine
% ax14 = subplot(6,3,14);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine + data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine - data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% title('Moderate movement Blood Volume (Z)')
% ylabel('Zscored \DeltaF/F Ach')
% ax14.YLim = [-2 6];
% 
% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine + data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine - data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F NE')
% ax14.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax14.YAxis(2).Color = [0.4660 0.6740 0.1880];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax14.TickLength = [0.03,0.03];
% ax14.YLim = [-2 6];
 %% [1-S3o] extended movements Rhodamine
% ax15 = subplot(6,3,15);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine + data.LongMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine - data.LongMovement.Z_AchstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F Ach')
% ax15.YLim = [-2 6];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine + data.LongMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine - data.LongMovement.Z_NEstdRhodamine,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F NE')
% ax15.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax15.YAxis(2).Color = [0.4660 0.6740 0.1880];
% title('Extended movement Blood Volume (Z)')
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax15.TickLength = [0.03,0.03];
% ax15.YLim = [-2 6];
%% [1-S3p] brief movements GFP
ax16 = subplot(6,3,16);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP + data.ShortMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP - data.ShortMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('Zscored \DeltaF/F Ach')
ax16.YLim = [-2 6];

yyaxis right
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP + data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP - data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('Zscored \DeltaF/F NE')
ax16.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax16.YAxis(2).Color = [0.4660 0.6740 0.1880];
title('Brief movement GFP')
xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
ax16.YLim = [-2 6];
%% [1-S3q] moderate movements GFP
% ax17 = subplot(6,3,17);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP + data.IntermediateMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP - data.IntermediateMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F Ach')
% ax17.YLim = [-2 6];
% 
% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP + data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP - data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F NE')
% ax17.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax17.YAxis(2).Color = [0.4660 0.6740 0.1880];
% title('Moderate movement GFP')
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax17.TickLength = [0.03,0.03];
% ax17.YLim = [-2 6];
 %% [1-S3r] extended movements GFP
% ax18 = subplot(6,3,18);
% 
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP + data.LongMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP - data.LongMovement.Z_AchstdGFP,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F Ach')
% ax18.YLim = [-2 6];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP + data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP - data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('Zscored \DeltaF/F NE')
% ax18.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax18.YAxis(2).Color = [0.4660 0.6740 0.1880];
% title('Extended movement GFP')
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax18.TickLength = [0.03,0.03];
% ax18.YLim = [-2 6];
%% axes positions
linkaxes([ax1],'xy') %,ax2,ax3
linkaxes([ax4],'xy') %,ax5,ax6
linkaxes([ax13],'xy') %,ax14,ax15
linkaxes([ax16],'xy') %,ax17,ax18
ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
% ax5Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3:4) = ax3Pos(3:4);
% ax10Pos(3:4) = ax1Pos(3:4);
% ax11Pos(3:4) = ax2Pos(3:4);
% ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
% set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
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
    savefig(summaryFigure,[dirpath 'Movement']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Movement'])
    close 
end
%% Only the fiber signals
summaryFigureN = figure('Name','Movement_fiber');
sgtitle('Movement evoked responses in fiber photometry signals')

% Ach ShortMovement
ax1 = subplot(3,3,1);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine + data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine - data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Brief movement Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax1.YLim = [-2 6];

yyaxis right
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP + data.ShortMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP - data.ShortMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')

ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = [0 0.4470 0.7410];

xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax1.YLim = [-2 6];
xlim([-5 15])

% Ach IntermediateMovement
% ax2 = subplot(3,3,2);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine + data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine - data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% title('Moderate movement Ach fiber')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax2.YLim = [-2 6];
% 
% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP + data.IntermediateMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP - data.IntermediateMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% ylabel('\DeltaF/F GRAB ACh (Z)')
% ax2.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax2.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ax2.YLim = [-2 6];
% xlim([-5 15])

% Ach LongMovement
% ax3 = subplot(3,3,3);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine + data.LongMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine - data.LongMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% title('Long movement Ach fiber')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax3.YLim = [-2 6];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP + data.LongMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP - data.LongMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% ylabel('\DeltaF/F GRAB ACh (Z)')
% ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax3.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ax3.YLim = [-2 6];
% xlim([-5 15])

% NE ShortMovement
ax4 = subplot(3,3,4);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine + data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine - data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
title('Brief movement NE fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax4.YLim = [-2 6];

yyaxis right
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP + data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP - data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax4.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax4.YAxis(2).Color = [0.4660 0.6740 0.1880];
xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
ax4.YLim = [-2 6];
xlim([-5 15])

% NE IntermediateMovement
% ax5 = subplot(3,3,5);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine + data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine - data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% title('Moderate movement NE fibers')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax5.YLim = [-2 6];
% 
% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP + data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP - data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('\DeltaF/F GRABNE (Z)')
% ax5.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax5.YAxis(2).Color = [0.4660 0.6740 0.1880];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax5.TickLength = [0.03,0.03];
% ax5.YLim = [-2 6];
% xlim([-5 15])

% NE LongMovement
% ax6 = subplot(3,3,6);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine + data.LongMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine - data.LongMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% title('Long movement NE fibers')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax6.YLim = [-2 6];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP + data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP - data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('\DeltaF/F GRABNE (Z)')
% ax6.YAxis(1).Color = [0.6350 0.0780 0.1840];
% ax6.YAxis(2).Color = [0.4660 0.6740 0.1880];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax6.TickLength = [0.03,0.03];
% ax6.YLim = [-2 6];
% xlim([-5 15])

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
    savefig(summaryFigureN,[dirpath 'Movement-FiberSignals']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Movement-FiberSignals'])
    close 
end
%}
%% Only the fiber signals consolidated
summaryFigureN = figure('Name','Movement_fiber_consolidated');
sgtitle('Movement evoked responses in fiber photometry signals')

% Short Movement
ax1 = subplot(1,3,1);
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine + data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanRhodamine - data.ShortMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
title('Brief movement')

plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine + data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanRhodamine - data.ShortMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('\DeltaF/F Blood Volume (Z)')
ax1.YLim = [-2 6];

yyaxis right
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP + data.ShortMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_AchmeanGFP - data.ShortMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
%
plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP + data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.Z_NEmeanGFP - data.ShortMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

ylabel('\DeltaF/F GRAB(Z)')

ax1.YAxis(1).Color = 'k';%[0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = 'k';%[0 0.4470 0.7410];
legend('CBV-LH','CBV-RH', 'ACh', 'NE')
xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax1.YLim = [-2 6];
xlim([-5 15])


% Intermediate Movement
% ax2 = subplot(1,3,2);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine + data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanRhodamine - data.IntermediateMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine + data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanRhodamine - data.IntermediateMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% 
% title('Moderate movement')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax2.YLim = [-2 6];
% 
% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP + data.IntermediateMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_AchmeanGFP - data.IntermediateMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP + data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.Z_NEmeanGFP - data.IntermediateMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('\DeltaF/F GRAB(Z)')
% ax2.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax2.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ax2.YLim = [-2 6];
% xlim([-5 15])

% Long Movement
% ax3 = subplot(1,3,3);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine + data.LongMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanRhodamine - data.LongMovement.Z_AchstdRhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine + data.LongMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanRhodamine - data.LongMovement.Z_NEstdRhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% 
% title('Long movement')
% ylabel('\DeltaF/F Blood Volume (Z)')
% ax3.YLim = [-2 6];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP + data.LongMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_AchmeanGFP - data.LongMovement.Z_AchstdGFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP + data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.Z_NEmeanGFP - data.LongMovement.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% 
% ylabel('\DeltaF/F GRAB(Z)')
% ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax3.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ax3.YLim = [-2 6];
% xlim([-5 15])

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
    savefig(summaryFigureN,[dirpath 'Movement-FiberSignals-consolidated']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Movement-FiberSignals-consolidated'])
    close 
end
end
