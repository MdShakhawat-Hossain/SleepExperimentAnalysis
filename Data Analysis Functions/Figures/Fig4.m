function [AnalysisResults] = Fig4(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T282'};
% TwoP_animalIDs = {'T115','T117','T118','T125','T126'};   % T116 has no REM events
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(transition).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
        data.(transition).Hip(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Hip;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Cort;
        data.(transition).TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).TRITC;
        data.(transition).GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GCaMP7s;

    end
end
% take average for each behavioral transition
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);
    data.(transition).meanTRITC = mean(data.(transition).TRITC,1);
    data.(transition).stdTRITC = std(data.(transition).TRITC,0,1);
    data.(transition).meanGCaMP7s = mean(data.(transition).GCaMP7s,1);
    data.(transition).stdGCaMP7s = std(data.(transition).GCaMP7s,0,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;

%% Fig. 4 (part two)
summaryFigure_B = figure('Name','Fig4 (f-g)');
sgtitle('Figure panel 4 - Turner et al. 2020')
%% [4f] NREM to REM transition
ax1 = subplot(2,2,2);
p1 = plot(T1,data.NREMtoREM.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p2 = plot(T1,data.NREMtoREM.meanGCaMP7s,'color',colors('rich black'),'LineWidth',2);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'TRITC','GCaMP7s','Location','NorthWest')
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('rich black');
%% [4g] REM to Awake transition
ax2 = subplot(2,2,4);
plot(T1,data.REMtoAWAKE.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.REMtoAWAKE.meanGCaMP7s,'color',colors('rich black'),'LineWidth',2)
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YAxis(1).Color = colors('dark candy apple red');
ax2.YAxis(2).Color = colors('rich black');
%%  NREM to Awake transition
ax3 = subplot(2,2,3);
plot(T1,data.NREMtoAWAKE.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'color',colors('rich black'),'LineWidth',2)
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YAxis(1).Color = colors('dark candy apple red');
ax3.YAxis(2).Color = colors('rich black');
%%   Awake to NREM transition
ax3 = subplot(2,2,1);
plot(T1,data.AWAKEtoNREM.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'color',colors('rich black'),'LineWidth',2)
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YAxis(1).Color = colors('dark candy apple red');
ax3.YAxis(2).Color = colors('rich black');
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_B,[dirpath 'Fig4_B']);
    set(summaryFigure_B,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_B'])
end
%% Fig. 4
summaryFigure_A = figure('Name','Fig4 (b-e)');
sgtitle('Figure 4 - Turner et al. 2020')
%% [4b] Awake to NREM
ax1 = subplot(6,2,1);
% TRITC and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanTRITC + data.AWAKEtoNREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanTRITC - data.AWAKEtoNREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaTRITC')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'TRITC','EMG')
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('rich black');
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(6,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [4c] NREM to Awake
ax4 = subplot(6,2,2);
% TRITC and EMG
plot(T1,data.NREMtoAWAKE.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanTRITC + data.NREMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanTRITC - data.NREMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaTRITC')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = colors('rich black');
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(6,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [4d] NREM to REM
ax7 = subplot(6,2,7);
% TRITC and EMG
plot(T1,data.NREMtoREM.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanTRITC + data.NREMtoREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanTRITC - data.NREMtoREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel(' \DeltaTRITC')
xlim([-30,30])
% ylim([35,90])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('dark candy apple red');
ax7.YAxis(2).Color = colors('rich black');
% ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
% hippocampal neural
ax9 = subplot(6,2,11);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [4e] REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanTRITC + data.REMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanTRITC - data.REMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('  \DeltaTRITC')
xlim([-30,30])
% ylim([0,90])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('dark candy apple red');
ax10.YAxis(2).Color = colors('rich black');
% ylim([-2,1])
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
% hippocampal neural
ax12 = subplot(6,2,12);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax4Pos(3:4);
ax6Pos(3:4) = ax4Pos(3:4);
ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax10Pos(3:4);
ax12Pos(3:4) = ax10Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A,[dirpath 'Fig4_A']);
    set(summaryFigure_A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_A'])
end
%% Fig. 4
summaryFigure_C = figure('Name','Fig4 (b-e)');
sgtitle('Figure 4 - Turner et al. 2020')
%% [4b] Awake to NREM
ax1 = subplot(6,2,1);
% GCaMP7s and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanGCaMP7s + data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanGCaMP7s - data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel(' \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'GCaMP7s','EMG')
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('rich black');
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(6,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [4c] NREM to Awake
ax4 = subplot(6,2,2);
% GCaMP7s and EMG
plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanGCaMP7s + data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanGCaMP7s - data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = colors('rich black');
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(6,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [4d] NREM to REM
ax7 = subplot(6,2,7);
% GCaMP7s and EMG
plot(T1,data.NREMtoREM.meanGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanGCaMP7s + data.NREMtoREM.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanGCaMP7s - data.NREMtoREM.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([35,90])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('dark candy apple red');
ax7.YAxis(2).Color = colors('rich black');
% ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
% hippocampal neural
ax9 = subplot(6,2,11);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [4e] REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanGCaMP7s + data.REMtoAWAKE.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanGCaMP7s - data.REMtoAWAKE.stdGCaMP7s,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([0,90])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('rich black'),'LineWidth',0.5);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('dark candy apple red');
ax10.YAxis(2).Color = colors('rich black');
% ylim([-2,1])
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
% hippocampal neural
ax12 = subplot(6,2,12);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,300])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax4Pos(3:4);
ax6Pos(3:4) = ax4Pos(3:4);
ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax10Pos(3:4);
ax12Pos(3:4) = ax10Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_C,[dirpath 'Fig4_C']);
    set(summaryFigure_C,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_C'])
end
end
