function [AnalysisResults] = Fig4_FP_Bilateral_No_REM(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T282','T285'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
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
        data.(transition).Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).Rhodamine;
        data.(transition).GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GCaMP7s;

        data.(transition).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;
        data.(transition).LH_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).LH_Rhodamine;
        data.(transition).LH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s;
        data.(transition).RH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).RH_Cort;
        data.(transition).RH_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_Rhodamine;
        data.(transition).RH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s;
    end
end
% take average for each behavioral transition _RH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;

    data.(transition).meanRhodamine = mean(data.(transition).Rhodamine,1);
    data.(transition).stdRhodamine = std(data.(transition).Rhodamine,0,1);
    data.(transition).meanGCaMP7s = mean(data.(transition).GCaMP7s,1);
    data.(transition).stdGCaMP7s = std(data.(transition).GCaMP7s,0,1);
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;

    data.(transition).mean_LH_Rhodamine = mean(data.(transition).LH_Rhodamine,1);
    data.(transition).std_LH_Rhodamine = std(data.(transition).LH_Rhodamine,0,1);
    data.(transition).mean_LH_GCaMP7s = mean(data.(transition).LH_GCaMP7s,1);
    data.(transition).std_LH_GCaMP7s = std(data.(transition).LH_GCaMP7s,0,1);
    data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;

    data.(transition).mean_RH_Rhodamine = mean(data.(transition).RH_Rhodamine,1);
    data.(transition).std_RH_Rhodamine = std(data.(transition).RH_Rhodamine,0,1);
    data.(transition).mean_RH_GCaMP7s = mean(data.(transition).RH_GCaMP7s,1);
    data.(transition).std_RH_GCaMP7s = std(data.(transition).RH_GCaMP7s,0,1);
    data.(transition).mean_RH_Cort = mean(data.(transition).RH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;

%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs
summaryFigure_A = figure('Name','Fig4 A');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(3,2,1);
% Rhodamine and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.meanRhodamine + data.AWAKEtoNREM.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanRhodamine - data.AWAKEtoNREM.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'Rhodamine','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(3,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(3,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
% NREM to Awake
ax4 = subplot(3,2,2);
% Rhodamine and EMG
plot(T1,data.NREMtoAWAKE.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanRhodamine + data.NREMtoAWAKE.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanRhodamine - data.NREMtoAWAKE.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(3,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(3,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
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
close

%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs Bilateral
summaryFigure_A_bilateral = figure('Name','Fig4 A');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(4,2,1);
% Rhodamine and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_LH_Rhodamine + data.AWAKEtoNREM.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_LH_Rhodamine - data.AWAKEtoNREM.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

p2 = plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine + data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine - data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort LH neural
ax2 = subplot(4,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(4,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% hippocampal neural
ax4 = subplot(4,2,7);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(4,2,2);

% Rhodamine and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_LH_Rhodamine + data.NREMtoAWAKE.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_LH_Rhodamine - data.NREMtoAWAKE.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

p2 = plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine + data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine - data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG','Location','best','FontSize',5)
ax5.YAxis(1).Color = colors('dark candy apple red');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% LH cort neural
ax6 = subplot(4,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(4,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];

% hippocampal neural
ax8 = subplot(4,2,8);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);

ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
ax8Pos(3:4) = ax5Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A_bilateral,[dirpath 'Fig4_A_Bilateral']);
    set(summaryFigure_A_bilateral,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_A_Bilateral'])
end
close
%% Fig. 4C
summaryFigure_C = figure('Name','Fig4 GCaMP7s');
sgtitle('GCaMP7s Transition with Cortical LFPs')
%  Awake to NREM
ax1 = subplot(3,2,1);
% GCaMP7s and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanGCaMP7s + data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanGCaMP7s - data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
ylabel(' \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'GCaMP7s','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('alizarin crimson');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(3,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(3,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')%caxis
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
% NREM to Awake
ax4 = subplot(3,2,2);
% GCaMP7s and EMG
plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanGCaMP7s + data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanGCaMP7s - data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('alizarin crimson');
ax4.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(3,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(3,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_C,[dirpath 'Fig4_C']);
    saveas(summaryFigure_C,[dirpath 'Fig4_C.tif']);
    set(summaryFigure_C,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_C'])
end
close
%% Fig. 4 Plot GCaMP7s Transition with Cortical LFPs Bilateral
summaryFigure_C_bilateral = figure('Name','Fig4 A');
sgtitle('GCaMP7s Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(4,2,1);
% GCaMP7s and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s + data.AWAKEtoNREM.std_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s - data.AWAKEtoNREM.std_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);

p2 = plot(T1,data.AWAKEtoNREM.mean_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_RH_GCaMP7s + data.AWAKEtoNREM.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_RH_GCaMP7s - data.AWAKEtoNREM.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
ylabel('\DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color','k','LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('alizarin crimson');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort LH neural
ax2 = subplot(4,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(4,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% hippocampal neural
ax4 = subplot(4,2,7);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(4,2,2);

% GCaMP7s and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s + data.NREMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s - data.NREMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('alizarin crimson'),'LineWidth',0.25);

p2 = plot(T1,data.NREMtoAWAKE.mean_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_RH_GCaMP7s + data.NREMtoAWAKE.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_RH_GCaMP7s - data.NREMtoAWAKE.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
ylabel('\DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color','k','LineWidth',0.25);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax5.YAxis(1).Color = colors('alizarin crimson');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% LH cort neural
ax6 = subplot(4,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(4,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];

% hippocampal neural
ax8 = subplot(4,2,8);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);

ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
ax8Pos(3:4) = ax5Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_C_bilateral,[dirpath 'Fig4_C_Bilateral']);
    set(summaryFigure_C_bilateral,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_C_Bilateral'])
end
close
end
