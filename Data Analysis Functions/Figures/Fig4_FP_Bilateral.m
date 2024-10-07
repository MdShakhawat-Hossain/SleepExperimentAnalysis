function [AnalysisResults] = Fig4_FP_Bilateral(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T282','T285'};
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
%% Fig. 4 Plot Rhodamine GCaMP7s transitions
summaryFigure_B = figure('Name','Fig4 GCaMP7-Rhodamine');
sgtitle('Rhodamine GCaMP7s transitions')
%NREM to REM transition
ax1 = subplot(2,2,2);
p1 = plot(T1,data.NREMtoREM.meanRhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p2 = plot(T1,data.NREMtoREM.meanGCaMP7s,'color',colors('deep jungle green'),'LineWidth',1);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'Rhodamine','GCaMP7s','Location','best','FontSize',5)
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('deep jungle green');
% REM to Awake transition
ax2 = subplot(2,2,4);
plot(T1,data.REMtoAWAKE.meanRhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.REMtoAWAKE.meanGCaMP7s,'color',colors('deep jungle green'),'LineWidth',1)
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YAxis(1).Color = colors('dark candy apple red');
ax2.YAxis(2).Color = colors('deep jungle green');
%  NREM to Awake transition
ax3 = subplot(2,2,3);
plot(T1,data.NREMtoAWAKE.meanRhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'color',colors('deep jungle green'),'LineWidth',1)
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YAxis(1).Color = colors('dark candy apple red');
ax3.YAxis(2).Color = colors('deep jungle green');
% Awake to NREM transition
ax3 = subplot(2,2,1);
plot(T1,data.AWAKEtoNREM.meanRhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([0,80])
yyaxis right
plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'color',colors('deep jungle green'),'LineWidth',1)
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
xlim([-30,30])
% ylim([0,42])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YAxis(1).Color = colors('dark candy apple red');
ax3.YAxis(2).Color = colors('deep jungle green');
% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_B,[dirpath 'Fig4_B']);
    set(summaryFigure_B,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_B'])
end
close
%% Fig. 4 Plot bilateral GCaMP7 and Rhodamine Transition
summaryFigure_B_bilateral = figure('Name','Fig4 GCaMP7-Rhodamine Bilateral');
sgtitle('Plot bilateral GCaMP7 and Rhodamine Transition')
% NREM to REM transition
ax1 = subplot(2,2,2);
p1 = plot(T1,data.NREMtoREM.mean_LH_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
hold on
p2 = plot(T1,data.NREMtoREM.mean_RH_Rhodamine,'color',colors('indian red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p3 = plot(T1,data.NREMtoREM.mean_LH_GCaMP7s,'color',colors('british racing green'),'LineWidth',1);
hold on
p4 = plot(T1,data.NREMtoREM.mean_RH_GCaMP7s,'color','g','LineWidth',1);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2,p3,p4],'LH Rhodamine','RH Rhodamine','LH GCaMP7s','RH GCaMP7s','Location','best','FontSize',5)
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = colors('deep jungle green');
% REM to Awake transition
ax2 = subplot(2,2,4);
p1 = plot(T1,data.REMtoAWAKE.mean_LH_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
hold on
p2 = plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine,'color',colors('indian red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p3 = plot(T1,data.REMtoAWAKE.mean_LH_GCaMP7s,'color',colors('british racing green'),'LineWidth',1);

hold on
p4 = plot(T1,data.REMtoAWAKE.mean_RH_GCaMP7s,'color','g','LineWidth',1);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% legend([p1,p2,p3,p4],'LH_Rhodamine','RH_Rhodamine','LH_GCaMP7s','RH_GCaMP7s','Location','best','FontSize',5)
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YAxis(1).Color = colors('dark candy apple red');
ax2.YAxis(2).Color = colors('deep jungle green');
%  NREM to Awake transition
ax3 = subplot(2,2,3);
p1 = plot(T1,data.NREMtoAWAKE.mean_LH_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
hold on
p2 = plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine,'color',colors('indian red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p3 = plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s,'color',colors('british racing green'),'LineWidth',1);
hold on
p4 = plot(T1,data.NREMtoAWAKE.mean_RH_GCaMP7s,'color','g','LineWidth',1);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% legend([p1,p2,p3,p4],'LH_Rhodamine','RH_Rhodamine','LH_GCaMP7s','RH_GCaMP7s','Location','NorthWest')
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YAxis(1).Color = colors('dark candy apple red');
ax3.YAxis(2).Color = colors('deep jungle green');
%   Awake to NREM transition
ax4 = subplot(2,2,1);
p1 = plot(T1,data.AWAKEtoNREM.mean_LH_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
hold on
p2 = plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine,'color',colors('indian red'),'LineWidth',1);
ylabel('\DeltaF/F')
% ylim([40,80])
yyaxis right
p3 = plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s,'color',colors('british racing green'),'LineWidth',1);
hold on
p4 = plot(T1,data.AWAKEtoNREM.mean_RH_GCaMP7s,'color','g','LineWidth',1);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% legend([p1,p2,p3,p4],'LH_Rhodamine','RH_Rhodamine','LH_GCaMP7s','RH_GCaMP7s','Location','best','FontSize',5)
xlim([-30,30])
% ylim([5,42])
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = colors('deep jungle green');
% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_B_bilateral,[dirpath 'Fig4_B_bilateral']);
    set(summaryFigure_B_bilateral,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_B_bilateral'])
end
close
%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs
summaryFigure_A = figure('Name','Fig4 A');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(6,2,1);
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

% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
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
ax2 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
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
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
% NREM to Awake
ax4 = subplot(6,2,2);
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
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('dark candy apple red');
ax4.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
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
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
% NREM to REM
ax7 = subplot(6,2,7);
% Rhodamine and EMG
plot(T1,data.NREMtoREM.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanRhodamine + data.NREMtoREM.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.meanRhodamine - data.NREMtoREM.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
ylabel(' \DeltaRhodamine')
xlim([-30,30])
% ylim([35,90])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('dark candy apple red');
ax7.YAxis(2).Color = 'k';
% ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-100,300])
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
%caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
% REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanRhodamine + data.REMtoAWAKE.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.meanRhodamine - data.REMtoAWAKE.stdRhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
ylabel('  \DeltaRhodamine')
xlim([-30,30])
% ylim([0,90])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('dark candy apple red');
ax10.YAxis(2).Color = 'k';
% ylim([-2,1])
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-100,300])
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
%caxis([-100,300])
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
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);


ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax10Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax7Pos(3:4);
ax12Pos(3:4) = ax7Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax10,'position',ax10Pos);
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
close

%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs Bilateral
summaryFigure_A = figure('Name','Fig4 A');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(8,2,1);
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
% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
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
ax2 = subplot(8,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(8,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% hippocampal neural
ax4 = subplot(8,2,7);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(8,2,2);

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
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% LH cort neural
ax6 = subplot(8,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(8,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];

% hippocampal neural
ax8 = subplot(8,2,8);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];

% NREM to REM
ax9 = subplot(8,2,9);
% Rhodamine and EMG
p1 = plot(T1,data.NREMtoREM.mean_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_LH_Rhodamine + data.NREMtoREM.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_LH_Rhodamine - data.NREMtoREM.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

p2 = plot(T1,data.NREMtoREM.mean_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_RH_Rhodamine + data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_RH_Rhodamine - data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG','Location','best','FontSize',5)
ax9.YAxis(1).Color = colors('dark candy apple red');
ax9.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax9.TickLength = [0.03,0.03];

% cort LH neural
ax10 = subplot(8,2,11);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];

% cort RH neural
ax11 = subplot(8,2,13);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];

% hippocampal neural
ax12 = subplot(8,2,15);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];

% REM to Awake
ax13 = subplot(8,2,10);

% Rhodamine and EMG
p1 = plot(T1,data.REMtoAWAKE.mean_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_LH_Rhodamine + data.REMtoAWAKE.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_LH_Rhodamine - data.REMtoAWAKE.std_LH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

p2 = plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine + data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine - data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.25);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.REMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('REM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH Rhodamine','RH Rhodamine','EMG','Location','best','FontSize',5)
ax13.YAxis(1).Color = colors('dark candy apple red');
ax13.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax13.TickLength = [0.03,0.03];

% LH cort neural
ax14 = subplot(8,2,12);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];

% RH cort neural
ax15 = subplot(8,2,14);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];

% hippocampal neural
ax16 = subplot(8,2,16);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
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
ax13Pos = get(ax13,'position');
ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
ax16Pos = get(ax16,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);

ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
ax8Pos(3:4) = ax5Pos(3:4);

ax10Pos(3:4) = ax9Pos(3:4);
ax11Pos(3:4) = ax9Pos(3:4);
ax12Pos(3:4) = ax9Pos(3:4);

ax14Pos(3:4) = ax13Pos(3:4);
ax15Pos(3:4) = ax13Pos(3:4);
ax16Pos(3:4) = ax13Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
set(ax16,'position',ax16Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A,[dirpath 'Fig4_A_Bilateral']);
    set(summaryFigure_A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_A_Bilateral'])
end
close
%% Fig. 4C
summaryFigure_C = figure('Name','Fig4 GCaMP7s');
sgtitle('GCaMP7s Transition with Cortical LFPs')
%  Awake to NREM
ax1 = subplot(6,2,1);
% GCaMP7s and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanGCaMP7s + data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.meanGCaMP7s - data.AWAKEtoNREM.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
ylabel(' \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'GCaMP7s','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('british racing green');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
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
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')%caxis
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
% NREM to Awake
ax4 = subplot(6,2,2);
% GCaMP7s and EMG
plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanGCaMP7s + data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.meanGCaMP7s - data.NREMtoAWAKE.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors('british racing green');
ax4.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
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
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
% NREM to REM
ax7 = subplot(6,2,7);
% GCaMP7s and EMG
plot(T1,data.NREMtoREM.meanGCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanGCaMP7s + data.NREMtoREM.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.meanGCaMP7s - data.NREMtoREM.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([35,90])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors('british racing green');
ax7.YAxis(2).Color = 'k';
% ylim([-2,-0.5])
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% %caxis([-100,300])
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
% %caxis([-100,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
% REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanGCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanGCaMP7s + data.REMtoAWAKE.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.meanGCaMP7s - data.REMtoAWAKE.stdGCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
ylabel('  \DeltaGCaMP7s')
xlim([-30,30])
% ylim([0,90])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors('british racing green');
ax10.YAxis(2).Color = 'k';
% ylim([-2,1])
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% %caxis([-100,300])
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
% %caxis([-100,300])
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
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);


ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax10Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax7Pos(3:4);
ax12Pos(3:4) = ax7Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax10,'position',ax10Pos);
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
close
%% Fig. 4 Plot GCaMP7s Transition with Cortical LFPs Bilateral
summaryFigure_C_bilateral = figure('Name','Fig4 A');
sgtitle('GCaMP7s Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(8,2,1);
% GCaMP7s and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s + data.AWAKEtoNREM.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_LH_GCaMP7s - data.AWAKEtoNREM.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);

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
% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('british racing green');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort LH neural
ax2 = subplot(8,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(8,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% hippocampal neural
ax4 = subplot(8,2,7);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(8,2,2);

% GCaMP7s and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s + data.NREMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_LH_GCaMP7s - data.NREMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);

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
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax5.YAxis(1).Color = colors('british racing green');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% LH cort neural
ax6 = subplot(8,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(8,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];

% hippocampal neural
ax8 = subplot(8,2,8);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];

% NREM to REM
ax9 = subplot(8,2,9);
% GCaMP7s and EMG
p1 = plot(T1,data.NREMtoREM.mean_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_LH_GCaMP7s + data.NREMtoREM.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_LH_GCaMP7s - data.NREMtoREM.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);

p2 = plot(T1,data.NREMtoREM.mean_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_RH_GCaMP7s + data.NREMtoREM.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_RH_GCaMP7s - data.NREMtoREM.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
ylabel('\DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax9.YAxis(1).Color = colors('british racing green');
ax9.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax9.TickLength = [0.03,0.03];

% cort LH neural
ax10 = subplot(8,2,11);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];

% cort RH neural
ax11 = subplot(8,2,13);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_RH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'RH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];

% hippocampal neural
ax12 = subplot(8,2,15);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];

% REM to Awake
ax13 = subplot(8,2,10);

% GCaMP7s and EMG
p1 = plot(T1,data.REMtoAWAKE.mean_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_LH_GCaMP7s + data.REMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_LH_GCaMP7s - data.REMtoAWAKE.std_LH_GCaMP7s,'-','color',colors('british racing green'),'LineWidth',0.25);

p2 = plot(T1,data.REMtoAWAKE.mean_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_RH_GCaMP7s + data.REMtoAWAKE.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_RH_GCaMP7s - data.REMtoAWAKE.std_RH_GCaMP7s,'-','color',colors('deep jungle green'),'LineWidth',0.25);
ylabel('\DeltaGCaMP7s')
xlim([-30,30])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.REMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
% plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.25);
title('REM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'LH GCaMP7s','RH GCaMP7s','EMG','Location','best','FontSize',5)
ax13.YAxis(1).Color = colors('british racing green');
ax13.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax13.TickLength = [0.03,0.03];

% LH cort neural
ax14 = subplot(8,2,12);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];

% RH cort neural
ax15 = subplot(8,2,14);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_RH_Cort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];

% hippocampal neural
ax16 = subplot(8,2,16);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-100,200])
xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
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
ax13Pos = get(ax13,'position');
ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
ax16Pos = get(ax16,'position');

ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax4Pos(3:4) = ax1Pos(3:4);

ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
ax8Pos(3:4) = ax5Pos(3:4);

ax10Pos(3:4) = ax9Pos(3:4);
ax11Pos(3:4) = ax9Pos(3:4);
ax12Pos(3:4) = ax9Pos(3:4);

ax14Pos(3:4) = ax13Pos(3:4);
ax15Pos(3:4) = ax13Pos(3:4);
ax16Pos(3:4) = ax13Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
set(ax16,'position',ax16Pos);
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
