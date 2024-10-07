function [AnalysisResults] = Fig4_FP_Transition_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'SHF025'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(transition).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).RH_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_Rhodamine;
        data.(transition).RH_GFP(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_GFP;
        data.(transition).LH_EEG(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).LH_EEG;
    end
end
% take average for each behavioral transition _RH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);

    data.(transition).mean_RH_Rhodamine = mean(data.(transition).RH_Rhodamine,1);
    data.(transition).std_RH_Rhodamine = std(data.(transition).RH_Rhodamine,0,1);
    data.(transition).mean_RH_GFP = mean(data.(transition).RH_GFP,1);
    data.(transition).std_RH_GFP = std(data.(transition).RH_GFP,0,1);

    data.(transition).mean_LH_EEG = mean(data.(transition).LH_EEG,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;

%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs Bilateral
summaryFigure_A = figure('Name','Fig4_Transition Rhodamine');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(6,2,1);
% Rhodamine and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine + data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine - data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

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
legend([p1,p3],'RH Rhodamine','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH  EEG LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];


% NREM to Awake
ax5 = subplot(6,2,2);

% Rhodamine and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine + data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine - data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

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
% legend([p1,p2,p3],'RH Rhodamine','NE Rhodamine','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_EEG,'y')
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


% NREM to REM
ax9 = subplot(6,2,7);
% Rhodamine and EMG
p1 = plot(T1,data.NREMtoREM.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_RH_Rhodamine + data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_RH_Rhodamine - data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

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
% legend([p1,p2,p3],'RH Rhodamine','NE Rhodamine','EMG','Location','best','FontSize',5)
ax9.YAxis(1).Color = colors('dark candy apple red');
ax9.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax9.TickLength = [0.03,0.03];

% cort RH neural
ax11 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH  EEG LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];


% REM to Awake
ax13 = subplot(6,2,8);

% Rhodamine and EMG
p1 = plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine + data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine - data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.25);

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
% legend([p1,p2,p3],'RH Rhodamine','NE Rhodamine','EMG','Location','best','FontSize',5)
ax13.YAxis(1).Color = colors('dark candy apple red');
ax13.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax13.TickLength = [0.03,0.03];

% RH cort neural
ax15 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_EEG,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];


%%%% axes positionns
ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax13Pos = get(ax13,'position');
% ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
% ax16Pos = get(ax16,'position');

% ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax1Pos(3:4);

% ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
% ax8Pos(3:4) = ax5Pos(3:4);

% ax10Pos(3:4) = ax9Pos(3:4);
ax11Pos(3:4) = ax9Pos(3:4);
% ax12Pos(3:4) = ax9Pos(3:4);

% ax14Pos(3:4) = ax13Pos(3:4);
ax15Pos(3:4) = ax13Pos(3:4);
% ax16Pos(3:4) = ax13Pos(3:4);

% set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
% set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
% set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
% set(ax16,'position',ax16Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A,[dirpath 'Fig4_Transition_Rhodamine']);
    set(summaryFigure_A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_Transition_Rhodamine'])
end
close

%% Fig. 4 Plot GFP Transition with Cortical LFPs Bilateral
summaryFigure_C = figure('Name','Fig4_Transition GFP');
sgtitle('GFP Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(6,2,1);
% GFP and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_RH_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_RH_GFP + data.AWAKEtoNREM.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.AWAKEtoNREM.mean_RH_GFP - data.AWAKEtoNREM.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);

ylabel('\DeltaGFP')
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
legend([p1,p3],'RH GFP','EMG','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('indian red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(6,2,3);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH  EEG LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];


% NREM to Awake
ax5 = subplot(6,2,2);

% GFP and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_RH_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_RH_GFP + data.NREMtoAWAKE.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.NREMtoAWAKE.mean_RH_GFP - data.NREMtoAWAKE.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);

ylabel('\DeltaGFP')
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
% legend([p1,p2,p3],'RH GFP','GRAB NE','EMG','Location','best','FontSize',5)
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% RH cort neural
ax7 = subplot(6,2,4);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_EEG,'y')
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

% NREM to REM
ax9 = subplot(6,2,7);
% GFP and EMG
p1 = plot(T1,data.NREMtoREM.mean_RH_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_RH_GFP + data.NREMtoREM.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.NREMtoREM.mean_RH_GFP - data.NREMtoREM.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);

ylabel('\DeltaGFP')
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
% legend([p1,p2,p3],'RH GFP','GRAB NE','EMG','Location','best','FontSize',5)
ax9.YAxis(1).Color = colors('indian red');
ax9.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax9.TickLength = [0.03,0.03];

% cort RH neural
ax11 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH  EEG LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];


% REM to Awake
ax13 = subplot(6,2,8);

% GFP and EMG
p1 = plot(T1,data.REMtoAWAKE.mean_RH_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_RH_GFP + data.REMtoAWAKE.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);
plot(T1,data.REMtoAWAKE.mean_RH_GFP - data.REMtoAWAKE.std_RH_GFP,'-','color',colors('indian red'),'LineWidth',0.25);

ylabel('\DeltaGFP')
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
% legend([p1,p2,p3],'RH GFP','GRAB NE','EMG','Location','best','FontSize',5)
ax13.YAxis(1).Color = colors('indian red');
ax13.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax13.TickLength = [0.03,0.03];

% RH cort neural
ax15 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_EEG,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];

%% axes positionns
ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax13Pos = get(ax13,'position');
% ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
% ax16Pos = get(ax16,'position');

% ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax1Pos(3:4);

% ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
% ax8Pos(3:4) = ax5Pos(3:4);

% ax10Pos(3:4) = ax9Pos(3:4);
ax11Pos(3:4) = ax9Pos(3:4);
% ax12Pos(3:4) = ax9Pos(3:4);

% ax14Pos(3:4) = ax13Pos(3:4);
ax15Pos(3:4) = ax13Pos(3:4);
% ax16Pos(3:4) = ax13Pos(3:4);

% set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
% set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
% set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
% set(ax16,'position',ax16Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_C,[dirpath 'Fig4_Transition_GFP']);
    set(summaryFigure_C,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_Transition_GFP'])
end
close

%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs Rhodamine');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(8,2,1);
% Rhodamine and GFP
p1 = plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine + data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_RH_Rhodamine - data.AWAKEtoNREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta RH Rhodamine')
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.mean_RH_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_RH_GFP + data.AWAKEtoNREM.std_RH_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_RH_GFP - data.AWAKEtoNREM.std_RH_GFP,'-','color','k','LineWidth',0.1);
ylabel('\Delta RH GFP','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'RH Rhodamine','RH GFP','Location','best','FontSize',5)

title('Awake to NREM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
ax1.TickLength = [0.03,0.03];

% cort RH neural
ax3 = subplot(8,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(8,2,2);
% Rhodamine and GFP
p1 = plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine + data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_RH_Rhodamine - data.NREMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta RH Rhodamine')
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.mean_RH_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_RH_GFP + data.NREMtoAWAKE.std_RH_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_RH_GFP - data.NREMtoAWAKE.std_RH_GFP,'-','color','k','LineWidth',0.1);
ylabel('\Delta RH GFP','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'RH Rhodamine','RH GFP','Location','best','FontSize',5)

title('NREM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax5.YAxis(1).Color = colors('dark candy apple red');
ax5.YAxis(2).Color = 'k';
ax5.TickLength = [0.03,0.03];

% cort RH neural
ax7 = subplot(8,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];


% NREM to REM
ax9 = subplot(8,2,9);
% Rhodamine and GFP
p1 = plot(T1,data.NREMtoREM.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.mean_RH_Rhodamine + data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_RH_Rhodamine - data.NREMtoREM.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta RH Rhodamine')
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.NREMtoREM.mean_RH_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_RH_GFP + data.NREMtoREM.std_RH_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_RH_GFP - data.NREMtoREM.std_RH_GFP,'-','color','k','LineWidth',0.1);
ylabel('\Delta RH GFP','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'RH Rhodamine','RH GFP','Location','best','FontSize',5)

title('NREM to REM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax9.YAxis(1).Color = colors('dark candy apple red');
ax9.YAxis(2).Color = 'k';
ax9.TickLength = [0.03,0.03];

% cort RH neural
ax11 = subplot(8,2,13);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];

% REM to AWAKE
ax13 = subplot(8,2,10);
% Rhodamine and GFP
p1 = plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine + data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_RH_Rhodamine - data.REMtoAWAKE.std_RH_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\Delta RH Rhodamine')
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.REMtoAWAKE.mean_RH_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_RH_GFP + data.REMtoAWAKE.std_RH_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_RH_GFP - data.REMtoAWAKE.std_RH_GFP,'-','color','k','LineWidth',0.1);
ylabel('\Delta RH GFP','rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'RH Rhodamine','RH GFP','Location','best','FontSize',5)

title('REM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax13.YAxis(1).Color = colors('dark candy apple red');
ax13.YAxis(2).Color = 'k';
ax13.TickLength = [0.03,0.03];

% cort RH neural
ax15 = subplot(8,2,14);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_EEG,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
% ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];

% axes positionns
ax1Pos = get(ax1,'position');
ax3Pos = get(ax3,'position');
ax5Pos = get(ax5,'position');
ax7Pos = get(ax7,'position');
ax9Pos = get(ax9,'position');
ax11Pos = get(ax11,'position');
ax13Pos = get(ax13,'position');
ax15Pos = get(ax15,'position');

ax3Pos(3:4) = ax1Pos(3:4);

ax7Pos(3:4) = ax5Pos(3:4);

ax11Pos(3:4) = ax9Pos(3:4);

ax15Pos(3:4) = ax13Pos(3:4);

set(ax3,'position',ax3Pos);
set(ax7,'position',ax7Pos);
set(ax11,'position',ax11Pos);
set(ax15,'position',ax15Pos);
% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_D,[dirpath 'Fig4_Transition_Rhodamine_vs_GFP']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_Transition_Rhodamine_vs_GFP'])
end
close

end
