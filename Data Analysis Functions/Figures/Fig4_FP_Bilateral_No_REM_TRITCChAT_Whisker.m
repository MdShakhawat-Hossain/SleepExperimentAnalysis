function [AnalysisResults] = Fig4_FP_Bilateral_No_REM_TRITCChAT_Whisker(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'T282','T285'};%{'T281','T282','T284','T285'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(transition).whisk(aa,:) = AnalysisResults.(animalID).Transitions.(transition).whisk;
        data.(transition).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
        data.(transition).Hip(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Hip;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Cort;
        data.(transition).TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).TRITC;
        data.(transition).GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GCaMP7s;

        data.(transition).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;
        data.(transition).LH_TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).LH_TRITC;
        data.(transition).LH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s;
        data.(transition).RH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).RH_Cort;
        data.(transition).RH_TRITC(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_TRITC;
        data.(transition).RH_GCaMP7s(aa,:) = AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s;
    end
end
% take average for each behavioral transition _RH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);
    data.(transition).meanwhisk = mean(data.(transition).whisk,1);
    data.(transition).stdwhisk = std(data.(transition).whisk,0,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;

    data.(transition).meanTRITC = mean(data.(transition).TRITC,1);
    data.(transition).stdTRITC = std(data.(transition).TRITC,0,1);
    data.(transition).meanGCaMP7s = mean(data.(transition).GCaMP7s,1);
    data.(transition).stdGCaMP7s = std(data.(transition).GCaMP7s,0,1);
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;

    data.(transition).mean_LH_TRITC = mean(data.(transition).LH_TRITC,1);
    data.(transition).std_LH_TRITC = std(data.(transition).LH_TRITC,0,1);
    data.(transition).mean_LH_GCaMP7s = mean(data.(transition).LH_GCaMP7s,1);
    data.(transition).std_LH_GCaMP7s = std(data.(transition).LH_GCaMP7s,0,1);
    data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;

    data.(transition).mean_RH_TRITC = mean(data.(transition).RH_TRITC,1);
    data.(transition).std_RH_TRITC = std(data.(transition).RH_TRITC,0,1);
    data.(transition).mean_RH_GCaMP7s = mean(data.(transition).RH_GCaMP7s,1);
    data.(transition).std_RH_GCaMP7s = std(data.(transition).RH_GCaMP7s,0,1);
    data.(transition).mean_RH_Cort = mean(data.(transition).RH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;
%% Fig. 4 Plot TRITC Transition with Cortical LFPs
summaryFigure_A = figure('Name','Fig4 A');
sgtitle('Transition with Cortical LFPs')
% Awake to NREM
% TRITC and GCaMP7
ax1 = subplot(4,2,1);
p1 = plot(T1,data.AWAKEtoNREM.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.meanTRITC + data.AWAKEtoNREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.meanTRITC - data.AWAKEtoNREM.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.meanGCaMP7s,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanGCaMP7s + data.AWAKEtoNREM.stdGCaMP7s,'-','color','k','LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.meanGCaMP7s - data.AWAKEtoNREM.stdGCaMP7s,'-','color','k','LineWidth',0.1);

title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('\DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')

set(gca,'box','off')
legend([p1,p2],'Rhodamine','GCaMP7s','Location','best','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% EMG and Whisk
ax2 = subplot(4,2,3);
p1 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('EMG a.u.')
xlim([-30,30])
% ylim([-5,50])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.meanwhisk,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanwhisk + data.AWAKEtoNREM.stdwhisk,'-','color','k','LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.meanwhisk - data.AWAKEtoNREM.stdwhisk,'-','color','k','LineWidth',0.1);

title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('Whisk angle','rotation',-90,'VerticalAlignment','bottom')

set(gca,'box','off')
legend([p1,p2],'EMG','whisk angle','Location','best','FontSize',5)
ax2.YAxis(1).Color = colors('dark candy apple red');
ax2.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax2.TickLength = [0.03,0.03];

% cort neural
ax3 = subplot(4,2,5);
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
% TRITC and GCaMP7s
plot(T1,data.NREMtoAWAKE.meanTRITC,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanTRITC + data.NREMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.meanTRITC - data.NREMtoAWAKE.stdTRITC,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('\DeltaRhodamine')
xlim([-30,30])
% ylim([-5,50])
yyaxis right

plot(T1,data.NREMtoAWAKE.meanGCaMP7s,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanGCaMP7s + data.NREMtoAWAKE.stdGCaMP7s,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.meanGCaMP7s - data.NREMtoAWAKE.stdGCaMP7s,'-','color','k','LineWidth',0.1);

title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('  \DeltaGCaMP7s','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax5.YAxis(1).Color = colors('dark candy apple red');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% EMG and Whisk
ax6 = subplot(4,2,4);
p1 = plot(T1,data.NREMtoAWAKE.meanEMG,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel('EMG a.u.')
xlim([-30,30])
% ylim([-5,50])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.meanwhisk,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanwhisk + data.NREMtoAWAKE.stdwhisk,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.meanwhisk - data.NREMtoAWAKE.stdwhisk,'-','color','k','LineWidth',0.1);

title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('Whisk angle','rotation',-90,'VerticalAlignment','bottom')

set(gca,'box','off')
legend([p1,p2],'EMG','whisk angle','Location','best','FontSize',5)
ax6.YAxis(1).Color = colors('dark candy apple red');
ax6.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax6.TickLength = [0.03,0.03];

% cort neural
ax7 = subplot(4,2,6);
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
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);
ax8Pos(3:4) = ax1Pos(3:4);

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure_A,[dirpath 'Fig4_A_TRITCChAT_whisk']);
    set(summaryFigure_A,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_A_TRITCChAT_whisk'])
end
close


end
