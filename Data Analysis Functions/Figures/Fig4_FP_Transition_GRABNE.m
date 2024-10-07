function [AnalysisResults] = Fig4_FP_Transition_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
IOSanimalIDs = {'NEACh001'};

    if firstHrs == "false"
         transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    elseif firstHrs == "true"
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};%,'NREMtoREM','REMtoAWAKE'};
    end
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(transition).EMG(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
%         data.(transition).Hip(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).Hip;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).Ach_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).Ach_Rhodamine;
        data.(transition).Ach_GFP(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_ACh;
        data.(transition).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;
        data.(transition).NE_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).NE_Rhodamine;
        data.(transition).NE_GFP(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NE;
    end
end
% take average for each behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).stdEMG = std(data.(transition).EMG,0,1);
%     data.(transition).meanHip = mean(data.(transition).Hip,3)*100;

    data.(transition).mean_Ach_Rhodamine = mean(data.(transition).Ach_Rhodamine,1);
    data.(transition).std_Ach_Rhodamine = std(data.(transition).Ach_Rhodamine,0,1);
    data.(transition).mean_Ach_GFP = mean(data.(transition).Ach_GFP,1);
    data.(transition).std_Ach_GFP = std(data.(transition).Ach_GFP,0,1);

    data.(transition).mean_NE_Rhodamine = mean(data.(transition).NE_Rhodamine,1);
    data.(transition).std_NE_Rhodamine = std(data.(transition).NE_Rhodamine,0,1);
    data.(transition).mean_NE_GFP = mean(data.(transition).NE_GFP,1);
    data.(transition).std_NE_GFP = std(data.(transition).NE_GFP,0,1);
    data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
T2 = -30 + (1/10):(1/10):30;

%% Fig. 4 Plot bilateral GCaMP7 and Rhodamine Transition
% summaryFigure_B_bilateral = figure('Name','Fig4_Transition GCaMP7-Rhodamine Bilateral');
% sgtitle('Plot bilateral GCaMP7 and Rhodamine Transition')
% % NREM to REM transition
% ax1 = subplot(2,2,2);
% p1 = plot(T1,data.NREMtoREM.mean_Ach_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
% hold on
% p2 = plot(T1,data.NREMtoREM.mean_NE_Rhodamine,'color',colors('army green'),'LineWidth',1);
% ylabel('\DeltaF/F')
% % ylim([40,80])
% yyaxis right
% p3 = plot(T1,data.NREMtoREM.mean_Ach_GFP,'color',colors('indian red'),'LineWidth',1);
% hold on
% p4 = plot(T1,data.NREMtoREM.mean_NE_GFP,'color','g','LineWidth',1);
% title('NREM to REM transition')
% xlabel('Time (s)')
% ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% legend([p1,p2,p3,p4],'Ach Rhodamine','NE Rhodamine','GRAB ACh','GRAB NE','Location','northeast','FontSize',5)
% xlim([-30,30])
% % ylim([5,42])
% axis square
% set(gca,'box','off')
% ax1.TickLength = [0.03,0.03];
% ax1.YAxis(1).Color = colors('dark candy apple red');
% ax1.YAxis(2).Color = colors('deep jungle green');
% % REM to Awake transition
% ax2 = subplot(2,2,4);
% p1 = plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
% hold on
% p2 = plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine,'color',colors('army green'),'LineWidth',1);
% ylabel('\DeltaF/F')
% % ylim([40,80])
% yyaxis right
% p3 = plot(T1,data.REMtoAWAKE.mean_Ach_GFP,'color',colors('indian red'),'LineWidth',1);
% 
% hold on
% p4 = plot(T1,data.REMtoAWAKE.mean_NE_GFP,'color','g','LineWidth',1);
% title('REM to Awake transition')
% xlabel('Time (s)')
% ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% % legend([p1,p2,p3,p4],'Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','Location','northeast','FontSize',5)
% xlim([-30,30])
% % ylim([5,42])
% axis square
% set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ax2.YAxis(1).Color = colors('dark candy apple red');
% ax2.YAxis(2).Color = colors('deep jungle green');
% %  NREM to Awake transition
% ax3 = subplot(2,2,3);
% p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
% hold on
% p2 = plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine,'color',colors('army green'),'LineWidth',1);
% ylabel('\DeltaF/F')
% % ylim([40,80])
% yyaxis right
% p3 = plot(T1,data.NREMtoAWAKE.mean_Ach_GFP,'color',colors('indian red'),'LineWidth',1);
% hold on
% p4 = plot(T1,data.NREMtoAWAKE.mean_NE_GFP,'color','g','LineWidth',1);
% title('NREM to Awake transition')
% xlabel('Time (s)')
% ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% % legend([p1,p2,p3,p4],'Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','Location','NorthWest')
% xlim([-30,30])
% % ylim([5,42])
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ax3.YAxis(1).Color = colors('dark candy apple red');
% ax3.YAxis(2).Color = colors('deep jungle green');
% %   Awake to NREM transition
% ax4 = subplot(2,2,1);
% p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
% hold on
% p2 = plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine,'color',colors('army green'),'LineWidth',1);
% ylabel('\DeltaF/F')
% % ylim([40,80])
% yyaxis right
% p3 = plot(T1,data.AWAKEtoNREM.mean_Ach_GFP,'color',colors('indian red'),'LineWidth',1);
% hold on
% p4 = plot(T1,data.AWAKEtoNREM.mean_NE_GFP,'color','g','LineWidth',1);
% title('Awake to NREM transition')
% xlabel('Time (s)')
% ylabel('\DeltaF/F','rotation',-90,'VerticalAlignment','bottom')
% % legend([p1,p2,p3,p4],'Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP','Location','northeast','FontSize',5)
% xlim([-30,30])
% % ylim([5,42])
% axis square
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];
% ax4.YAxis(1).Color = colors('dark candy apple red');
% ax4.YAxis(2).Color = colors('deep jungle green');
% % save figure(s)
% if strcmp(saveFigs,'y') == true
%     dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(summaryFigure_B_bilateral,[dirpath 'Fig4_Transition_GFP_Rhodamine']);
%     set(summaryFigure_B_bilateral,'PaperPositionMode','auto');
%     print('-painters','-dpdf','-bestfit',[dirpath 'Fig4_Transition_GFP_Rhodamine'])
% end
% close
%% Fig. 4 Plot Rhodamine Transition with Cortical LFPs Bilateral
summaryFigure_A = figure('Name','Fig4_Transition Rhodamine');
sgtitle('Rhodamine Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(6,2,1);
% Rhodamine and EMG
p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine + data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine - data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);

p2 = plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine + data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine - data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
ylabel('\DeltaF/F Rhodamine (Z)')
xlim([-30,30])
ylim([-1.75,1.75])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'Ach Rhodamine','NE Rhodamine','EMG','Location','northeast','FontSize',5)
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

axis square
% cort LH neural
ax3 = subplot(6,2,3);
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
ax3.TickLength = [0.03,0.03];
axis square
% hippocampal neural
% ax4 = subplot(6,2,5);
% Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(6,2,2);

% Rhodamine and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine + data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine - data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);

p2 = plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine + data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine - data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
ylabel('\DeltaF/F Rhodamine (Z)')
xlim([-30,30])
ylim([-1.75,1.75])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'Ach Rhodamine','NE Rhodamine','EMG','Location','northeast','FontSize',5)
ax5.YAxis(1).Color = colors('dark candy apple red');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];
axis square
% LH cort neural
ax7 = subplot(6,2,4);
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
ax7.TickLength = [0.03,0.03];
axis square
% hippocampal neural
% ax8 = subplot(6,2,6);
% Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];


% NREM to REM
if firstHrs == "false"
    ax9 = subplot(6,2,7);
    % Rhodamine and EMG
    p1 = plot(T1,data.NREMtoREM.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.NREMtoREM.mean_Ach_Rhodamine + data.NREMtoREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.NREMtoREM.mean_Ach_Rhodamine - data.NREMtoREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    
    p2 = plot(T1,data.NREMtoREM.mean_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',2);
    hold on
    plot(T1,data.NREMtoREM.mean_NE_Rhodamine + data.NREMtoREM.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
    plot(T1,data.NREMtoREM.mean_NE_Rhodamine - data.NREMtoREM.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
    ylabel('\DeltaF/F Rhodamine (Z)')
    xlim([-30,30])
    ylim([-3,3])
    % ylim([-5,50])
    yyaxis right
    p3 = plot(T1,data.NREMtoREM.meanEMG,'-','color','k','LineWidth',2);
    hold on
    % plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
    % plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
    title('NREM to REM transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    % legend([p1,p2,p3],'Ach Rhodamine','NE Rhodamine','EMG','Location','northeast','FontSize',5)
    ax9.YAxis(1).Color = colors('dark candy apple red');
    ax9.YAxis(2).Color = 'k';
    % ylim([-1,0.5])
    ax9.TickLength = [0.03,0.03];
    axis square
    % cort LH neural
    ax11 = subplot(6,2,9);
    Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_Cort,'y')
    axis xy
    c1 = colorbar;
    ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    xlabel('Time (s)')
    ylabel({'LH Cortical LFP';'Frequency (Hz)'})
    set(gca,'Yticklabel','10^1')
    xlim([-30,30])
    set(gca,'box','off')
    ax11.TickLength = [0.03,0.03];
    axis square
    % hippocampal neural
%     ax12 = subplot(6,2,11);
%     Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
%     c2 = colorbar;
%     ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%     caxis([-100,100])
%     xlabel('Time (s)')
%     ylabel({'Hippocampal LFP';'Frequency (Hz)'})
%     set(gca,'Yticklabel','10^1')
%     xlim([-30,30])
%     set(gca,'box','off')
%     ax12.TickLength = [0.03,0.03];
    
    % REM to Awake
    ax13 = subplot(6,2,8);
    
    % Rhodamine and EMG
    p1 = plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine + data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine - data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.5);
    
    p2 = plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine + data.REMtoAWAKE.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
    plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine - data.REMtoAWAKE.std_NE_Rhodamine,'-','color',colors('army green'),'LineWidth',0.5);
    ylabel('\DeltaF/F Rhodamine (Z)')
    xlim([-30,30])
    ylim([-3,3])
    % ylim([-5,50])
    yyaxis right
    p3 = plot(T1,data.REMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
    hold on
    % plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
    % plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
    title('REM to AWAKE transition')
    xlabel('Time (s)')
    ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
    set(gca,'box','off')
    % legend([p1,p2,p3],'Ach Rhodamine','NE Rhodamine','EMG','Location','northeast','FontSize',5)
    ax13.YAxis(1).Color = colors('dark candy apple red');
    ax13.YAxis(2).Color = 'k';
    % ylim([-1,0.5])
    ax13.TickLength = [0.03,0.03];
    axis square
    % LH cort neural
    ax15 = subplot(6,2,10);
    Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_Cort,'y')
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
    axis square
    % hippocampal neural
%     ax16 = subplot(6,2,12);
%     Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
%     c4 = colorbar;
%     ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%     caxis([-100,100])
%     xlabel('Time (s)')
%     % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
%     set(gca,'Yticklabel','10^1')
%     xlim([-30,30])
%     set(gca,'box','off')
%     ax16.TickLength = [0.03,0.03];
end

%%%% axes positionns
ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
if firstHrs == "false"
ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax13Pos = get(ax13,'position');
% ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
% ax16Pos = get(ax16,'position');
end

% ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax1Pos(3:4);

% ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
% ax8Pos(3:4) = ax5Pos(3:4);
if firstHrs == "false"
    % ax10Pos(3:4) = ax9Pos(3:4);
    ax11Pos(3:4) = ax9Pos(3:4);
%     ax12Pos(3:4) = ax9Pos(3:4);
    
    % ax14Pos(3:4) = ax13Pos(3:4);
    ax15Pos(3:4) = ax13Pos(3:4);
%     ax16Pos(3:4) = ax13Pos(3:4);
end
% set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
if firstHrs == "false"
    % set(ax10,'position',ax10Pos);
    set(ax11,'position',ax11Pos);
%     set(ax12,'position',ax12Pos);
    % set(ax14,'position',ax14Pos);
    set(ax15,'position',ax15Pos);
%     set(ax16,'position',ax16Pos);
end
%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
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
p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_Ach_GFP + data.AWAKEtoNREM.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.mean_Ach_GFP - data.AWAKEtoNREM.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);

p2 = plot(T1,data.AWAKEtoNREM.mean_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_NE_GFP + data.AWAKEtoNREM.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.mean_NE_GFP - data.AWAKEtoNREM.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
ylabel('\DeltaGFP')
xlim([-30,30])
ylim([-1.75,1.75])
yyaxis right
p3 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2,p3],'GRAB ACh','GRAB NE','EMG','Location','northeast','FontSize',5)
ax1.YAxis(1).Color = colors('indian red');
ax1.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax1.TickLength = [0.03,0.03];

% cort LH neural
ax3 = subplot(6,2,3);
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
ax3.TickLength = [0.03,0.03];

% hippocampal neural
% ax4 = subplot(6,2,5);
% Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(6,2,2);

% GFP and EMG
p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_Ach_GFP + data.NREMtoAWAKE.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.mean_Ach_GFP - data.NREMtoAWAKE.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);

p2 = plot(T1,data.NREMtoAWAKE.mean_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_NE_GFP + data.NREMtoAWAKE.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.mean_NE_GFP - data.NREMtoAWAKE.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
ylabel('\DeltaGFP')
xlim([-30,30])
ylim([-1.75,1.75])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('NREM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'GRAB ACh','GRAB NE','EMG','Location','northeast','FontSize',5)
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax5.TickLength = [0.03,0.03];

% LH cort neural
ax7 = subplot(6,2,4);
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
ax7.TickLength = [0.03,0.03];

% hippocampal neural
% ax8 = subplot(6,2,6);
% Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];

% NREM to REM
if firstHrs == "false"
ax9 = subplot(6,2,7);
% GFP and EMG
p1 = plot(T1,data.NREMtoREM.mean_Ach_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_Ach_GFP + data.NREMtoREM.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.mean_Ach_GFP - data.NREMtoREM.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);

p2 = plot(T1,data.NREMtoREM.mean_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_NE_GFP + data.NREMtoREM.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.mean_NE_GFP - data.NREMtoREM.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
ylabel('\DeltaGFP')
xlim([-30,30])
ylim([-3,3])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.NREMtoREM.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'GRAB ACh','GRAB NE','EMG','Location','northeast','FontSize',5)
ax9.YAxis(1).Color = colors('indian red');
ax9.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax9.TickLength = [0.03,0.03];

% cort LH neural
ax11 = subplot(6,2,9);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_Cort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (s)')
ylabel({'LH Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];

% hippocampal neural
% ax12 = subplot(6,2,11);
% Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];

% REM to Awake
ax13 = subplot(6,2,8);

% GFP and EMG
p1 = plot(T1,data.REMtoAWAKE.mean_Ach_GFP,'-','color',colors('indian red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_Ach_GFP + data.REMtoAWAKE.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.mean_Ach_GFP - data.REMtoAWAKE.std_Ach_GFP,'-','color',colors('indian red'),'LineWidth',0.5);

p2 = plot(T1,data.REMtoAWAKE.mean_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_NE_GFP + data.REMtoAWAKE.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.mean_NE_GFP - data.REMtoAWAKE.std_NE_GFP,'-','color',colors('deep jungle green'),'LineWidth',0.5);
ylabel('\DeltaGFP')
xlim([-30,30])
ylim([-3,3])
% ylim([-5,50])
yyaxis right
p3 = plot(T1,data.REMtoAWAKE.meanEMG,'-','color','k','LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
% plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors('deep jungle green'),'LineWidth',0.5);
title('REM to AWAKE transition')
xlabel('Time (s)')
ylabel('EMG power (a.u.)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
% legend([p1,p2,p3],'GRAB ACh','GRAB NE','EMG','Location','northeast','FontSize',5)
ax13.YAxis(1).Color = colors('indian red');
ax13.YAxis(2).Color = 'k';
% ylim([-1,0.5])
ax13.TickLength = [0.03,0.03];

% LH cort neural
ax15 = subplot(6,2,10);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_Cort,'y')
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

% hippocampal neural
% ax16 = subplot(6,2,12);
% Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
% c4 = colorbar;
% ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax16.TickLength = [0.03,0.03];
end
%% axes positionns
ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
if firstHrs == "false"
ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax13Pos = get(ax13,'position');
% ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
% ax16Pos = get(ax16,'position');
end

% ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax1Pos(3:4);

% ax6Pos(3:4) = ax5Pos(3:4);
ax7Pos(3:4) = ax5Pos(3:4);
% ax8Pos(3:4) = ax5Pos(3:4);
if firstHrs == "false"

% ax10Pos(3:4) = ax9Pos(3:4);
ax11Pos(3:4) = ax9Pos(3:4);
% ax12Pos(3:4) = ax9Pos(3:4);

% ax14Pos(3:4) = ax13Pos(3:4);
ax15Pos(3:4) = ax13Pos(3:4);
% ax16Pos(3:4) = ax13Pos(3:4);
end

% set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
if firstHrs == "false"
% set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
% set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
% set(ax16,'position',ax16Pos);
end
%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
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
sgtitle('Rhodamine and GFP Transition with Cortical LFPs')
% Awake to NREM
ax1 = subplot(8,2,1);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine + data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine - data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';' Ach (Z)'})
xlim([-30,30])
ylim([-1,1])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.mean_Ach_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_Ach_GFP + data.AWAKEtoNREM.std_Ach_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_Ach_GFP - data.AWAKEtoNREM.std_Ach_GFP,'-','color','k','LineWidth',0.1);
ylim([-1,1])
ylabel({'\DeltaF/F';'Ach (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('Awake to NREM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax1.YAxis(1).Color = colors('dark candy apple red');
ax1.YAxis(2).Color = 'k';
ax1.TickLength = [0.03,0.03];
% Rhodamine and GRAB NE
ax2 = subplot(8,2,3);
% Rhodamine and 
p1 = plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine + data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine - data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'NE (Z)'})
xlim([-30,30])
ylim([-1,1.5])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.mean_NE_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.mean_NE_GFP + data.AWAKEtoNREM.std_NE_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.AWAKEtoNREM.mean_NE_GFP - data.AWAKEtoNREM.std_NE_GFP,'-','color','k','LineWidth',0.1);
ylabel({'\DeltaF/F GRABNE';'NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
ylim([-1,1.5])
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

title('Awake to NREM transition')
set(gca,'box','off')
ax2.YAxis(1).Color = colors('dark candy apple red');
ax2.YAxis(2).Color = 'k';
ax2.TickLength = [0.03,0.03];

% cort LH neural
ax3 = subplot(8,2,5);
Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.mean_LH_Cort,'y')
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
% ax4 = subplot(8,2,7);
% Semilog_ImageSC(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax4.TickLength = [0.03,0.03];

% NREM to Awake
ax5 = subplot(8,2,2);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine + data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine - data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'Ach (Z)'})
xlim([-30,30])
ylim([-1,1])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.mean_Ach_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_Ach_GFP + data.NREMtoAWAKE.std_Ach_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_Ach_GFP - data.NREMtoAWAKE.std_Ach_GFP,'-','color','k','LineWidth',0.1);
ylim([-1,1])
ylabel({'\DeltaF/F';'Ach (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('NREM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax5.YAxis(1).Color = colors('dark candy apple red');
ax5.YAxis(2).Color = 'k';
ax5.TickLength = [0.03,0.03];
% Rhodamine and GRAB NE
ax6 = subplot(8,2,4);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine + data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine - data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'NE (Z)'})
xlim([-30,30])
ylim([-1,1.5])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.mean_NE_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.mean_NE_GFP + data.NREMtoAWAKE.std_NE_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoAWAKE.mean_NE_GFP - data.NREMtoAWAKE.std_NE_GFP,'-','color','k','LineWidth',0.1);
ylim([-1,1.5])
ylabel({'\DeltaF/F GRABNE';'NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax6.YAxis(1).Color = colors('dark candy apple red');
ax6.YAxis(2).Color = 'k';
ax6.TickLength = [0.03,0.03];

% cort LH neural
ax7 = subplot(8,2,6);
Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.mean_LH_Cort,'y')
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
% hippocampal neural
% ax8 = subplot(8,2,8);
% Semilog_ImageSC(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];

% NREM to REM
if firstHrs == "false"
ax9 = subplot(8,2,9);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.NREMtoREM.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.mean_Ach_Rhodamine + data.NREMtoREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_Ach_Rhodamine - data.NREMtoREM.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'Ach (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoREM.mean_Ach_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_Ach_GFP + data.NREMtoREM.std_Ach_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_Ach_GFP - data.NREMtoREM.std_Ach_GFP,'-','color','k','LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F';'Ach (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('NREM to REM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax9.YAxis(1).Color = colors('dark candy apple red');
ax9.YAxis(2).Color = 'k';
ax9.TickLength = [0.03,0.03];
% Rhodamine and GRAB NE
ax10 = subplot(8,2,11);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.NREMtoREM.mean_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.mean_NE_Rhodamine + data.NREMtoREM.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_NE_Rhodamine - data.NREMtoREM.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'NE (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoREM.mean_NE_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.NREMtoREM.mean_NE_GFP + data.NREMtoREM.std_NE_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.NREMtoREM.mean_NE_GFP - data.NREMtoREM.std_NE_GFP,'-','color','k','LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F GRABNE';'NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax10.YAxis(1).Color = colors('dark candy apple red');
ax10.YAxis(2).Color = 'k';
ax10.TickLength = [0.03,0.03];

% cort LH neural
ax11 = subplot(8,2,13);
Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.mean_LH_Cort,'y')
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
% hippocampal neural
% ax12 = subplot(8,2,15);
% Semilog_ImageSC(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];

% REM to AWAKE

ax13 = subplot(8,2,10);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine + data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine - data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F Rhodamine';'Ach (Z)'})
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.REMtoAWAKE.mean_Ach_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_Ach_GFP + data.REMtoAWAKE.std_Ach_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_Ach_GFP - data.REMtoAWAKE.std_Ach_GFP,'-','color','k','LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F';'Ach (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('REM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax13.YAxis(1).Color = colors('dark candy apple red');
ax13.YAxis(2).Color = 'k';
ax13.TickLength = [0.03,0.03];
% Rhodamine and GRAB NE
ax14 = subplot(8,2,12);
% Rhodamine and GCaMP7s
p1 = plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',2);
hold on

plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine + data.REMtoAWAKE.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine - data.REMtoAWAKE.std_NE_Rhodamine,'-','color',colors('dark candy apple red'),'LineWidth',0.1);
ylabel({'\DeltaF/F Rhodamine';'NE (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.REMtoAWAKE.mean_NE_GFP,'-','color','k','LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.mean_NE_GFP + data.REMtoAWAKE.std_NE_GFP,'-','color','k','LineWidth',0.1);
plot(T1,data.REMtoAWAKE.mean_NE_GFP - data.REMtoAWAKE.std_NE_GFP,'-','color','k','LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F GRABNE';'NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax14.YAxis(1).Color = colors('dark candy apple red');
ax14.YAxis(2).Color = 'k';
ax14.TickLength = [0.03,0.03];

% cort LH neural
ax15 = subplot(8,2,14);
Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.mean_LH_Cort,'y')
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
% % hippocampal neural
% ax16 = subplot(8,2,16);
% Semilog_ImageSC(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% xlabel('Time (s)')
% % ylabel({'Hippocampal LFP';'Frequency (Hz)'})
% set(gca,'Yticklabel','10^1')
% xlim([-30,30])
% set(gca,'box','off')
% ax16 .TickLength = [0.03,0.03];
end
% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');
if firstHrs == "false"
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax13Pos = get(ax13,'position');
ax14Pos = get(ax14,'position');
ax15Pos = get(ax15,'position');
% ax16Pos = get(ax16,'position');
end

ax3Pos(3:4) = ax1Pos(3:4);
% ax4Pos(3:4) = ax1Pos(3:4);

ax7Pos(3:4) = ax5Pos(3:4);
% ax8Pos(3:4) = ax5Pos(3:4);
if firstHrs == "false"
ax11Pos(3:4) = ax9Pos(3:4);
% ax12Pos(3:4) = ax9Pos(3:4);

ax15Pos(3:4) = ax13Pos(3:4);
% ax16Pos(3:4) = ax13Pos(3:4);
end

set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
if firstHrs == "false"
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
set(ax14,'position',ax14Pos);
set(ax15,'position',ax15Pos);
% set(ax16,'position',ax16Pos);
end
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
    savefig(summaryFigure_D,[dirpath 'Fig4_Transition_Rhodamine_vs_GFP']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4_Transition_Rhodamine_vs_GFP'])
end
close

end
