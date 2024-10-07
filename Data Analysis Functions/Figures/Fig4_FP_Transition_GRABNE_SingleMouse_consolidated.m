function [AnalysisResults] = Fig4_FP_Transition_GRABNE_SingleMouse_consolidated(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,IOSanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
% IOSanimalIDs = {'NEACh008'};

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
        % data.(transition).EMG_std(aa,:) = AnalysisResults.(animalID).Transitions.(transition).EMG_std;

        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).LH_Cort(:,:,aa) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;

        data.(transition).Ach_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).Ach_Rhodamine;
        data.(transition).Ach_GFP(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_ACh;
        
        data.(transition).NE_Rhodamine(aa,:) = AnalysisResults.(animalID).Transitions.(transition).NE_Rhodamine;
        data.(transition).NE_GFP(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NE;

        % data.(transition).Ach_Rhodamine_std(aa,:) = AnalysisResults.(animalID).Transitions.(transition).Ach_Rhodamine_std;
        % data.(transition).Ach_GFP_std(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_ACh_std;
        
        % data.(transition).NE_Rhodamine_std(aa,:) = AnalysisResults.(animalID).Transitions.(transition).NE_Rhodamine_std;
        % data.(transition).NE_GFP_std(aa,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NE_std;
    end
end
% take average for each behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    % data.(transition).stdEMG = mean(data.(transition).EMG_std,1);

    data.(transition).mean_Ach_Rhodamine = mean(data.(transition).Ach_Rhodamine,1);
    % data.(transition).std_Ach_Rhodamine = mean(data.(transition).Ach_Rhodamine_std,1);

    data.(transition).mean_Ach_GFP = mean(data.(transition).Ach_GFP,1);
    % data.(transition).std_Ach_GFP = mean(data.(transition).Ach_GFP_std,1);

    data.(transition).mean_NE_Rhodamine = mean(data.(transition).NE_Rhodamine,1);
    % data.(transition).std_NE_Rhodamine = mean(data.(transition).NE_Rhodamine_std,1);

    data.(transition).mean_NE_GFP = mean(data.(transition).NE_GFP,1);
    % data.(transition).std_NE_GFP = mean(data.(transition).NE_GFP_std,1);

    data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
% T2 = -30 + (1/10):(1/10):30;
%% Fig. 4 Plot Rhodamine and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs Rhodamine');
sgtitle('Blood Volume and Neuromodulators Transition')

% Awake to NREM

% Blood Volume and NE, ACh
ax1 = subplot(4,2,1);
p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,mean([data.AWAKEtoNREM.mean_NE_Rhodamine;data.AWAKEtoNREM.mean_Ach_Rhodamine]) + mean([data.AWAKEtoNREM.NE_Rhodamine_std;data.AWAKEtoNREM.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,mean([data.AWAKEtoNREM.mean_NE_Rhodamine;data.AWAKEtoNREM.mean_Ach_Rhodamine]) - mean([data.AWAKEtoNREM.NE_Rhodamine_std;data.AWAKEtoNREM.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine + data.AWAKEtoNREM.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine - data.AWAKEtoNREM.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);


ylabel({'\DeltaF/F (Z)'})

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.AWAKEtoNREM.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.mean_NE_GFP + data.AWAKEtoNREM.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_NE_GFP - data.AWAKEtoNREM.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
p3 =  plot(T1,data.AWAKEtoNREM.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.mean_Ach_GFP + data.AWAKEtoNREM.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_Ach_GFP - data.AWAKEtoNREM.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-2,2])
legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)
xlim([-30,30])
title('Awake to NREM')
set(gca,'box','off')
ax1.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
ax1.TickLength = [0.03,0.03];
axis square

% NREM to Awake
ax2 = subplot(4,2,2);
p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,mean([data.NREMtoAWAKE.mean_NE_Rhodamine;data.NREMtoAWAKE.mean_Ach_Rhodamine]) + mean([data.NREMtoAWAKE.NE_Rhodamine_std;data.NREMtoAWAKE.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,mean([data.NREMtoAWAKE.mean_NE_Rhodamine;data.NREMtoAWAKE.mean_Ach_Rhodamine]) - mean([data.NREMtoAWAKE.NE_Rhodamine_std;data.NREMtoAWAKE.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine + data.NREMtoAWAKE.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine - data.NREMtoAWAKE.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

ylabel({'\DeltaF/F (Z)'})

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.NREMtoAWAKE.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.mean_NE_GFP + data.NREMtoAWAKE.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_NE_GFP - data.NREMtoAWAKE.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
p3 =  plot(T1,data.NREMtoAWAKE.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.mean_Ach_GFP + data.NREMtoAWAKE.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_Ach_GFP - data.NREMtoAWAKE.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-2,2])
legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)

title('NREM to Awake')
set(gca,'box','off')
ax2.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
ax2.TickLength = [0.03,0.03];
xlim([-30,30])
axis square
%%
% NREM to REM

% Blood Volume and NE, ACh
ax3 = subplot(4,2,3);
p1 = plot(T1,data.NREMtoREM.mean_Ach_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,mean([data.NREMtoREM.mean_NE_Rhodamine;data.NREMtoREM.mean_Ach_Rhodamine]) + mean([data.NREMtoREM.NE_Rhodamine_std;data.NREMtoREM.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,mean([data.NREMtoREM.mean_NE_Rhodamine;data.NREMtoREM.mean_Ach_Rhodamine]) - mean([data.NREMtoREM.NE_Rhodamine_std;data.NREMtoREM.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_Ach_Rhodamine + data.NREMtoREM.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_Ach_Rhodamine - data.NREMtoREM.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.NREMtoREM.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.mean_NE_GFP + data.NREMtoREM.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_NE_GFP - data.NREMtoREM.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
p3 =  plot(T1,data.NREMtoREM.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.mean_Ach_GFP + data.NREMtoREM.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_Ach_GFP - data.NREMtoREM.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-2.5,2.5])
% legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)
xlim([-30,30])
title('NREM to REM')
set(gca,'box','off')
ax3.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
ax3.TickLength = [0.03,0.03];
axis square
% REM to Awake
ax4 = subplot(4,2,4);
p1 = plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,mean([data.REMtoAWAKE.mean_NE_Rhodamine;data.REMtoAWAKE.mean_Ach_Rhodamine]) + mean([data.REMtoAWAKE.NE_Rhodamine_std;data.REMtoAWAKE.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,mean([data.REMtoAWAKE.mean_NE_Rhodamine;data.REMtoAWAKE.mean_Ach_Rhodamine]) - mean([data.REMtoAWAKE.NE_Rhodamine_std;data.REMtoAWAKE.Ach_Rhodamine_std]),'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine + data.REMtoAWAKE.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine - data.REMtoAWAKE.Ach_Rhodamine_std,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.REMtoAWAKE.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.mean_NE_GFP + data.REMtoAWAKE.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_NE_GFP - data.REMtoAWAKE.NE_GFP_std,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
p3 =  plot(T1,data.REMtoAWAKE.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.mean_Ach_GFP + data.REMtoAWAKE.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_Ach_GFP - data.REMtoAWAKE.Ach_GFP_std,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-2.5,2.5])
% legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)

title('REM to Awake')
set(gca,'box','off')
ax4.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
ax4.TickLength = [0.03,0.03];
xlim([-30,30])
axis square
% % axes positionns
% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% % ax3Pos = get(ax3,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% % ax7Pos = get(ax7,'position');
% if firstHrs == "false"
% ax9Pos = get(ax9,'position');
% ax10Pos = get(ax10,'position');
% % ax11Pos = get(ax11,'position');
% ax13Pos = get(ax13,'position');
% ax14Pos = get(ax14,'position');
% % ax15Pos = get(ax15,'position');
% end

% ax3Pos(3:4) = ax1Pos(3:4);

% % ax7Pos(3:4) = ax5Pos(3:4);
% if firstHrs == "false"
% % ax11Pos(3:4) = ax9Pos(3:4);
% 
% % ax15Pos(3:4) = ax13Pos(3:4);
% end
% 
% set(ax2,'position',ax2Pos);
% % set(ax3,'position',ax3Pos);
% set(ax6,'position',ax6Pos);
% % set(ax7,'position',ax7Pos);
% if firstHrs == "false"
% set(ax10,'position',ax10Pos);
% % set(ax11,'position',ax11Pos);
% set(ax14,'position',ax14Pos);
% % set(ax15,'position',ax15Pos);
% end
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
    savefig(summaryFigure_D,[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse_Sleep_consolidated']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse_Sleep_consolidated'])
end
close

end
