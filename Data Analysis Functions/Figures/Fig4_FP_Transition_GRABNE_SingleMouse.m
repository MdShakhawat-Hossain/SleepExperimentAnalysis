function [AnalysisResults] = Fig4_FP_Transition_GRABNE_SingleMouse(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,IOSanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
% IOSanimalIDs = {'NEACh001'};

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
ax1 = subplot(4,2,1);
% Blood Volume and GRABACh
p1 = plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

% plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine + data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_Ach_Rhodamine - data.AWAKEtoNREM.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
ylabel({'\DeltaF/F LH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.mean_Ach_GFP + data.AWAKEtoNREM.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_Ach_GFP - data.AWAKEtoNREM.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F ACh (Z)'},'rotation',-90,'VerticalAlignment','bottom')

title('Awake to NREM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax1.YAxis(2).Color = [0 0.4470 0.7410];
ax1.TickLength = [0.03,0.03];
axis square
% Blood Volume and GRAB NE
ax2 = subplot(4,2,3);

p1 = plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine + data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_NE_Rhodamine - data.AWAKEtoNREM.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F RH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.AWAKEtoNREM.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.mean_NE_GFP + data.AWAKEtoNREM.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.AWAKEtoNREM.mean_NE_GFP - data.AWAKEtoNREM.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
ylim([-3,3])
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

% title('Awake to NREM transition')
set(gca,'box','off')
ax2.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax2.YAxis(2).Color = [0.4660 0.6740 0.1880];
ax2.TickLength = [0.03,0.03];
axis square
% NREM to Awake
ax5 = subplot(4,2,2);
% Blood Volume and GRABACh
p1 = plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

% plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine + data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_Ach_Rhodamine - data.NREMtoAWAKE.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
ylabel({'\DeltaF/F LH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.mean_Ach_GFP + data.NREMtoAWAKE.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_Ach_GFP - data.NREMtoAWAKE.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F ACh (Z)'},'rotation',-90,'VerticalAlignment','bottom')

title('NREM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax5.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax5.YAxis(2).Color = [0 0.4470 0.7410];
ax5.TickLength = [0.03,0.03];
axis square
% Blood Volume and GRAB NE
ax6 = subplot(4,2,4);
% Blood Volume and GRABACh
p1 = plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine + data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_NE_Rhodamine - data.NREMtoAWAKE.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F RH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoAWAKE.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.NREMtoAWAKE.mean_NE_GFP + data.NREMtoAWAKE.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.NREMtoAWAKE.mean_NE_GFP - data.NREMtoAWAKE.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax6.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax6.YAxis(2).Color = [0.4660 0.6740 0.1880];
ax6.TickLength = [0.03,0.03];
axis square
% NREM to REM
if firstHrs == "false"
ax9 = subplot(4,2,5);
% Blood Volume and GRABACh
p1 = plot(T1,data.NREMtoREM.mean_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

% plot(T1,data.NREMtoREM.mean_Ach_Rhodamine + data.NREMtoREM.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_Ach_Rhodamine - data.NREMtoREM.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
ylabel({'\DeltaF/F LH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoREM.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% % plot(T1,data.NREMtoREM.mean_Ach_GFP + data.NREMtoREM.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_Ach_GFP - data.NREMtoREM.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F ACh (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('NREM to REM transition')
xlabel('Time (s)')
set(gca,'box','off')
ax9.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax9.YAxis(2).Color = [0 0.4470 0.7410];
ax9.TickLength = [0.03,0.03];
axis square

% Blood Volume and GRAB NE
ax10 = subplot(4,2,7);

% Blood Volume and GRABACh
p1 = plot(T1,data.NREMtoREM.mean_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,data.NREMtoREM.mean_NE_Rhodamine + data.NREMtoREM.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_NE_Rhodamine - data.NREMtoREM.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F RH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.NREMtoREM.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.NREMtoREM.mean_NE_GFP + data.NREMtoREM.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.NREMtoREM.mean_NE_GFP - data.NREMtoREM.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax10.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax10.YAxis(2).Color = [0.4660 0.6740 0.1880];
ax10.TickLength = [0.03,0.03];
axis square
% REM to AWAKE

ax13 = subplot(4,2,6);
% Blood Volume and GRABACh
p1 = plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

% plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine + data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_Ach_Rhodamine - data.REMtoAWAKE.std_Ach_Rhodamine,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F LH CBV (Z)'})
xlim([-30,30])
yyaxis right

p2 =  plot(T1,data.REMtoAWAKE.mean_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.mean_Ach_GFP + data.REMtoAWAKE.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_Ach_GFP - data.REMtoAWAKE.std_Ach_GFP,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F ACh (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'Ach Rhodamine','GRAB ACh','Location','northeast','FontSize',5)

title('REM to Awake transition')
xlabel('Time (s)')
set(gca,'box','off')
ax13.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax13.YAxis(2).Color = [0 0.4470 0.7410];
ax13.TickLength = [0.03,0.03];
axis square
% Blood Volume and GRAB NE
ax14 = subplot(4,2,8);
% Blood Volume and GRABACh
p1 = plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

% plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine + data.REMtoAWAKE.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_NE_Rhodamine - data.REMtoAWAKE.std_NE_Rhodamine,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F RH CBV (Z)'})
xlim([-30,30])
ylim([-3,3])
yyaxis right

p2 =  plot(T1,data.REMtoAWAKE.mean_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
% plot(T1,data.REMtoAWAKE.mean_NE_GFP + data.REMtoAWAKE.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% plot(T1,data.REMtoAWAKE.mean_NE_GFP - data.REMtoAWAKE.std_NE_GFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
ylim([-3,3])
ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
%legend([p1,p2],'NE Rhodamine','GRABNE','Location','northeast','FontSize',5)

set(gca,'box','off')
ax14.YAxis(1).Color = [0.6350 0.0780 0.1840];
ax14.YAxis(2).Color = [0.4660 0.6740 0.1880];
ax14.TickLength = [0.03,0.03];
axis square
end
% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
% ax7Pos = get(ax7,'position');
if firstHrs == "false"
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
ax13Pos = get(ax13,'position');
ax14Pos = get(ax14,'position');
% ax15Pos = get(ax15,'position');
end

% ax3Pos(3:4) = ax1Pos(3:4);

% ax7Pos(3:4) = ax5Pos(3:4);
if firstHrs == "false"
% ax11Pos(3:4) = ax9Pos(3:4);

% ax15Pos(3:4) = ax13Pos(3:4);
end

set(ax2,'position',ax2Pos);
% set(ax3,'position',ax3Pos);
set(ax6,'position',ax6Pos);
% set(ax7,'position',ax7Pos);
if firstHrs == "false"
set(ax10,'position',ax10Pos);
% set(ax11,'position',ax11Pos);
set(ax14,'position',ax14Pos);
% set(ax15,'position',ax15Pos);
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
    savefig(summaryFigure_D,[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse'])
end
close

end
