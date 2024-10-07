function [AnalysisResults] = Fig7_FP_Transition_Pupil(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data

    if firstHrs == "false"
         transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    elseif firstHrs == "true"
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};%,'NREMtoREM','REMtoAWAKE'};
    end
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};

        % the size of the matrix
        MatLength = size(AnalysisResults.(animalID).Transitions.(transition).Ach_RhodamineRaw,1);
        if aa == 1
            DataLength = 0;
        else
            DataLength = size(data.(transition).AchRhodamineRaw,1);
        end           

        data.(transition).PupilDiameterRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).PupilDiameterRaw;

        % data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        % data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        % data.(transition).LH_Cort(:,:,DataLength+1:DataLength+MatLength) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;

        data.(transition).AchRhodamineRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).Ach_RhodamineRaw;
        data.(transition).AchGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_AChRaw;
        
        data.(transition).NERhodamineRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).NE_RhodamineRaw;
        data.(transition).NEGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NERaw;
    end
end
% take average for each behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    data.(transition).PupilDiameter_N_Exp = size(data.(transition).PupilDiameterRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).PupilDiameter_Mean = mean(data.(transition).PupilDiameterRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).PupilDiameter_SEM = std(data.(transition).PupilDiameterRaw,1)/sqrt(data.(transition).PupilDiameter_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).PupilDiameter_CI95 = tinv([0.025 0.975], data.(transition).PupilDiameter_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).PupilDiameter_yCI95 = bsxfun(@times, data.(transition).PupilDiameter_SEM, data.(transition).PupilDiameter_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
    

    data.(transition).AchRhodamine_N_Exp = size(data.(transition).AchRhodamineRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AchRhodamine_Mean = mean(data.(transition).AchRhodamineRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AchRhodamine_SEM = std(data.(transition).AchRhodamineRaw,1)/sqrt(data.(transition).AchRhodamine_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AchRhodamine_CI95 = tinv([0.025 0.975], data.(transition).AchRhodamine_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AchRhodamine_yCI95 = bsxfun(@times, data.(transition).AchRhodamine_SEM, data.(transition).AchRhodamine_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
    
    data.(transition).NERhodamine_N_Exp = size(data.(transition).NERhodamineRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NERhodamine_Mean = mean(data.(transition).NERhodamineRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NERhodamine_SEM = std(data.(transition).NERhodamineRaw,1)/sqrt(data.(transition).NERhodamine_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NERhodamine_CI95 = tinv([0.025 0.975], data.(transition).NERhodamine_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NERhodamine_yCI95 = bsxfun(@times, data.(transition).NERhodamine_SEM, data.(transition).NERhodamine_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).AchGFP_N_Exp = size(data.(transition).AchGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AchGFP_Mean = mean(data.(transition).AchGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AchGFP_SEM = std(data.(transition).AchGFPRaw,1)/sqrt(data.(transition).AchGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AchGFP_CI95 = tinv([0.025 0.975], data.(transition).AchGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AchGFP_yCI95 = bsxfun(@times, data.(transition).AchGFP_SEM, data.(transition).AchGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).NEGFP_N_Exp = size(data.(transition).NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NEGFP_Mean = mean(data.(transition).NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_SEM = std(data.(transition).NEGFPRaw,1)/sqrt(data.(transition).NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_CI95 = tinv([0.025 0.975], data.(transition).NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NEGFP_yCI95 = bsxfun(@times, data.(transition).NEGFP_SEM, data.(transition).NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

end
T1 = -30 + (1/30):(1/30):30;
%% Fig. 4 Plot Rhodamine and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Pupil vs GFP vs Rhodamine');
sgtitle('Pupil Blood Volume and Neuromodulators Transition')

% Awake to NREM
scatter3(data.AWAKEtoNREM.PupilDiameter_Mean(1:end/2),data.AWAKEtoNREM.NEGFP_Mean(1:end/2),data.AWAKEtoNREM.AchGFP_Mean(1:end/2),'MarkerFaceColor','k','MarkerEdgeColor','k')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')
hold on
% scatter3(data.AWAKEtoNREM.PupilDiameter_Mean(1+end/2:end),data.AWAKEtoNREM.NEGFP_Mean(1+end/2:end),data.AWAKEtoNREM.AchGFP_Mean(1+end/2:end),'MarkerFaceColor','m','MarkerEdgeColor','m')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')

scatter3(data.NREMtoAWAKE.PupilDiameter_Mean(1:end/2),data.NREMtoAWAKE.NEGFP_Mean(1:end/2),data.NREMtoAWAKE.AchGFP_Mean(1:end/2),'MarkerFaceColor','m','MarkerEdgeColor','m')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')
scatter3(data.NREMtoAWAKE.PupilDiameter_Mean(1+end/2:end),data.NREMtoAWAKE.NEGFP_Mean(1+end/2:end),data.NREMtoAWAKE.AchGFP_Mean(1+end/2:end),'MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('Pupil'); ylabel('NE');zlabel('ACh')
scatter3(data.NREMtoREM.PupilDiameter_Mean(1:end/2),data.NREMtoREM.NEGFP_Mean(1:end/2),data.NREMtoREM.AchGFP_Mean(1:end/2),'MarkerFaceColor','m','MarkerEdgeColor','m')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')
hold on
scatter3(data.NREMtoREM.PupilDiameter_Mean(1+end/2:end),data.NREMtoREM.NEGFP_Mean(1+end/2:end),data.NREMtoREM.AchGFP_Mean(1+end/2:end),'MarkerFaceColor','g','MarkerEdgeColor','g')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')

scatter3(data.REMtoAWAKE.PupilDiameter_Mean(1:end/2),data.REMtoAWAKE.NEGFP_Mean(1:end/2),data.REMtoAWAKE.AchGFP_Mean(1:end/2),'MarkerFaceColor','g','MarkerEdgeColor','g')
% xlabel('Pupil'); ylabel('NE');zlabel('ACh')
% scatter3(data.REMtoAWAKE.PupilDiameter_Mean(1+end/2:end),data.REMtoAWAKE.NEGFP_Mean(1+end/2:end),data.REMtoAWAKE.AchGFP_Mean(1+end/2:end),'MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('Pupil'); ylabel('NE');zlabel('ACh')

% Blood Volume and NE, ACh
ax1 = subplot(2,2,1);
p1 = plot(T1,data.AWAKEtoNREM.AchRhodamine_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.AchRhodamine_Mean + data.AWAKEtoNREM.AchRhodamine_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F (Z)'})

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.AWAKEtoNREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.NEGFP_Mean + data.AWAKEtoNREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);

p3 =  plot(T1,data.AWAKEtoNREM.AchGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.AchGFP_Mean + data.AWAKEtoNREM.AchGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-1.5 1.5])
legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)
xlim([-30,30])
title('Awake to NREM')
set(gca,'box','off')
ax1.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
ax1.TickLength = [0.03,0.03];
axis square

% NREM to Awake
ax2 = subplot(2,2,2);
p1 = plot(T1,data.NREMtoAWAKE.AchRhodamine_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.AchRhodamine_Mean + data.NREMtoAWAKE.AchRhodamine_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

ylabel({'\DeltaF/F (Z)'})

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.NREMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.NEGFP_Mean + data.NREMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
p3 =  plot(T1,data.NREMtoAWAKE.AchGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.AchGFP_Mean + data.NREMtoAWAKE.AchGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-1.5 1.5])
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
ax3 = subplot(2,2,3);
p1 = plot(T1,data.NREMtoREM.AchRhodamine_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.AchRhodamine_Mean + data.NREMtoREM.AchRhodamine_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.NREMtoREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.NEGFP_Mean + data.NREMtoREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE (Z)'},'rotation',-90,'VerticalAlignment','bottom')
p3 =  plot(T1,data.NREMtoREM.AchGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.AchGFP_Mean + data.NREMtoREM.AchGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3 3])
% legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)
xlim([-30,30])
title('NREM to REM')
set(gca,'box','off')
ax3.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
ax3.TickLength = [0.03,0.03];
axis square
% REM to Awake
ax4 = subplot(2,2,4);
p1 = plot(T1,data.REMtoAWAKE.AchRhodamine_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

plot(T1,data.REMtoAWAKE.AchRhodamine_Mean + data.REMtoAWAKE.AchRhodamine_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.REMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.NEGFP_Mean + data.REMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
p3 =  plot(T1,data.REMtoAWAKE.AchGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.AchGFP_Mean + data.REMtoAWAKE.AchGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
ylim([-3 3])
% legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)

title('REM to Awake')
set(gca,'box','off')
ax4.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
ax4.TickLength = [0.03,0.03];
xlim([-30,30])
axis square
% figure;
% subplot(2,2,4)
% plot(T1,data.REMtoAWAKE.AchGFPRaw,'-','color',[0 0.4470 0.7410],'LineWidth',2);

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
    savefig(summaryFigure_D,[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse_Sleep_consolidated_CI']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Transition_Rhodamine_vs_GFP_singleMouse_Sleep_consolidated_CI'])
end
close

end
