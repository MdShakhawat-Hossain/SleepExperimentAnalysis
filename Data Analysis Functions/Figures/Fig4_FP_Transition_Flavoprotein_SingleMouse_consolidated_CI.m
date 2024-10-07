function [AnalysisResults] = Fig4_FP_Transition_Flavoprotein_SingleMouse_consolidated_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,IOSanimalIDs)
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
%% IOS mean transitions between eACh arousal-state
% cd through eACh animal's directory and extract the appropriate analysis results
for aa = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};

        % the size of the matrix
        MatLength = size(AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw,1);
        if aa == 1
            DataLength = 0;
        else
            DataLength = size(data.(transition).AChCBVRaw,1);
        end           

        % data.(transition).EMG(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;

        % data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        % data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        % data.(transition).LH_Cort(:,:,DataLength+1:DataLength+MatLength) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;

        data.(transition).AChCBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw;
        data.(transition).AChGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_AChRaw;
        
        data.(transition).NECBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).NE_CBVRaw;
        data.(transition).NEGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NERaw;
    end
end
% take average for eACh behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    % data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).AChCBV_N_Exp = size(data.(transition).AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChCBV_Mean = mean(data.(transition).AChCBVRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(transition).AChCBV_SEM = std(data.(transition).AChCBVRaw,1)/sqrt(data.(transition).AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(transition).AChCBV_CI95 = tinv([0.025 0.975], data.(transition).AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChCBV_yCI95 = bsxfun(@times, data.(transition).AChCBV_SEM, data.(transition).AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’
    
    data.(transition).NECBV_N_Exp = size(data.(transition).NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NECBV_Mean = mean(data.(transition).NECBVRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(transition).NECBV_SEM = std(data.(transition).NECBVRaw,1)/sqrt(data.(transition).NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(transition).NECBV_CI95 = tinv([0.025 0.975], data.(transition).NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NECBV_yCI95 = bsxfun(@times, data.(transition).NECBV_SEM, data.(transition).NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    data.(transition).AChGFP_N_Exp = size(data.(transition).AChGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChGFP_Mean = mean(data.(transition).AChGFPRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(transition).AChGFP_SEM = std(data.(transition).AChGFPRaw,1)/sqrt(data.(transition).AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(transition).AChGFP_CI95 = tinv([0.025 0.975], data.(transition).AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChGFP_yCI95 = bsxfun(@times, data.(transition).AChGFP_SEM, data.(transition).AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    data.(transition).NEGFP_N_Exp = size(data.(transition).NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NEGFP_Mean = mean(data.(transition).NEGFPRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(transition).NEGFP_SEM = std(data.(transition).NEGFPRaw,1)/sqrt(data.(transition).NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(transition).NEGFP_CI95 = tinv([0.025 0.975], data.(transition).NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NEGFP_yCI95 = bsxfun(@times, data.(transition).NEGFP_SEM, data.(transition).NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    % data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
%% Fig. 4 Plot CBV and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs CBV');
sgtitle('Blood Volume and Neuromodulators Transition')

% Awake to NREM

% Blood Volume and NE, ACh
ax1 = subplot(2,2,1);
% p1 = plot(T1,data.AWAKEtoNREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
% plot(T1,data.AWAKEtoNREM.AChCBV_Mean + data.AWAKEtoNREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
% ylabel({'\DeltaF/F (%)'})

p4 = plot(T1,data.AWAKEtoNREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.NECBV_Mean + data.AWAKEtoNREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F (%)'})

ylim([-1,1])

% yyaxis right
p2 =  plot(T1,data.AWAKEtoNREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.NEGFP_Mean + data.AWAKEtoNREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);

% p3 =  plot(T1,data.AWAKEtoNREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% plot(T1,data.AWAKEtoNREM.AChGFP_Mean + data.AWAKEtoNREM.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% ylim([-1.5 1.5])
% legend([p1,p4,p2,p3],'CBV LH','CBV RH','NE','ACh','Location','northeast','FontSize',5)
xlim([-30,30])
title('Awake to NREM')
set(gca,'box','off')
ax1.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
ax1.TickLength = [0.03,0.03];
axis square

% NREM to Awake
ax2 = subplot(2,2,2);
% p1 = plot(T1,data.NREMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on
% plot(T1,data.NREMtoAWAKE.AChCBV_Mean + data.NREMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

p4 = plot(T1,data.NREMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.NECBV_Mean + data.NREMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F (%)'})

ylim([-1,1])

% yyaxis right
p2 =  plot(T1,data.NREMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.NEGFP_Mean + data.NREMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% p3 =  plot(T1,data.NREMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% plot(T1,data.NREMtoAWAKE.AChGFP_Mean + data.NREMtoAWAKE.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% ylim([-1.5 1.5])
% legend([p1,p4,p2,p3],'CBV LH','CBV RH','NE','ACh','Location','northeast','FontSize',5)
% 
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
% p1 = plot(T1,data.NREMtoREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on

% plot(T1,data.NREMtoREM.AChCBV_Mean + data.NREMtoREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
p4 = plot(T1,data.NREMtoREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
ylim([-1,1])

plot(T1,data.NREMtoREM.NECBV_Mean + data.NREMtoREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
% ylim([-4,4])

% yyaxis right
p2 =  plot(T1,data.NREMtoREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.NEGFP_Mean + data.NREMtoREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% ylabel({'\DeltaF/F NE Percentage'},'rotation',-90,'VerticalAlignment','bottom')
% p3 =  plot(T1,data.NREMtoREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% plot(T1,data.NREMtoREM.AChGFP_Mean + data.NREMtoREM.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% ylim([-3 3])
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
% % p1 = plot(T1,data.REMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on
% plot(T1,data.REMtoAWAKE.AChCBV_Mean + data.REMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

p4 = plot(T1,data.REMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.NECBV_Mean + data.REMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylim([-1,1])

% yyaxis right
p2 =  plot(T1,data.REMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.NEGFP_Mean + data.REMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
% p3 =  plot(T1,data.REMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% plot(T1,data.REMtoAWAKE.AChGFP_Mean + data.REMtoAWAKE.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
% ylim([-3 3])
% legend([p1,p2,p3],'CBV','NE','ACh','Location','northeast','FontSize',5)

title('REM to Awake')
set(gca,'box','off')
ax4.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
ax4.TickLength = [0.03,0.03];
xlim([-30,30])
axis square

%}
% figure;
% subplot(2,2,4)
% plot(T1,data.REMtoAWAKE.AChGFPRaw,'-','color',[0 0.4470 0.7410],'LineWidth',2);

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
    savefig(summaryFigure_D,[dirpath 'Transition_CBV_vs_GFP_singleMouse_Sleep_consolidated_CI']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Transition_CBV_vs_GFP_singleMouse_Sleep_consolidated_CI'])
end
close

end
