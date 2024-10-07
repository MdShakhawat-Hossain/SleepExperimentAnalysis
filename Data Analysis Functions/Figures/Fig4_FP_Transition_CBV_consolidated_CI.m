function [AnalysisResults] = Fig4_FP_Transition_CBV_consolidated_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,IOSanimalIDs)
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
% take average for each behavioral transition _LH_
for cc = 1:length(transitions)
    transition = transitions{1,cc};
    % data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).AChCBV_N_Exp = size(data.(transition).AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChCBV_Mean = mean(data.(transition).AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AChCBV_SEM = std(data.(transition).AChCBVRaw,1)/sqrt(data.(transition).AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AChCBV_CI95 = tinv([0.025 0.975], data.(transition).AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChCBV_yCI95 = bsxfun(@times, data.(transition).AChCBV_SEM, data.(transition).AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
    
    data.(transition).NECBV_N_Exp = size(data.(transition).NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NECBV_Mean = mean(data.(transition).NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NECBV_SEM = std(data.(transition).NECBVRaw,1)/sqrt(data.(transition).NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NECBV_CI95 = tinv([0.025 0.975], data.(transition).NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NECBV_yCI95 = bsxfun(@times, data.(transition).NECBV_SEM, data.(transition).NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).AChGFP_N_Exp = size(data.(transition).AChGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).AChGFP_Mean = mean(data.(transition).AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).AChGFP_SEM = std(data.(transition).AChGFPRaw,1)/sqrt(data.(transition).AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).AChGFP_CI95 = tinv([0.025 0.975], data.(transition).AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).AChGFP_yCI95 = bsxfun(@times, data.(transition).AChGFP_SEM, data.(transition).AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(transition).NEGFP_N_Exp = size(data.(transition).NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).NEGFP_Mean = mean(data.(transition).NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_SEM = std(data.(transition).NEGFPRaw,1)/sqrt(data.(transition).NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).NEGFP_CI95 = tinv([0.025 0.975], data.(transition).NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).NEGFP_yCI95 = bsxfun(@times, data.(transition).NEGFP_SEM, data.(transition).NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
%% Fig. 4 Plot CBV and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs CBV');
sgtitle('Blood Volume and Neuromodulators Transition')

% Awake to NREM

% Blood Volume and NE, ACh
ax1 = subplot(2,2,1);
p1 = plot(T1,data.AWAKEtoNREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.AChCBV_Mean + data.AWAKEtoNREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
ylabel({'\DeltaF/F (%)'})

p4 = plot(T1,data.AWAKEtoNREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.NECBV_Mean + data.AWAKEtoNREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F (%)'})
ylim([-.5 .5])
yticks(-.5:.25:.5)

legend([p1,p4],'CBV LH','CBV RH','Location','northeast','FontSize',5)
xlim([-30,30])
title('Awake to NREM')
set(gca,'box','off')
% ax1.YAxis(1).Color = 'k';
% ax1.YAxis(2).Color = 'k';
% set(gca,'TickLength',[0,0])
axis square

% NREM to Awake
ax2 = subplot(2,2,2);
p1 = plot(T1,data.NREMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.AChCBV_Mean + data.NREMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

p4 = plot(T1,data.NREMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.NECBV_Mean + data.NREMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylabel({'\DeltaF/F (%)'})
ylim([-.5 .5])
yticks(-.5:.25:.5)

title('NREM to Awake')
set(gca,'box','off')
% ax2.YAxis(1).Color = 'k';
% ax2.YAxis(2).Color = 'k';
% set(gca,'TickLength',[0,0])
xlim([-30,30])
axis square
%%

% NREM to REM

% Blood Volume and NE, ACh
ax3 = subplot(2,2,3);
p1 = plot(T1,data.NREMtoREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.AChCBV_Mean + data.NREMtoREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
p4 = plot(T1,data.NREMtoREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on

plot(T1,data.NREMtoREM.NECBV_Mean + data.NREMtoREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylim([-1.2,1.2])
% set(gca,'TickLength',[0,0])

xlim([-30,30])
title('NREM to REM')
set(gca,'box','off')
% ax3.YAxis(1).Color = 'k';
% ax3.YAxis(2).Color = 'k';
% set(gca,'TickLength',[0,0])
axis square


% REM to Awake
ax4 = subplot(2,2,4);
p1 = plot(T1,data.REMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.AChCBV_Mean + data.REMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

p4 = plot(T1,data.REMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.NECBV_Mean + data.REMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
ylim([-1.2 1.2])

title('REM to Awake')
xlabel('time(s)')
set(gca,'box','off')
% ax4.YAxis(1).Color = 'k';
% ax4.YAxis(2).Color = 'k';
xlim([-30,30])
axis square


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
    savefig(summaryFigure_D,[dirpath 'Transition_CBV_Sleep_consolidated_CI']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Transition_CBV_Sleep_consolidated_CI'])
end
close

end
