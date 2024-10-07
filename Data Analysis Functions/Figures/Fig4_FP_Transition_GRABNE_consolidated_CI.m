function [AnalysisResults] = Fig4_FP_Transition_GRABNE_consolidated_CI(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);
%% check if there is any REM transitions
    transitionsCount = 0;
    
    for aa = 1:length(FPanimalIDs)
        animalID = FPanimalIDs{1,aa};
    
            if isfield(AnalysisResults.(animalID).Transitions,'NREMtoREM') || isfield(AnalysisResults.(animalID).Transitions,'REMtoAwake') % there is REM transitions
                transitionsCount = transitionsCount + 1;
            end
    end
%% decide what transitions are available
    if transitionsCount > 0 % there is REM transitions
        transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    end
    
    if transitionsCount == 0 % there is no REM transitions
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
    end
%% IOS mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
    
    if isfield(AnalysisResults.(animalID).Transitions,'NREMtoREM')
        transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    end

    if ~isfield(AnalysisResults.(animalID).Transitions,'NREMtoREM')
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
    end

    for bb = 1:length(transitions)
        transition = transitions{1,bb};

        % the size of the matrix
        MatLength = size(AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw,1);
        MatLength2 = 2*size(AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw,1);
        if aa == 1
            DataLength = 0;
            DataLength2 = 0;
        else
            DataLength = size(data.(transition).AChCBVRaw,1);
            DataLength2 = 2*size(data.(transition).AChCBVRaw,1);
        end           

        % data.(transition).EMG(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;

        % data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        % data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        % data.(transition).LH_Cort(:,:,DataLength+1:DataLength+MatLength) = AnalysisResults.(animalID).Transitions.(transition).LH_Cort;

        data.(transition).AChCBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw;
        data.(transition).AChGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_AChRaw;
        
        data.(transition).NECBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).NE_CBVRaw;
        data.(transition).NEGFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Transitions.(transition).GRAB_NERaw;

        data.(transition).CBVRaw(DataLength2+1:DataLength2+MatLength2,:) = cat(1,AnalysisResults.(animalID).Transitions.(transition).NE_CBVRaw,AnalysisResults.(animalID).Transitions.(transition).ACh_CBVRaw);

    end
end
% take average for each behavioral transition
    if isfield(data,'NREMtoREM')
        transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
    end

    if ~isfield(data,'NREMtoREM')
        transitions = {'AWAKEtoNREM','NREMtoAWAKE'};
    end

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

    data.(transition).CBV_N_Exp = size(data.(transition).CBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(transition).CBV_Mean = mean(data.(transition).CBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(transition).CBV_SEM = std(data.(transition).CBVRaw,1)/sqrt(data.(transition).CBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(transition).CBV_CI95 = tinv([0.025 0.975], data.(transition).CBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(transition).CBV_yCI95 = bsxfun(@times, data.(transition).CBV_SEM, data.(transition).CBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(transition).mean_LH_Cort = mean(data.(transition).LH_Cort,3)*100;
end
T1 = -30 + (1/30):(1/30):30;
%% Fig. 4 Plot CBV and GFP Transition with Cortical LFPs
summaryFigure_D = figure('Name','Fig4 GFP vs CBV');
sgtitle('Blood Volume and Neuromodulators Transition')
%% Awake to NREM

% Blood Volume and NE, ACh
    ax1 = subplot(2,2,1);
    p1 = plot(T1,data.AWAKEtoNREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
    hold on
    plot(T1,data.AWAKEtoNREM.AChCBV_Mean + data.AWAKEtoNREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
    % ylabel({'\DeltaF/F (%)'})
    % 
    p4 = plot(T1,data.AWAKEtoNREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
    hold on
    plot(T1,data.AWAKEtoNREM.NECBV_Mean + data.AWAKEtoNREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
    % ylabel({'\DeltaF/F (%)'})
    
    % pC = plot(T1,data.AWAKEtoNREM.CBV_Mean,'-','color','red','LineWidth',2);
    % hold on
    % plot(T1,data.AWAKEtoNREM.CBV_Mean + data.AWAKEtoNREM.CBV_yCI95,'-','color','red','LineWidth',0.1);
    
    ylabel({'\DeltaF/F (%)'})
    ylim([-1 1])
    yticks(-1:.25:1)
        
    yyaxis right
    p2 =  plot(T1,data.AWAKEtoNREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    plot(T1,data.AWAKEtoNREM.NEGFP_Mean + data.AWAKEtoNREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
    
    p3 =  plot(T1,data.AWAKEtoNREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
    hold on
    plot(T1,data.AWAKEtoNREM.AChGFP_Mean + data.AWAKEtoNREM.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
    ylim([-4 4])
    yticks(-4:1:4)
    
    % pq = plot(T1,data.AWAKEtoNREM.AChGFP_Mean - data.AWAKEtoNREM.NEGFP_Mean,'-','color',[0 0.2 0.1],'LineWidth',2);
        
    legend([p2,p3,p1,p4],'NE','ACh','LH CBV','RH CBV','Location','northeast','FontSize',5)
    xlim([-30,30])
    title('Awake to NREM')
    set(gca,'box','off')
    ax1.YAxis(1).Color = 'k';
    ax1.YAxis(2).Color = 'k';
    axis square
%% NREM to Awake
    ax2 = subplot(2,2,2);
    p1 = plot(T1,data.NREMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoAWAKE.AChCBV_Mean + data.NREMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

    p4 = plot(T1,data.NREMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoAWAKE.NECBV_Mean + data.NREMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
    
    % pC = plot(T1,data.NREMtoAWAKE.CBV_Mean,'-','color','red','LineWidth',2);
    % hold on
    % plot(T1,data.NREMtoAWAKE.CBV_Mean + data.NREMtoAWAKE.CBV_yCI95,'-','color','red','LineWidth',0.1);
    
    ylabel({'\DeltaF/F (%)'})
    ylim([-1 1])
    yticks(-1:.25:1)
    
    yyaxis right
    p2 =  plot(T1,data.NREMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoAWAKE.NEGFP_Mean + data.NREMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
    p3 =  plot(T1,data.NREMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoAWAKE.AChGFP_Mean + data.NREMtoAWAKE.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);

    % pq = plot(T1,data.NREMtoAWAKE.AChGFP_Mean - data.NREMtoAWAKE.NEGFP_Mean,'-','color',[0 0.2 0.1],'LineWidth',2);
    
    ylim([-4 4])
    yticks(-4:1:4)
    ylabel('\Delta F/F (%)')

    title('NREM to Awake')
    set(gca,'box','off')
    ax2.YAxis(1).Color = 'k';
    ax2.YAxis(2).Color = 'k';
    xlim([-30,30])
    axis square
%% if REM transitions are available plot those data
if isfield(data,'NREMtoREM')
    % NREM to REM
    
    % Blood Volume and NE, ACh
    ax3 = subplot(2,2,3);
    p1 = plot(T1,data.NREMtoREM.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
    hold on 
    plot(T1,data.NREMtoREM.AChCBV_Mean + data.NREMtoREM.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);

    p4 = plot(T1,data.NREMtoREM.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoREM.NECBV_Mean + data.NREMtoREM.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);

    % pC = plot(T1,data.NREMtoREM.CBV_Mean,'-','color','red','LineWidth',2);
    % hold on
    % plot(T1,data.NREMtoREM.CBV_Mean + data.NREMtoREM.CBV_yCI95,'-','color','red','LineWidth',0.1);
    
    ylabel('\Delta F/F (%)')
    ylim([-5 5])
    yticks(-5:1:5)
    
    yyaxis right
    p2 =  plot(T1,data.NREMtoREM.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoREM.NEGFP_Mean + data.NREMtoREM.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
    p3 =  plot(T1,data.NREMtoREM.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
    hold on
    plot(T1,data.NREMtoREM.AChGFP_Mean + data.NREMtoREM.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
    ylim([-11 11])
    yticks(-11:2:11)
    ylabel('\Delta F/F (%)')

    xlim([-30,30])
    title('NREM to REM')
    set(gca,'box','off')
    ax3.YAxis(1).Color = 'k';
    ax3.YAxis(2).Color = 'k';
    axis square
    % px =  plot(T1,data.NREMtoREM.AChGFP_Mean - data.NREMtoREM.NEGFP_Mean,'-','color',[0 0.2 0.1],'LineWidth',2);
end
if isfield(data,'REMtoAWAKE')
    %% REM to Awake
    ax4 = subplot(2,2,4);
    p1 = plot(T1,data.REMtoAWAKE.AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.AChCBV_Mean + data.REMtoAWAKE.AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.1);
    % 
    p4 = plot(T1,data.REMtoAWAKE.NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.NECBV_Mean + data.REMtoAWAKE.NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.1);
    
    % pC = plot(T1,data.REMtoAWAKE.CBV_Mean,'-','color','red','LineWidth',2);
    % hold on
    % plot(T1,data.REMtoAWAKE.CBV_Mean + data.REMtoAWAKE.CBV_yCI95,'-','color','red','LineWidth',0.1);
    
    ylabel('\Delta F/F (%)')
    ylim([-5 5])
    yticks(-5:1:5)
    
    yyaxis right
    p2 =  plot(T1,data.REMtoAWAKE.NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.NEGFP_Mean + data.REMtoAWAKE.NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
    p3 =  plot(T1,data.REMtoAWAKE.AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
    hold on
    plot(T1,data.REMtoAWAKE.AChGFP_Mean + data.REMtoAWAKE.AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.1);
    ylim([-11 11])
    yticks(-11:2:11)
    ylabel('\Delta F/F (%)')

    % px =  plot(T1,data.REMtoAWAKE.AChGFP_Mean - data.REMtoAWAKE.NEGFP_Mean,'-','color',[0 0.2 0.1],'LineWidth',2);
        
    title('REM to Awake')
    xlabel('time(s)')
    set(gca,'box','off')
    ax4.YAxis(1).Color = 'k';
    ax4.YAxis(2).Color = 'k';
    xlim([-30,30])
    axis square
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
    savefig(summaryFigure_D,[dirpath ManipulationType '_Bilateral_Transition_CBV_vs_GFP_Sleep_CI_ACh_NE']);
    set(summaryFigure_D,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_Bilateral_Transition_CBV_vs_GFP_Sleep_CI_ACh_NE'])
end
close

end
