function [AnalysisResults] = Fig1_S5_Movement_GRABNE_Response_CI_OnlyCBV(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%________________________________________________________________________________________________________________________
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);
%% set-up and process data
movementDataTypes = {'ShortMovement','IntermediateMovement'};%,'LongMovement'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(movementDataTypes)
        movementDataType = movementDataTypes{1,cc};
         % the size of the matrix
        MatLength = size(AnalysisResults.(animalID).Movement.P_NE.(movementDataType).CBV.CBVRaw,1);
        if aa == 1
           DataLength = 0;
        else
           DataLength = size(data.(movementDataType).P_NE.CBVRaw,1);
        end

        data.(movementDataType).P_ACh.CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Movement.P_ACh.(movementDataType).CBV.CBVRaw;
        % data.(movementDataType).P_ACh.GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Movement.P_ACh.(movementDataType).GFP.GFPRaw;
        data.(movementDataType).P_NE.CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Movement.P_NE.(movementDataType).CBV.CBVRaw;
        % data.(movementDataType).P_NE.GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Movement.P_NE.(movementDataType).GFP.GFPRaw;
        % time vector
        try
            data.(movementDataType).timeVector(:,aa) = AnalysisResults.(animalID).Movement.cortical_LH.(movementDataType).timeVector;
        catch
            data.(movementDataType).timeVector(:,aa) = AnalysisResults.(animalID).Movement.cortical.(movementDataType).timeVector;
        end

    end
end
% 
for ee = 1:length(movementDataTypes)
    movementDataType = movementDataTypes{1,ee};
    data.(movementDataType).P_AChCBVRaw = data.(movementDataType).P_ACh.CBVRaw;
    % data.(movementDataType).P_AChGFPRaw = data.(movementDataType).P_ACh.GFPRaw;
    data.(movementDataType).P_NECBVRaw = data.(movementDataType).P_NE.CBVRaw;
    % data.(movementDataType).P_NEGFPRaw = data.(movementDataType).P_NE.GFPRaw;
    data.(movementDataType).P_CBVRaw = cat(1,data.(movementDataType).P_NE.CBVRaw,data.(movementDataType).P_ACh.CBVRaw);
end
% mean
for ee = 1:length(movementDataTypes)
    movementDataType = movementDataTypes{1,ee};

    data.(movementDataType).P_AChCBV_N_Exp = size(data.(movementDataType).P_AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(movementDataType).P_AChCBV_Mean = mean(data.(movementDataType).P_AChCBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_AChCBV_SEM = std(data.(movementDataType).P_AChCBVRaw,1)/sqrt(data.(movementDataType).P_AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_AChCBV_CI95 = tinv([0.025 0.975], data.(movementDataType).P_AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(movementDataType).P_AChCBV_yCI95 = bsxfun(@times, data.(movementDataType).P_AChCBV_SEM, data.(movementDataType).P_AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(movementDataType).P_NECBV_N_Exp = size(data.(movementDataType).P_NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(movementDataType).P_NECBV_Mean = mean(data.(movementDataType).P_NECBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_NECBV_SEM = std(data.(movementDataType).P_NECBVRaw,1)/sqrt(data.(movementDataType).P_NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_NECBV_CI95 = tinv([0.025 0.975], data.(movementDataType).P_NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(movementDataType).P_NECBV_yCI95 = bsxfun(@times, data.(movementDataType).P_NECBV_SEM, data.(movementDataType).P_NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    % data.(movementDataType).P_AChGFP_N_Exp = size(data.(movementDataType).P_AChGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    % data.(movementDataType).P_AChGFP_Mean = mean(data.(movementDataType).P_AChGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    % data.(movementDataType).P_AChGFP_SEM = std(data.(movementDataType).P_AChGFPRaw,1)/sqrt(data.(movementDataType).P_AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    % data.(movementDataType).P_AChGFP_CI95 = tinv([0.025 0.975], data.(movementDataType).P_AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    % data.(movementDataType).P_AChGFP_yCI95 = bsxfun(@times, data.(movementDataType).P_AChGFP_SEM, data.(movementDataType).P_AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’
    % 
    % data.(movementDataType).P_NEGFP_N_Exp = size(data.(movementDataType).P_NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    % data.(movementDataType).P_NEGFP_Mean = mean(data.(movementDataType).P_NEGFPRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    % data.(movementDataType).P_NEGFP_SEM = std(data.(movementDataType).P_NEGFPRaw,1)/sqrt(data.(movementDataType).P_NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    % data.(movementDataType).P_NEGFP_CI95 = tinv([0.025 0.975], data.(movementDataType).P_NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    % data.(movementDataType).P_NEGFP_yCI95 = bsxfun(@times, data.(movementDataType).P_NEGFP_SEM, data.(movementDataType).P_NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(movementDataType).P_CBV_N_Exp = size(data.(movementDataType).P_CBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(movementDataType).P_CBV_Mean = mean(data.(movementDataType).P_CBVRaw,1);                                   % Mean Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_CBV_SEM = std(data.(movementDataType).P_CBVRaw,1)/sqrt(data.(movementDataType).P_CBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At Each Value Of _x’
    data.(movementDataType).P_CBV_CI95 = tinv([0.025 0.975], data.(movementDataType).P_CBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(movementDataType).P_CBV_yCI95 = bsxfun(@times, data.(movementDataType).P_CBV_SEM, data.(movementDataType).P_CBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of _x’

    data.(movementDataType).meanTimeVector = mean(data.(movementDataType).timeVector(:,aa),2);
end

%% Only the fiber signals consolidated
summaryFigureN = figure('Name','Movement_fiber_consolidated');
sgtitle('Movement evoked responses in fiber photometry signals')

% Short Movement
ax1 = subplot(1,2,1);
% p1 = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_AChCBV_Mean + data.ShortMovement.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% 
% p2 = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
% hold on
% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_NECBV_Mean + data.ShortMovement.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

pC = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_CBV_Mean,'-','color','red','LineWidth',2);
hold on
NEn =  data.ShortMovement.P_CBV_Mean + data.ShortMovement.P_CBV_yCI95;
patch([data.ShortMovement.meanTimeVector' fliplr(data.ShortMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],'red','FaceColor','red','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_CBV_Mean + data.ShortMovement.P_CBV_yCI95,'-','color','red','LineWidth',0.10)

title('Brief movement')
ylabel('\DeltaF/F CBV (%)')
ylim([-0.5 3])

% yyaxis right
% p3 = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% NEn =  data.ShortMovement.P_AChGFP_Mean + data.ShortMovement.P_AChGFP_yCI95;
% patch([data.ShortMovement.meanTimeVector' fliplr(data.ShortMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0 0.4470 0.7410],'FaceColor',[0 0.4470 0.7410],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);
% 
% % plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_AChGFP_Mean + data.ShortMovement.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% p4 = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
% hold on
% NEn =  data.ShortMovement.P_NEGFP_Mean + data.ShortMovement.P_NEGFP_yCI95;
% patch([data.ShortMovement.meanTimeVector' fliplr(data.ShortMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.4660 0.6740 0.1880],'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_NEGFP_Mean + data.ShortMovement.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

% pq = plot(data.ShortMovement.meanTimeVector,data.ShortMovement.P_AChGFP_Mean - data.ShortMovement.P_NEGFP_Mean,'-','color','black','LineWidth',2);


% ylabel('\DeltaF/F (%)')

% ax1.YAxis(1).Color = 'k';%[0.8500 0.3250 0.0980];
% ax1.YAxis(2).Color = 'k';%[0 0.4470 0.7410];
% legend([p1, p2, p3, p4, pq, pC],'CBV-LH','CBV-RH', 'ACh', 'NE', 'ACh-NE', 'CBV')
legend([pC],'CBV')

xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ylim([-0.5 3])
xlim([-5 15])


% Intermediate Movement
ax2 = subplot(1,2,2);
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_AChCBV_Mean + data.IntermediateMovement.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_NECBV_Mean + data.IntermediateMovement.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_CBV_Mean,'-','color','red','LineWidth',2)
hold on
NEn =  data.IntermediateMovement.P_CBV_Mean + data.IntermediateMovement.P_CBV_yCI95;
patch([data.IntermediateMovement.meanTimeVector' fliplr(data.IntermediateMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],'red','FaceColor','red','EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);

% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_CBV_Mean + data.IntermediateMovement.P_CBV_yCI95,'-','color','red','LineWidth',0.10)
title('Moderate movement')
ylabel('\DeltaF/F CBV (%)')
ax2.YLim = [-0.5 3];

% yyaxis right
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% NEn =  data.IntermediateMovement.P_AChGFP_Mean + data.IntermediateMovement.P_AChGFP_yCI95;
% patch([data.IntermediateMovement.meanTimeVector' fliplr(data.IntermediateMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0 0.4470 0.7410],'FaceColor',[0 0.4470 0.7410],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);
% 
% % plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_AChGFP_Mean + data.IntermediateMovement.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% NEn =  data.IntermediateMovement.P_NEGFP_Mean + data.IntermediateMovement.P_NEGFP_yCI95;
% patch([data.IntermediateMovement.meanTimeVector' fliplr(data.IntermediateMovement.meanTimeVector')], [NEn(1,:) fliplr(NEn(2,:))],[0.4660 0.6740 0.1880],'FaceColor',[0.4660 0.6740 0.1880],'EdgeColor','none','LineWidth',0.1,'FaceAlpha',0.2);
% 
% % plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_NEGFP_Mean + data.IntermediateMovement.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% ylabel('\DeltaF/F (%)')

% plot(data.IntermediateMovement.meanTimeVector,data.IntermediateMovement.P_AChGFP_Mean - data.IntermediateMovement.P_NEGFP_Mean,'-','color','black','LineWidth',2)

xlabel('Peri-movement time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-0.5 3];
xlim([-5 15])


% Long Movement
% ax3 = subplot(1,3,3);
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChCBV_Mean + data.LongMovement.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChCBV_Mean - data.LongMovement.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NECBV_Mean + data.LongMovement.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NECBV_Mean - data.LongMovement.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% 
% title('Long movement')
% ylabel('\DeltaF/F Blood Volume')
% ax3.YLim = [-0.5 3];
% 
% yyaxis right
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChGFP_Mean + data.LongMovement.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_AChGFP_Mean - data.LongMovement.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NEGFP_N_Exp,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NEGFP_N_Exp + data.LongMovement.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongMovement.meanTimeVector,data.LongMovement.P_NEGFP_N_Exp - data.LongMovement.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% 
% ylabel('\DeltaF/F GRAB(Z)')
% ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax3.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-movement time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ax3.YLim = [-0.5 3];
% xlim([-5 15])

% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath ManipulationType '_Movement_Response_CI_ACh_NE']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_Movement_Response_CI_ACh_NE'])
    close 
end
end
