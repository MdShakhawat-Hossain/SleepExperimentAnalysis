function [AnalysisResults] = Fig1_S3_Whisk_Flavoprotein_Response_CI_Flavoprotein(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%________________________________________________________________________________________________________________________

%% set-up and process data
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through eACh animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
        for cc = 1:length(whiskDataTypes)
            whiskDataType = whiskDataTypes{1,cc};
            % the size of the matrix
            MatLength = size(AnalysisResults.(animalID).Whisk.P_NE.(whiskDataType).CBV.CBVRaw,1);
            if aa == 1
                DataLength = 0;
            else
                DataLength = size(data.(whiskDataType).P_NE.CBVRaw,1);
            end
            
            data.(whiskDataType).P_NE.CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Whisk.P_NE.(whiskDataType).CBV.CBVRaw;
            data.(whiskDataType).P_ACh.CBVRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Whisk.P_ACh.(whiskDataType).CBV.CBVRaw;
            data.(whiskDataType).P_NE.GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Whisk.P_NE.(whiskDataType).GFP.GFPRaw;
            data.(whiskDataType).P_ACh.GFPRaw(DataLength+1:DataLength+MatLength,:) = AnalysisResults.(animalID).Whisk.P_ACh.(whiskDataType).GFP.GFPRaw;
            data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).timeVector;

        end
end

%% raw data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).P_AChCBVRaw = data.(whiskDataType).P_ACh.CBVRaw;
    data.(whiskDataType).P_AChGFPRaw = data.(whiskDataType).P_ACh.GFPRaw;
    data.(whiskDataType).P_NECBVRaw = data.(whiskDataType).P_NE.CBVRaw;
    data.(whiskDataType).P_NEGFPRaw = data.(whiskDataType).P_NE.GFPRaw;

    % data.(whiskDataType).cortMUA = data.(whiskDataType).cortical.cortMUA;
    % data.(whiskDataType).cortGam = data.(whiskDataType).cortical.cortGam;
    % data.(whiskDataType).cortS = data.(whiskDataType).cortical.cortS;
    % data.(whiskDataType).cortS_Gam = data.(whiskDataType).cortical.cortS_Gam;
    % data.(whiskDataType).cortT = data.(whiskDataType).cortical.cortT;
    % data.(whiskDataType).cortF = data.(whiskDataType).cortical.cortF;
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
      
    data.(whiskDataType).P_AChCBV_N_Exp = size(data.(whiskDataType).P_AChCBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(whiskDataType).P_AChCBV_Mean = mean(data.(whiskDataType).P_AChCBVRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_AChCBV_SEM = std(data.(whiskDataType).P_AChCBVRaw,1)/sqrt(data.(whiskDataType).P_AChCBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_AChCBV_CI95 = tinv([0.025 0.975], data.(whiskDataType).P_AChCBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(whiskDataType).P_AChCBV_yCI95 = bsxfun(@times, data.(whiskDataType).P_AChCBV_SEM, data.(whiskDataType).P_AChCBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    data.(whiskDataType).P_NECBV_N_Exp = size(data.(whiskDataType).P_NECBVRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(whiskDataType).P_NECBV_Mean = mean(data.(whiskDataType).P_NECBVRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_NECBV_SEM = std(data.(whiskDataType).P_NECBVRaw,1)/sqrt(data.(whiskDataType).P_NECBV_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_NECBV_CI95 = tinv([0.025 0.975], data.(whiskDataType).P_NECBV_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(whiskDataType).P_NECBV_yCI95 = bsxfun(@times, data.(whiskDataType).P_NECBV_SEM, data.(whiskDataType).P_NECBV_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    data.(whiskDataType).P_AChGFP_N_Exp = size(data.(whiskDataType).P_AChGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(whiskDataType).P_AChGFP_Mean = mean(data.(whiskDataType).P_AChGFPRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_AChGFP_SEM = std(data.(whiskDataType).P_AChGFPRaw,1)/sqrt(data.(whiskDataType).P_AChGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_AChGFP_CI95 = tinv([0.025 0.975], data.(whiskDataType).P_AChGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(whiskDataType).P_AChGFP_yCI95 = bsxfun(@times, data.(whiskDataType).P_AChGFP_SEM, data.(whiskDataType).P_AChGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’

    data.(whiskDataType).P_NEGFP_N_Exp = size(data.(whiskDataType).P_NEGFPRaw,1);                                      % Number of _ExperimentsE In Data Set
    data.(whiskDataType).P_NEGFP_Mean = mean(data.(whiskDataType).P_NEGFPRaw,1);                                   % Mean Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_NEGFP_SEM = std(data.(whiskDataType).P_NEGFPRaw,1)/sqrt(data.(whiskDataType).P_NEGFP_N_Exp);                              % Compute tStandard Error Of The Meand Of All Experiments At EACh Value Of _x’
    data.(whiskDataType).P_NEGFP_CI95 = tinv([0.025 0.975], data.(whiskDataType).P_NEGFP_N_Exp-1);                    % Calculate 95% Probability Intervals Of t-Distribution
    data.(whiskDataType).P_NEGFP_yCI95 = bsxfun(@times, data.(whiskDataType).P_NEGFP_SEM, data.(whiskDataType).P_NEGFP_CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At EACh Value Of _x’
    
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end
%% Only the fiber signals consolidated
summaryFigureN = figure('Name','Whisk_fiber_consolidated');
sgtitle('Whisking evoked responses in fiber photometry signals')

% Short Movement
ax1 = subplot(1,2,1);
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_AChCBV_Mean + data.ShortWhisks.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)


plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_NECBV_Mean + data.ShortWhisks.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
ylabel('\DeltaF/F (%)')
% ax1.YLim = [-.5 2];
title('Brief Whisks')
% yyaxis right
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_AChGFP_Mean + data.ShortWhisks.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
%
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.P_NEGFP_Mean + data.ShortWhisks.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)

ylabel('\DeltaF/F (%)')

% ax1.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax1.YAxis(2).Color = [0 0.4470 0.7410];
xlabel('Peri-Whisks time (s)')
axis square
set(gca,'box','off')
ylim([-.3 .3])
xlim([-5 15])


% Intermediate Movement
ax2 = subplot(1,2,2);
% p1 = plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2);
% hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_AChCBV_Mean + data.IntermediateWhisks.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)

p2 = plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_NECBV_Mean + data.IntermediateWhisks.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)

title('Moderate Whisks')
ylabel('\DeltaF/F (%)')
% ax2.YLim = [-.5 2];

% yyaxis right
% p3 = plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2);
% hold on
% plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_AChGFP_Mean + data.IntermediateWhisks.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)

p4 = plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.P_NEGFP_Mean + data.IntermediateWhisks.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F (%)')
% ax2.YAxis(1).Color = 'k';%[0.8500 0.3250 0.0980];
% ax2.YAxis(2).Color = 'k';%[;%0 0.4470 0.7410];
xlabel('Peri-Whisks time (s)')
axis square
set(gca,'box','off')
% ax2.TickLength = [0.03,0.03];
% ax2.YLim = [-.5 2];
ylim([-.3 .3])
xlim([-5 15])
legend([ p2, p4],'CBV-RH','NE')

% % Long Whisk
% ax3 = subplot(2,2,3);
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChCBV_Mean,'-','color',[0.8500 0.3250 0.0980],'LineWidth',2)
% hold on
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChCBV_Mean + data.LongWhisks.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChCBV_Mean - data.LongWhisks.P_AChCBV_yCI95,'-','color',[0.8500 0.3250 0.0980],'LineWidth',0.10)
% 
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NECBV_Mean,'-','color',[0.6350 0.0780 0.1840],'LineWidth',2)
% hold on
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NECBV_Mean + data.LongWhisks.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NECBV_Mean - data.LongWhisks.P_NECBV_yCI95,'-','color',[0.6350 0.0780 0.1840],'LineWidth',0.10)
% 
% title('Long Whisks')
% ylabel('\DeltaF/F Blood Volume (P)')
% % ax3.YLim = [-.5 2];
% 
% % yyaxis right
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChGFP_Mean,'-','color',[0 0.4470 0.7410],'LineWidth',2)
% hold on
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChGFP_Mean + data.LongWhisks.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_AChGFP_Mean - data.LongWhisks.P_AChGFP_yCI95,'-','color',[0 0.4470 0.7410],'LineWidth',0.10)
% 
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NEGFP_Mean,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2)
% hold on
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NEGFP_Mean + data.LongWhisks.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% plot(data.LongWhisks.meanTimeVector,data.LongWhisks.P_NEGFP_Mean - data.LongWhisks.P_NEGFP_yCI95,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
% 
% ylabel('\DeltaF/F(P)')
% % ax3.YAxis(1).Color = [0.8500 0.3250 0.0980];
% % ax3.YAxis(2).Color = [0 0.4470 0.7410];
% xlabel('Peri-Whisks time (s)')
% axis square
% set(gca,'box','off')
% ax3.TickLength = [0.03,0.03];
% ylim([-.5 2])
% xlim([-5 15])

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
    savefig(summaryFigureN,[dirpath 'Whisk-FiberSignals-consolidated_CI']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Whisk-FiberSignals-consolidated_CI'])
    close 
end
end
