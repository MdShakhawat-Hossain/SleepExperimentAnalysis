function [AnalysisResults] = Whisk_GRABNE_Response_2Photon(rootFolder,saveFigs,delim,AnalysisResults,animalIDs)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
%________________________________________________________________________________________________________________________

%% set-up and process data
% animalIDs = {'NEACh008'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        data.(whiskDataType).Z_NE.GFP(:,aa) = AnalysisResults.(animalID).Whisk.Z_NE.(whiskDataType).GFP.GFP;

        % data.(whiskDataType).cortical.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).MUA.corticalData;
        % data.(whiskDataType).cortical.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).Gam.corticalData;
        % data.(whiskDataType).cortical.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS;
        % data.(whiskDataType).cortical.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.corticalS(49:end,20:23);
        % data.(whiskDataType).cortical.cortT(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.T;
        % data.(whiskDataType).cortical.cortF(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.F;
        % hippocampal
%         data.(whiskDataType).Hip.hipMUA(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).MUA.hippocampalData;
%         data.(whiskDataType).Hip.hipGam(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).Gam.hippocampalData;
%         data.(whiskDataType).Hip.hipS(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.hippocampalS;
%         data.(whiskDataType).Hip.hipS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.hippocampalS(49:end,20:23);
%         data.(whiskDataType).Hip.hipT(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.T;
%         data.(whiskDataType).Hip.hipF(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).LFP.F;
        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.cortical.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).Z_NEGFP = data.(whiskDataType).Z_NE.GFP;

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

    data.(whiskDataType).Z_NEmeanGFP = mean(data.(whiskDataType).Z_NEGFP,2);
    data.(whiskDataType).Z_NEstdGFP = std(data.(whiskDataType).Z_NEGFP,0,2);

    % data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    % data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    % data.(whiskDataType).meanCortGam = mean(data.(whiskDataType).cortGam,2);
    % data.(whiskDataType).stdCortGam = std(data.(whiskDataType).cortGam,0,2);
    % data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    % data.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),3);
    % data.(whiskDataType).std_CortS_Gam = std(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),0,3);
    % data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    % data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
%     data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
%     data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
%     data.(whiskDataType).meanHipGam = mean(data.(whiskDataType).Hip.hipGam,2);
%     data.(whiskDataType).stdHipGam = std(data.(whiskDataType).Hip.hipGam,0,2);
%     data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3).*100;
%     data.(whiskDataType).mean_HipS_Gam = mean(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),3);
%     data.(whiskDataType).std_HipS_Gam = std(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),0,3);
%     data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
%     data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end




%% Fig. 1-S3
summaryFigure = figure('Name','Whisk-GRABNE');
sgtitle('Whisk evoked responses')
%% brief whisks GFP
ax16 = subplot(1,3,1);

plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP + data.ShortWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.Z_NEmeanGFP - data.ShortWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')

title('Brief whisk GRABNE')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
% ax16.YLim = [-1 1];
%% moderate whisks GFP
ax17 = subplot(1,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP + data.IntermediateWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.Z_NEmeanGFP - data.IntermediateWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
title('Moderate whisk GRABNE')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
% ax17.YLim = [-1 1];
%%  extended whisks GFP
ax18 = subplot(1,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP + data.LongWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.Z_NEmeanGFP - data.LongWhisks.Z_NEstdGFP,'-','color',[0.4660 0.6740 0.1880],'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
title('Extended whisk GRABNE')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
% ax18.YLim = [-1 1];

%% save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim ];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim ];
    end
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Whisk-GRABNE']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Whisk-GRABNE'])
    close 
end

end
