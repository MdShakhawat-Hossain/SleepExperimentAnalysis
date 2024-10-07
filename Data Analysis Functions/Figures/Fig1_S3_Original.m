function [AnalysisResults] = Fig1_S3_Original(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs =  {'T282'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        % left cortical
        data.(whiskDataType).LH.HbT(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).HbT.HbT;
        data.(whiskDataType).LH.GCaMP7s(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).GCaMP7s.GCaMP7s;

        data.(whiskDataType).LH.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).LH.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).Gam.corticalData;
        data.(whiskDataType).LH.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).LH.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.(whiskDataType).LH.cortT(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.T;
        data.(whiskDataType).LH.cortF(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.F;
        % right cortical
        data.(whiskDataType).RH.HbT(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).HbT.HbT;
%         data.(whiskDataType).RH.CBV(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).HbT.CBV;
        data.(whiskDataType).RH.GCaMP7s(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).GCaMP7s.GCaMP7s;

        data.(whiskDataType).RH.cortMUA(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).RH.cortGam(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).Gam.corticalData;
        data.(whiskDataType).RH.cortS(:,:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).RH.cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.corticalS(49:end,20:23);
        data.(whiskDataType).RH.cortT(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.T;
        data.(whiskDataType).RH.cortF(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.F;
        % hippocampal
        data.(whiskDataType).Hip.hipMUA(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).MUA.hippocampalData;
        data.(whiskDataType).Hip.hipGam(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).Gam.hippocampalData;
        data.(whiskDataType).Hip.hipS(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.hippocampalS;
        data.(whiskDataType).Hip.hipS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.hippocampalS(49:end,20:23);
        data.(whiskDataType).Hip.hipT(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.T;
        data.(whiskDataType).Hip.hipF(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.F;
        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).HbT = cat(2,data.(whiskDataType).LH.HbT,data.(whiskDataType).RH.HbT);
    data.(whiskDataType).GCaMP7s = cat(2,data.(whiskDataType).LH.GCaMP7s,data.(whiskDataType).RH.GCaMP7s);

    data.(whiskDataType).cortMUA = cat(2,data.(whiskDataType).LH.cortMUA,data.(whiskDataType).RH.cortMUA);
    data.(whiskDataType).cortGam = cat(2,data.(whiskDataType).LH.cortGam,data.(whiskDataType).RH.cortGam);
    data.(whiskDataType).cortS = cat(3,data.(whiskDataType).LH.cortS,data.(whiskDataType).RH.cortS);
    data.(whiskDataType).cortS_Gam = cat(3,data.(whiskDataType).LH.cortS_Gam,data.(whiskDataType).RH.cortS_Gam);
    data.(whiskDataType).cortT = cat(2,data.(whiskDataType).LH.cortT,data.(whiskDataType).RH.cortT);
    data.(whiskDataType).cortF = cat(2,data.(whiskDataType).LH.cortF,data.(whiskDataType).RH.cortF);
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).meanHbT = mean(data.(whiskDataType).HbT,2);
    data.(whiskDataType).stdHbT = std(data.(whiskDataType).HbT,0,2);
    data.(whiskDataType).meanGCaMP7s = mean(data.(whiskDataType).GCaMP7s,2);
    data.(whiskDataType).stdGCaMP7s = std(data.(whiskDataType).GCaMP7s,0,2);

    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortGam = mean(data.(whiskDataType).cortGam,2);
    data.(whiskDataType).stdCortGam = std(data.(whiskDataType).cortGam,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    data.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),3);
    data.(whiskDataType).std_CortS_Gam = std(mean(mean(data.(whiskDataType).cortS_Gam.*100,2),1),0,3);
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
    data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
    data.(whiskDataType).meanHipGam = mean(data.(whiskDataType).Hip.hipGam,2);
    data.(whiskDataType).stdHipGam = std(data.(whiskDataType).Hip.hipGam,0,2);
    data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3).*100;
    data.(whiskDataType).mean_HipS_Gam = mean(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),3);
    data.(whiskDataType).std_HipS_Gam = std(mean(mean(data.(whiskDataType).Hip.hipS_Gam.*100,2),1),0,3);
    data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
    data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end
%% Fig. 1-S3
summaryFigure = figure('Name','Fig1-S3 (a-r)');
sgtitle('Figure 1-S3 - Turner et al. 2020')
%% [1-S3a] brief whisks cortical MUA
ax1 = subplot(6,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA + data.ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA - data.ShortWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title(' Brief whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S3b] moderate whisks cortical MUA
ax2 = subplot(6,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA + data.IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA - data.IntermediateWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title('Moderate whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S3c] extended whisks cortical MUA
ax3 = subplot(6,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA + data.LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA - data.LongWhisks.stdCortMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title('Extended whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S3d] brief whisks cortical LFP
ax4 = subplot(6,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title(' Brief whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S3e] moderate whisks cortical LFP
ax5 = subplot(6,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title(' Moderate whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S3f] extended whisks cortical LFP
ax6 = subplot(6,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('Extended whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S3g] brief whisks hippocampal MUA
ax7 = subplot(6,3,7);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA + data.ShortWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA - data.ShortWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title('Brief whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [1-S3h] moderate whisks hippocampal MUA
ax8 = subplot(6,3,8);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA + data.IntermediateWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA - data.IntermediateWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title('Moderate whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [1-S3i] extended whisks hippocampal MUA
ax9 = subplot(6,3,9);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA + data.LongWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA - data.LongWhisks.stdHipMUA,'color',colors('battleship grey'),'LineWidth',0.5)
title('Extended whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [1-S3j] brief whisks hippocampal LFP
ax10 = subplot(6,3,10);
imagesc(data.ShortWhisks.meanHipT,data.ShortWhisks.meanHipF,data.ShortWhisks.meanHipS)
title('Brief whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [1-S3k] moderate whisks hippocampal LFP
ax11 = subplot(6,3,11);
imagesc(data.IntermediateWhisks.meanHipT,data.IntermediateWhisks.meanHipF,data.IntermediateWhisks.meanHipS)
title(' Moderate whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [1-S3l] extended whisks hippocampal LFP
ax12 = subplot(6,3,12);
imagesc(data.LongWhisks.meanHipT,data.LongWhisks.meanHipF,data.LongWhisks.meanHipS)
title('Extended whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [1-S3m] brief whisks HbT
ax13 = subplot(6,3,13);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT + data.ShortWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT - data.ShortWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Brief whisk CBV')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S3n] moderate whisks HbT
ax14 = subplot(6,3,14);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT + data.IntermediateWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT - data.IntermediateWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Moderate whisk CBV')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S3o] extended whisks HbT
ax15 = subplot(6,3,15);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT + data.LongWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT - data.LongWhisks.stdHbT,'color',colors('battleship grey'),'LineWidth',0.5)
title('Extended whisk CBV')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S3p] brief whisks GCaMP7s
ax16 = subplot(6,3,16);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGCaMP7s,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGCaMP7s + data.ShortWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGCaMP7s - data.ShortWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
title('Brief whisk GCaMP7s')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S3q] moderate whisks GCaMP7s
ax17 = subplot(6,3,17);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGCaMP7s,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGCaMP7s + data.IntermediateWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGCaMP7s - data.IntermediateWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
title('Moderate whisk GCaMP7s')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S3r] extended whisks GCaMP7s
ax18 = subplot(6,3,18);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGCaMP7s,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGCaMP7s + data.LongWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGCaMP7s - data.LongWhisks.stdGCaMP7s,'color',colors('battleship grey'),'LineWidth',0.5)
title('Extended whisk GCaMP7s')
ylabel('\DeltaF/F (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% axes positions
linkaxes([ax1,ax2,ax3,ax7,ax8,ax9],'xy')
linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
linkaxes([ax13,ax14,ax15],'xy')
linkaxes([ax16,ax17,ax18],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
ax10Pos(3:4) = ax1Pos(3:4);
ax11Pos(3:4) = ax2Pos(3:4);
ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S3']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S3'])
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S3_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % text values
%     disp('======================================================================================================================')
%     disp('[1-S3] Text values for gamma/HbT changes')
%     disp('======================================================================================================================')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%      % cortical MUA/LFP
%     [~,index] = max(data.ShortWhisks.meanCortMUA);
%     disp(['Brief whisk Cort gamma MUA P/P (%): ' num2str(round(data.ShortWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdCortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanCortMUA);
%     disp(['Moderate whisk Cort gamma MUA P/P (%): ' num2str(round(data.IntermediateWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdCortMUA(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanCortMUA);
%     disp(['Extended whisk Cort gamma MUA P/P (%): ' num2str(round(data.LongWhisks.meanCortMUA(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdCortMUA(index),1))]); disp(' ')
%     % cortical LFP
%     disp(['Brief whisk Cort gamma LFP P/P (%): ' num2str(round(data.ShortWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.ShortWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Cort gamma LFP P/P (%): ' num2str(round(data.IntermediateWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.IntermediateWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Cort gamma LFP P/P (%): ' num2str(round(data.LongWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.LongWhisks.std_CortS_Gam,1))]); disp(' ')
%     % hippocampal MUA
%     [~,index] = max(data.ShortWhisks.meanHipMUA);
%     disp(['Brief whisk Hip gamma MUA P/P (%): ' num2str(round(data.ShortWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdHipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanHipMUA);
%     disp(['Moderate whisk Hip gamma MUA P/P (%): ' num2str(round(data.IntermediateWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdHipMUA(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanHipMUA);
%     disp(['Extended whisk Hip gamma MUA P/P (%): ' num2str(round(data.LongWhisks.meanHipMUA(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdHipMUA(index),1))]); disp(' ')
%     % hippocampal LFP
%     disp(['Brief whisk Hip gamma LFP P/P (%): ' num2str(round(data.ShortWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.ShortWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Hip gamma LFP P/P (%): ' num2str(round(data.IntermediateWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.IntermediateWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Hip gamma LFP P/P (%): ' num2str(round(data.LongWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.LongWhisks.std_HipS_Gam,1))]); disp(' ')
%     % HbT
%     [~,index] = max(data.ShortWhisks.meanHbT);
%     disp(['Brief whisk [HbT]: ' num2str(round(data.ShortWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdHbT(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanHbT);
%     disp(['Moderate whisk [HbT]: ' num2str(round(data.IntermediateWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdHbT(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanHbT);
%     disp(['Extended whisk [HbT]: ' num2str(round(data.LongWhisks.meanHbT(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdHbT(index),1))]); disp(' ')
%     % R/R
%     [~,index] = min(data.ShortWhisks.meanGCaMP7s);
%     disp(['Brief whisk GCaMP7s: ' num2str(round(data.ShortWhisks.meanGCaMP7s(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdGCaMP7s(index),1))]); disp(' ')
%     [~,index] = min(data.IntermediateWhisks.meanGCaMP7s);
%     disp(['Moderate whisk GCaMP7s: ' num2str(round(data.IntermediateWhisks.meanGCaMP7s(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdGCaMP7s(index),1))]); disp(' ')
%     [~,index] = min(data.LongWhisks.meanGCaMP7s);
%     disp(['Extended whisk GCaMP7s: ' num2str(round(data.LongWhisks.meanGCaMP7s(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdGCaMP7s(index),1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
end

end
