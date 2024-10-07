function [AnalysisResults] = Fig1_S3_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomed Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S3 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
animalIDs = {'SHF025','SHF026','SHF027'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for cc = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,cc};
        % left EEG
        data.(whiskDataType).LH.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).LH.GFP(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).GFP.GFP;

        data.(whiskDataType).LH.EEGEEG(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).EEG.EEGData;
        data.(whiskDataType).LH.EEGGam(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).Gam.EEGData;
        data.(whiskDataType).LH.EEGS(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.EEGS;
        data.(whiskDataType).LH.EEGS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.EEGS(49:end,20:23);
        data.(whiskDataType).LH.EEGT(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.T;
        data.(whiskDataType).LH.EEGF(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).LFP.F;
        % right EEG
        data.(whiskDataType).RH.Rhodamine(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).Rhodamine.Rhodamine;
        data.(whiskDataType).RH.GFP(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).GFP.GFP;

        data.(whiskDataType).RH.EEGEEG(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).EEG.EEGData;
        data.(whiskDataType).RH.EEGGam(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).Gam.EEGData;
        data.(whiskDataType).RH.EEGS(:,:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.EEGS;
        data.(whiskDataType).RH.EEGS_Gam(:,:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.EEGS(49:end,20:23);
        data.(whiskDataType).RH.EEGT(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.T;
        data.(whiskDataType).RH.EEGF(:,aa) = AnalysisResults.(animalID).Whisk.RH.(whiskDataType).LFP.F;
  
        % time vector
        data.(whiskDataType).timeVector(:,aa) = AnalysisResults.(animalID).Whisk.LH.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
%     data.(whiskDataType).Rhodamine = cat(2,data.(whiskDataType).LH.Rhodamine,data.(whiskDataType).RH.Rhodamine);
%     data.(whiskDataType).GFP = cat(2,data.(whiskDataType).LH.GFP,data.(whiskDataType).RH.GFP);
% 
%     data.(whiskDataType).EEGEEG = cat(2,data.(whiskDataType).LH.EEGEEG,data.(whiskDataType).RH.EEGEEG);
%     data.(whiskDataType).EEGGam = cat(2,data.(whiskDataType).LH.EEGGam,data.(whiskDataType).RH.EEGGam);
%     data.(whiskDataType).EEGS = cat(3,data.(whiskDataType).LH.EEGS,data.(whiskDataType).RH.EEGS);
%     data.(whiskDataType).EEGS_Gam = cat(3,data.(whiskDataType).LH.EEGS_Gam,data.(whiskDataType).RH.EEGS_Gam);
%     data.(whiskDataType).EEGT = cat(2,data.(whiskDataType).LH.EEGT,data.(whiskDataType).RH.EEGT);
%     data.(whiskDataType).EEGF = cat(2,data.(whiskDataType).LH.EEGF,data.(whiskDataType).RH.EEGF);
    data.(whiskDataType).Rhodamine = data.(whiskDataType).RH.Rhodamine;
    data.(whiskDataType).GFP = data.(whiskDataType).RH.GFP;

    data.(whiskDataType).EEGEEG = data.(whiskDataType).LH.EEGEEG;
    data.(whiskDataType).EEGGam = data.(whiskDataType).LH.EEGGam;
    data.(whiskDataType).EEGS = data.(whiskDataType).LH.EEGS;
    data.(whiskDataType).EEGS_Gam = data.(whiskDataType).LH.EEGS_Gam;
    data.(whiskDataType).EEGT = data.(whiskDataType).LH.EEGT;
    data.(whiskDataType).EEGF = data.(whiskDataType).LH.EEGF;
end
% concatenate the data from the contra and ipsi data
for ee = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,ee};
    data.(whiskDataType).meanRhodamine = mean(data.(whiskDataType).Rhodamine,2);
    data.(whiskDataType).stdRhodamine = std(data.(whiskDataType).Rhodamine,0,2);
    data.(whiskDataType).meanGFP = mean(data.(whiskDataType).GFP,2);
    data.(whiskDataType).stdGFP = std(data.(whiskDataType).GFP,0,2);

    data.(whiskDataType).meanCortEEG = mean(data.(whiskDataType).EEGEEG,2);
    data.(whiskDataType).stdCortEEG = std(data.(whiskDataType).EEGEEG,0,2);
    data.(whiskDataType).meanCortGam = mean(data.(whiskDataType).EEGGam,2);
    data.(whiskDataType).stdCortGam = std(data.(whiskDataType).EEGGam,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).EEGS,3).*100;
    data.(whiskDataType).mean_CortS_Gam = mean(mean(mean(data.(whiskDataType).EEGS_Gam.*100,2),1),3);
    data.(whiskDataType).std_CortS_Gam = std(mean(mean(data.(whiskDataType).EEGS_Gam.*100,2),1),0,3);
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).EEGT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).EEGF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,aa),2);
end
%% Fig. 1-S3
summaryFigure = figure('Name','Fig1-S3 (a-r)');
sgtitle('Whisk evoked responses')
%% [1-S3a] brief whisks EEG EEG
ax1 = subplot(4,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortEEG,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortEEG + data.ShortWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortEEG - data.ShortWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title(' Brief whisk EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S3b] moderate whisks EEG EEG
ax2 = subplot(4,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortEEG,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortEEG + data.IntermediateWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortEEG - data.IntermediateWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Moderate whisk EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S3c] extended whisks EEG EEG
ax3 = subplot(4,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortEEG,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortEEG + data.LongWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortEEG - data.LongWhisks.stdCortEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Extended whisk EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S3d] brief whisks EEG LFP
ax4 = subplot(4,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title(' Brief whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S3e] moderate whisks EEG LFP
ax5 = subplot(4,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title(' Moderate whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S3f] extended whisks EEG LFP
ax6 = subplot(4,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('Extended whisk EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S3m] brief whisks Rhodamine
ax13 = subplot(4,3,7);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanRhodamine + data.ShortWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanRhodamine - data.ShortWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Brief whisk Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S3n] moderate whisks Rhodamine
ax14 = subplot(4,3,8);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanRhodamine + data.IntermediateWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanRhodamine - data.IntermediateWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Moderate whisk Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S3o] extended whisks Rhodamine
ax15 = subplot(4,3,9;
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanRhodamine,'color',colors('indian red'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanRhodamine + data.LongWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanRhodamine - data.LongWhisks.stdRhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Extended whisk Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S3p] brief whisks GFP
ax16 = subplot(4,3,10);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGFP,'color',colors('army green'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGFP + data.ShortWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanGFP - data.ShortWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
title('Brief whisk GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S3q] moderate whisks GFP
ax17 = subplot(4,3,11);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGFP,'color',colors('army green'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGFP + data.IntermediateWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanGFP - data.IntermediateWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
title('Moderate whisk GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S3r] extended whisks GFP
ax18 = subplot(4,3,12);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGFP,'color',colors('army green'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGFP + data.LongWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanGFP - data.LongWhisks.stdGFP,'color',colors('army green'),'LineWidth',0.25)
title('Extended whisk GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-whisk time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% axes positions
linkaxes([ax1,ax2,ax3],'xy')
linkaxes([ax4,ax5,ax6],'xy')
linkaxes([ax13,ax14,ax15],'xy')
linkaxes([ax16,ax17,ax18],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
% ax10Pos = get(ax10,'position');
% ax11Pos = get(ax11,'position');
% ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
% ax10Pos(3:4) = ax1Pos(3:4);
% ax11Pos(3:4) = ax2Pos(3:4);
% ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
% set(ax10,'position',ax10Pos);
% set(ax11,'position',ax11Pos);
% set(ax12,'position',ax12Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S3']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S3'])
    close 
    %% Text diary
%     diaryFile = [dirpath 'Fig1-S3_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % text values
%     disp('======================================================================================================================')
%     disp('[1-S3] Text values for gamma/Rhodamine changes')
%     disp('======================================================================================================================')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%      % EEG EEG/LFP
%     [~,index] = max(data.ShortWhisks.meanCortEEG);
%     disp(['Brief whisk Cort gamma EEG P/P (%): ' num2str(round(data.ShortWhisks.meanCortEEG(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdCortEEG(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanCortEEG);
%     disp(['Moderate whisk Cort gamma EEG P/P (%): ' num2str(round(data.IntermediateWhisks.meanCortEEG(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdCortEEG(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanCortEEG);
%     disp(['Extended whisk Cort gamma EEG P/P (%): ' num2str(round(data.LongWhisks.meanCortEEG(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdCortEEG(index),1))]); disp(' ')
%     % EEG LFP
%     disp(['Brief whisk Cort gamma LFP P/P (%): ' num2str(round(data.ShortWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.ShortWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Cort gamma LFP P/P (%): ' num2str(round(data.IntermediateWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.IntermediateWhisks.std_CortS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Cort gamma LFP P/P (%): ' num2str(round(data.LongWhisks.mean_CortS_Gam,1)) ' +/- ' num2str(round(data.LongWhisks.std_CortS_Gam,1))]); disp(' ')
%     % hippocampal EEG
%     [~,index] = max(data.ShortWhisks.meanHipEEG);
%     disp(['Brief whisk Hip gamma EEG P/P (%): ' num2str(round(data.ShortWhisks.meanHipEEG(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdHipEEG(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanHipEEG);
%     disp(['Moderate whisk Hip gamma EEG P/P (%): ' num2str(round(data.IntermediateWhisks.meanHipEEG(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdHipEEG(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanHipEEG);
%     disp(['Extended whisk Hip gamma EEG P/P (%): ' num2str(round(data.LongWhisks.meanHipEEG(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdHipEEG(index),1))]); disp(' ')
%     % hippocampal LFP
%     disp(['Brief whisk Hip gamma LFP P/P (%): ' num2str(round(data.ShortWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.ShortWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Moderate whisk Hip gamma LFP P/P (%): ' num2str(round(data.IntermediateWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.IntermediateWhisks.std_HipS_Gam,1))]); disp(' ')
%     disp(['Extended whisk Hip gamma LFP P/P (%): ' num2str(round(data.LongWhisks.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.LongWhisks.std_HipS_Gam,1))]); disp(' ')
%     % Rhodamine
%     [~,index] = max(data.ShortWhisks.meanRhodamine);
%     disp(['Brief whisk [Rhodamine]: ' num2str(round(data.ShortWhisks.meanRhodamine(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdRhodamine(index),1))]); disp(' ')
%     [~,index] = max(data.IntermediateWhisks.meanRhodamine);
%     disp(['Moderate whisk [Rhodamine]: ' num2str(round(data.IntermediateWhisks.meanRhodamine(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdRhodamine(index),1))]); disp(' ')
%     [~,index] = max(data.LongWhisks.meanRhodamine);
%     disp(['Extended whisk [Rhodamine]: ' num2str(round(data.LongWhisks.meanRhodamine(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdRhodamine(index),1))]); disp(' ')
%     % R/R
%     [~,index] = min(data.ShortWhisks.meanGFP);
%     disp(['Brief whisk GFP: ' num2str(round(data.ShortWhisks.meanGFP(index),1)) ' +/- ' num2str(round(data.ShortWhisks.stdGFP(index),1))]); disp(' ')
%     [~,index] = min(data.IntermediateWhisks.meanGFP);
%     disp(['Moderate whisk GFP: ' num2str(round(data.IntermediateWhisks.meanGFP(index),1)) ' +/- ' num2str(round(data.IntermediateWhisks.stdGFP(index),1))]); disp(' ')
%     [~,index] = min(data.LongWhisks.meanGFP);
%     disp(['Extended whisk GFP: ' num2str(round(data.LongWhisks.meanGFP(index),1)) ' +/- ' num2str(round(data.LongWhisks.stdGFP(index),1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
end

end
