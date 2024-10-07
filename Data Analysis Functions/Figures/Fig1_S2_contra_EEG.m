function [AnalysisResults] = Fig1_S2_contra_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomed Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
FPanimalIDs = {'SHF025','SHF026','SHF027'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
dataTypes = {'LH','RH'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.(dataType).(solenoidName).count(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).count;
            data.(dataType).(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).Rhodamine.Rhodamine;
            data.(dataType).(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).GFP.GFP;
            data.(dataType).(solenoidName).EEGEEG(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).EEG.EEGData;
            data.(dataType).(solenoidName).EEGGam(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).Gam.EEGData;
            data.(dataType).(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).timeVector;
            data.(dataType).(solenoidName).EEGS(:,:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).LFP.EEGS;
            data.(dataType).(solenoidName).EEGS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).LFP.EEGS(49:end,20:23);
            data.(dataType).(solenoidName).T(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).LFP.T;
            data.(dataType).(solenoidName).F(:,aa) = AnalysisResults.(animalID).Stim.(dataType).(solenoidName).LFP.F;
        end
    end
end
% concatenate the data from the contra and ipsi data
% data.Contra.count = cat(2,data.LH.RPadSol.count,data.RH.LPadSol.count);
% data.Contra.Rhodamine = cat(2,data.LH.RPadSol.Rhodamine,data.RH.LPadSol.Rhodamine);
% data.Contra.GFP = cat(2,data.LH.RPadSol.GFP,data.RH.LPadSol.GFP);
% data.Contra.EEGEEG = cat(2,data.LH.RPadSol.EEGEEG,data.RH.LPadSol.EEGEEG);
% data.Contra.EEGGam = cat(2,data.LH.RPadSol.EEGGam,data.RH.LPadSol.EEGGam);
% data.Contra.timeVector = cat(2,data.LH.RPadSol.timeVector,data.RH.LPadSol.timeVector);
% data.Contra.EEGS = cat(3,data.LH.RPadSol.EEGS,data.RH.LPadSol.EEGS);
% data.Contra.EEGS_Gam = cat(3,data.LH.RPadSol.EEGS_Gam,data.RH.LPadSol.EEGS_Gam);
% data.Contra.T = cat(2,data.LH.RPadSol.T,data.RH.LPadSol.T);
% data.Contra.F = cat(2,data.LH.RPadSol.F,data.RH.LPadSol.F);
% data.Ipsi.count = cat(2,data.LH.LPadSol.count,data.RH.RPadSol.count);
% data.Ipsi.Rhodamine = cat(2,data.LH.LPadSol.Rhodamine,data.RH.RPadSol.Rhodamine);
% data.Ipsi.GFP = cat(2,data.LH.LPadSol.GFP,data.RH.RPadSol.GFP);
% data.Ipsi.EEGEEG = cat(2,data.LH.LPadSol.EEGEEG,data.RH.RPadSol.EEGEEG);
% data.Ipsi.EEGGam = cat(2,data.LH.LPadSol.EEGGam,data.RH.RPadSol.EEGGam);
% data.Ipsi.timeVector = cat(2,data.LH.LPadSol.timeVector,data.RH.RPadSol.timeVector);
% data.Ipsi.EEGS = cat(3,data.LH.LPadSol.EEGS,data.RH.RPadSol.EEGS);
% data.Ipsi.EEGS_Gam = cat(3,data.LH.LPadSol.EEGS_Gam,data.RH.RPadSol.EEGS_Gam);
% data.Ipsi.T = cat(2,data.LH.LPadSol.T,data.RH.RPadSol.T);
% data.Ipsi.F = cat(2,data.LH.LPadSol.F,data.RH.RPadSol.F);
% data.Auditory.count = cat(2,data.LH.AudSol.count,data.RH.AudSol.count);
% data.Auditory.Rhodamine = cat(2,data.LH.AudSol.Rhodamine,data.RH.AudSol.Rhodamine);
% data.Auditory.GFP = cat(2,data.LH.AudSol.GFP,data.RH.AudSol.GFP);
% data.Auditory.EEGEEG = cat(2,data.LH.AudSol.EEGEEG,data.RH.AudSol.EEGEEG);
% data.Auditory.EEGGam = cat(2,data.LH.AudSol.EEGGam,data.RH.AudSol.EEGGam);
% data.Auditory.timeVector = cat(2,data.LH.AudSol.timeVector,data.RH.AudSol.timeVector);
% data.Auditory.EEGS = cat(3,data.LH.AudSol.EEGS,data.RH.AudSol.EEGS);
% data.Auditory.EEGS_Gam = cat(3,data.LH.AudSol.EEGS_Gam,data.RH.AudSol.EEGS_Gam);
% data.Auditory.T = cat(2,data.LH.AudSol.T,data.RH.AudSol.T);
% data.Auditory.F = cat(2,data.LH.AudSol.F,data.RH.AudSol.F);

data.Contra.count = data.RH.LPadSol.count;
data.Contra.Rhodamine = data.RH.LPadSol.Rhodamine;
data.Contra.GFP = data.RH.LPadSol.GFP;
data.Contra.EEGEEG = data.LH.RPadSol.EEGEEG;
data.Contra.EEGGam = data.LH.RPadSol.EEGGam;
data.Contra.timeVector = data.LH.RPadSol.timeVector;
data.Contra.EEGS = data.LH.RPadSol.EEGS;
data.Contra.EEGS_Gam = data.LH.RPadSol.EEGS_Gam;
data.Contra.T = data.LH.RPadSol.T;
data.Contra.F = data.LH.RPadSol.F;
data.Ipsi.count = data.RH.RPadSol.count;
data.Ipsi.Rhodamine = data.RH.RPadSol.Rhodamine;
data.Ipsi.GFP = data.RH.RPadSol.GFP;
data.Ipsi.EEGEEG = data.LH.LPadSol.EEGEEG;
data.Ipsi.EEGGam = data.LH.LPadSol.EEGGam;
data.Ipsi.timeVector = data.LH.LPadSol.timeVector;
data.Ipsi.EEGS = data.LH.LPadSol.EEGS;
data.Ipsi.EEGS_Gam = data.LH.LPadSol.EEGS_Gam;
data.Ipsi.T = data.LH.LPadSol.T;
data.Ipsi.F = data.LH.LPadSol.F;
data.Auditory.count = data.RH.AudSol.count;
data.Auditory.Rhodamine = data.RH.AudSol.Rhodamine;
data.Auditory.GFP = data.RH.AudSol.GFP;
data.Auditory.EEGEEG = data.LH.AudSol.EEGEEG;
data.Auditory.EEGGam = data.LH.AudSol.EEGGam;
data.Auditory.timeVector = data.LH.AudSol.timeVector;
data.Auditory.EEGS = cdata.LH.AudSol.EEGS;
data.Auditory.EEGS_Gam = data.LH.AudSol.EEGS_Gam;
data.Auditory.T = data.LH.AudSol.T;
data.Auditory.F = data.LH.AudSol.F;
% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    data.(compDataType).mean_Rhodamine = mean(data.(compDataType).Rhodamine,2);
    data.(compDataType).std_Rhodamine = std(data.(compDataType).Rhodamine,0,2);
    data.(compDataType).mean_GFP = mean(data.(compDataType).GFP,2);
    data.(compDataType).std_GFP = std(data.(compDataType).GFP,0,2);
    data.(compDataType).mean_EEGEEG = mean(data.(compDataType).EEGEEG,2);
    data.(compDataType).std_EEGEEG = std(data.(compDataType).EEGEEG,0,2);
    data.(compDataType).mean_EEGGam = mean(data.(compDataType).EEGGam,2);
    data.(compDataType).std_EEGGam = std(data.(compDataType).EEGGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_EEGS = mean(data.(compDataType).EEGS,3).*100;
    data.(compDataType).mean_EEGS_Gam = mean(mean(mean(data.(compDataType).EEGS_Gam.*100,2),1),3);
    data.(compDataType).std_EEGS_Gam = std(mean(mean(data.(compDataType).EEGS_Gam.*100,2),1),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end
%% Fig. 1-S2
summaryFigure = figure('Name','Fig1-S2 (a-r)');
sgtitle('Stimulus evoked responses')
%% [1-S2a] EEG EEG contra stim
ax1 = subplot(6,3,1);
plot(data.Contra.mean_timeVector,data.Contra.mean_EEGEEG,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_EEGEEG + data.Contra.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.mean_EEGEEG - data.Contra.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Contra stim EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S2b] EEG EEG ispi stim
ax2 = subplot(6,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_EEGEEG,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_EEGEEG + data.Ipsi.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_EEGEEG - data.Ipsi.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Ipsi stim EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S2c] EEG EEG auditory stim
ax3 = subplot(6,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_EEGEEG,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_EEGEEG + data.Auditory.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_EEGEEG - data.Auditory.std_EEGEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Aud stim EEG EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S2d] EEG LFP contra stim
ax4 = subplot(6,3,4);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_EEGS)
title('Contra stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S2e] EEG LFP ispi stim
ax5 = subplot(6,3,5);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_EEGS)
title('Ipsi stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S2f] EEG LFP auditory stim
ax6 = subplot(6,3,6);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_EEGS)
title('Aud stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S2g] hippocampal EEG contra stim
ax7 = subplot(6,3,7);
plot(data.Contra.mean_timeVector,data.Contra.mean_HipEEG,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_HipEEG + data.Contra.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.mean_HipEEG - data.Contra.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Contra stim hippocampal EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [1-S2h] hippocampal EEG ispi stim
ax8 = subplot(6,3,8);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipEEG,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipEEG + data.Ipsi.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipEEG - data.Ipsi.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Ipsi stim hippocampal EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [1-S2i] hippocampal EEG auditory stim
ax9 = subplot(6,3,9);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipEEG,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipEEG + data.Auditory.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipEEG - data.Auditory.std_HipEEG,'color',colors('battleship grey'),'LineWidth',0.25)
title('Aud stim hippocampal EEG')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP contra stim
ax10 = subplot(6,3,10);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_HipS)
title('Contra stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP ispi stim
ax11 = subplot(6,3,11);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_HipS)
title('Ipsi stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [1-S2l] hippocampal LFP auditory stim
ax12 = subplot(6,3,12);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_HipS)
title('Aud stim hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VertAlignment','bottom')
caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [1-S2m] Rhodamine contra stim
ax13 = subplot(6,3,13);
plot(data.Contra.mean_timeVector,data.Contra.mean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_Rhodamine + data.Contra.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.mean_Rhodamine - data.Contra.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Contra stim Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S2n] Rhodamine ispi stim
ax14 = subplot(6,3,14);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_Rhodamine + data.Ipsi.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_Rhodamine - data.Ipsi.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Ipsi stim Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S2o] Rhodamine auditory stim
ax15 = subplot(6,3,15);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_Rhodamine + data.Auditory.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_Rhodamine - data.Auditory.std_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Aud stim Rhodamine')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S2p] GFP contra stim
ax16 = subplot(6,3,16);
plot(data.Contra.mean_timeVector,data.Contra.mean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_GFP + data.Contra.std_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.mean_GFP - data.Contra.std_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Contra stim GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S2q] GFP ispi stim
ax17 = subplot(6,3,17);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_GFP + data.Ipsi.std_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_GFP - data.Ipsi.std_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Ipsi stim GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S2r] GFP auditory stim
ax18 = subplot(6,3,18);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_GFP + data.Auditory.std_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_GFP - data.Auditory.std_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Aud stim GFP')
ylabel('Zscored \DeltaF/F')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% adjust and link axes
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
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S2']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2'])
    close 
    %% text diary
%     diaryFile = [dirpath 'Fig1-S2_Statistics.txt'];
%     if exist(diaryFile,'file') == 2
%         delete(diaryFile)
%     end
%     diary(diaryFile)
%     diary on
%     % text values
%     disp('======================================================================================================================')
%     disp('[1-S2] Text values for gamma/Rhodamine/reflectance changes')
%     disp('======================================================================================================================')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     % EEG EEG/LFP
%     [~,index] = max(data.Contra.mean_EEGEEG);
%     disp(['Contra stim EEG gamma EEG P/P (%): ' num2str(round(data.Contra.mean_EEGEEG(index),1)) ' +/- ' num2str(round(data.Contra.std_EEGEEG(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_EEGEEG);
%     disp(['Ipsil stim EEG gamma EEG P/P (%): ' num2str(round(data.Ipsi.mean_EEGEEG(index),1)) ' +/- ' num2str(round(data.Ipsi.std_EEGEEG(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_EEGEEG);
%     disp(['Audit stim EEG gamma EEG P/P (%): ' num2str(round(data.Auditory.mean_EEGEEG(index),1)) ' +/- ' num2str(round(data.Auditory.std_EEGEEG(index),1))]); disp(' ')
%     % EEG LFP
%     disp(['Contra stim EEG gamma LFP P/P (%): ' num2str(round(data.Contra.mean_EEGS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_EEGS_Gam,1))]); disp(' ')
%     disp(['Ipsil stim EEG gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_EEGS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_EEGS_Gam,1))]); disp(' ')
%     disp(['Audit stim EEG gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_EEGS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_EEGS_Gam,1))]); disp(' ')
%     % hippocampal EEG
%     [~,index] = max(data.Contra.mean_HipEEG);
%     disp(['Contra stim Hip gamma EEG P/P (%): ' num2str(round(data.Contra.mean_HipEEG(index),1)) ' +/- ' num2str(round(data.Contra.std_HipEEG(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_HipEEG);
%     disp(['Ipsil stim Hip gamma EEG P/P (%): ' num2str(round(data.Ipsi.mean_HipEEG(index),1)) ' +/- ' num2str(round(data.Ipsi.std_HipEEG(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_HipEEG);
%     disp(['Audit stim Hip gamma EEG P/P (%): ' num2str(round(data.Auditory.mean_HipEEG(index),1)) ' +/- ' num2str(round(data.Auditory.std_HipEEG(index),1))]); disp(' ')
%     % hippocampal LFP
%     disp(['Contra stim Hip gamma LFP P/P (%): ' num2str(round(data.Contra.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Contra.std_HipS_Gam,1))]); disp(' ')
%     disp(['Ipsil stim Hip gamma LFP P/P (%): ' num2str(round(data.Ipsi.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Ipsi.std_HipS_Gam,1))]); disp(' ')
%     disp(['Auditory stim Hip gamma LFP P/P (%): ' num2str(round(data.Auditory.mean_HipS_Gam,1)) ' +/- ' num2str(round(data.Auditory.std_HipS_Gam,1))]); disp(' ')
%     % Rhodamine
%     [~,index] = max(data.Contra.mean_Rhodamine);
%     disp(['Contra stim [Rhodamine]: ' num2str(round(data.Contra.mean_Rhodamine(index),1)) ' +/- ' num2str(round(data.Contra.std_Rhodamine(index),1))]); disp(' ')
%     [~,index] = max(data.Ipsi.mean_Rhodamine);
%     disp(['Ipsil stim [Rhodamine]: ' num2str(round(data.Ipsi.mean_Rhodamine(index),1)) ' +/- ' num2str(round(data.Ipsi.std_Rhodamine(index),1))]); disp(' ')
%     [~,index] = max(data.Auditory.mean_Rhodamine);
%     disp(['Audit stim [Rhodamine]: ' num2str(round(data.Auditory.mean_Rhodamine(index),1)) ' +/- ' num2str(round(data.Auditory.std_Rhodamine(index),1))]); disp(' ')
%     % GFP
%     [~,index] = min(data.Contra.mean_GFP);
%     disp(['Contra stim GFP: ' num2str(round(data.Contra.mean_GFP(index),1)) ' +/- ' num2str(round(data.Contra.std_GFP(index),1))]); disp(' ')
%     [~,index] = min(data.Ipsi.mean_GFP);
%     disp(['Ipsil stim GFP: ' num2str(round(data.Ipsi.mean_GFP(index),1)) ' +/- ' num2str(round(data.Ipsi.std_GFP(index),1))]); disp(' ')
%     [~,index] = min(data.Auditory.mean_GFP);
%     disp(['Audit stim GFP: ' num2str(round(data.Auditory.mean_GFP(index),1)) ' +/- ' num2str(round(data.Auditory.std_GFP(index),1))]); disp(' ')
%     disp('----------------------------------------------------------------------------------------------------------------------')
%     diary off
end

end
