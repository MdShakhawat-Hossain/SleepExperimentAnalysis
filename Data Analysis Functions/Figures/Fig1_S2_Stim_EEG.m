function [AnalysisResults] = Fig1_S2_Stim_EEG(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
FPanimalIDs = {'SHF025'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
% dataTypes = {'LH','RH'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.RH.(solenoidName).count(:,aa) = AnalysisResults.(animalID).Stim.RH.(solenoidName).count;
            data.RH.(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.RH.(solenoidName).Rhodamine.Rhodamine;
            data.RH.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.RH.(solenoidName).GFP.GFP;
            data.EEG.(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).MUA.EEGData;
            data.EEG.(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).Gam.EEGData;
            data.EEG.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).timeVector;
            data.EEG.(solenoidName).cortS(:,:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).LFP.EEGS;
            data.EEG.(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).LFP.EEGS(49:end,20:23);
            data.EEG.(solenoidName).T(:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).LFP.T;
            data.EEG.(solenoidName).F(:,aa) = AnalysisResults.(animalID).Stim.EEG.(solenoidName).LFP.F;
        end
end
% concatenate the data from the contra and ipsi data
data.Contra.count = data.RH.LPadSol.count;
data.Contra.RHRhodamine = data.RH.LPadSol.Rhodamine;
data.Contra.RHGFP = data.RH.LPadSol.GFP;
data.Contra.cortMUA = data.EEG.LPadSol.cortMUA;
data.Contra.cortGam = data.EEG.LPadSol.cortGam;
data.Contra.timeVector = data.EEG.LPadSol.timeVector;
data.Contra.cortS = data.EEG.LPadSol.cortS;
data.Contra.cortS_Gam = data.EEG.LPadSol.cortS_Gam;
data.Contra.T = data.EEG.LPadSol.T;
data.Contra.F = data.EEG.LPadSol.F;


data.Ipsi.count = data.RH.RPadSol.count;
data.Ipsi.RHRhodamine = data.RH.RPadSol.Rhodamine;
data.Ipsi.RHGFP = data.RH.RPadSol.GFP; 
data.Ipsi.cortMUA = data.EEG.RPadSol.cortMUA; 
data.Ipsi.cortGam = data.EEG.RPadSol.cortGam; 
data.Ipsi.timeVector = data.EEG.RPadSol.timeVector;
data.Ipsi.cortS = data.EEG.RPadSol.cortS; 
data.Ipsi.cortS_Gam = data.EEG.RPadSol.cortS_Gam;
data.Ipsi.T = data.EEG.RPadSol.T;
data.Ipsi.F = data.EEG.RPadSol.F;


data.Auditory.count = data.RH.AudSol.count;
data.Auditory.RHGFP = data.RH.AudSol.GFP;
data.Auditory.RHRhodamine = data.RH.AudSol.Rhodamine;
data.Auditory.cortMUA = data.EEG.AudSol.cortMUA;
data.Auditory.cortGam = data.EEG.AudSol.cortGam;
data.Auditory.timeVector = data.EEG.AudSol.timeVector;
data.Auditory.cortS = data.EEG.AudSol.cortS;
data.Auditory.cortS_Gam = data.EEG.AudSol.cortS_Gam;
data.Auditory.T = data.EEG.AudSol.T;
data.Auditory.F = data.EEG.AudSol.F;


% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    data.(compDataType).RHmean_Rhodamine = mean(data.(compDataType).RHRhodamine,2);
    data.(compDataType).RHstd_Rhodamine = std(data.(compDataType).RHRhodamine,0,2);
    data.(compDataType).RHmean_GFP = mean(data.(compDataType).RHGFP,2);
    data.(compDataType).RHstd_GFP = std(data.(compDataType).RHGFP,0,2);
    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),1),3);
    data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),1),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end
%% Fig. 1-S2
summaryFigure = figure('Name','Fig1-S2 Stim');
sgtitle('Stimulus evoked responses')
%% [1-S2a] EEG MUA contra stim
ax1 = subplot(4,3,1);
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA,'color',colors('rich black'),'LineWidth',1);
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA + data.Contra.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA - data.Contra.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title('Contra stim EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [1-S2b] EEG MUA ispi stim
ax2 = subplot(4,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA + data.Ipsi.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA - data.Ipsi.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title('Ipsi stim EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S2c] EEG MUA auditory stim
ax3 = subplot(4,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA,'color',colors('rich black'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA + data.Auditory.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA - data.Auditory.std_CortMUA,'color',colors('battleship grey'),'LineWidth',0.25)
title('Aud stim EEG MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S2d] EEG LFP contra stim
ax4 = subplot(4,3,4);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_CortS)
title('Contra stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%%caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S2e] EEG LFP ispi stim
ax5 = subplot(4,3,5);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_CortS)
title('Ipsi stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S2f] EEG LFP auditory stim
ax6 = subplot(4,3,6);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_CortS)
title('Aud stim EEG LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-50,75])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];

%% [1-S2m] Rhodamine contra stim
ax13 = subplot(4,3,7);
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine + data.Contra.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine - data.Contra.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F (RH)')


xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [1-S2n] Rhodamine ispi stim
ax14 = subplot(4,3,8);
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine + data.Ipsi.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine - data.Ipsi.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F (RH)')

xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [1-S2o] Rhodamine auditory stim
ax15 = subplot(4,3,9);
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine + data.Auditory.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine - data.Auditory.RHstd_Rhodamine,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F (RH)')

xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [1-S2p] GFP contra stim
ax16 = subplot(4,3,10);

plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP + data.Contra.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP - data.Contra.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Contra stim GFP')
ylabel('Z \DeltaF/F GFP RH')

xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [1-S2q] GFP ispi stim
ax17 = subplot(4,3,11);

plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP + data.Ipsi.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP - data.Ipsi.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Ipsi stim GFP')
ylabel('Z \DeltaF/F GFP RH')

xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [1-S2r] GFP auditory stim
ax18 = subplot(4,3,12);

plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP + data.Auditory.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP - data.Auditory.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
title('Auditory stim GFP')
ylabel('Z \DeltaF/F GFP RH')

xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% adjust and link axes
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
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S2-Stim']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2-Stim'])
    close 
end
%% plot the comparison of GCaMP and Rhodamine with GFPRH and Rhodamine
summaryFigureN = figure('Name','Fig1-S2 Stim');
sgtitle('Stimulus evoked responses in fiber photometry signals')


% RH contra stim
ax4 = subplot(1,3,1);
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine + data.Contra.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.RHmean_Rhodamine - data.Contra.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Contra stim RH fiber')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP + data.Contra.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Contra.mean_timeVector,data.Contra.RHmean_GFP - data.Contra.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GFPRH')
ax4.YAxis(1).Color = colors('indian red');
ax4.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
% RH ispi stim
ax5 = subplot(1,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine + data.Ipsi.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_Rhodamine - data.Ipsi.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Ipsi stim RH fibers')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP + data.Ipsi.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Ipsi.mean_timeVector,data.Ipsi.RHmean_GFP - data.Ipsi.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GFPRH')
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% RH auditory stim
ax6 = subplot(1,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine,'color',colors('indian red'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine + data.Auditory.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_Rhodamine - data.Auditory.RHstd_Rhodamine,'color',colors('indian red'),'LineWidth',0.25)
title('Auditory stim RH fibers')
ylabel('Z \DeltaF/F Rhodamine')

yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP,'color',colors('army green'),'LineWidth',1)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP + data.Auditory.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
plot(data.Auditory.mean_timeVector,data.Auditory.RHmean_GFP - data.Auditory.RHstd_GFP,'color',colors('army green'),'LineWidth',0.25)
ylabel('Z \DeltaF/F GFPRH')
ax6.YAxis(1).Color = colors('indian red');
ax6.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath 'Fig1-S2-Stim-FiberSignals']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig1-S2-Stim-FiberSignals'])
    close 
end

end
