function [AnalysisResults] = Fig1_S2_Stim_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,FPanimalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S2 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________
% AnalysisResults =  AnalysisResults_firstHrs; 
% saveFigs = 'y';
%% set-up and process data
% FPanimalIDs = {'NEACh005'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};
% dataTypes = {'LH','RH'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
        for dd = 1:length(solenoidNames)
            solenoidName = solenoidNames{1,dd};
            data.Z_NE.(solenoidName).count(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).count;
            data.Z_NE.(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).Rhodamine.Rhodamine;
            data.Z_Ach.(solenoidName).Rhodamine(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).Rhodamine.Rhodamine;
            data.Z_NE.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.Z_NE.(solenoidName).GFP.GFP;
            data.Z_Ach.(solenoidName).GFP(:,aa) = AnalysisResults.(animalID).Stim.Z_Ach.(solenoidName).GFP.GFP;
            data.cortical.(solenoidName).cortMUA(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).MUA.corticalData;
%             data.cortical.(solenoidName).hipMUA(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).MUA.hippocampalData;
            data.cortical.(solenoidName).cortGam(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).Gam.corticalData;
%             data.cortical.(solenoidName).hipGam(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).Gam.hippocampalData;
            data.cortical.(solenoidName).timeVector(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).timeVector;
            data.cortical.(solenoidName).cortS(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS;
            data.cortical.(solenoidName).cortS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.corticalS(49:end,20:23);
%             data.cortical.(solenoidName).hipS(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.hippocampalS;
%             data.cortical.(solenoidName).hipS_Gam(:,:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.hippocampalS(49:end,20:23);
            data.cortical.(solenoidName).T(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.T;
            data.cortical.(solenoidName).F(:,aa) = AnalysisResults.(animalID).Stim.cortical.(solenoidName).LFP.F;
        end
end
% concatenate the data from the contra and ipsi data
data.Contra.count = data.Z_NE.LPadSol.count;
data.Contra.Z_AchRhodamine = data.Z_Ach.RPadSol.Rhodamine;
data.Contra.Z_AchGFP = data.Z_Ach.RPadSol.GFP;

data.Contra.Z_NERhodamine = data.Z_NE.LPadSol.Rhodamine;
data.Contra.Z_NEGFP = data.Z_NE.LPadSol.GFP;

data.Contra.cortMUA = data.cortical.LPadSol.cortMUA;
% data.Contra.hipMUA = data.cortical.RPadSol.hipMUA;
data.Contra.cortGam = data.cortical.LPadSol.cortGam;
% data.Contra.hipGam = data.cortical.RPadSol.hipGam;
data.Contra.timeVector = data.cortical.LPadSol.timeVector;
data.Contra.cortS = data.cortical.LPadSol.cortS;
data.Contra.cortS_Gam = data.cortical.LPadSol.cortS_Gam;
% data.Contra.hipS = data.cortical.RPadSol.hipS;
% data.Contra.hipS_Gam = data.cortical.RPadSol.hipS_Gam;
data.Contra.T = data.cortical.LPadSol.T;
data.Contra.F = data.cortical.LPadSol.F;


data.Ipsi.count = data.Z_NE.RPadSol.count;
data.Ipsi.Z_AchRhodamine = data.Z_Ach.LPadSol.Rhodamine; 
data.Ipsi.Z_AchGFP = data.Z_Ach.LPadSol.GFP;

data.Ipsi.Z_NERhodamine = data.Z_NE.RPadSol.Rhodamine;
data.Ipsi.Z_NEGFP = data.Z_NE.RPadSol.GFP; 

data.Ipsi.cortMUA = data.cortical.RPadSol.cortMUA; 
% data.Ipsi.hipMUA = data.cortical.LPadSol.hipMUA;
data.Ipsi.cortGam = data.cortical.RPadSol.cortGam; 
% data.Ipsi.hipGam = data.cortical.LPadSol.hipGam;
data.Ipsi.timeVector = data.cortical.RPadSol.timeVector;
data.Ipsi.cortS = data.cortical.RPadSol.cortS; 
data.Ipsi.cortS_Gam = data.cortical.RPadSol.cortS_Gam;
% data.Ipsi.hipS = data.cortical.LPadSol.hipS;
% data.Ipsi.hipS_Gam = data.cortical.LPadSol.hipS_Gam;
data.Ipsi.T = data.cortical.RPadSol.T;
data.Ipsi.F = data.cortical.RPadSol.F;


data.Auditory.count = data.Z_NE.AudSol.count;
data.Auditory.Z_NERhodamine = data.Z_NE.AudSol.Rhodamine;
data.Auditory.Z_AchRhodamine = data.Z_Ach.AudSol.Rhodamine;

data.Auditory.Z_NEGFP = data.Z_NE.AudSol.GFP;
data.Auditory.Z_AchGFP = data.Z_Ach.AudSol.GFP;

data.Auditory.cortMUA = data.cortical.AudSol.cortMUA;
% data.Auditory.hipMUA = data.cortical.AudSol.hipMUA;
data.Auditory.cortGam = data.cortical.AudSol.cortGam;
% data.Auditory.hipGam = data.cortical.AudSol.hipGam;
data.Auditory.timeVector = data.cortical.AudSol.timeVector;
data.Auditory.cortS = data.cortical.AudSol.cortS;
data.Auditory.cortS_Gam = data.cortical.AudSol.cortS_Gam;
% data.Auditory.hipS = data.cortical.AudSol.hipS;
% data.Auditory.hipS_Gam = data.cortical.AudSol.hipS_Gam;
data.Auditory.T = data.cortical.AudSol.T;
data.Auditory.F = data.cortical.AudSol.F;


% take the averages of each field through the proper dimension
for ff = 1:length(compDataTypes)
    compDataType = compDataTypes{1,ff};
    data.(compDataType).mean_Count = mean(data.(compDataType).count,2);
    data.(compDataType).std_Count = std(data.(compDataType).count,0,2);   
    
    data.(compDataType).Z_Achmean_Rhodamine = mean(data.(compDataType).Z_AchRhodamine,2);
    data.(compDataType).Z_Achstd_Rhodamine = std(data.(compDataType).Z_AchRhodamine,0,2);

    data.(compDataType).Z_NEmean_Rhodamine = mean(data.(compDataType).Z_NERhodamine,2);
    data.(compDataType).Z_NEstd_Rhodamine = std(data.(compDataType).Z_NERhodamine,0,2);

    data.(compDataType).Z_Achmean_GFP = mean(data.(compDataType).Z_AchGFP,2);
    data.(compDataType).Z_Achstd_GFP = std(data.(compDataType).Z_AchGFP,0,2);

    data.(compDataType).Z_NEmean_GFP = mean(data.(compDataType).Z_NEGFP,2);
    data.(compDataType).Z_NEstd_GFP = std(data.(compDataType).Z_NEGFP,0,2);

    data.(compDataType).mean_CortMUA = mean(data.(compDataType).cortMUA,2);
    data.(compDataType).std_CortMUA = std(data.(compDataType).cortMUA,0,2);
%     data.(compDataType).mean_HipMUA = mean(data.(compDataType).hipMUA,2);
%     data.(compDataType).std_HipMUA = std(data.(compDataType).hipMUA,0,2);
    data.(compDataType).mean_CortGam = mean(data.(compDataType).cortGam,2);
    data.(compDataType).std_CortGam = std(data.(compDataType).cortGam,0,2);
%     data.(compDataType).mean_HipGam = mean(data.(compDataType).hipGam,2);
%     data.(compDataType).std_HipGam = std(data.(compDataType).hipGam,0,2);
    data.(compDataType).mean_timeVector = mean(data.(compDataType).timeVector,2);
    data.(compDataType).mean_CortS = mean(data.(compDataType).cortS,3).*100;
    data.(compDataType).mean_CortS_Gam = mean(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),3);
    data.(compDataType).std_CortS_Gam = std(mean(mean(data.(compDataType).cortS_Gam.*100,2),2),0,3);
%     data.(compDataType).mean_HipS = mean(data.(compDataType).hipS,3).*100;
%     data.(compDataType).mean_HipS_Gam = mean(mean(mean(data.(compDataType).hipS_Gam.*100,2),2),3);
%     data.(compDataType).std_HipS_Gam = std(mean(mean(data.(compDataType).hipS_Gam.*100,2),2),0,3);
    data.(compDataType).mean_T = mean(data.(compDataType).T,2);
    data.(compDataType).mean_F = mean(data.(compDataType).F,2);
end
%% Fig. 1-S2
summaryFigure = figure('Name','Fig1-S2 Stim');
sgtitle('Stimulus evoked responses')
%% [1-S2a] cortical MUA contra stim
ax1 = subplot(6,3,1);
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA,'-','color',colors('rich black'),'LineWidth',2);
hold on
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA + data.Contra.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.mean_CortMUA - data.Contra.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title('Contra stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];

%% [1-S2b] cortical MUA ispi stim
ax2 = subplot(6,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA,'-','color',colors('rich black'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA + data.Ipsi.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_CortMUA - data.Ipsi.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title('Ipsi stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [1-S2c] cortical MUA auditory stim
ax3 = subplot(6,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA,'-','color',colors('rich black'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA + data.Auditory.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.mean_CortMUA - data.Auditory.std_CortMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
title('Aud stim cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [1-S2d] cortical LFP contra stim
ax4 = subplot(6,3,4);
imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_CortS)
title('Contra stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
%caxis([-100,100])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [1-S2e] cortical LFP ispi stim
ax5 = subplot(6,3,5);
imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_CortS)
title('Ipsi stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [1-S2f] cortical LFP auditory stim
ax6 = subplot(6,3,6);
imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_CortS)
title('Aud stim cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-stimulus time (s)')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [1-S2g] hippocampal MUA contra stim
% ax7 = subplot(6,3,7);
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA,'-','color',colors('rich black'),'LineWidth',2)
% hold on
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA + data.Contra.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.Contra.mean_timeVector,data.Contra.mean_HipMUA - data.Contra.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Contra stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax7.TickLength = [0.03,0.03];
%% [1-S2h] hippocampal MUA ispi stim
% ax8 = subplot(6,3,8);
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA,'-','color',colors('rich black'),'LineWidth',2)
% hold on
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA + data.Ipsi.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.Ipsi.mean_timeVector,data.Ipsi.mean_HipMUA - data.Ipsi.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Ipsi stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax8.TickLength = [0.03,0.03];
%% [1-S2i] hippocampal MUA auditory stim
% ax9 = subplot(6,3,9);
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA,'-','color',colors('rich black'),'LineWidth',2)
% hold on
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA + data.Auditory.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% plot(data.Auditory.mean_timeVector,data.Auditory.mean_HipMUA - data.Auditory.std_HipMUA,'-','color',colors('battleship grey'),'LineWidth',0.10)
% title('Aud stim hippocampal MUA')
% ylabel('\DeltaP/P (%)')
% xlabel('Peri-stimulus time (s)')
% axis square
% set(gca,'box','off')
% ax9.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP contra stim
% ax10 = subplot(6,3,10);
% imagesc(data.Contra.mean_T,data.Contra.mean_F,data.Contra.mean_HipS)
% title('Contra stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c10 = colorbar;
% ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% axis square
% axis xy
% set(gca,'box','off')
% ax10.TickLength = [0.03,0.03];
%% [1-S2j] hippocampal LFP ispi stim
% ax11 = subplot(6,3,11);
% imagesc(data.Ipsi.mean_T,data.Ipsi.mean_F,data.Ipsi.mean_HipS)
% title('Ipsi stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c11 = colorbar;
% ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% axis square
% axis xy
% set(gca,'box','off')
% ax11.TickLength = [0.03,0.03];
%% [1-S2l] hippocampal LFP auditory stim
% ax12 = subplot(6,3,12);
% imagesc(data.Auditory.mean_T,data.Auditory.mean_F,data.Auditory.mean_HipS)
% title('Aud stim hippocampal LFP')
% ylabel('Freq (Hz)')
% xlabel('Peri-stimulus time (s)')
% c12 = colorbar;
% ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,100])
% axis square
% axis xy
% set(gca,'box','off')
% ax12.TickLength = [0.03,0.03];
%% [1-S2m] Rhodamine contra stim
ax13 = subplot(6,3,13);
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine + data.Contra.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine - data.Contra.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Contra stim Blood Volume')
ylabel('Z \DeltaF/F (Ach)')
ax13.YLim = [-3 6];

yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine + data.Contra.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine - data.Contra.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Z \DeltaF/F (NE)')

ax13.YAxis(1).Color = colors('indian red');
ax13.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
ax13.YLim = [-3 6];
%% [1-S2n] Rhodamine ispi stim
ax14 = subplot(6,3,14);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine + data.Ipsi.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine - data.Ipsi.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Ipsi stim Blood Volume')
ylabel('Z \DeltaF/F (Ach)')
ax14.YLim = [-3 6];

yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine + data.Ipsi.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine - data.Ipsi.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Z \DeltaF/F (NE)')
ax14.YAxis(1).Color = colors('indian red');
ax14.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
ax14.YLim = [-3 6];
%% [1-S2o] Rhodamine auditory stim
ax15 = subplot(6,3,15);
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine + data.Auditory.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine - data.Auditory.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Auditory stim Blood Volume')
ylabel('Z \DeltaF/F (Ach)')
ax15.YLim = [-3 6];

yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine + data.Auditory.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine - data.Auditory.Z_NEstd_Rhodamine,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('Z \DeltaF/F (NE)')
ax15.YAxis(1).Color = colors('indian red');
ax15.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
ax15.YLim = [-3 6];
%% [1-S2p] GFP contra stim
ax16 = subplot(6,3,16);
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP + data.Contra.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP - data.Contra.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
title('Contra Green')
ylabel('Z \DeltaF/F ACh')
ax16.YLim = [-3 6];

yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP + data.Contra.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP - data.Contra.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
title('Contra Green')
ylabel('Z \DeltaF/F GRAB NE')
ax16.YAxis(1).Color = colors('indian red');
ax16.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
ax16.YLim = [-3 6];
%% [1-S2q] GFP ispi stim
ax17 = subplot(6,3,17);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP + data.Ipsi.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP - data.Ipsi.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
title('Ipsi Green')
ylabel('Z \DeltaF/F ACh')
ax17.YLim = [-3 6];

yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP + data.Ipsi.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP - data.Ipsi.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
title('Ipsi Green')
ylabel('Z \DeltaF/F GRAB NE')
ax17.YAxis(1).Color = colors('indian red');
ax17.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
ax17.YLim = [-3 6];
%% [1-S2r] GFP auditory stim
ax18 = subplot(6,3,18);
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP + data.Auditory.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP - data.Auditory.Z_Achstd_GFP,'-','color',colors('indian red'),'LineWidth',0.10)
title('Auditory Green')
ylabel('Z \DeltaF/F ACh')
ax18.YLim = [-3 6];

yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP + data.Auditory.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP - data.Auditory.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
title('Auditory Green')
ylabel('Z \DeltaF/F GRAB NE')
ax18.YAxis(1).Color = colors('indian red');
ax18.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
ax18.YLim = [-3 6];
%% adjust and link axes
% linkaxes([ax1,ax2,ax3,ax7,ax8,ax9],'xy')
% linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
% linkaxes([ax13,ax14,ax15],'xy')
% linkaxes([ax16,ax17,ax18],'xy')
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
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath animalID 'Stim']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath animalID 'Stim'])
    close 
end
%% plot the comparison of GCaMP and Rhodamine with GRABNE and Rhodamine
summaryFigureN = figure('Name','Fig1-S2 Stim');
sgtitle('Stimulus evoked responses in fiber photometry signals')

% Ach contra stim
ax1 = subplot(3,3,1);
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine + data.Contra.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_Rhodamine - data.Contra.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Contra stim Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax1.YLim = [-3 6];

yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP + data.Contra.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_Achmean_GFP - data.Contra.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')

ax1.YAxis(1).Color = colors('indian red');
ax1.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
ax1.YLim = [-3 6];
xlim([-5 10])

% Ach ispi stim
ax2 = subplot(3,3,2);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine + data.Ipsi.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_Rhodamine - data.Ipsi.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Ipsi stim Ach fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax2.YLim = [-3 6];

yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP + data.Ipsi.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_Achmean_GFP - data.Ipsi.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')
ax2.YAxis(1).Color = colors('indian red');
ax2.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
ax2.YLim = [-3 6];
xlim([-5 10])

% Ach auditory stim
ax3 = subplot(3,3,3);
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine + data.Auditory.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_Rhodamine - data.Auditory.Z_Achstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Auditory stim Ach fiber ')
ylabel('\DeltaF/F Blood Volume (Z)')
ax3.YLim = [-3 6];

yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP + data.Auditory.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_Achmean_GFP - data.Auditory.Z_Achstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRAB ACh (Z)')
ax3.YAxis(1).Color = colors('indian red');
ax3.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
ax3.YLim = [-3 6];
xlim([-5 10])

% NE contra stim
ax4 = subplot(3,3,4);
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine + data.Contra.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_Rhodamine - data.Contra.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Contra stim NE fiber')
ylabel('\DeltaF/F Blood Volume (Z)')
ax4.YLim = [-3 6];

yyaxis right
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP + data.Contra.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Contra.mean_timeVector,data.Contra.Z_NEmean_GFP - data.Contra.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax4.YAxis(1).Color = colors('indian red');
ax4.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
ax4.YLim = [-3 6];
xlim([-5 10])

% NE ispi stim
ax5 = subplot(3,3,5);
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine + data.Ipsi.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_Rhodamine - data.Ipsi.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Ipsi stim NE fibers')
ylabel('\DeltaF/F Blood Volume (Z)')
ax5.YLim = [-3 6];

yyaxis right
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP + data.Ipsi.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Ipsi.mean_timeVector,data.Ipsi.Z_NEmean_GFP - data.Ipsi.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax5.YAxis(1).Color = colors('indian red');
ax5.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
ax5.YLim = [-3 6];
xlim([-5 10])

% NE auditory stim
ax6 = subplot(3,3,6);
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine,'-','color',colors('indian red'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine + data.Auditory.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_Rhodamine - data.Auditory.Z_NEstd_Rhodamine,'-','color',colors('indian red'),'LineWidth',0.10)
title('Auditory stim NE fibers')
ylabel('\DeltaF/F Blood Volume (Z)')
ax6.YLim = [-3 6];

yyaxis right
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP,'-','color',colors('army green'),'LineWidth',2)
hold on
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP + data.Auditory.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
plot(data.Auditory.mean_timeVector,data.Auditory.Z_NEmean_GFP - data.Auditory.Z_NEstd_GFP,'-','color',colors('army green'),'LineWidth',0.10)
ylabel('\DeltaF/F GRABNE (Z)')
ax6.YAxis(1).Color = colors('indian red');
ax6.YAxis(2).Color = colors('army green');
xlabel('Peri-stimulus time (s)')
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
ax6.YLim = [-3 6];
xlim([-5 10])

%save figure(s)
if strcmp(saveFigs,'y') == true
    if firstHrs == "false"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'lastHrs' delim];
    elseif firstHrs == "true"
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim 'firstHrs' delim];
    end    
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigureN,[dirpath animalID 'Stim-FiberSignals']);
    set(summaryFigureN,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath animalID 'Stim-FiberSignals'])
    close 
end

% end
