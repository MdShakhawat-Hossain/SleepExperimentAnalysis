function [figHandle,ax1,ax2,ax3,ax4,ax5,ax6] = GenerateSingleFigures_FP_GRABNE_SingleFiber(procDataFileID,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure for a single n minute FP trial
%________________________________________________________________________________________________________________________

% load file and gather information
load(procDataFileID)
[animalID,fileDate,fileID] = GetFileInfo_FP(procDataFileID);
strDay = ConvertDate_FP(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,1/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;
% pupil area
filteredpupilarea = filtfilt(sos1,g1,ProcData.data.Pupil.pupilArea);
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.emg;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% CBV and GFP data
% no LH
% Ach_CBV = ProcData.data.CBV.P_Ach;
% filtAch_CBV = filtfilt(sos2,g2,Ach_CBV);
NE_CBV = ProcData.data.CBV.P_NE;
filtNE_CBV = filtfilt(sos2,g2,NE_CBV);
% Ach_GFP = ProcData.data.GFP.P_Ach;
% filtAch_GFP = filtfilt(sos2,g2,Ach_GFP);
NE_GFP = ProcData.data.GFP.P_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
% cortical_RHnormS = SpecData.cortical_RH.normS.*100;
% hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;

EEG_LH = ProcData.data.cortical_LH.corticalSignal;
% remove some extra data
EEG_LH(1:ProcData.notes.dsFs) = EEG_LH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
EEG_LH = medfilt1(EEG_LH,3);

% Yvals for behavior Indices
indecesMax = max([filtNE_CBV;filtNE_GFP]);
% indecesMax = max([filtAch_CBV;filtAch_GFP]);
whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));
forceInds = binForce.*force_Yvals;
whiskInds = binWhiskers.*whisking_Yvals;
% set force indeces
for x = 1:length(forceInds)
    if forceInds(1,x) == 0
        forceInds(1,x) = NaN;
    end
end
% set whisk indeces
for x = 1:length(whiskInds)
    if whiskInds(1,x) == 0
        whiskInds(1,x) = NaN;
    end
end
%% Figure
figHandle = figure;
% force sensor and EMG
ax1 = subplot(6,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' FP behavioral characterization and CBV dynamics for ' fileID2])
ylabel('Force Sensor (V)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
ylabel('EMG (V^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'force sensor','EMG')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Whisker angle and pupil area
ax2 = subplot(6,1,2);
p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
ylabel('Angle (deg)')
xlim([0,ProcData.notes.trialDuration_sec])
% ylim([-20,60])

yyaxis right
p2 = plot((1:length(filteredpupilarea))/ProcData.notes.dsFs,filteredpupilarea,'color',colors('north texas green'),'LineWidth',1);
ylabel('Pupil Area (pixels)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Whisker Angle','Pupil Area')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% CBV and behavioral indeces
ax3 = subplot(6,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
% p5 = plot((1:length(filtAch_CBV))/ProcData.notes.dsFs,filtAch_CBV,'color',colors('dark candy apple red'),'LineWidth',1);
p6 = plot((1:length(filtNE_CBV))/ProcData.notes.dsFs,filtNE_CBV,'color',colors('sapphire'),'LineWidth',1);
ylabel('\Delta F/F')
legend([p6,s1,s2,s3,s4,s5],'NE CBV','movement','whisking','LSol','RSol','AudSol')
% legend([p6,s1,s2],'RH CBV','movement','whisking')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Ach Ca2+ activity
% ax4 = subplot(5,1,4);
% % p7 = plot((1:length(filtAch_GFP))/ProcData.notes.dsFs,filtAch_GFP,'color',colors('electric purple'),'LineWidth',1);
% % hold on
% % p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('vegas gold'),'LineWidth',1);
% ylabel('\delta F/F')
% legend(p7,'Ach GCaMP7s')
% xlim([0,ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight

% GRAB NE
ax4 = subplot(6,1,4);
p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('\Delta F/F')
legend(p8,'GRAB NE')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% EEG
ax5 = subplot(6,1,5);
plot((1:length(EEG_LH))/ProcData.notes.dsFs,EEG_LH,'color',colors('deep carrot orange'),'LineWidth',1);
xlim([0,ProcData.notes.trialDuration_sec])
ylabel('ECoG (V)')
% ylim(ax5,[min_EEG_LH,max_EEG_LH])
legend('ECoG')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% Left cortical electrode spectrogram
ax6 = subplot(6,1,6);
Semilog_ImageSC(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
xlabel('Time (sec)')
set(gca,'box','off')
yyaxis right
ylabel('Left cort LFP')
set(gca,'Yticklabel', [])
% Right cortical electrode spectrogram
% ax7 = subplot(5,1,7);
% Semilog_ImageSC(T,F,cortical_RHnormS,'y')
% axis xy
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)')
% caxis([-100,100])
% ylabel('Frequency (Hz)')
% set(gca,'Yticklabel','10^1')
% xlim([0,ProcData.notes.trialDuration_sec])
% xlabel('Time (sec)')
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% yyaxis right
% ylabel('Right cort LFP')
% set(gca,'Yticklabel',[])
%Hippocampal electrode spectrogram
% ax8 = subplot(5,1,8);
% Semilog_ImageSC(T,F,hippocampusNormS,'y')
% c7 = colorbar;
% ylabel(c7,'\DeltaP/P (%)')
% caxis([-100,100])
% xlabel('Time (sec)')
% ylabel('Frequency (Hz)')
% xlim([0,ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0,0])
% set(gca,'box','off')
% yyaxis right
% ylabel('Hippocamp LFP')
% set(gca,'Yticklabel',[])
% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
% ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');

% ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
% ax7Pos(3:4) = ax1Pos(3:4);
% ax8Pos(3:4) = ax1Pos(3:4);


ax1Pos(4) = ax1Pos(4)+0.03;
ax2Pos(4) = ax2Pos(4)+0.03;
ax3Pos(4) = ax3Pos(4)+0.03;
ax4Pos(4) = ax4Pos(4)+0.03;
ax5Pos(4) = ax5Pos(4)+0.03;
ax6Pos(4) = ax6Pos(4)+0.03;
% ax7Pos(4) = ax7Pos(4)+0.03;
% ax8Pos(4) = ax8Pos(4)+0.03;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
ax4Pos(3) = ax4Pos(3)+0.05;
ax5Pos(3) = ax5Pos(3)+0.05;
ax6Pos(3) = ax6Pos(3)+0.05;
% ax7Pos(3) = ax7Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
ax4Pos(1) = ax4Pos(1)-0.05;
ax5Pos(1) = ax5Pos(1)-0.05;
ax6Pos(1) = ax6Pos(1)-0.05;
% ax7Pos(1) = ax7Pos(1)-0.05;
% ax8Pos(1) = ax8Pos(1)-0.05;

ax1Pos(2) = ax1Pos(2)-0.05;
ax2Pos(2) = ax2Pos(2)-0.05;
ax3Pos(2) = ax3Pos(2)-0.05;
ax4Pos(2) = ax4Pos(2)-0.05;
ax5Pos(2) = ax5Pos(2)-0.05;
ax6Pos(2) = ax6Pos(2)-0.05;
% ax7Pos(2) = ax7Pos(2)-0.05;
% ax8Pos(2) = ax8Pos(2)-0.05;

ax6Pos(3:4) = ax1Pos(3:4);
% ax7Pos(3:4) = ax1Pos(3:4);
% ax8Pos(3:4) = ax1Pos(3:4);

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
% set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
%% save the file to directory.
if strcmp(saveFigs,'y') == true
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig']);
    saveas(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig'],'tiff')

end

end
