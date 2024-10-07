function [figHandle,ax1,ax2,ax3,ax4,ax5,ax6] = GenerateSingleFigures_FP_EEG(procDataFileID,saveFigs)
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
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.power;
% eeg
EEG_RH = ProcData.data.EEG_RH.EEGSignal;
EEG_LH = ProcData.data.EEG_LH.EEGSignal;
% pupil area
filteredpupilarea = filtfilt(sos1,g1,ProcData.data.Pupil.pupilArea);
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% Rhodamine and GFP data
% no LH
% LH_Rhodamine = ProcData.data.Rhodamine.LH;
% filtLH_Rhodamine = filtfilt(sos2,g2,LH_Rhodamine);
RH_Rhodamine = ProcData.data.Rhodamine.RH;
filtRH_Rhodamine = filtfilt(sos2,g2,RH_Rhodamine);
% LH_GFP = ProcData.data.GFP.LH;
% filtLH_GFP = filtfilt(sos2,g2,LH_GFP);
RH_GFP = ProcData.data.GFP.RH;
filtRH_GFP = filtfilt(sos2,g2,RH_GFP);
% EEG and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
EEG_LHnormS = SpecData.EEG_LH.normS.*100;
EEG_RHnormS = SpecData.EEG_RH.normS.*100;
% hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.EEG_LH.T;
F = SpecData.EEG_LH.F;
% Yvals for behavior Indices
% indecesMax = max([filtLH_Rhodamine;filtRH_Rhodamine;filtLH_GFP;filtRH_GFP]);
indecesMax = max([filtRH_Rhodamine;filtRH_GFP]);
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
title([animalID ' FP behavioral characterization and Rhodamine dynamics for ' fileID2])
ylabel('Force Sensor (Volts)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color',colors('deep carrot orange'),'LineWidth',1);
ylabel('EMG (Volts^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'force sensor','EMG')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Whisker angle and pupil area
ax2 = subplot(6,1,2);
p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
ylabel('Angle (deg)')
xlim([0,ProcData.notes.trialDuration_sec])
ylim([-20,60])

yyaxis right
p2 = plot((1:length(filteredpupilarea))/ProcData.notes.dsFs,filteredpupilarea,'color',colors('north texas green'),'LineWidth',1);
ylabel('Pupil Area (pixels)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Whisker Angle','Pupil Area')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% EEG 
ax3 = subplot(6,1,3);
p1 = plot((1:length(EEG_LH))/ProcData.notes.dsFs,EEG_LH,'color',colors('north texas green'),'LineWidth',1);
ylabel('ECOG LH(Volts)')
xlim([0,ProcData.notes.trialDuration_sec])
legend('ECOG LH')
% yyaxis right
% p2 = plot((1:length(EEG_RH))/ProcData.notes.dsFs,EEG_RH,'color',colors('deep carrot orange'),'LineWidth',1);
% ylabel('EEG RH')
% xlim([0,ProcData.notes.trialDuration_sec])
% legend([p1,p2],'EEG LH','EEG RH')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% Rhodamine and behavioral indeces
ax4 = subplot(6,1,4);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
% p5 = plot((1:length(filtLH_Rhodamine))/ProcData.notes.dsFs,filtLH_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
p6 = plot((1:length(filtRH_Rhodamine))/ProcData.notes.dsFs,filtRH_Rhodamine,'color',colors('sapphire'),'LineWidth',1);
ylabel('Rhodamine')
% legend([p5,p6,s1,s2],'LH Rhodamine','RH Rhodamine','movement','whisking')
legend([p6,s1,s2],'RH Rhodamine','movement','whisking')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% GFP and behavioral indeces
ax5 = subplot(6,1,5);
% p7 = plot((1:length(filtLH_GFP))/ProcData.notes.dsFs,filtLH_GFP,'color',colors('electric purple'),'LineWidth',1);
% hold on
p8 = plot((1:length(filtRH_GFP))/ProcData.notes.dsFs,filtRH_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('GFP')
% legend([p7,p8],'LH GFP','RH GFP')
legend('RH GFP')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Left EEG electrode spectrogram
ax7 = subplot(6,1,6);
Semilog_ImageSC(T,F,EEG_LHnormS,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('LH ECOG LFP')
set(gca,'Yticklabel', [])
xlabel('Time (sec)')
% Right EEG electrode spectrogram
% ax7 = subplot(6,1,6);
% Semilog_ImageSC(T,F,EEG_RHnormS,'y')
% axis xy
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)')
% caxis([-100,100])
% ylabel('Frequency (Hz)')
% set(gca,'Yticklabel','10^1')
% xlim([0,ProcData.notes.trialDuration_sec])
% xlabel('Time (sec)')
% set(gca,'TickLength',[0,0])
% % set(gca,'Xticklabel',[])
% set(gca,'box','off')
% yyaxis right
% ylabel('RH ECOG LFP')
% set(gca,'Yticklabel',[])

% Axes properties
% linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'x')
linkaxes([ax1,ax2,ax3,ax4,ax5,ax7],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');

% ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);

ax1Pos(4) = ax1Pos(4)+0.05;
ax2Pos(4) = ax2Pos(4)+0.05;
ax3Pos(4) = ax3Pos(4)+0.05;
ax4Pos(4) = ax4Pos(4)+0.05;
ax5Pos(4) = ax5Pos(4)+0.05;
% ax6Pos(4) = ax6Pos(4)+0.05;
ax7Pos(4) = ax7Pos(4)+0.05;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
ax4Pos(3) = ax4Pos(3)+0.05;
ax5Pos(3) = ax5Pos(3)+0.05;
% ax6Pos(3) = ax6Pos(3)+0.05;
ax7Pos(3) = ax7Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
ax4Pos(1) = ax4Pos(1)-0.05;
ax5Pos(1) = ax5Pos(1)-0.05;
% ax6Pos(1) = ax6Pos(1)-0.05;
ax7Pos(1) = ax7Pos(1)-0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);

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
