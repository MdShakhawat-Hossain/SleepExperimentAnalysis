function [figHandle,ax1,ax2,ax3,ax4,ax5,min_EEG_LH,max_EEG_LH] = GenerateSingleFigures_Sleep_FP_EEG(procDataFileID,saveFigs)
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
%pupil 
filteredpupildiameter = filtfilt(sos1,g1,ProcData.data.Pupil.zDiameter);
% eeg
% EEG_RH = ProcData.data.EEG_RH.EEGSignal;
EEG_LH = ProcData.data.EEG_LH.EEGSignal;
% remove some extra data
% EEG_RH(1:ProcData.notes.dsFs) = EEG_RH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
EEG_LH(1:ProcData.notes.dsFs) = EEG_LH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
EEG_LH = medfilt1(EEG_LH,3);

max_EEG_LH = max(EEG_LH);
min_EEG_LH = min(EEG_LH);
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% Rhodamine and GFP data
% no LH
% LH_Rhodamine = ProcData.data.Rhodamine.LH;
% filtLH_Rhodamine = filtfilt(sos2,g2,LH_Rhodamine);
% RH_Rhodamine = ProcData.data.Rhodamine.RH;
% filtRH_Rhodamine = filtfilt(sos2,g2,RH_Rhodamine);
% LH_GFP = ProcData.data.GFP.LH;
% filtLH_GFP = filtfilt(sos2,g2,LH_GFP);
% RH_GFP = ProcData.data.GFP.RH;
% filtRH_GFP = filtfilt(sos2,g2,RH_GFP);
% EEG and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
EEG_LHnormS = SpecData.EEG_LH.normS.*100;
% EEG_RHnormS = SpecData.EEG_RH.normS.*100;
% hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.EEG_LH.T;
F = SpecData.EEG_LH.F;
% Yvals for behavior Indices
% indecesMax = max([filtLH_Rhodamine;filtRH_Rhodamine;filtLH_GFP;filtRH_GFP]);
indecesMax = max(filteredpupildiameter);% change this to match the range of the data
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
ax1 = subplot(5,1,1);
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
% Whisker angle
ax2 = subplot(5,1,2);
plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
ylabel('Angle (deg)')
xlim([0,ProcData.notes.trialDuration_sec])
ylim([-20,60])
legend('Whisker Angle')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% pupil and behavioral indeces
ax3 = subplot(5,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p5 = plot((1:length(filteredpupildiameter))/ProcData.notes.dsFs,filteredpupildiameter,'color',colors('dark candy apple red'),'LineWidth',1);
legend([p5,s1,s2,s3,s4,s5],'Pupil Diameter','movement','whisking','LSol','RSol','AudSol')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% EEG
ax4 = subplot(5,1,4);
p5 = plot((1:length(EEG_LH))/ProcData.notes.dsFs,EEG_LH,'color',colors('deep carrot orange'),'LineWidth',1);
xlim([0,ProcData.notes.trialDuration_sec])
ylabel('EEG LH (Volts)')
ylim(ax4,[min_EEG_LH,max_EEG_LH])
legend([p5],'EEG LH')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Left EEG electrode spectrogram
ax5 = subplot(5,1,5);
Semilog_ImageSC(T,F,EEG_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('Left EEG LFP')
set(gca,'Yticklabel', [])
xlabel('Time(S)')


linkaxes([ax1,ax2,ax3,ax4,ax5],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax5Pos(3:4) = ax1Pos(3:4);

ax1Pos(4) = ax1Pos(4)+0.05;
ax2Pos(4) = ax2Pos(4)+0.05;
ax3Pos(4) = ax3Pos(4)+0.05;
ax4Pos(4) = ax4Pos(4)+0.05;
ax5Pos(4) = ax5Pos(4)+0.05;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
ax4Pos(3) = ax4Pos(3)+0.05;
ax5Pos(3) = ax5Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
ax4Pos(1) = ax4Pos(1)-0.05;
ax5Pos(1) = ax5Pos(1)-0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
%% save the file to directory.
if strcmp(saveFigs,'y') == true
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/Single Trial Figures/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig_Sleep']);
    saveas(figHandle,[dirpath animalID '_' fileID '_SingleTrialFig_Sleep'],'tiff')

end

end
