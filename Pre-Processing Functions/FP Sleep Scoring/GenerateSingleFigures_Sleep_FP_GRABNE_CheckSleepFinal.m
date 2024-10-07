function [figHandle,ax1,ax2,ax3,ax4,ax6,ax7] = GenerateSingleFigures_Sleep_FP_GRABNE_CheckSleepFinal(procDataFileID,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
% Adopted from Kevin Turner
%________________________________________________________________________________________________________________________
% Purpose: Create a summary figure for a single n minute FP trial and the
% sleep scores
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
%pupil 
filteredpupildiameter = filtfilt(sos1,g1,ProcData.data.Pupil.zDiameter);
% force sensor
filtForceSensor = filtfilt(sos1,g1,ProcData.data.forceSensor);
binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.emg;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% NE data
Ach_GFP = ProcData.data.GFP.Z_Ach;
filtAch_GFP = filtfilt(sos2,g2,Ach_GFP);
% NE data
NE_GFP = ProcData.data.GFP.Z_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);
% NE data
Ach_Rhodamine = ProcData.data.Rhodamine.Z_Ach;
filtAch_Rhodamine = filtfilt(sos2,g2,Ach_Rhodamine);
% NE data
NE_Rhodamine = ProcData.data.Rhodamine.Z_NE;
filtNE_Rhodamine = filtfilt(sos2,g2,NE_Rhodamine);
% cortical and hippocampal spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_RH.T;
F = SpecData.cortical_RH.F;
% Yvals for behavior Indices
indecesMax = max(filteredpupildiameter);
whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
forceInds = binForce.*force_Yvals;
whiskInds = binWhiskers.*whisking_Yvals;
LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));

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
%% get the sleep score
AwakeStage = 1*(ProcData.sleep.logicals.Manual.awakeLogical); 
NREMStage = 1*ProcData.sleep.logicals.Manual.nremLogical;
REMStage = 1*ProcData.sleep.logicals.Manual.remLogical;
AwakeStage (AwakeStage==0) = nan ;
NREMStage (NREMStage==0) = nan ;
REMStage (REMStage==0) = nan ;
indecesMax = max(filtNE_GFP);
AwakeYvals = 1.20 * indecesMax * AwakeStage; 
NREMYvals = 1.20 * indecesMax * NREMStage; 
REMYvals = 1.20 * indecesMax * REMStage; 
SleepDummy = 1:5:ProcData.notes.trialDuration_sec;
%% Figure
figHandle = figure;
% force sensor and EMG
ax1 = subplot(8,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' FP behavioral characterization and sleep scoring ' fileID2])
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
ax2 = subplot(8,1,2);
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
ax3 = subplot(8,1,3);
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
% GRABNE and SleepScore data
ax4 = subplot(8,1,4);
s1 = scatter(SleepDummy',AwakeYvals,'v','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
s2 = scatter(SleepDummy',NREMYvals,'v','MarkerFaceColor','b','MarkerEdgeColor','b');
s3 = scatter(SleepDummy',REMYvals,'v','MarkerFaceColor','r','MarkerEdgeColor','r');
p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('army green'),'LineWidth',1);
ylabel('\DeltaF/F  (Z)')
legend([s1,s2,s3,p8],'Awake','NREM','REM','GRAB NE')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
hold off

% Ach GCaMP7s
ax5 = subplot(8,1,5);
p8 = plot((1:length(filtAch_GFP))/ProcData.notes.dsFs,filtAch_GFP,'color',colors('north texas green'),'LineWidth',1);
ylabel('\DeltaF/F  (Z)')
legend('Ach GCaMP7s')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
hold off

% Ach GCaMP7s
ax8 = subplot(8,1,6);
p9 = plot((1:length(filtAch_Rhodamine))/ProcData.notes.dsFs,filtAch_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
hold on
p10 = plot((1:length(filtNE_Rhodamine))/ProcData.notes.dsFs,filtNE_Rhodamine,'color',colors('deep carrot orange'),'LineWidth',1);
ylabel('\DeltaF/F Rhodamine (Z)')
legend([p9,p10],'Ach','NE')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
hold off

% Right cortical electrode spectrogram
ax6 = subplot(8,1,7);
Semilog_ImageSC(T,F,cortical_RHnormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-100,100])
ylabel('Frequency (Hz)')
set(gca,'Yticklabel','10^1')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
yyaxis right
ylabel('Right cortical LFP')
set(gca,'Yticklabel',[])
%Hippocampal electrode spectrogram
ax7 = subplot(8,1,8);
Semilog_ImageSC(T,F,hippocampusNormS,'y')
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)')
caxis([-100,100])
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('Hippocampal LFP')
set(gca,'Yticklabel',[])
% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax6,ax7,ax8],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');

% ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);

ax1Pos(4) = ax1Pos(4)+0.03;
ax2Pos(4) = ax2Pos(4)+0.03;
ax3Pos(4) = ax3Pos(4)+0.03;
ax4Pos(4) = ax4Pos(4)+0.03;
ax5Pos(4) = ax5Pos(4)+0.03;
ax6Pos(4) = ax6Pos(4)+0.03;
ax7Pos(4) = ax7Pos(4)+0.03;
ax8Pos(4) = ax8Pos(4)+0.03;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
ax4Pos(3) = ax4Pos(3)+0.05;
ax5Pos(3) = ax5Pos(3)+0.05;
ax6Pos(3) = ax6Pos(3)+0.05;
ax7Pos(3) = ax7Pos(3)+0.05;
ax8Pos(3) = ax8Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
ax4Pos(1) = ax4Pos(1)-0.05;
ax5Pos(1) = ax5Pos(1)-0.05;
ax6Pos(1) = ax6Pos(1)-0.05;
ax7Pos(1) = ax7Pos(1)-0.05;
ax8Pos(1) = ax8Pos(1)-0.05;

ax1Pos(2) = ax1Pos(2)-0.05;
ax2Pos(2) = ax2Pos(2)-0.05;
ax3Pos(2) = ax3Pos(2)-0.05;
ax4Pos(2) = ax4Pos(2)-0.05;
ax5Pos(2) = ax5Pos(2)-0.05;
ax6Pos(2) = ax6Pos(2)-0.05;
ax7Pos(2) = ax7Pos(2)-0.05;
ax8Pos(2) = ax8Pos(2)-0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
%% save the file to directory.
if strcmp(saveFigs,'y') == true
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/SleepScoreCheck/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(figHandle,[dirpath animalID '_' fileID '_SleepScoreCheck']);
    saveas(figHandle,[dirpath animalID '_' fileID '_SleepScoreCheck'],'tiff')

end

end
