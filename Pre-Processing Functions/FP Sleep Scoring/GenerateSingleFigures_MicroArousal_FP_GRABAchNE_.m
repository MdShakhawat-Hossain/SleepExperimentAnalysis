function [figHandle,ax1,ax2,ax3,ax4,ax5,ax7] = GenerateSingleFigures_MicroArousal_FP_GRABAchNE_(procDataFileID,saveFigs,MALabels,EMGArousalLabels,Notes,MicroLabels)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
% Adopted from Kevin Turner
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
[z3,p3,k3] = butter(4,0.1/(ProcData.notes.dsFs/2),'low');
[sos3,g3] = zp2sos(z3,p3,k3);

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
EMGSignal = ProcData.data.EMG.emgSignal;
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% NE data
NE_GFP = ProcData.data.GFP.Z_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);

% Ach data
Ach_GFP = ProcData.data.GFP.Z_Ach;
filtAch_GFP = filtfilt(sos2,g2,Ach_GFP);
% Rhodamine data
NE_Rhodamine = ProcData.data.Rhodamine.Z_NE;
filtNE_Rhodamine = filtfilt(sos2,g2,NE_Rhodamine);
% cortical and hippocampal spectrograms
% specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
% load(specDataFile,'-mat');
% cortical_LHnormS = SpecData.cortical_LH.normS.*100;
% T = SpecData.cortical_LH.T;
% F = SpecData.cortical_LH.F;

% EEG_LH = ProcData.data.cortical_LH.corticalSignal;
% % remove some extra data
% EEG_LH(1:ProcData.notes.dsFs) = EEG_LH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
% EEG_LH = medfilt1(EEG_LH,3);
% Yvals for behavior Indices
indecesMax = max(filteredpupildiameter);
% indecesMax = max([filtAch_Rhodamine;filtAch_GFP]);
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
AwakeStage = zeros(3120,1);
NREMStage = zeros(3120,1);
REMStage = zeros(3120,1);

AwakeStage(MALabels==1) = 1;
NREMStage(MALabels==2) = 1;
REMStage(MALabels==3) = 1;

AwakeStage (AwakeStage==0) = nan ;
NREMStage (NREMStage==0) = nan ;
REMStage (REMStage==0) = nan ;
indecesMax = max(filtNE_GFP);
AwakeYvals = 1.20 * indecesMax * AwakeStage; 
NREMYvals = 1.20 * indecesMax * NREMStage; 
REMYvals = 1.20 * indecesMax * REMStage; 
SleepDummy = 1:1:ProcData.notes.trialDuration_sec;


indecesMax = max(filtAch_GFP);
MicroLabels(MicroLabels == 0) = nan;
MicroLabelYvals = 1.20 * indecesMax * MicroLabels; 
%% Detect the microarousals based on EMG activation
EMGD = abs(ProcData.data.EMG.emgSignal);
EMGThreshold = prctile(EMGD,Notes.EMGThresholdPerc)*ones(size(EMGD));
emgBinarized = EMGD;
emgBinarized(emgBinarized > EMGThreshold) = 1;
emgBinarized(emgBinarized < EMGThreshold) = 0;
emgBinarized(emgBinarized ==0 ) = nan;
indecesMax = max(filteredpupildiameter);
emg_Yvals = 1.30*max(indecesMax)*ones(size(emgBinarized));
emg_Idx = emg_Yvals.*emgBinarized;

indecesMax = max(filtNE_Rhodamine);
EMGArousalLabels(EMGArousalLabels==0) = nan;
EMGArousalLabels = 1.30*max(indecesMax)*EMGArousalLabels;
%% gamma power
% gammaPower = ProcData.data.cortical_LH.gammaBandPower;
% filtgammaPower = filtfilt(sos3,g3,gammaPower);
% %derivative 
% gammaPower_Diff = diff(gammaPower);
% filtgammaPower_Diff = filtfilt(sos3,g3,gammaPower_Diff);
% filtgammaPower_Diff(1:90) = filtgammaPower_Diff(91:180);
% % Threshold
% GammaThreshold = prctile(filtgammaPower_Diff,Notes.GammaThresholdPerc)*ones(size(filtgammaPower_Diff));
% 
% indecesMax = max(filtNE_Rhodamine);
% GammaArousalLabels(GammaArousalLabels==0) = nan;
% GammaArousalLabels = 1.30*max(indecesMax)*GammaArousalLabels;
%% Figure
figHandle = figure;
% force sensor and EMG
ax1 = subplot(6,1,1);
fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' FP behavioral characterization and sleep scoring ' fileID2])
ylabel('Force Sensor (V)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/ProcData.notes.dsFs,EMG,'color','k','LineWidth',1);
ylabel('EMG Power (V^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Force Sensor','EMG Power')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% Whisker angle
ax2 = subplot(6,1,2);
p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
hold on
ylabel('Angle (deg)')
yyaxis right
p2 = plot((1:length(EMGSignal))/ProcData.notes.dsFs,EMGSignal,'color','k','LineWidth',1);
p3 = plot((1:length(EMGD))/ProcData.notes.dsFs,EMGThreshold,'-r');
ylabel('EMG (V)')

xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2,p3],'Whisker Angle','EMG','EMGThreshold')
% ylim([-20,60])

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% pupil and behavioral indeces
ax3 = subplot(6,1,3);
s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
s6 = scatter((1:length(emg_Idx))/ProcData.notes.dsFs,emg_Idx,'.','MarkerEdgeColor','k');
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p5 = plot((1:length(filteredpupildiameter))/ProcData.notes.dsFs,filteredpupildiameter,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('Diameter (Z)')
legend([p5,s1,s2,s6,s3,s4,s5],'Pupil Diameter','movement','whisking','emg','LSol','RSol','AudSol')

% yyaxis right
% p9 = plot((1:length(filtgammaPower))/ProcData.notes.dsFs,filtgammaPower,'color','b','LineWidth',1);
% ylabel('sigmaPower')
% legend([p5,s1,s2,s6,s3,s4,s5,p9],'Pupil Diameter','movement','whisking','emg','LSol','RSol','AudSol','gammaPower')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% GRABNE and Sleep Score data
ax4 = subplot(6,1,4);
s1 = scatter(SleepDummy',AwakeYvals,'v','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
s2 = scatter(SleepDummy',NREMYvals,'v','MarkerFaceColor','b','MarkerEdgeColor','b');
s3 = scatter(SleepDummy',REMYvals,'v','MarkerFaceColor','r','MarkerEdgeColor','r');
% sn = scatter(1:length(EMGArousalLabels),EMGArousalLabels,'v','MarkerEdgeColor','m','MarkerFaceColor','m');
% sg = scatter(1:length(GammaArousalLabels),GammaArousalLabels,'v','MarkerEdgeColor','g','MarkerFaceColor','g');
p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('\Delta F/F (Z)')
legend([s1,s2,s3,p8],'Awake','NREM','REM','GRAB NE')

% yyaxis right

% p9 = plot((1:length(filtgammaPower_Diff))/ProcData.notes.dsFs,abs(filtgammaPower_Diff),'color','b','LineWidth',1);
% p10 = plot((1:length(GammaThreshold))/ProcData.notes.dsFs,GammaThreshold,'-k','LineWidth',1);

% ylabel('Power Derivative')
% legend([s1,s2,s3,p8,p9,p10],'Awake','NREM','REM','GRAB NE','Gamma Power Der','Gamma Threshold')
xlim([0,ProcData.notes.trialDuration_sec])

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% Rhodamine and MicroArousals Score data
ax7 = subplot(6,1,5);
% sn = scatter(1:length(EMGArousalLabels),EMGArousalLabels,'v','MarkerEdgeColor','m','MarkerFaceColor','m');
% hold on
% sg = scatter(1:length(GammaArousalLabels),GammaArousalLabels,'v','MarkerEdgeColor','g','MarkerFaceColor','g');
p18 = plot((1:length(filtNE_Rhodamine))/ProcData.notes.dsFs,filtNE_Rhodamine,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\Delta F/F (Z)')
legend('Blood Volume')
% legend([sn,p18],'EMG Arousals','Rhodamine')
% legend([sn,sg,p18],'EMG Arousals','Gamma Arousals','Rhodamine')
xlim([0,ProcData.notes.trialDuration_sec])

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% Ach
ax5 = subplot(6,1,6);
% GRABNE and Sleep Score data
plot((1:length(filtAch_GFP))/ProcData.notes.dsFs,filtAch_GFP,'color',colors('vegas gold'),'LineWidth',1);
hold on
s1 = scatter(SleepDummy',MicroLabelYvals,'v','MarkerFaceColor','k','MarkerEdgeColor','k');
ylabel('\Delta F/F (Z)')
legend('GRAB Ach','MicroArousals')

xlim([0,ProcData.notes.trialDuration_sec])

set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

% plot((1:length(EEG_LH))/ProcData.notes.dsFs,EEG_LH,'color',colors('deep carrot orange'),'LineWidth',1);
% xlim([0,ProcData.notes.trialDuration_sec])
% ylabel('ECoG (V)')
% legend('ECoG')
% xlim([0,ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight

% % Left cortical electrode spectrogram
% ax6 = subplot(6,1,[7,8]);
% Semilog_ImageSC(T,F,cortical_LHnormS,'y')
% axis xy
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)')
% caxis([-100,100])
% ylabel('Frequency (Hz)')
% % set(gca,'Yticklabel','10^1')
% yticks([1 4 8 15 30 100])
% xlim([0,ProcData.notes.trialDuration_sec])
% xlabel('Time (sec)')
% set(gca,'TickLength',[0,0])
% set(gca,'box','off')
% yyaxis right
% ylabel('Right cortical LFP')
% set(gca,'Yticklabel',[])

% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax5,ax7],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax6Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3) = ax2Pos(3);

ax1Pos(4) = ax1Pos(4)+0.05;
ax2Pos(4) = ax2Pos(4)+0.05;
ax3Pos(4) = ax3Pos(4)+0.05;
ax4Pos(4) = ax4Pos(4)+0.05;
ax5Pos(4) = ax5Pos(4)+0.05;
ax7Pos(4) = ax7Pos(4)+0.05;
% ax6Pos(4) = ax6Pos(4)+0.03;

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

ax1Pos(2) = ax1Pos(2)-0.05;
ax2Pos(2) = ax2Pos(2)-0.05;
ax3Pos(2) = ax3Pos(2)-0.05;
ax4Pos(2) = ax4Pos(2)-0.05;
ax5Pos(2) = ax5Pos(2)-0.05;
% ax6Pos(2) = ax6Pos(2)-0.05;
ax7Pos(2) = ax7Pos(2)-0.05;

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
    savefig(figHandle,[dirpath animalID '_' fileID '_MicroArousalSleep']);
%     saveas(figHandle,[dirpath animalID '_' fileID '_MicroArousalSleep'],'tiff')

end

end
