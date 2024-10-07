function [figHandle,ax1,ax2,ax3,ax4,ax5,ax7] = Generate_Sleep_Figures(procDataFileID,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure with sleep scores
%________________________________________________________________________________________________________________________
% [figHandle,ax1,ax2,ax3,ax4,ax5,ax7] = Generate_Sleep_Figures(procDataFileID,saveFigs)

%% load file and gather information
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
EMGSignal = ProcData.data.EMG.emgSignal;



% NE data
NE_GFP = ProcData.data.GFP.P_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);

% ACh data
ACh_GFP = ProcData.data.GFP.P_ACh;
filtACh_GFP = filtfilt(sos2,g2,ACh_GFP);

% CBV data
NE_CBV = ProcData.data.CBV.P_NE;
filtNE_CBV = filtfilt(sos2,g2,NE_CBV);
ACh_CBV = ProcData.data.CBV.P_ACh;
filtACh_CBV = filtfilt(sos2,g2,ACh_CBV);

% ECoG spectrograms
specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
load(specDataFile,'-mat');
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;

% ECoG
EEG_LH = ProcData.data.cortical_LH.corticalSignal;
EEG_LH = medfilt1(EEG_LH,3);

% Yvals for behavior Indices
indecesMax = max(filteredpupildiameter);
whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
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

% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
OptoStim = ProcData.data.stimulations.OptoStim;

LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));
OptoStim_Yvals = 1.30*max(indecesMax)*ones(size(OptoStim));
%% get the sleep score
AwakeStage = double(ProcData.sleep.logicals.Manual.awakeLogical);
NREMStage = double(ProcData.sleep.logicals.Manual.nremLogical);
REMStage = double(ProcData.sleep.logicals.Manual.remLogical);

AwakeStage (AwakeStage==0) = nan ;
NREMStage (NREMStage==0) = nan ;
REMStage (REMStage==0) = nan ;

% assign new values based on the NE data
indecesMax = max(filtNE_GFP);
AwakeYvals = 1.20 * indecesMax * AwakeStage; 
NREMYvals = 1.20 * indecesMax * NREMStage; 
REMYvals = 1.20 * indecesMax * REMStage; 
SleepDummy = 1:5:ProcData.notes.trialDuration_sec;
%% Figure
figHandle = figure;
fileID2 = strrep(fileID,'_',' ');
figHandle.WindowState = 'minimized';
%% force sensor and EMG
ax1 = subplot(7,1,1);

p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' Manual sleep scores: ' fileID2])

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
%% Whisker angle
ax2 = subplot(7,1,2);

p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
hold on
ylabel('Angle (deg)')

yyaxis right
p2 = plot((1:length(EMGSignal))/ProcData.notes.dsFs,EMGSignal,'color','k','LineWidth',1);
ylabel('EMG (V)')

xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Whisker Angle','EMG')

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% pupil and behavioral indeces
ax3 = subplot(7,1,3);

s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
s6 = scatter(OptoStim,OptoStim_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p5 = plot((1:length(filteredpupildiameter))/ProcData.notes.dsFs,filteredpupildiameter,'color',colors('dark candy apple red'),'LineWidth',1);

ylabel('Diameter (Z)')

if (isempty(LPad_Yvals) == 1) && (isempty(RPad_Yvals) == 1) && (isempty(Aud_Yvals) == 1)
    if isempty(OptoStim_Yvals) == 1
            legend([p5,s1,s2],'Pupil Diameter','movement','whisking')
    else
            legend([p5,s1,s2,s6],'Pupil Diameter','movement','whisking','OptoStim')
    end
else
    legend([p5,s1,s2,s3,s4,s5,s6],'Pupil Diameter','movement','whisking','LSol','RSol','AudSol','OptoStim')
end

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% GRABNE and Sleep Score data
ax4 = subplot(7,1,4);

s1 = scatter(SleepDummy',AwakeYvals,'v','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
s2 = scatter(SleepDummy',NREMYvals,'v','MarkerFaceColor','b','MarkerEdgeColor','b');
s3 = scatter(SleepDummy',REMYvals,'v','MarkerFaceColor','r','MarkerEdgeColor','r');
p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('vegas gold'),'LineWidth',1);
p9 = plot((1:length(filtNE_CBV))/ProcData.notes.dsFs,filtNE_CBV,'color','red','LineWidth',1);

ylabel('\Delta F/F (%)')
legend([s1,s2,s3,p8,p9],'Awake','NREM','REM','GRAB NE','CBV')

xlim([0,ProcData.notes.trialDuration_sec])

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% ACh
ax5 = subplot(7,1,5);

plot((1:length(filtACh_GFP))/ProcData.notes.dsFs,filtACh_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('\Delta F/F (%) ACh')
hold on

% CBV 
yyaxis right
p18 = plot((1:length(filtACh_CBV))/ProcData.notes.dsFs,filtACh_CBV,'color',colors('dark candy apple red'),'LineWidth',1);
ylabel('\Delta F/F (%) CBV')

legend('GRAB ACh','Blood Volume')
xlim([0,ProcData.notes.trialDuration_sec])

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% EEG
ax6 = subplot(7,1,6);

plot((1:length(EEG_LH))/ProcData.notes.dsFs,EEG_LH,'color',colors('deep carrot orange'),'LineWidth',1);
xlim([0,ProcData.notes.trialDuration_sec])
ylabel('ECoG (V)')
legend('ECoG')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% ECoG spectrogram
ax7 = subplot(7,1,7);

Semilog_ImageSC(T,F,cortical_LHnormS,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)')
clim([-100,100])
ylabel('Frequency (Hz)')
yticks([1 4 8 15 30 100])
xlim([0,ProcData.notes.trialDuration_sec])
xlabel('Time (sec)')
set(gca,'TickLength',[0,0])
set(gca,'box','off')
yyaxis right
ylabel('Right cortical LFP')
%% Axes properties
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');

ax7Pos(3:4) = ax6Pos(3:4);

ax1Pos(4) = ax1Pos(4)+0.05;
ax2Pos(4) = ax2Pos(4)+0.05;
ax3Pos(4) = ax3Pos(4)+0.05;
ax4Pos(4) = ax4Pos(4)+0.05;
ax5Pos(4) = ax5Pos(4)+0.05;
ax6Pos(4) = ax6Pos(4)+0.03;
ax7Pos(4) = ax7Pos(4)+0.05;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
ax4Pos(3) = ax4Pos(3)+0.05;
ax5Pos(3) = ax5Pos(3)+0.05;
ax6Pos(3) = ax6Pos(3)+0.05;
ax7Pos(3) = ax7Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
ax4Pos(1) = ax4Pos(1)-0.05;
ax5Pos(1) = ax5Pos(1)-0.05;
ax6Pos(1) = ax6Pos(1)-0.05;
ax7Pos(1) = ax7Pos(1)-0.05;

ax1Pos(2) = ax1Pos(2)-0.05;
ax2Pos(2) = ax2Pos(2)-0.05;
ax3Pos(2) = ax3Pos(2)-0.05;
ax4Pos(2) = ax4Pos(2)-0.05;
ax5Pos(2) = ax5Pos(2)-0.05;
ax6Pos(2) = ax6Pos(2)-0.05;
ax7Pos(2) = ax7Pos(2)-0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
    %% save the file to directory.
    if strcmp(saveFigs,'y') == true
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Single Trial Figures/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(figHandle,[dirpath animalID '_' fileID '_SleepScore']);
        saveas(figHandle,[dirpath animalID '_' fileID '_SleepScore'],'tiff')
    end
    close all
end
