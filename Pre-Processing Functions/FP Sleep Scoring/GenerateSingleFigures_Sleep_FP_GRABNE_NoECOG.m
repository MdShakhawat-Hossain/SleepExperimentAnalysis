function [figHandle,ax1,ax2,ax3,ax4,ax5,ax6,ax7] = GenerateSingleFigures_Sleep_FP_GRABNE_NoECOG(procDataFileID,saveFigs)
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
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
binWhiskers = ProcData.data.binWhiskerAngle;
%pupil 
filteredpupildiameter = filtfilt(sos1,g1,ProcData.data.Pupil.zDiameter);
% filtereddistanceTraveled = ProcData.data.Pupil.distanceTraveled;
% filteredCentroidX = ProcData.data.Pupil.CentroidX;
% filteredCentroidY = ProcData.data.Pupil.CentroidY;

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
NE_GFP = ProcData.data.GFP.P_NE;
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);
% ACh data
ACh_GFP = ProcData.data.GFP.P_ACh;
filtACh_GFP = filtfilt(sos2,g2,ACh_GFP);
% CBV data
ACh_CBV = ProcData.data.Rhodamine.Z_ACh;
filtACh_CBV = filtfilt(sos2,g2,ACh_CBV);

% Yvals for behavior Indices
indecesMax = max(filteredpupildiameter);
% indecesMax = max([filtACh_Rhodamine;filtACh_GFP]);
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
%% Figure
figHandle = figure;
fileID2 = strrep(fileID,'_',' ');
%% force sensor and EMG
ax1 = subplot(6,1,1);

p1 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
title([animalID ' FP behavioral characterization and sleep scoring ' fileID2])
ylabel('Force Sensor (V)')

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
ax2 = subplot(6,1,2);

p1 = plot((1:length(filteredWhiskerAngle))/ProcData.notes.dsFs,filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
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
ax3 = subplot(6,1,3);

s1 = scatter((1:length(binForce))/ProcData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors('north texas green'));
hold on
s2 = scatter((1:length(binWhiskers))/ProcData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors('dark pink'));
s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
p5 = plot((1:length(filteredpupildiameter))/ProcData.notes.dsFs,filteredpupildiameter,'color',colors('dark candy apple red'),'LineWidth',1);

legend([p5,s1,s2,s3,s4,s5],'Pupil Diameter','movement','whisking','LSol','RSol','AudSol')
ylabel('Diameter (Z)')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight

%% pupil centriod locations
% ax4 = subplot(6,1,4);
% 
% p1 = plot((1:length(filteredCentroidX))/ProcData.notes.dsFs,filteredCentroidX,'color','b','LineWidth',1);
% ylabel('Centriod X')
% 
% yyaxis right
% p2 = plot((1:length(filteredCentroidY))/ProcData.notes.dsFs,filteredCentroidY,'color','k','LineWidth',1);
% ylabel('Centriod Y')
% 
% xlim([0,ProcData.notes.trialDuration_sec])
% legend([p1,p2],'centriodX','centriodY')
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight
%% pupil motion
% ax5 = subplot(6,1,5);
% 
% p1 = plot((1:length(filtereddistanceTraveled))/ProcData.notes.dsFs,filtereddistanceTraveled,'color',colors('dark candy apple red'),'LineWidth',1);
% ylabel('Distance Travelled')
% 
% xlim([0,ProcData.notes.trialDuration_sec])
% legend(p1,'Distance Travelled')
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight
%% GRABNE
ax6 = subplot(6,1,4);

p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('\Delta F/F')
legend(p8,'GRAB NE')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
%% GRABACh
ax7 = subplot(6,1,5);

p8 = plot((1:length(filtACh_GFP))/ProcData.notes.dsFs,filtACh_GFP,'color',colors('vegas gold'),'LineWidth',1);
ylabel('\Delta F/F')
legend(p8,'GRAB ACh')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
axis tight
%% Blood Volume
ax8 = subplot(6,1,6);

p8 = plot((1:length(filtACh_CBV))/ProcData.notes.dsFs,filtACh_CBV,'color','r','LineWidth',1);
ylabel('\Delta F/F')
legend(p8,'CBV ACh')

xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'TickLength',[0,0])
set(gca,'box','off')
axis tight
%% Axes properties
linkaxes([ax1,ax2,ax3,ax6,ax7,ax8],'x')

ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');

ax1Pos(4) = ax1Pos(4)+0.05;
ax2Pos(4) = ax2Pos(4)+0.05;
ax3Pos(4) = ax3Pos(4)+0.05;
% ax4Pos(4) = ax4Pos(4)+0.05;
% ax5Pos(4) = ax5Pos(4)+0.05;
ax6Pos(4) = ax6Pos(4)+0.05;
ax7Pos(4) = ax7Pos(4)+0.05;
ax8Pos(4) = ax8Pos(4)+0.05;

ax1Pos(3) = ax1Pos(3)+0.05;
ax2Pos(3) = ax2Pos(3)+0.05;
ax3Pos(3) = ax3Pos(3)+0.05;
% ax4Pos(3) = ax4Pos(3)+0.05;
% ax5Pos(3) = ax5Pos(3)+0.05;
ax6Pos(3) = ax6Pos(3)+0.05;
ax7Pos(3) = ax7Pos(3)+0.05;
ax8Pos(3) = ax8Pos(3)+0.05;

ax1Pos(1) = ax1Pos(1)-0.05;
ax2Pos(1) = ax2Pos(1)-0.05;
ax3Pos(1) = ax3Pos(1)-0.05;
% ax4Pos(1) = ax4Pos(1)-0.05;
% ax5Pos(1) = ax5Pos(1)-0.05;
ax6Pos(1) = ax6Pos(1)-0.05;
ax7Pos(1) = ax7Pos(1)-0.05;
ax8Pos(1) = ax8Pos(1)-0.05;

ax1Pos(2) = ax1Pos(2)-0.05;
ax2Pos(2) = ax2Pos(2)-0.05;
ax3Pos(2) = ax3Pos(2)-0.05;
% ax4Pos(2) = ax4Pos(2)-0.05;
% ax5Pos(2) = ax5Pos(2)-0.05;
ax6Pos(2) = ax6Pos(2)-0.05;
ax7Pos(2) = ax7Pos(2)-0.05;
ax8Pos(2) = ax8Pos(2)-0.05;

set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
set(ax8,'position',ax8Pos);
    %% save the file to directory.
    if strcmp(saveFigs,'y') == true
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Single Trial Figures/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(figHandle,[dirpath animalID '_' fileID '_SingleTrialFigSleep']);
        saveas(figHandle,[dirpath animalID '_' fileID '_SingleTrialFigSleep'],'tiff')
    end

end
