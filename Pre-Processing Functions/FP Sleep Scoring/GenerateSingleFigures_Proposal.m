% function [figHandle,ax2,ax4,ax7] = GenerateSingleFigures_Proposal(procDataFileID,saveFigs,MALabels,EMGArousalLabels,Notes,MicroLabels)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
% Adopted from Kevin Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Create a summary figure for a single n minute FP trial
%________________________________________________________________________________________________________________________

% load file and gather information
% load(procDataFileID)
% [animalID,fileDate,fileID] = GetFileInfo_FP(procDataFileID);
% strDay = ConvertDate_FP(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.25/(ProcData.notes.dsFs/2),'low');
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
% binForce = ProcData.data.binForceSensor;
% emg
EMG = ProcData.data.EMG.emg;
% U_Time = 2200;
StartTime = 1;
EndTime = 3120;

EMGSignal = ProcData.data.EMG.emgSignal;%((StartTime*ProcData.notes.dsFs) + 1 : EndTime*ProcData.notes.dsFs);
% EMGSignal = sgolayfilt(EMGSignal,3,17);
% stimulations
% LPadSol = ProcData.data.stimulations.LPadSol;
% RPadSol = ProcData.data.stimulations.RPadSol;
% AudSol = ProcData.data.stimulations.AudSol;
% NE data
NE_GFP = ProcData.data.GFP.Z_NE;%((StartTime*ProcData.notes.dsFs) + 1 : EndTime*ProcData.notes.dsFs);
filtNE_GFP = filtfilt(sos2,g2,NE_GFP);
% filtNE_GFP = sgolayfilt(filtNE_GFP,3,17);

% Ach data
Ach_GFP = ProcData.data.GFP.Z_Ach;%((StartTime*ProcData.notes.dsFs) + 1 : EndTime*ProcData.notes.dsFs);
filtAch_GFP = filtfilt(sos2,g2,Ach_GFP);
% filtAch_GFP = sgolayfilt(filtAch_GFP,3,17);
% Rhodamine data
NE_Rhodamine = ProcData.data.Rhodamine.Z_NE;%((StartTime*ProcData.notes.dsFs) + 1 : EndTime*ProcData.notes.dsFs);
filtNE_Rhodamine = filtfilt(sos2,g2,NE_Rhodamine);
% filtNE_Rhodamine = sgolayfilt(filtNE_Rhodamine,3,17);

% Rhodamine data
Ach_Rhodamine = ProcData.data.Rhodamine.Z_Ach;%((StartTime*ProcData.notes.dsFs) + 1 : EndTime*ProcData.notes.dsFs);
filtAch_Rhodamine = filtfilt(sos2,g2,Ach_Rhodamine);

EEG_LH = ProcData.data.cortical_LH.corticalSignal;
% remove some extra data
EEG_LH(1:ProcData.notes.dsFs) = EEG_LH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
EEG_LH = medfilt1(EEG_LH,3);
% filtAch_Rhodamine = sgolayfilt(filtAch_Rhodamine,3,17);
% Yvals for behavior Indices
% indecesMax = max(filteredpupildiameter);
% % indecesMax = max([filtAch_Rhodamine;filtAch_GFP]);
% whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
% force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
% forceInds = binForce.*force_Yvals;
% whiskInds = binWhiskers.*whisking_Yvals;
% LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
% RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
% Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));

% set force indeces
% for x = 1:length(forceInds)
%     if forceInds(1,x) == 0
%         forceInds(1,x) = NaN;
%     end
% end
% % set whisk indeces
% for x = 1:length(whiskInds)
%     if whiskInds(1,x) == 0
%         whiskInds(1,x) = NaN;
%     end
% end
%% get the sleep score
% AwakeStage = zeros(3120,1);
% NREMStage = zeros(3120,1);
% REMStage = zeros(3120,1);
% 
% % AwakeStage(MALabels==1) = 1;
% % NREMStage(MALabels==2) = 1;
% % REMStage(MALabels==3) = 1;
% 
% AwakeStage (AwakeStage==0) = nan ;
% NREMStage (NREMStage==0) = nan ;
% REMStage (REMStage==0) = nan ;
% % indecesMax = max(max(filtNE_GFP),max(filtAch_GFP));
% indecesMax = max(EMGSignal(60000:80000));
% 
% AwakeYvals = 1.3 * indecesMax * AwakeStage;%(StartTime : EndTime); 
% NREMYvals = 1.3 * indecesMax * NREMStage;%(StartTime : EndTime); 
% REMYvals = 1.3 * indecesMax * REMStage;%(StartTime : EndTime);
% 
% SleepDummy = 1:1:length(REMYvals);%ProcData.notes.trialDuration_sec;


% indecesMax = max(filtAch_GFP);
% MicroLabels(MicroLabels == 0) = nan;
% MicroLabelYvals = 1.20 * indecesMax * MicroLabels; 
% %% Detect the microarousals based on EMG activation
% EMGD = abs(ProcData.data.EMG.emgSignal);
% EMGThreshold = prctile(EMGD,Notes.EMGThresholdPerc)*ones(size(EMGD));
% emgBinarized = EMGD;
% emgBinarized(emgBinarized > EMGThreshold) = 1;
% emgBinarized(emgBinarized < EMGThreshold) = 0;
% emgBinarized(emgBinarized ==0 ) = nan;
% indecesMax = max(filteredpupildiameter);
% emg_Yvals = 1.30*max(indecesMax)*ones(size(emgBinarized));
% emg_Idx = emg_Yvals.*emgBinarized;
% 
% indecesMax = max(filtNE_Rhodamine);
% EMGArousalLabels(EMGArousalLabels==0) = nan;
% EMGArousalLabels = 1.30*max(indecesMax)*EMGArousalLabels;

%%
%% get the sleep score
% TableSize = length(ProcData.sleep.logicals.Manual.awakeLogical);
            
% DataSize = TableSize*5;
% AwakeStage = zeros(DataSize,1);
% NREMStage = zeros(DataSize,1);
% REMStage = zeros(DataSize,1);

% OGLabels_awake = ProcData.sleep.logicals.Manual.awakeLogical;
% OGLabels_nrem = ProcData.sleep.logicals.Manual.nremLogical;
% OGLabels_rem = ProcData.sleep.logicals.Manual.remLogical;
AwakeStage = double(ProcData.sleep.logicals.Manual.awakeLogical);
NREMStage = double(ProcData.sleep.logicals.Manual.nremLogical);
REMStage = double(ProcData.sleep.logicals.Manual.remLogical);
% for SL = 1:1:TableSize
%     if OGLabels_awake(SL) == 1
%         for ML = 1:1:5
%             AwakeStage(((SL-1)*5)+ML) = 1;
%         end
%     elseif OGLabels_nrem(SL) == 1
%         for ML = 1:1:5
%             NREMStage(((SL-1)*5)+ML) = 1;
%         end
%     elseif OGLabels_rem(SL) == 1
%         for ML = 1:1:5
%             REMStage(((SL-1)*5)+ML) = 1;
%         end
%     end
% end

AwakeStage (AwakeStage==0) = nan ;
NREMStage (NREMStage==0) = nan ;
REMStage (REMStage==0) = nan ;

% indecesMax = max(filtNE_GFP);
% AwakeYvals = 1.20 * indecesMax * AwakeStage; 
% NREMYvals = 1.20 * indecesMax * NREMStage; 
% REMYvals = 1.20 * indecesMax * REMStage; 
SleepDummy = 1:5:ProcData.notes.trialDuration_sec;
%% Figure
figHandle = figure;
% force sensor and EMG
ax2 = subplot(6,1,1);
s1 = scatter(SleepDummy',AwakeYvals,'s','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',0.025);
hold on
s2 = scatter(SleepDummy',NREMYvals,'s','MarkerFaceColor','c','MarkerEdgeColor','c','LineWidth',0.025);
s3 = scatter(SleepDummy',REMYvals,'s','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',0.025);
plot((1:length(EMGSignal))/ProcData.notes.dsFs,EMGSignal,'color',[0.9290 0.6940 0.1250],'LineWidth',1);
ylabel('EMG (V)')
xlim([0,ProcData.notes.trialDuration_sec])
ylim([2e-5 2e5]);
ax1.YAxis(1).Limits = [-2e-5 2e-5];
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% GRABNE and Sleep Score data
ax4 = subplot(6,1,[2,3]);
hold on
p8 = plot((1:length(filtNE_GFP))/ProcData.notes.dsFs,filtNE_GFP,'Color',[0.4660 0.6740 0.1880],'LineWidth',1);
ylabel('\Delta F/F (Z)')

% yyaxis right
p18 = plot((1:length(filtAch_GFP))/ProcData.notes.dsFs,filtAch_GFP,'Color', [0 0.4470 0.7410], 'LineWidth',1);
ylim([-4 4])
% ylabel('\Delta F/F (Z)')
% ax4.YAxis(1).Limits = [-4 4];
ax4.YAxis(1).Color = 'k';% [0.4660 0.6740 0.1880];
% ax4.YAxis(2).Color = [0 0.4470 0.7410];
legend([p8,p18],'NE','ACh')
xlim([0,ProcData.notes.trialDuration_sec])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% ax4.YAxis(2).Limits = [-3 3];

axis tight

% Rhodamine and MicroArousals Score data
% ax7 = subplot(6,1,[4,5]);
% ylim([-3 3])
plot((1:length(filtNE_Rhodamine))/ProcData.notes.dsFs,filtNE_Rhodamine,'color',[0.6350 0.0780 0.1840],'LineWidth',1);

ylabel('\Delta F/F (Z)')
 
hold on
% yyaxis right
% plot((1:length(filtAch_Rhodamine))/ProcData.notes.dsFs,filtAch_Rhodamine,'color',[0.8500 0.3250 0.0980],'LineWidth',1);
ylim([-4 4])
% ylabel('\Delta F/F CBV (Z)')
% legend('RH CBV','LH CBV')
legend('NE','ACh','CBV')
xlim([0,ProcData.notes.trialDuration_sec])
% ax7.YAxis(1).Limits = [-4 4];
ax7.YAxis(1).Color = 'k';% [0.6350 0.0780 0.1840];
% ax7.YAxis(2).Color = [0.8500 0.3250 0.0980];
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])

set(gca,'box','off')
axis tight
% ax7.YAxis(2).Limits = [-3 3];


% Axes properties
% linkaxes([ax2,ax4,ax7],'x')
linkaxes([ax2,ax4],'x')
% ax4.YAxis(2).Limits = [-3 3];

% ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
% ax6Pos(3:4) = ax2Pos(3:4);
% ax6Pos(3) = ax2Pos(3);

% ax1Pos(4) = ax1Pos(4)+0.03;
ax2Pos(4) = ax2Pos(4)+0.03;
% ax3Pos(4) = ax3Pos(4)+0.03;
ax4Pos(4) = ax4Pos(4)+0.03;
% ax5Pos(4) = ax5Pos(4)+0.03;
ax7Pos(4) = ax7Pos(4)+0.03;
% ax6Pos(4) = ax6Pos(4)+0.03;

% ax1Pos(3) = ax1Pos(3)+0.03;
ax2Pos(3) = ax2Pos(3)+0.03;
% ax3Pos(3) = ax3Pos(3)+0.03;
ax4Pos(3) = ax4Pos(3)+0.03;
% ax5Pos(3) = ax5Pos(3)+0.03;
% ax6Pos(3) = ax6Pos(3)+0.03;
ax7Pos(3) = ax7Pos(3)+0.03;

% ax1Pos(1) = ax1Pos(1)-0.03;
ax2Pos(1) = ax2Pos(1)-0.03;
% ax3Pos(1) = ax3Pos(1)-0.03;
ax4Pos(1) = ax4Pos(1)-0.03;
% ax5Pos(1) = ax5Pos(1)-0.03;
% ax6Pos(1) = ax6Pos(1)-0.03;
ax7Pos(1) = ax7Pos(1)-0.03;

% ax1Pos(2) = ax1Pos(2)-0.03;
ax2Pos(2) = ax2Pos(2)-0.03;
% ax3Pos(2) = ax3Pos(2)-0.03;
ax4Pos(2) = ax4Pos(2)-0.03;
% ax5Pos(2) = ax5Pos(2)-0.03;
% ax6Pos(2) = ax6Pos(2)-0.03;
ax7Pos(2) = ax7Pos(2)-0.03;

% set(ax1,'position',ax1Pos);
set(ax2,'position',ax2Pos);
% set(ax3,'position',ax3Pos);
set(ax4,'position',ax4Pos);
% set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
% xlim([1350 2050])
%% save the file to directory.
% if strcmp(saveFigs,'y') == true
%     [pathstr,~,~] = fileparts(cd);
%     dirpath = [pathstr '/Figures/Single Trial Figures/'];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(figHandle,[dirpath animalID '_' fileID '_Proposal']);
%     saveas(figHandle,[dirpath animalID '_' fileID '_Proposal'],'tiff')
% 
% end

% end
