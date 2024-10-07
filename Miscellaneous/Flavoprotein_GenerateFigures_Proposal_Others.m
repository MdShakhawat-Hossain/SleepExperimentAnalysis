% function [figHandle,ax1,ax2,ax3,ax4,ax5,ax7] = 
% function Flavoprotein_GenerateFigures_Proposal_Others()
% _(procDataFileID,saveFigs,MALabels,EMGArousalLabels,Notes,MicroLabels)
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
[z2,p2,k2] = butter(4,1/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
[z3,p3,k3] = butter(4,0.1/(ProcData.notes.dsFs/2),'low');
[sos3,g3] = zp2sos(z3,p3,k3);
resample_rate = 10;
% whisker angle
filteredWhiskerAngle = filtfilt(sos1,g1,resample(-ProcData.data.whiskerAngle,resample_rate,ProcData.notes.dsFs));
% binWhiskers = ProcData.data.binWhiskerAngle;
%pupil 
% filteredpupildiameter = filtfilt(sos1,g1,resample(ProcData.data.Pupil.zDiameter,resample_rate,ProcData.notes.dsFs));
% force sensor
filtForceSensor = filtfilt(sos1,g1,resample(ProcData.data.forceSensor,resample_rate,ProcData.notes.dsFs));
% binForce = ProcData.data.binForceSensor;
% emg
EMG = filtfilt(sos1,g1,resample(ProcData.data.EMG.emg,resample_rate,ProcData.notes.dsFs));% ProcData.data.EMG.emg; 
EMGSignal = resample(ProcData.data.EMG.emgSignal,resample_rate,ProcData.notes.dsFs);
% stimulations
LPadSol = ProcData.data.stimulations.LPadSol;
RPadSol = ProcData.data.stimulations.RPadSol;
AudSol = ProcData.data.stimulations.AudSol;
% NE data

NE_GFP = resample(ProcData.data.GFP.Z_NE,resample_rate,ProcData.notes.dsFs);
NE_GFP1 = medfilt1(NE_GFP,20);
filtNE_GFP = filtfilt(sos2,g2,NE_GFP1);
filtNE_GFP = sgolayfilt(filtNE_GFP,3,15);

% Ach data
Ach_GFP = resample(ProcData.data.GFP.Z_Ach,resample_rate,ProcData.notes.dsFs);
Ach_GFP1 = medfilt1(Ach_GFP,20);
filtAch_GFP = filtfilt(sos2,g2,Ach_GFP1);
filtAch_GFP = sgolayfilt(filtAch_GFP,3,15);
% Rhodamine data
NE_Rhodamine = resample(ProcData.data.Rhodamine.Z_NE,resample_rate,ProcData.notes.dsFs);
NE_Rhodamine1 = medfilt1(NE_Rhodamine,20);
filtNE_Rhodamine = filtfilt(sos2,g2,NE_Rhodamine1); 
filtNE_Rhodamine = sgolayfilt(filtNE_Rhodamine,3,15);

Ach_Rhodamine = resample(ProcData.data.Rhodamine.Z_Ach,resample_rate,ProcData.notes.dsFs);
Ach_Rhodamine1 = medfilt1(Ach_Rhodamine,20);
filtAch_Rhodamine = filtfilt(sos2,g2,Ach_Rhodamine1); 
filtAch_Rhodamine = sgolayfilt(filtAch_Rhodamine,3,15);
% cortical and hippocampal spectrograms
% specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
% load(specDataFile,'-mat');
% cortical_LHnormS = SpecData.cortical_LH.normS.*100;
% T = SpecData.cortical_LH.T;
% F = SpecData.cortical_LH.F;

EEG_LH = resample(ProcData.data.cortical_LH.corticalSignal,resample_rate,ProcData.notes.dsFs);
% remove some extra data
EEG_LH(1:resample_rate) = EEG_LH(resample_rate+1:resample_rate+resample_rate);
EEG_LH = medfilt1(EEG_LH,3);
% cortical and hippocampal spectrograms
% specDataFile = [animalID '_' fileID '_SpecDataA.mat'];
% load(specDataFile,'-mat');
% cortical_LHnormS = SpecData.cortical_LH.normS.*100;
% T = SpecData.cortical_LH.T;
% F = SpecData.cortical_LH.F;
% % remove some extra data
% EEG_LH(1:ProcData.notes.dsFs) = EEG_LH(ProcData.notes.dsFs+1:ProcData.notes.dsFs+ProcData.notes.dsFs);
% EEG_LH = medfilt1(EEG_LH,3);
%% Yvals for behavior Indices
% indecesMax = max(filteredpupildiameter);
indecesMax = max([filtAch_Rhodamine;filtAch_GFP]);
% whisking_Yvals = 1.10*max(indecesMax)*ones(size(binWhiskers));
% force_Yvals = 1.20*max(indecesMax)*ones(size(binForce));
% forceInds = binForce.*force_Yvals;
% whiskInds = binWhiskers.*whisking_Yvals;
LPad_Yvals = 1.30*max(indecesMax)*ones(size(LPadSol));
RPad_Yvals = 1.30*max(indecesMax)*ones(size(RPadSol));
Aud_Yvals = 1.30*max(indecesMax)*ones(size(AudSol));

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
%% Detect the microarousals based on EMG activation
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
ax1 = subplot(7,1,1);
% fileID2 = strrep(fileID,'_',' ');
p1 = plot((1:length(filtForceSensor))/resample_rate,filtForceSensor,'color',colors('north texas green'),'LineWidth',1);
% title([animalID ' FP behavioral characterization and sleep scoring ' fileID2])
ylabel('Force Sensor (V)')
xlim([0,ProcData.notes.trialDuration_sec])
yyaxis right
p2 = plot((1:length(EMG))/resample_rate,EMG,'color','k','LineWidth',1);
ylabel('EMG Power (V^2)')
xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Force Sensor','EMG Power')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% axis tight
% Whisker angle
ax2 = subplot(7,1,2);
p1 = plot((1:length(filteredWhiskerAngle))/resample_rate,-filteredWhiskerAngle,'color',colors('dark pink'),'LineWidth',1);
hold on
ylabel('Angle (deg)')
yyaxis right
p2 = plot((1:length(EMGSignal))/resample_rate,EMGSignal,'color','k','LineWidth',1);
ylabel('EMG (V)')

xlim([0,ProcData.notes.trialDuration_sec])
legend([p1,p2],'Whisker Angle','EMG')

set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
% % pupil and behavioral indeces
% ax3 = subplot(9,1,3);
% s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c'); hold on;
% s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
% s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
% p5 = plot((1:length(filteredpupildiameter))/ProcData.notes.dsFs,filteredpupildiameter,'color',colors('dark candy apple red'),'LineWidth',1);
% ylabel('Diameter Percentage')
% legend([p5,s3,s4,s5],'Pupil Diameter','LSol','RSol','AudSol')
% 
% 
% xlim([0,ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight
% Sleep Score data
ax4 = subplot(7,1,3);
s1 = scatter(SleepDummy',AwakeStage,'s','MarkerFaceColor','k','MarkerEdgeColor','k');
hold on
s2 = scatter(SleepDummy',NREMStage,'s','MarkerFaceColor',[0.3010 0.7450 0.9330],'MarkerEdgeColor',[0.3010 0.7450 0.9330]);
s3 = scatter(SleepDummy',REMStage,'s','MarkerFaceColor','r','MarkerEdgeColor','r');
legend([s1,s2,s3],'Awake','NREM','REM')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
xlim([0,(ProcData.notes.trialDuration_sec)])

% GRABNE and GRAB ACh
% ax5 = subplot(9,1,[3:5]);
% 
% plot((1:length(filtNE_GFP))/(resample_rate),filtNE_GFP,'color',[0.4660 0.6740 0.1880],'LineWidth',1);
% % ylim([-3.5 4])
% ylabel('\Delta F/F Percentage')
% hold on
% % yyaxis right
% % plot((1:length(filtAch_GFP))/(resample_rate),filtAch_GFP,'color',[0 0.4470 0.7410],'LineWidth',1);
% % ylim([-3.5 4])
% legend('Flavoprotein')
% % ylabel('\Delta F/F Percentage')
% xlim([0,(ProcData.notes.trialDuration_sec)])

% ax6 = subplot(9,1,[7:9]);
% % % yyaxis right
% % plot((1:length(filtAch_Rhodamine))/(resample_rate),filtAch_Rhodamine,'color',[0.1350 0.0780 0.6840],'LineWidth',1);
% % ylim([-3.5 4])
% % ylabel('\Delta F/F Percentage')
% % legend('CBV')
% hold on
% xlim([0,(ProcData.notes.trialDuration_sec)]) 
% % yyaxis right
% 
% plot((1:length(filtNE_Rhodamine))/(resample_rate),filtNE_Rhodamine,'color',[ 0.6350 0.0780 0.1840],'LineWidth',1);
% % ylim([-3.5 4])
% ylabel('\Delta F/F Percentage')
% legend('CBV-LH','CBV-RH')
% 
% % set(gca,'TickLength',[0,0])
% % set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight

%ECoG
% ax7 = subplot(7,1,6:7);
% plot((1:length(EEG_LH))/(resample_rate*60),EEG_LH,'color',colors('deep carrot orange'),'LineWidth',1);
% % xlim([0,ProcData.notes.trialDuration_sec/60])
% ylabel('ECoG (V)')
% legend('ECoG')
% % xlim([0,ProcData.notes.trialDuration_sec])
% set(gca,'TickLength',[0,0])
% % set(gca,'Xticklabel',[])
% set(gca,'box','off')
% axis tight

% Left cortical electrode spectrogram
% ax8 = subplot(8,1,8);
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
% linkaxes([ax4,ax5,ax6,ax7],'x') %,ax8 

% ax1Pos = get(ax1,'position');
% ax2Pos = get(ax2,'position');
% ax3Pos = get(ax3,'position');
% ax4Pos = get(ax4,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% ax7Pos = get(ax7,'position');
% ax8Pos = get(ax8,'position');

% ax8Pos(3) = ax3Pos(3);

% ax1Pos(4) = ax1Pos(4)+0.08;
% ax2Pos(4) = ax2Pos(4)+0.08;
% ax3Pos(4) = ax3Pos(4)+0.08;
% ax4Pos(4) = ax4Pos(4)+0.07;
% ax5Pos(4) = ax5Pos(4)+0.07;
% ax6Pos(4) = ax6Pos(4)+0.07;
% ax7Pos(4) = ax7Pos(4)+0.07;
% ax8Pos(4) = ax8Pos(4)+0.05;

% ax8Pos(3) = ax3Pos(3);

% ax1Pos(3) = ax1Pos(3)+0.05;
% ax2Pos(3) = ax2Pos(3)+0.05;
% ax3Pos(3) = ax3Pos(3)+0.05;
% ax4Pos(3) = ax4Pos(3)+0.05;
% ax5Pos(3) = ax5Pos(3)+0.05;
% ax6Pos(3) = ax6Pos(3)+0.05;
% ax7Pos(3) = ax7Pos(3)+0.05;
% ax8Pos(3) = ax8Pos(3)+0.05;

% ax1Pos(1) = ax1Pos(1)-0.05;
% ax2Pos(1) = ax2Pos(1)-0.05;
% ax3Pos(1) = ax3Pos(1)-0.05;
% ax4Pos(1) = ax4Pos(1)-0.05;
% ax5Pos(1) = ax5Pos(1)-0.05;
% ax6Pos(1) = ax6Pos(1)-0.05;
% ax7Pos(1) = ax5Pos(1)-0.05;
% ax8Pos(1) = ax6Pos(1)-0.05;

% ax1Pos(2) = ax1Pos(2)-0.05;
% ax2Pos(2) = ax2Pos(2)-0.05;
% ax3Pos(2) = ax3Pos(2)-0.05;
% ax4Pos(2) = ax4Pos(2)-0.05;
% ax5Pos(2) = ax5Pos(2)-0.05;
% ax6Pos(2) = ax6Pos(2)-0.05;
% ax7Pos(2) = ax7Pos(2)-0.05;
% ax8Pos(2) = ax8Pos(2)-0.05;

% set(ax1,'position',ax1Pos);
% set(ax2,'position',ax2Pos);
% set(ax3,'position',ax3Pos);
% set(ax4,'position',ax4Pos);
% % set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
% set(ax7,'position',ax7Pos);
% set(ax8,'position',ax8Pos);
%% save the file to directory.
% if strcmp(saveFigs,'y') == true
%     [pathstr,~,~] = fileparts(cd);
%     dirpath = [pathstr '/Figures/Single Trial Figures/'];
%     if ~exist(dirpath,'dir')
%         mkdir(dirpath);
%     end
%     savefig(figHandle,[dirpath animalID '_' fileID '_MicroArousalSleep']);
% %     saveas(figHandle,[dirpath animalID '_' fileID '_MicroArousalSleep'],'tiff')
% 
% end

% end
