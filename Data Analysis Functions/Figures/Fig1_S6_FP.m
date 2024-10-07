function [AnalysisResults] = Fig1_S6_FP(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 1-S6 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
% colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
% colorNREM = [(191/256),(0/256),(255/256)];
% colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
% information and data for first example
if isfield(AnalysisResults,'ExampleTrials') == false
    AnalysisResults.ExampleTrials = [];
end
if isfield(AnalysisResults.ExampleTrials,'T123B') == true
    dsFs = AnalysisResults.ExampleTrials.T123B.dsFs;
    filtEMG = AnalysisResults.ExampleTrials.T123B.filtEMG;
    filtForceSensor = AnalysisResults.ExampleTrials.T123B.filtForceSensor;
    filtWhiskerAngle = AnalysisResults.ExampleTrials.T123B.filtWhiskerAngle;
    heartRate = AnalysisResults.ExampleTrials.T123B.heartRate;
    filtLH_HbT = AnalysisResults.ExampleTrials.T123B.filtLH_HbT;
    filtRH_HbT = AnalysisResults.ExampleTrials.T123B.filtRH_HbT;
    T = AnalysisResults.ExampleTrials.T123B.T;
    F = AnalysisResults.ExampleTrials.T123B.F;
    cortical_LHnormS = AnalysisResults.ExampleTrials.T123B.cortical_LHnormS;
    cortical_RHnormS = AnalysisResults.ExampleTrials.T123B.cortical_RHnormS;
    hippocampusNormS = AnalysisResults.ExampleTrials.T123B.hippocampusNormS;
else
    animalID = 'T123';
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
    exampleProcDataFileID = 'T123_200224_16_27_59_ProcData.mat';
    load(exampleProcDataFileID,'-mat')
    exampleSpecDataFileID = 'T123_200224_16_27_59_SpecDataA.mat';
    load(exampleSpecDataFileID,'-mat')
    exampleBaselineFileID = 'T123_RestingBaselines.mat';
    load(exampleBaselineFileID,'-mat')
    [~,fileDate,~] = GetFileInfo_IOS(exampleProcDataFileID);
    strDay = ConvertDate_IOS(fileDate);
    dsFs = ProcData.notes.dsFs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z1,p1_A,k1] = butter(4,10/(dsFs/2),'low');
    [sos1,g1] = zp2sos(z1,p1_A,k1);
    [z2,p2_A,k2] = butter(4,0.5/(dsFs/2),'low');
    [sos2,g2] = zp2sos(z2,p2_A,k2);
    % whisker angle
    filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
    % force sensor
    filtForceSensor = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
    % EMG
    EMG = ProcData.data.EMG.emg;
    normEMG = EMG - RestingBaselines.manualSelection.EMG.emg.(strDay);
    filtEMG = filtfilt(sos1,g1,normEMG);
    % heart rate
    heartRate = ProcData.data.heartRate;
    % HbT data
    LH_HbT = ProcData.data.CBV_HbT.adjLH;
    filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
    RH_HbT = ProcData.data.CBV_HbT.adjRH;
    filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
    % cortical and hippocampal spectrograms
    cortical_LHnormS = SpecData.cortical_LH.normS.*100;
    cortical_RHnormS = SpecData.cortical_RH.normS.*100;
    hippocampusNormS = SpecData.hippocampus.normS.*100;
    T = SpecData.cortical_LH.T;
    F = SpecData.cortical_LH.F;
    % update analysis structure
    AnalysisResults.ExampleTrials.T123B.dsFs = dsFs;
    AnalysisResults.ExampleTrials.T123B.filtEMG = filtEMG;
    AnalysisResults.ExampleTrials.T123B.filtForceSensor = filtForceSensor;
    AnalysisResults.ExampleTrials.T123B.filtWhiskerAngle = filtWhiskerAngle;
    AnalysisResults.ExampleTrials.T123B.heartRate = heartRate;
    AnalysisResults.ExampleTrials.T123B.filtLH_HbT = filtLH_HbT;
    AnalysisResults.ExampleTrials.T123B.filtRH_HbT = filtRH_HbT;
    AnalysisResults.ExampleTrials.T123B.T = T;
    AnalysisResults.ExampleTrials.T123B.F = F;
    AnalysisResults.ExampleTrials.T123B.cortical_LHnormS = cortical_LHnormS;
    AnalysisResults.ExampleTrials.T123B.cortical_RHnormS = cortical_RHnormS;
    AnalysisResults.ExampleTrials.T123B.hippocampusNormS = hippocampusNormS;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. 1-S6
summaryFigure = figure('Name','Fig1-S6 (a-f)');
sgtitle('Figure 1-S6 - Turner et al. 2020')
%% EMG and force sensor
ax1 = subplot(7,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors('rich black'),'LineWidth',0.5);
ylabel({'EMG','power (a.u.)'})  
ylim([-2.5,3])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
%% whisker angle and heart rate
ax2 = subplot(7,1,2);
p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4 = plot((1:length(heartRate)),heartRate,'color',colors('deep carrot orange'),'LineWidth',0.5);
ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p3,p4],'Whisker angle','Heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors('rich black');
ax2.YAxis(2).Color = colors('deep carrot orange');
%% HbT and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6 = plot((1:length(filtRH_HbT))/dsFs,filtRH_HbT,'color',colors('sapphire'),'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/dsFs,filtLH_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
x1 = xline(0,'color',colorRfcNREM,'LineWidth',2);
ylabel('\Delta[HbT] (\muM)')
legend([p5,p6,x1],'Left hem','Right hem','NREM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
axis tight
xlim([0,600])
ax34.TickLength = [0.01,0.01];
%% left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ax5.TickLength = [0.01,0.01];
%% right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc(T,F,cortical_RHnormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ax6.TickLength = [0.01,0.01];
%% hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
axis xy
semilog_imagesc(T,F,hippocampusNormS,'y')
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([0,600])
ax7.TickLength = [0.01,0.01];
%% axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig1-S6']);
    % remove surface subplots because they take forever to render
    cla(ax5);
    set(ax5,'YLim',[1,99]);
    cla(ax6);
    set(ax6,'YLim',[1,99]);
    cla(ax7);
    set(ax7,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S6'])
    close(summaryFigure)
    %% subplot figures
    summaryFigure_imgs = figure;
    % example 3 LH cortical LFP
    subplot(3,1,1);
    semilog_imagesc(T,F,cortical_LHnormS,'y')
    caxis([-100,200])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([0,600])
    % example 3 RH cortical LFP
    subplot(3,1,2);
    semilog_imagesc(T,F,cortical_RHnormS,'y')
    caxis([-100,200])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([0,600])
    % example 3 hippocampal LFP
    subplot(3,1,3);
    semilog_imagesc(T,F,hippocampusNormS,'y')
    caxis([-100,200])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([0,600])
    print('-painters','-dtiffn',[dirpath 'Fig1-S6_SpecImages'])
    close(summaryFigure_imgs)
    %% Fig. 1-S6
    figure('Name','Fig1-S6 (a-f)');
    sgtitle('Figure 1-S6 - Turner et al. 2020')
    %% EMG and force sensor
    ax1 = subplot(7,1,1);
    p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors('rich black'),'LineWidth',0.5);
    ylabel({'EMG','power (a.u.)'})  
    ylim([-2.5,3])
    yyaxis right
    p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
    ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'EMG','Pressure')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xlim([0,600])
    ylim([-0.1,2.5])
    ax1.TickLength = [0.01,0.01];
    ax1.YAxis(1).Color = colors('rich black');
    ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
    %% whisker angle and heart rate
    ax2 = subplot(7,1,2);
    p3 = plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors('rich black'),'LineWidth',0.5);
    ylabel({'Whisker','angle (deg)'})
    ylim([-20,60])
    yyaxis right
    p4 = plot((1:length(heartRate)),heartRate,'color',colors('deep carrot orange'),'LineWidth',0.5);
    ylabel({'Heart rate','Freq (Hz)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p3,p4],'Whisker angle','Heart rate')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xlim([0,600])
    ylim([5,15])
    ax2.TickLength = [0.01,0.01];
    ax2.YAxis(1).Color = colors('rich black');
    ax2.YAxis(2).Color = colors('deep carrot orange');
    %% HbT and behavioral indeces
    ax34 =subplot(7,1,[3,4]);
    p6 = plot((1:length(filtRH_HbT))/dsFs,filtRH_HbT,'color',colors('sapphire'),'LineWidth',1);
    hold on
    p5 = plot((1:length(filtLH_HbT))/dsFs,filtLH_HbT,'color',colors('dark candy apple red'),'LineWidth',1);
    x1 = xline(0,'color',colorRfcNREM,'LineWidth',2);
    ylabel('\Delta[HbT] (\muM)')
    legend([p5,p6,x1],'Left hem','Right hem','NREM')
    set(gca,'TickLength',[0,0])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    axis tight
    xlim([0,600])
    ax34.TickLength = [0.01,0.01];
    %% left cortical electrode spectrogram
    ax5 = subplot(7,1,5);
    semilog_imagesc(T,F,cortical_LHnormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,200])
    ylabel({'LH cortical LFP','Freq (Hz)'})
    set(gca,'Yticklabel','10^1')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xlim([0,600])
    ax5.TickLength = [0.01,0.01];
    %% right cortical electrode spectrogram
    ax6 = subplot(7,1,6);
    semilog_imagesc(T,F,cortical_RHnormS,'y')
    axis xy
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,200])
    ylabel({'RH cortical LFP','Freq (Hz)'})
    set(gca,'Yticklabel','10^1')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xlim([0,600])
    ax6.TickLength = [0.01,0.01];
    %% hippocampal electrode spectrogram
    ax7 = subplot(7,1,7);
    semilog_imagesc(T,F,hippocampusNormS,'y')
    axis xy
    c7 = colorbar;
    ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,200])
    xlabel('Time (min)')
    ylabel({'Hippocampal LFP','Freq (Hz)'})
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    xlim([0,600])
    ax7.TickLength = [0.01,0.01];
    %% axes properties
    ax1Pos = get(ax1,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax7Pos = get(ax7,'position');
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    ax7Pos(3:4) = ax1Pos(3:4);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
    set(ax7,'position',ax7Pos);
end

end
