function [AnalysisResults] = Fig4_S1_GCaMP7s(rootFolder,saveFigs,delim,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel 4-S1 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

%% set-up and process data
FPanimalIDs = {'T282','T285'};%{'T281','T282','T284','T285'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
%% mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,aa};
    for bb = 1:length(transitions)
        transition = transitions{1,bb};
        data.(animalID).(transition) = [];
        uniqueFileDates = AnalysisResults.(animalID).Transitions.(transition).fileDates;
        for cc = 1:length(uniqueFileDates)
            data.(animalID).(transition).(uniqueFileDates{cc,1}).LH_GCaMP7s = [];
            data.(animalID).(transition).(uniqueFileDates{cc,1}).RH_GCaMP7s = [];
        end
        for dd = 1:length(AnalysisResults.(animalID).Transitions.(transition).indFileDate)
            strDay = AnalysisResults.(animalID).Transitions.(transition).indFileDate{dd,1};
            data.(animalID).(transition).(strDay).LH_GCaMP7s = cat(1,data.(animalID).(transition).(strDay).LH_GCaMP7s,AnalysisResults.(animalID).Transitions.(transition).LH_GCaMP7s);%(dd,:));
            data.(animalID).(transition).(strDay).RH_GCaMP7s = cat(1,data.(animalID).(transition).(strDay).RH_GCaMP7s,AnalysisResults.(animalID).Transitions.(transition).RH_GCaMP7s);%(dd,:));
        end
    end
end
% put together L/R for each day
for ee = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,ee};
    for ff = 1:length(transitions)
        transition = transitions{1,ff};
        uniqueFileDates = fieldnames(data.(animalID).(transition));
        for gg = 1:length(uniqueFileDates)
            strDay = uniqueFileDates{gg,1};
            procData.(animalID).(transition).(strDay).GCaMP7s = cat(1,data.(animalID).(transition).(strDay).LH_GCaMP7s,data.(animalID).(transition).(strDay).RH_GCaMP7s);
        end
    end
end
% take average for each animal's behavioral transition per day
for ee = 1:length(FPanimalIDs)
    animalID = FPanimalIDs{1,ee};
    for ff = 1:length(transitions)
        transition = transitions{1,ff};
        uniqueFileDates = fieldnames(data.(animalID).(transition));
        for gg = 1:length(uniqueFileDates)
            strDay = uniqueFileDates{gg,1};
            if isempty(procData.(animalID).(transition).(strDay).GCaMP7s) == false
                finData.(transition){gg,1}(ee,:) = mean(procData.(animalID).(transition).(strDay).GCaMP7s,1);
            else
                finData.(transition){gg,1}(ee,:) = NaN(1,1800);
            end
        end
    end
end
% patch Animal T110 who only had 5 days
for hh = 1:length(transitions)
    transition = transitions{1,hh};
    for kk = 2:3
        for nnn = 3:4
        finData.(transition){kk,1}(nnn,:) = NaN(1,1800);
        end
    end
end
% take average across animals
for ii = 1:length(transitions)
    transition = transitions{1,ii};
    for jj = 1:3
        meanData.(transition){jj,1} = nanmean(finData.(transition){jj,1},1);
        % check nans
        ll = 1;
        for kk = 1:size(finData.(transition){jj,1},1)
            if isnan(finData.(transition){jj,1}(kk,1)) == false
                nCount.(transition){jj,1} = ll; %#ok<STRNU>
                ll = ll + 1;
            end
        end
    end
end
T1 = -30:(1/30):30;
T1 = T1(1:end - 1);
%% Fig. 4-S1
summaryFigure = figure('Name','Fig4-S1');
sgtitle('﻿Transitional changes in GCaMP7s across days')
%% [4-S1a] Awake to NREM
ax1 = subplot(2,2,1);
for kk = 1:3
    p(kk) = plot(T1,meanData.AWAKEtoNREM{kk,1},'-','LineWidth',2); %#ok<*AGROW>
    hold on
end
title('Awake to NREM transition')
xlabel('Time (s)')
ylabel('Zscored \Delta GCaMP7s')
legend([p(1),p(2),p(3)],'Day 1','Day 2','Day 3')
% xlim([-30,30])
% ylim([-5,45])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [4-S1b] NREM to Awake
ax2 = subplot(2,2,2);
for kk = 1:3
    plot(T1,meanData.NREMtoAWAKE{kk,1},'-','LineWidth',2);
    hold on
end
title('NREM to Awake transition')
xlabel('Time (s)')
ylabel('Zscored \Delta GCaMP7s')
% xlim([-30,30])
% ylim([-5,45])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [4-S1c] NREM to REM
ax3 = subplot(2,2,3);
for kk = 1:3
    plot(T1,meanData.NREMtoREM{kk,1},'-','LineWidth',2);
    hold on
end
title('NREM to REM transition')
xlabel('Time (s)')
ylabel('Zscored \Delta GCaMP7s')
% xlim([-30,30])
% ylim([35,80])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [4-S1d] REM to Awake
ax4 = subplot(2,2,4);
for kk = 1:3
    plot(T1,meanData.REMtoAWAKE{kk,1},'-','LineWidth',2);
    hold on
end
title('REM to Awake transition')
xlabel('Time (s)')
ylabel('Zscored \Delta GCaMP7s')
% xlim([-30,30])
% ylim([0,100])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig4-S1-GCaMP7s']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig4-S1-GCaMP7s'])
end

end
