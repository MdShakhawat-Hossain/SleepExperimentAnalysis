function [AnalysisResults] = PlotpowerSpectrum_GRABNE_AAA(rootFolder,saveFigs,delim,AnalysisResults,animalIDs)


behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
dataTypes = {'Ach_Rhodamine','NE_Rhodamine','Ach_GFP','NE_GFP'};
% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.PowerSpec.(behavField).(dataType).S = [];
        data.PowerSpec.(behavField).(dataType).f = [];
    end
end
% concatenate power spectra during different arousal states for each animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).S) == false
                data.PowerSpec.(behavField).(dataType).S = cat(2,data.PowerSpec.(behavField).(dataType).S,AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).S);
                data.PowerSpec.(behavField).(dataType).f = cat(1,data.PowerSpec.(behavField).(dataType).f,AnalysisResults.(animalID).PowerSpectrum.(behavField).(dataType).f);
            end
        end
    end
end
% mean and standard error for arousal state power
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        data.PowerSpec.(behavField).(dataType).meanS = mean(data.PowerSpec.(behavField).(dataType).S,2);
        data.PowerSpec.(behavField).(dataType).semS = std(data.PowerSpec.(behavField).(dataType).S,0,2)./sqrt(size(data.PowerSpec.(behavField).(dataType).S,2));
        data.PowerSpec.(behavField).(dataType).meanf = mean(data.PowerSpec.(behavField).(dataType).f,1);
    end
end
%% power spectrum during arousal states
    powerSpec_Fig = figure('Name','PowerSpectrum_Figure');
    sgtitle('Power Spectrum in different arousal stages')
       SubplotNeed = ceil(length(dataTypes)/3);
       splotNo = 1;
for PPS = 1:length(dataTypes)
    dataType = dataTypes{1,PPS};

    subplot(3,SubplotNeed,splotNo);
    splotNo = splotNo + 1;

    % hold on
    
    L4 = loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS,'color',[0.9290 0.6940 0.1250],'LineWidth',2);
    hold on
    loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS + data.PowerSpec.Awake.(dataType).semS,'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
    loglog(data.PowerSpec.Awake.(dataType).meanf,data.PowerSpec.Awake.(dataType).meanS - data.PowerSpec.Awake.(dataType).semS,'color',[0.4940 0.1840 0.5560],'LineWidth',0.5);
    
    L5 = loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS,'color',[1 0 1],'LineWidth',2);
    loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS + data.PowerSpec.Asleep.(dataType).semS,'color',[1 0 1],'LineWidth',0.5);
    loglog(data.PowerSpec.Asleep.(dataType).meanf,data.PowerSpec.Asleep.(dataType).meanS - data.PowerSpec.Asleep.(dataType).semS,'color',[1 0 1],'LineWidth',0.5);
    
    L6 = loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS,'color',[0.4660 0.6740 0.1880],'LineWidth',2);
    loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS + data.PowerSpec.All.(dataType).semS,'color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
    loglog(data.PowerSpec.All.(dataType).meanf,data.PowerSpec.All.(dataType).meanS - data.PowerSpec.All.(dataType).semS,'color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
    xline(1/10)
    xline(1/30)
    xline(0.35)
    xline(0.02)
    title([dataType])
    ylabel('power (a.u.)')
    xlabel('Freq (Hz)')
    if splotNo == length(dataTypes) + 1
    legend([ L4, L5, L6],'Alert','Asleep','All','Location','best')
    end
    axis square
    xlim([0.002,0.9])
    % ylim([0.01,100])
    set(gca,'box','off')
%     ax4.TickLength = [0.03,0.03];
end
    %% save the figure
    if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(powerSpec_Fig,[dirpath 'Fig_PowerSpectrum_AwakeAsleepAlert']);
        set(powerSpec_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_PowerSpectrum_AwakeAsleepAlert'])
    end
    close
