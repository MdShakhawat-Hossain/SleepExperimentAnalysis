function [AnalysisResults] = PlotCoherence_GRABNE_AAA(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)


%%
FigInital = strfind(rootFolder,'\');
ManipulationType = rootFolder(FigInital(end)+1:end);
%% coherence
    if firstHrs == "false"
         behavFields = {'Rest','NREM','Awake','Asleep','All'}; %'REM',
    elseif firstHrs == "true"
        behavFields = {'Rest','NREM','Awake'};
    end

dataTypes = {'NE_GFP','ACh_GFP'};
SdataTypes = {'ACh_CBV','NE_CBV','ACh_GFP','NE_GFP'};

% pre-allocate data structure
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
                data.Coherr.(behavField).(dataType).(CName) = [];
                data.Coherr.(behavField).(dataType).(fName) = [];
%                 data.Coherr.(behavField).(dataType).animalID = {};
                data.Coherr.(behavField).(dataType).behavField = {};
        end
    end
end
% concatenate coherence during different arousal states for each animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
            % don't concatenate empty arrays where there was no data for this behavior
                if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).C) == false
                    data.Coherr.(behavField).(dataType).(CName) = cat(2,data.Coherr.(behavField).(dataType).(CName),AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).C);
                    data.Coherr.(behavField).(dataType).(fName) = cat(1,data.Coherr.(behavField).(dataType).(fName),AnalysisResults.(animalID).Coherence.(behavField).(dataType).(SdataType).f);
%                     data.Coherr.(behavField).(dataType).animalID = cat(1,data.Coherr.(behavField).(dataType).animalID,animalID,animalID);
                    data.Coherr.(behavField).(dataType).behavField = cat(1,data.Coherr.(behavField).(dataType).behavField,behavField,behavField);
                end
            end
        end
    end
end
% mean and standard error for arousal state coherence
for aa = 1:length(behavFields)
    behavField = behavFields{1,aa};
    for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                CName = [SdataType 'C'];
                fName = [SdataType 'f'];
                meanCName = ['mean' SdataType 'C'];
                semCName = ['sem' SdataType 'C'];
                meanfName = ['mean' SdataType 'f'];
                data.Coherr.(behavField).(dataType).(meanCName) = mean(data.Coherr.(behavField).(dataType).(CName),2);
                data.Coherr.(behavField).(dataType).(semCName) = std(data.Coherr.(behavField).(dataType).(CName),0,2)./sqrt(size(data.Coherr.(behavField).(dataType).(CName),2));
                data.Coherr.(behavField).(dataType).(meanfName) = mean(data.Coherr.(behavField).(dataType).(fName),1);
        end
    end
end
%% plot coherence during arousal states

 for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        Coherence_Fig = figure('Name','Coherence_Figure');
        sgtitle([(dataType) 'Coherence in different arousal stages'])
        SubplotNeed = ceil(length(SdataTypes)/2);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
            meanCName = ['mean' SdataType 'C'];
            semCName = ['sem' SdataType 'C'];
            meanfName = ['mean' SdataType 'f'];
         subplot(2,SubplotNeed,splotNo);
            L4 =semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName),'color',[0.9290 0.6940 0.1250],'LineWidth',2);
            hold on

            semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName) + data.Coherr.Awake.(dataType).(semCName),'color',[0.9290 0.6940 0.1250],'LineWidth',0.5);
            semilogx(data.Coherr.Awake.(dataType).(meanfName),data.Coherr.Awake.(dataType).(meanCName) - data.Coherr.Awake.(dataType).(semCName),'color',[0.9290 0.6940 0.1250],'LineWidth',0.5);
            
            L5 = semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName),'color',[1 0 1],'LineWidth',2);
            semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName) + data.Coherr.Asleep.(dataType).(semCName),'color',[1 0 1],'LineWidth',0.5);
            semilogx(data.Coherr.Asleep.(dataType).(meanfName),data.Coherr.Asleep.(dataType).(meanCName) - data.Coherr.Asleep.(dataType).(semCName),'color',[1 0 1],'LineWidth',0.5);
            
            L6= semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName),'color',[0.4660 0.6740 0.1880],'LineWidth',2);
            semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName) + data.Coherr.All.(dataType).(semCName),'color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
            semilogx(data.Coherr.All.(dataType).(meanfName),data.Coherr.All.(dataType).(meanCName) - data.Coherr.All.(dataType).(semCName),'color',[0.4660 0.6740 0.1880],'LineWidth',0.5);
            xline(0.02,'color','b');
            xline(0.35,'color','r');
            % xline(1/3,'color','k');
            title({ (dataType) },...
                {[(SdataType) ' coherence']})
            ylabel('Coherence')
            xlabel('Freq (Hz)')
            if strcmp(dataType,SdataType) == 1
                if firstHrs == "false"
                    legend([L4,L5,L6],'Alert','Asleep','All','Location','best')
                elseif firstHrs == "true"
                    legend('Rest','NREM','Location','best')
                end
            end
            axis square
            xlim([0.004,1])
            xticks([0.01 0.1 1])
            ylim([0,1])
            set(gca,'box','off')
            splotNo = splotNo +1;
        end
        %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Coherence_Fig,[dirpath ManipulationType '_Coherence_logAxis_AwakeAsleepAll' (dataType)]);
        set(Coherence_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath ManipulationType '_Coherence_logAxis_AwakeAsleepAll' (dataType)])
    end
    close
 end

