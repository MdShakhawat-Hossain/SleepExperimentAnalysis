function [AnalysisResults] = PlotCoherence_GRABNE_RNR(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%% coherence
    if firstHrs == "false"
         behavFields = {'Rest','NREM','REM','Awake','Asleep','All'};
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
% concatenate coherence during different arousal states for eACh animal
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
        SubplotNeed = ceil(length(SdataTypes)/3);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
            meanCName = ['mean' SdataType 'C'];
            semCName = ['sem' SdataType 'C'];
            meanfName = ['mean' SdataType 'f'];
         subplot(3,SubplotNeed,splotNo);

            L1 = semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName),'color',[0 0 0],'LineWidth',2);
                    hold on
                    semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName) + data.Coherr.Rest.(dataType).(semCName),'color',[0 0 0],'LineWidth',0.5);
                    semilogx(data.Coherr.Rest.(dataType).(meanfName),data.Coherr.Rest.(dataType).(meanCName) - data.Coherr.Rest.(dataType).(semCName),'color',[0 0 0],'LineWidth',0.5);
                    % rectangle('Position',[0.005,0.1,0.1 - 0.005,1],'FaceColor','w','EdgeColor','w')
            L2 = semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName),'color',[0 0.4470 0.7410],'LineWidth',2);
                    hold on
                    semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName) + data.Coherr.NREM.(dataType).(semCName),'color',[0 0.4470 0.7410],'LineWidth',0.5);
                    semilogx(data.Coherr.NREM.(dataType).(meanfName),data.Coherr.NREM.(dataType).(meanCName) - data.Coherr.NREM.(dataType).(semCName),'color',[0 0.4470 0.7410],'LineWidth',0.5);
                    % rectangle('Position',[0.005,0.1,1/30 - 0.005,1],'FaceColor','w','EdgeColor','w')
           if firstHrs == "false"
            L3 = semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName),'color',[1 0 0],'LineWidth',2);
                hold on        
                semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName) + data.Coherr.REM.(dataType).(semCName),'color',[1 0 0],'LineWidth',0.5);
                semilogx(data.Coherr.REM.(dataType).(meanfName),data.Coherr.REM.(dataType).(meanCName) - data.Coherr.REM.(dataType).(semCName),'color',[1 0 0],'LineWidth',0.5);
                % rectangle('Position',[0.005,0.1,1/60 - 0.005,1],'FaceColor','w','EdgeColor','w')
           elseif firstHrs == "true"
           end
           xline(1/30,'color','b');
            xline(0.1,'color','r');
            xline(0.2,'color','k');
            title({ (dataType) },...
                {[(SdataType) ' coherence']})
            ylabel('Coherence')
            xlabel('Freq (Hz)')
            if strcmp(dataType,SdataType) == 1
                if firstHrs == "false"
                    legend([L1,L2,L3],'Rest','NREM','REM','Location','best')
                elseif firstHrs == "true"
                    legend('Rest','NREM','Location','best')
                end
            end
            axis square
            xlim([0.01,0.9])
            xticks([0.01 0.03 0.1 0.4])
            set(gca,'box','off')
            splotNo = splotNo +1;
        end
        %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Coherence_Fig,[dirpath 'Fig_Coherence_RestNREMREM' (dataType)]);
        set(Coherence_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_Coherence_RestNREMREM' (dataType)])
    end
    close
 end

