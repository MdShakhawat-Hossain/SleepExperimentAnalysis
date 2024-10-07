function [AnalysisResults] = PlotCrossCorrelation_GRABNE(rootFolder,saveFigs,delim,AnalysisResults,firstHrs,animalIDs)
%% cross correlation



    if firstHrs == "false"
        behavFields = {'NREM','REM','Rest','Asleep','All','Alert'}; % 'Rest',
    elseif firstHrs == "true"
        behavFields = {'Rest','NREM','Asleep','All'};
    end

SdataTypes ={'NE_GFP','ACh_GFP'}; %  {'AChNE'}; %
dataTypes = {'ACh_CBV','NE_CBV'}; %,'ACh_GFP','NE_GFP'

% concatenate the cross-correlation during different arousal states for eACh animal
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(dataTypes)
            dataType = dataTypes{1,cc};
            % pre-allocate necessary variable fields
            for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                lagsName = [SdataType 'lags'];
                xcValsName = [SdataType 'xcVals'];

                        data.XCorr.(behavField).(dataType).dummCheck = 1;
                        if isfield(data.XCorr.(behavField).(dataType),(xcValsName)) == false
                            data.XCorr.(behavField).(dataType).(xcValsName) = [];
                            % lags and stats fields
                            data.XCorr.(behavField).(dataType).(lagsName) = [];
                            data.XCorr.(behavField).(dataType).animalID = {};
                            data.XCorr.(behavField).(dataType).behavField = {};
                        end
                        % concatenate cross correlation during eACh arousal state, find peak + lag time
                        if isfield(AnalysisResults.(animalID).CrossCorrelation,(behavField)) == true
                            if isempty(AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).xcVals) == false
                                % peak + lag time
                                data.XCorr.(behavField).(dataType).(xcValsName) = cat(1,data.XCorr.(behavField).(dataType).(xcValsName),AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).xcVals);
                                % lags and stats fields
                                data.XCorr.(behavField).(dataType).(lagsName) = cat(1,data.XCorr.(behavField).(dataType).(lagsName),AnalysisResults.(animalID).CrossCorrelation.(behavField).(SdataType).(dataType).lags);
%                                 data.XCorr.(behavField).(dataType).animalID = cat(1,data.XCorr.(behavField).(dataType).animalID,animalID);
                                data.XCorr.(behavField).(dataType).behavField = cat(1,data.XCorr.(behavField).(dataType).behavField,behavField);
                            end
                        end
            end
        end
    end
end
samplingRate = 30;
% mean and standard error/standard deviation of cross correlation values
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for ff = 1:length(dataTypes)
        dataType = dataTypes{1,ff};
         for nn = 1:length(SdataTypes)
                SdataType = SdataTypes{1,nn};
                lagsName = [SdataType 'lags'];
                xcValsName = [SdataType 'xcVals'];
                meanlagsName = ['mean' SdataType 'lags'];
                meanxcValsName = ['mean' SdataType 'xcVals'];
                stdxcValsName = ['std' SdataType 'xcVals'];
                semxcValsName = ['sem' SdataType 'xcVals'];
                peakxcValsNameVal = ['peak' SdataType 'Val'];
                peakxcValsNameMean = ['peak' SdataType 'Val'];
                peakxcValsNameStd = ['peak' SdataType 'Std'];
                peakxcValsNameLag = ['peak' SdataType 'Lag'];

                % Lags time vector
                data.XCorr.(behavField).(dataType).(meanlagsName) = mean(data.XCorr.(behavField).(dataType).(lagsName),1);
                % XC mean/sem
                data.XCorr.(behavField).(dataType).(meanxcValsName) = mean(data.XCorr.(behavField).(dataType).(xcValsName),1);
                data.XCorr.(behavField).(dataType).(stdxcValsName) = std(data.XCorr.(behavField).(dataType).(xcValsName),0,1);
                data.XCorr.(behavField).(dataType).(semxcValsName) = std(data.XCorr.(behavField).(dataType).(xcValsName),0,1)./sqrt(size(data.XCorr.(behavField).(dataType).(xcValsName),1));
                % find peak lag time and value/std at that time
                [~,idx] = max(abs(data.XCorr.(behavField).(dataType).(meanxcValsName)));
                data.XCorr.(behavField).(dataType).(peakxcValsNameVal) = data.XCorr.(behavField).(dataType).((xcValsName))(1,idx);
                data.XCorr.(behavField).(dataType).(peakxcValsNameMean) = mean(data.XCorr.(behavField).(dataType).(peakxcValsNameVal),1);
                data.XCorr.(behavField).(dataType).(peakxcValsNameStd) = std(data.XCorr.(behavField).(dataType).(peakxcValsNameVal),0,1);
                data.XCorr.(behavField).(dataType).(peakxcValsNameLag) = data.XCorr.(behavField).(dataType).((meanlagsName))(1,idx)/samplingRate;
         end
    end
end
%% cross correlation [rest, NREM, REM]

for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        Crosscorrelation_Fig = figure('Name','CrossCorrelation_Figure');
        sgtitle([(dataType) 'CrossCorrelation in different arousal stages'])
        SubplotNeed = ceil(length(SdataTypes)/3);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
                meanlagsName = ['mean' SdataType 'lags'];
                meanxcValsName = ['mean' SdataType 'xcVals'];
%                 stdxcValsName = ['std' SdataType 'xcVals'];
                semxcValsName = ['sem' SdataType 'xcVals'];
%                 peakxcValsNameVal = ['peak' SdataType 'Val'];
%                 peakxcValsNameStd = ['peak' SdataType 'Std'];
%                 peakxcValsNameLag = ['peak' SdataType 'Lag'];
            
            subplot(3,SubplotNeed,splotNo);
            hold on
        
            freq = 30;
            lagSec = 30;
            L1 = plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName),'color',[0 0 0],'LineWidth',2);
            hold on
            plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName) + data.XCorr.Rest.(dataType).(semxcValsName),'color',[0 0 0],'LineWidth',0.1);
            plot(data.XCorr.Rest.(dataType).(meanlagsName),data.XCorr.Rest.(dataType).(meanxcValsName) - data.XCorr.Rest.(dataType).(semxcValsName),'color',[0 0 0],'LineWidth',0.1);

            L2 = plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName),'color',[0 0.4470 0.7410],'LineWidth',2);
            plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName) + data.XCorr.NREM.(dataType).(semxcValsName),'color',[0 0.4470 0.7410],'LineWidth',0.1);
            plot(data.XCorr.NREM.(dataType).(meanlagsName),data.XCorr.NREM.(dataType).(meanxcValsName) - data.XCorr.NREM.(dataType).(semxcValsName),'color',[0 0.4470 0.7410],'LineWidth',0.1);
            % if firstHrs == "false"
            hold on
            L3= plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName),'color',[1 0 0],'LineWidth',2);
                plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName) + data.XCorr.REM.(dataType).(semxcValsName),'color',[1 0 0],'LineWidth',0.1);
                plot(data.XCorr.REM.(dataType).(meanlagsName),data.XCorr.REM.(dataType).(meanxcValsName) - data.XCorr.REM.(dataType).(semxcValsName),'color',[1 0 0],'LineWidth',0.1);
            % elseif firstHrs == "true"
            
            % end
            L4 = plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName),'color',[1 0 1],'LineWidth',2);
            plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName) + data.XCorr.Asleep.(dataType).(semxcValsName),'color',[1 0 1],'LineWidth',0.1);
            plot(data.XCorr.Asleep.(dataType).(meanlagsName),data.XCorr.Asleep.(dataType).(meanxcValsName) - data.XCorr.Asleep.(dataType).(semxcValsName),'color',[1 0 1],'LineWidth',0.1);
            hold on
            L5= plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName),'color',[0.4660 0.6740 0.1880],'LineWidth',2);
            plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName) + data.XCorr.All.(dataType).(semxcValsName),'color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
            plot(data.XCorr.All.(dataType).(meanlagsName),data.XCorr.All.(dataType).(meanxcValsName) - data.XCorr.All.(dataType).(semxcValsName),'color',[0.4660 0.6740 0.1880],'LineWidth',0.1);
            hold on
            L6= plot(data.XCorr.Alert.(dataType).(meanlagsName),data.XCorr.Alert.(dataType).(meanxcValsName),'color',[0.9290 0.6940 0.1250],'LineWidth',2);
            plot(data.XCorr.Alert.(dataType).(meanlagsName),data.XCorr.Alert.(dataType).(meanxcValsName) + data.XCorr.Alert.(dataType).(semxcValsName),'color',[0.9290 0.6940 0.1250],'LineWidth',0.1);
            plot(data.XCorr.Alert.(dataType).(meanlagsName),data.XCorr.Alert.(dataType).(meanxcValsName) - data.XCorr.Alert.(dataType).(semxcValsName),'color',[0.9290 0.6940 0.1250],'LineWidth',0.1);

            xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
            xticklabels({'-30','-15','0','15','30'})
            xlim([-lagSec*freq,lagSec*freq])
%             ylim([-0.35,0.15])
            xlabel('Lags (s)')
            ylabel('Correlation')
            title({[(SdataType) ' '],...
                [(dataType) ' XCorr']})
            % if strcmp(SdataType,dataType) == 1
                 % if firstHrs == "false"
            legend([L1, L2, L3, L4, L5, L6],'Rest','NREM','REM','Asleep','All','Alert','Location','best');
            % legend('Rest','NREM','REM','All','Alert','Location','best');


                % elseif firstHrs == "true"
                %       legend([L1,L2,L5],'Rest','NREM','All','Location','best')
                % end   
            % end
            axis square
            set(gca,'box','off')

            splotNo = splotNo +1;
        end
        %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Crosscorrelation_Fig,[dirpath 'Fig_Crosscorrelation_' (SdataType) '_' (dataType)]);
        set(Crosscorrelation_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_Crosscorrelation_' (SdataType) '_' (dataType)])
    end
    close
end
    %%
    
for bb = 1:length(dataTypes)
        dataType = dataTypes{1,bb};
        Crosscorrelation_Fig = figure('Name','Peak CrossCorrelation_Figure');
        sgtitle([(dataType) 'Peak CrossCorrelation in different arousal stages'])
        SubplotNeed = ceil(length(SdataTypes)/3);
        splotNo = 1;
        for nn = 1:length(SdataTypes)
            SdataType = SdataTypes{1,nn};
                peakxcValsNameVal = ['peak' SdataType 'Val'];
                peakxcValsNameMean = ['peak' SdataType 'Val'];
                peakxcValsNameStd = ['peak' SdataType 'Std'];
                peakxcValsNameLag = ['peak' SdataType 'Lag'];
            
            subplot(3,SubplotNeed,splotNo);
            hold on

            scatter(ones(1,length(data.XCorr.Rest.(dataType).(peakxcValsNameVal)))*1,data.XCorr.Rest.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'jitter','on','jitterAmount',0.25);
            hold on
            e1 = errorbar(1,data.XCorr.Rest.(dataType).(peakxcValsNameMean),data.XCorr.Rest.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e1.Color = 'black';
            e1.MarkerSize = 10;
            e1.CapSize = 10;
            scatter(ones(1,length(data.XCorr.NREM.(dataType).(peakxcValsNameVal)))*2,data.XCorr.NREM.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[0 0.4470 0.7410],'jitter','on','jitterAmount',0.25);
            hold on
            e2 = errorbar(2,data.XCorr.NREM.(dataType).(peakxcValsNameMean),data.XCorr.NREM.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e2.Color = 'black';
            e2.MarkerSize = 10;
            e2.CapSize = 10;
            scatter(ones(1,length(data.XCorr.REM.(dataType).(peakxcValsNameVal)))*3,data.XCorr.REM.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'jitter','on','jitterAmount',0.25);
            hold on
            e3 = errorbar(3,data.XCorr.REM.(dataType).(peakxcValsNameMean),data.XCorr.REM.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e3.Color = 'black';
            e3.MarkerSize = 10;
            e3.CapSize = 10;
            scatter(ones(1,length(data.XCorr.Alert.(dataType).(peakxcValsNameVal)))*6,data.XCorr.Alert.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[0.9290 0.6940 0.1250],'jitter','on','jitterAmount',0.25);
            hold on
            e4 = errorbar(6,data.XCorr.Alert.(dataType).(peakxcValsNameVal),data.XCorr.Alert.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e4.Color = 'black';
            e4.MarkerSize = 10;
            e4.CapSize = 10;
            scatter(ones(1,length(data.XCorr.Asleep.(dataType).(peakxcValsNameVal)))*4,data.XCorr.Asleep.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'jitter','on','jitterAmount',0.25);
            hold on
            e5 = errorbar(4,data.XCorr.Asleep.(dataType).(peakxcValsNameMean),data.XCorr.Asleep.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e5.Color = 'black';
            e5.MarkerSize = 10;
            e5.CapSize = 10;
            scatter(ones(1,length(data.XCorr.All.(dataType).(peakxcValsNameVal)))*5,data.XCorr.All.(dataType).(peakxcValsNameVal),75,'MarkerEdgeColor','k','MarkerFaceColor',[0.4660 0.6740 0.1880],'jitter','on','jitterAmount',0.25);
            hold on
            e6 = errorbar(5,data.XCorr.All.(dataType).(peakxcValsNameMean),data.XCorr.All.(dataType).(peakxcValsNameStd),'d','MarkerEdgeColor','k','MarkerFaceColor','k');
            e6.Color = 'black';
            e6.MarkerSize = 10;
            e6.CapSize = 10;
            title({[(SdataType) ' '],...
                            [(dataType) 'Peak XCorr']})
            ylabel('Corr. Coef')
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            axis square
            xlim([0,7])
            set(gca,'box','off')
            splotNo = splotNo +1;
        end
% ax2.TickLength = [0.03,0.03];
  %% save the figure
     if strcmp(saveFigs,'y') == true
        dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(Crosscorrelation_Fig,[dirpath 'Fig_Crosscorrelation_Peak_' (SdataType) '_' (dataType) ]);
        set(Crosscorrelation_Fig,'PaperPositionMode','auto');
        print('-painters','-dpdf','-fillpage',[dirpath 'Fig_Crosscorrelation_Peak_' (SdataType) '_' (dataType)])
    end
    close
end
 


