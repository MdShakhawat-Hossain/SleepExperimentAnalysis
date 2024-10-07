clear all; close all; clc
% function GRABNE_Pupil_CrossCorrelation()
dataTypes = {'NE_Rhodamine'};
pupilDataTypes = {'zDiameter'};
animalIDs = {'GRABNE001','GRABNE002'};
rootFolder = pwd();
    for AID = 1:length(animalIDs)
        animalID = animalIDs{AID};
        % go to animal's data location
        delim = '\';

        dataLocation = [rootFolder delim animalID delim 'CombinedImaging'];
        cd(dataLocation)
        
        % character list of all ProcData file IDs
        procDataFileStruct = dir('*_ProcData.mat');
        procDataFiles = {procDataFileStruct.name}';
        procDataFiles = procDataFiles(3:end);
        procDataFileIDs = char(procDataFiles);
    
        for FLId = 1:length(procDataFiles)
            procDataFileId = procDataFileIDs(FLId,:);
            load(procDataFileId);
            Flix = strfind(procDataFileId,'_');
            fileName = procDataFileId(Flix(1)+1:Flix(4)+1);
            Analysis.GRABNEData.(animalID)(:,FLId) = ProcData.data.GFP.NE(1:93600)';
            Analysis.PupilData.(animalID)(:,FLId) = ProcData.data.Pupil.zDiameter(1:93600);
            Analysis.FileIDInfo.(animalID){:,FLId} = fileName;
        end
    end

    LagTime = 10; % seconds
    Frequency = 30; % Hz
    MaxLag = LagTime*Frequency;

    for AID = 1:length(animalIDs)
        animalID = animalIDs{AID};
        for FLId = 1:length(procDataFiles)
                    NEArray = Analysis.GRABNEData.(animalID)(:,FLId);
                    Pupilarray = Analysis.PupilData.(animalID)(:,FLId);
                    [XcVals(FLId,:),Lags(FLId,:)] = xcorr(NEArray,Pupilarray,MaxLag,'coeff');   
        end
        Analysis.CrossCor.(animalID).XVals = XcVals ;
        Analysis.CrossCor.(animalID).Lags = Lags ;
        clear XcVals Lags
    end

    freq = 30;
    lagSec = 10;
 for AID = 1:length(animalIDs)
        animalID = animalIDs{AID};
        SubplotNeed = ceil(length(procDataFiles)/4);
        splotNo = 1;
        Crosscorrelation_Fig = figure('Name','CrossCorrelation_Figure');
        sgtitle([(animalID) 'CrossCorrelation between GRABNE and Pupil Diameter'])
        for FLId = 1:length(procDataFiles)   
            subplot(SubplotNeed,4,splotNo);
            L1 = plot(Analysis.CrossCor.(animalID).Lags(FLId,:),Analysis.CrossCor.(animalID).XVals(FLId,:),'color','green','LineWidth',2);
            xticks([-lagSec*freq,-lagSec*freq/2,0,lagSec*freq/2,lagSec*freq])
            xticklabels({'-5','-2.5','0','2.5','5'})
            xlim([-lagSec*freq,lagSec*freq])
        %             ylim([-0.35,0.15])
            xlabel('Lags (s)')
            ylabel('Correlation')
            title([Analysis.FileIDInfo.(animalID){1,FLId} 'XCorr'])
            axis square
            set(gca,'box','off')
            r = Analysis.CrossCor.(animalID).XVals(FLId,:);
            lags = Analysis.CrossCor.(animalID).Lags(FLId,:);
            [~,index] = max(r(1:(length(r) - 1)/2) - freq/lagSec);
            offset = lags(index);
            offsetTime = abs(offset/freq);
            legend(['NE trailed pupil by' num2str(offsetTime, '%.1f') 's'],'Location','best')
            disp(['NE trailed pupil by' num2str(offsetTime, '%.1f') 's'])
            splotNo = splotNo +1;
        end

        %% save the figure
            dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(Crosscorrelation_Fig,[dirpath  'GRABNE_Pupil_Crosscorrelation_' (animalID)]);
            set(Crosscorrelation_Fig,'PaperPositionMode','auto');
            print('-painters','-dpdf','-fillpage',[dirpath 'GRABNE_Pupil_Fig_Crosscorrelation_' (animalID)])
        close

 end
 cd(rootFolder)
 save('GRABNE_Pupil_CrossCorrelation.mat','Analysis','-v7.3')