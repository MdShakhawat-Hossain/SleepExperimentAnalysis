function [AnalysisResults] = Fig_S8_PupilNEACHCBV(rootFolder,saveFigs,delim,AnalysisResults,FP_animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Generate figure panel 1-S4 for Turner_Gheres_Proctor_Drew
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
% colorRfcAwake = [(0/256),(64/256),(64/256)];
% colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
colorWhisk = [(0/256),(166/256),(81/256)];
colorRest = [(31/256),(120/256),(179/256)];
colorStim = [(255/256),(28/256),(206/256)];
colorNREM = [(191/256),(0/256),(255/256)];
colorREM = [(254/256),(139/256),(0/256)];
% colorRest = [(255/256),(191/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
behavFields = {'Rest','NREM','REM'}; 
%% Rhodamine comparison between behaviors
% pre-allocate the date for each day
for aa = 1:length(FP_animalIDs)
    animalID = FP_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};  
            CellSize = length(AnalysisResults.(animalID).MeanPupil.(behavField).Pupil.IndzDiameter);
            for ll = 1:CellSize
            % CBVAch = mean(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.IndAch{ll,1});        
            % CBVNE = mean(AnalysisResults.(animalID).MeanRhodamine.(behavField).Rhodamine.IndnNE{ll,1});
            GRABAch = (AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndAch{ll,1})';
            GRABNE = (AnalysisResults.(animalID).MeanGFP.(behavField).GFP.IndNE{ll,1})';   
            PupilzDiameter = (AnalysisResults.(animalID).MeanPupil.(behavField).Pupil.IndzDiameter{ll,1})';
                figure(100)
                hold on
                subplot(1,3,bb)
                if bb == 1
                    % scatter3(GRABAch,PupilzDiameter,GRABNE,'MarkerEdgeColor','k','MarkerFaceColor',colorRest);axis equal
                    % xlabel('ACh');ylabel('Pupil');zlabel('GRABNE')
                    hist3([PupilzDiameter,GRABNE]);
                    xlabel('Pupil');ylabel('NE')
                    hold on
                elseif bb==2
                    hist3([PupilzDiameter,GRABNE]);
                    xlabel('Pupil');ylabel('NE')
                    % scatter3(GRABAch,PupilzDiameter,GRABNE,'MarkerEdgeColor','k','MarkerFaceColor',colorNREM);axis equal
                    hold on
                elseif bb==3
                    hist3([PupilzDiameter,GRABNE]);
                    xlabel('Pupil');ylabel('NE')
                    % scatter3(GRABAch,PupilzDiameter,GRABNE,'MarkerEdgeColor','k','MarkerFaceColor',colorREM);axis equal
                    hold on
                end
            end
    end
end

%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder delim 'Summary Figures and Structures' delim 'MATLAB Analysis Figures' delim];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    set(summaryFigure,'PaperPositionMode','auto');
    savefig(summaryFigure,[dirpath 'Fig1-S8-PupilCBV']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'Fig1-S8-PupilCBV'])
    close 
end

end
