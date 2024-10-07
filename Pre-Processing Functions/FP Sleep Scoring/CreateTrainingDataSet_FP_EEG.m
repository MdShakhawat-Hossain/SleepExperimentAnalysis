function [] = CreateTrainingDataSet_FP_EEG(procDataFileIDs,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%CreateTrainingDataSet_FP(procDataFileIDs,RestingBaselines,baselineType,NBins)
%   Purpose: Go through each file and train a data set for the model or for model validation
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataFileID = [procDataFileID(1:end-12) 'ModelData.mat'];
    trainingDataFileID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    if ~exist(trainingDataFileID,'file') %|| exist(trainingDataFileID,'file') % just for these case. since I need to redo all the scoring
        disp(['Loading ' procDataFileID ' for manual sleep scoring.' ]); disp(' ')
        load(procDataFileID)
        load(modelDataFileID)
        saveFigs = 'n';
        [figHandle,ax1,ax2,ax3,ax4,ax5,min_EEG_LH,max_EEG_LH] = GenerateSingleFigures_Sleep_FP_EEG(procDataFileID,saveFigs);
        numBins = NBins;
        behavioralState = cell(numBins,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            xStartVal = (b*5) - 5;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            figHandle = gcf;

            subplot(ax3)
            hold on
            leftEdge3 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge3 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);

            subplot(ax4)
            hold on
            leftEdge6 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge6 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            
            if b <= 600
            xaxis_1st = (ceil(b/30)-1)*150;
            xaxis_2nd = (ceil(b/30))*150;
            elseif b >= 601
            xaxis_1st = 3000;
            xaxis_2nd = 3120;
            end
            xlim([xaxis_1st xaxis_2nd])
            ylim(ax4,[min_EEG_LH,max_EEG_LH])


            [updatedGUI] = SelectBehavioralStateGUI_IOS;
            while buttonState == 0
                drawnow()
                if buttonState == 1
                    guiResults = guidata(updatedGUI);
                    if guiResults.togglebutton1.Value == true
                        behavioralState{b,1} = 'Not Sleep';
                    elseif guiResults.togglebutton2.Value == true
                        behavioralState{b,1} = 'NREM Sleep';
                    elseif guiResults.togglebutton3.Value == true
                        behavioralState{b,1} = 'REM Sleep';
                    else
                        disp('No button pressed'); disp(' ')
                        keyboard
                    end
                    close(updatedGUI)
                    break;
                end
                ...
            end

            delete(leftEdge3)
            delete(leftEdge6)
            delete(rightEdge3)
            delete(rightEdge6)
        end
        close(figHandle)
        paramsTable.behavState = behavioralState;
        trainingTable = paramsTable;
        save(trainingDataFileID, 'trainingTable')
    else
        disp([trainingDataFileID ' already exists. Continuing...']); disp(' ')
    end
end

end
