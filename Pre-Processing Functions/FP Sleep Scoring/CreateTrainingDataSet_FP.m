function [] = CreateTrainingDataSet_FP(procDataFileIDs,NBins)
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
%         imagingType = 'bilateral';
%         hemoType = 'CBV';
        [figHandle,ax1,ax2,ax3,ax5,ax6] = GenerateSingleFigures_Sleep_FP(procDataFileID,saveFigs);
%         trialDuration = ProcData.notes.trialDuration_sec;
        numBins = NBins;
        behavioralState = cell(numBins,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            global ButtonValue %#ok<TLEV>
            ButtonValue = 0;
            xStartVal = (b*5) - 4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            figHandle = gcf;
%             subplot(ax1)
%             hold on
%             leftEdge1 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
%             rightEdge1 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
%             subplot(ax2)
%             hold on
%             leftEdge2 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
%             rightEdge2 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            subplot(ax3)
            hold on
            leftEdge3 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge3 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
%             subplot(ax3)
%             hold on
% %             leftEdge4 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
%             rightEdge4 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
%             subplot(ax5)
%             hold on
%             leftEdge5 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
%             rightEdge5 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            subplot(ax6)
            hold on
            leftEdge6 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge6 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            if b <= 600
            xaxis_1st = (ceil(b/60)-1)*300;
            xaxis_2nd = (ceil(b/60))*300;
            elseif b >= 601
            xaxis_1st = 3000;
            xaxis_2nd = 3120;
            end
            xlim([xaxis_1st xaxis_2nd])
%             if b == 1 % b <= 60
%                 xlim([1,300])
%             elseif b == 61 % b >= 61 && b <= 120
%                 xlim([300,600])
%             elseif b == 121 % b >= 121 && b <= 180
%                 xlim([600,900])
%             elseif b == 181 % b >= 181 && b <= 240
%                 xlim([900,1200])
%             elseif b == 241 % b >= 241 && b <= 300
%                 xlim([1200,1500])
%             elseif b == 301 % b >= 301 && b <= 360
%                 xlim([1500,1800])
%             elseif b == 361 % b >= 361 && b <= 420
%                 xlim([1800,2100])
%             elseif b == 421 % b >= 421 && b <= 480
%                 xlim([2100,2400])
%             elseif b == 481 % b >= 481 && b <= 540
%                 xlim([2400,2700])
%             elseif b == 541 % b >= 541 && b <= 600
%                 xlim([2700,3000])
%             elseif b == 601 % b >= 601 && b <= 660
%                 xlim([3000,3300])
%             elseif b == 661 % b >= 661 && b <= 720
%                 xlim([3300,3600])
%             end
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
%             delete(leftEdge1)
%             delete(leftEdge2)
            delete(leftEdge3)
%             delete(leftEdge4)
%             delete(leftEdge5)
            delete(leftEdge6)
%             delete(rightEdge1)
%             delete(rightEdge2)
            delete(rightEdge3)
%             delete(rightEdge4)
%             delete(rightEdge5)
            delete(rightEdge6)
        end
        close(figHandle)
        paramsTable.behavState = behavioralState;
        trainingTable = paramsTable;
        save(trainingDataFileID, 'trainingTable')
    else
        disp([trainingDataFileID ' already exists. Continuing...']); disp(' ')
%         disp(['Loading ' procDataFileID ' to adjust parameters.' ]); disp(' ')
%         load(modelDataFileID)
%         load(trainingDataFileID)
%         paramsTable.behavState = trainingTable.behavState;
%         trainingTable = paramsTable;
%         save(trainingDataFileID, 'trainingTable')
    end
end

end
