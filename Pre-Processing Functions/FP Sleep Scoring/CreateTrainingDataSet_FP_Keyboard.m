function [] = CreateTrainingDataSet_FP_Keyboard(procDataFileIDs,NBins)

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    modelDataFileID = [procDataFileID(1:end-12) 'ModelData.mat'];
    trainingDataFileID = [procDataFileID(1:end-12) 'TrainingData.mat'];
    if ~exist(trainingDataFileID,'file') % || exist(trainingDataFileID,'file') % just for these case. since I need to redo all the scoring
        disp(['Loading ' procDataFileID ' for manual sleep scoring.' ]); disp(' ')
        load(procDataFileID)
        load(modelDataFileID)
        saveFigs = 'n';

        %% not all animal has good ECoG recording
        ECoGChoice = input('Do you want to use ECoG (y/n)?','s'); disp(" ");
        ECoGLoop = 1;
        while ECoGLoop == 1
            if ECoGChoice == "y"
                ECoGLoop = 2;
                [figHandle,ax1,ax2,ax3,ax4,ax6] = GenerateSingleFigures_Sleep_FP_GRABNE(procDataFileID,saveFigs);
                
            elseif ECoGChoice == "n"
                ECoGLoop = 2;
                [figHandle,ax1,ax2,ax3,ax4,ax5,ax6,ax7] = GenerateSingleFigures_Sleep_FP_GRABNE_NoECOG(procDataFileID,saveFigs);
            else 
                ECoGLoop = 2;
            end 
        end
        %%
        numBins = NBins;
        behavioralState = cell(numBins,1);
        for b = 1:numBins
            global buttonState %#ok<TLEV>
            buttonState = 0;
            global ButtonValue %#ok<TLEV>
            ButtonValue = 0;
            global closeButtonState %#ok<TLEV>
            closeButtonState = 0;

            xStartVal = (b*5) - 4;
            xEndVal = b*5;
            xInds = xStartVal:1:xEndVal;
            figHandle = gcf;

            subplot(ax4)
            hold on
            leftEdge3 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge3 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);

            subplot(ax6)
            hold on
            leftEdge6 = xline(xInds(1),'color',colors('electric purple'),'LineWidth',2);
            rightEdge6 = xline(xInds(5),'color',colors('electric purple'),'LineWidth',2);
            if b <= 600
            xaxis_1st = (ceil(b/100)-1)*500;
            xaxis_2nd = (ceil(b/100))*500;
            elseif b >= 601
            xaxis_1st = 3000;
            xaxis_2nd = 3120;
            end
            xlim([xaxis_1st xaxis_2nd])
            % condition to close the GUI
            if mod(b,100) == 0 || b == numBins
                closeButtonState = 1;
            end
            %% GUI for stage selection
            SelectSleepState_GUI; %The GUI to select the sleep stages
            while buttonState == 0
                drawnow()
                if buttonState == 1
                    if ButtonValue == 1
                        behavioralState{b,1} = 'Not Sleep';
                    elseif ButtonValue == 2
                        behavioralState{b,1} = 'NREM Sleep';
                    elseif ButtonValue == 3
                        behavioralState{b,1} = 'REM Sleep';
                    elseif ButtonValue == 0
                        warndlg('Please select one of the buttons or type an appropriate letter','Warning');
                        SelectSleepState_GUI;
                    end
                    break;
                end
                ...
            end
            
            delete(leftEdge3)
            delete(leftEdge6)
            delete(rightEdge3)
            delete(rightEdge6)
            closeButtonState = 0;
            pause(300/1000);% pause for 300 milliseconds
        end
        close(figHandle)
        paramsTable.behavState = behavioralState;
        trainingTable = paramsTable;
        save(trainingDataFileID, 'trainingTable')
    else
        disp([trainingDataFileID ' already exists. Continuing...']); disp(' ')
    end
end
close all
end
