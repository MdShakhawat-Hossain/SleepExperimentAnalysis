function CheckSleepScore_EEG()

% written by Md Shakhawat Hossain
% purpose
% the function plot the manual sleep scoring along with other behavioral,
% neural, and hemodynamic activity
% if you are not happy with the score and want to make any changes to the
% sleep score you can do that.
clear all;close all;clc
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    TDataIx = strfind(procDataFileID,'_');
    TrainingDataFileID = [procDataFileID(1:TDataIx(end)) 'TrainingData.mat'];
    disp(['Checking sleep scoring for file ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    saveFig = 'y';%if you want to save the sleep score checked figure
    [figHandle,ax1,ax2,ax3,ax4,ax5,min_EEG_LH,max_EEG_LH] = GenerateSingleFigures_Sleep_FP_EEG_CheckSleep(procDataFileID,saveFig); % plot the sleep scores.
    SaveChoice = 'n'; % are we going to save the modified data. The answer is no since we haven't decided whether we want to chane the scores.
    commandwindow;
    ChangeChoice = input('Do you want to change the manual scores for this hour (y/n)? : ','s');disp(' ')
    
    if strcmp(ChangeChoice,'y') == true % load the training data only if we want to make any changes
        load(TrainingDataFileID);
        behavioralState = trainingTable.behavState;
    end

    while strcmp(ChangeChoice,'y') == true % loop through until you decide to change your mind and not change sleep scores
        commandwindow;
        Timechoice = 'y'; % this is used to make sure we are selecting a correct time range. if any ranges are needed the loop will run
        while strcmp(Timechoice,'y') == true % get the time range.time is in seconds 
            startTime = input('Input the start time for resting data: '); disp(' ')
            endTime = input('Input the end time for resting data: '); disp(' ')
            disp(['The time range is  ' num2str(startTime) ':' num2str(endTime)]); disp(' ')
            Timechoice = input('Do you want to change the time range (y/n)? : ','s');disp(' ') % answer n if you are happy with your time range
        end
        startBin = round(startTime/5);
        endBin = ceil(endTime/5);
            
        for b = startBin:endBin
            global buttonState %#ok<TLEV>
            buttonState = 0;
            global ButtonValue %#ok<TLEV>
            ButtonValue = 0;
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
            xaxis_1st = (ceil(b/60)-1)*300;
            xaxis_2nd = (ceil(b/60))*300;
            elseif b >= 601
            xaxis_1st = 3000;
            xaxis_2nd = 3120;
            end
            xlim([xaxis_1st xaxis_2nd])
%             ylim(ax4,[min_EEG_LH,max_EEG_LH])
            %%
            % if you want to use the sleep scoring GUI designed by
            % kevin, uncomment this part

            
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
            end
            
            %% GUI to perform manual stage selection
            % wrote by Shakhawat. This GUI allows you to type in keyboard
            % along with selecting a button
            % if you want to use the other sleep scoring GUI designed by
            % kevin comment this part
%             Select_SleepStage;  % GUI to select a sleep state
%             while buttonState == 0                
%                 drawnow()
%                 if buttonState == 1
%                          if ButtonValue == 1
%                                 behavioralState{b,1} = 'Not Sleep';
%                          elseif ButtonValue == 2
%                                 behavioralState{b,1} = 'NREM Sleep';
%                          elseif ButtonValue == 3
%                                 behavioralState{b,1} = 'REM Sleep';
%                          elseif ButtonValue == 0
%                             warndlg('Please select any of the buttons or type an appropriate letter','Warning');
%                             buttonState = 0;
%                             Select_SleepStage;
%                          end
%                 end
%             end
            %%
            delete(leftEdge3)
            delete(leftEdge6)
            delete(rightEdge3)
            delete(rightEdge6)
        end
        xlim([0 length(behavioralState)*5]);
        commandwindow;
        ChangeChoice = input('Do you want to change the manual scores for any other time range during this hour (y/n): ','s'); % do you want to change some other time points? 
        SaveChoice = 'y'; % do we want to save our modified data. The answer is yes. 
    end
        close(figHandle) % close the figure and save the modified data
        if strcmp(SaveChoice,'y') == true
            paramsTable.behavState = behavioralState;
            trainingTable = paramsTable;
            save(TrainingDataFileID, 'trainingTable')
        end
end
