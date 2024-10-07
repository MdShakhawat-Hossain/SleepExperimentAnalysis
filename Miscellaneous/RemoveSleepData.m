function [decData,decFileIDs,decDurations,decEventTimes] = RemoveSleepData(animalID,data,fileIDs,durations,eventTimes,SleepFileIDs,SleepBinTimes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Remove resting events from the various fields that aren't in the sleep selection
%________________________________________________________________________________________________________________________

% trialDuration_sec = 3120; % sec 3120
StimSize = 1; %
x = 1;
    for a = 1:size(fileIDs,1)
        fileID = fileIDs{a,1};
        startTime = floor(eventTimes(a,1))-10; 
        endTime = startTime + StimSize + 15;
        StimTime = startTime:1:endTime;

        trainingDataFileID = [animalID '_' fileID '_TrainingData.mat'];
        load(trainingDataFileID)


        OGLabels = trainingTable.behavState;
        TableSize = size(OGLabels,1);
            for SL = 1:1:TableSize
                if OGLabels(SL) == "Not Sleep"
                    for ML = 1:1:5
                        SleepLabels(((SL-1)*5)+ML) = 0;
                    end
                elseif OGLabels(SL) == "NREM Sleep"
                    for ML = 1:1:5
                        SleepLabels(((SL-1)*5)+ML) = 1;
                    end
                elseif OGLabels(SL) == "REM Sleep"
                    for ML = 1:1:5
                        SleepLabels(((SL-1)*5)+ML) = 1;
                    end
                end
            end

            StimLabels = zeros(1,size(SleepLabels,2));
            StimLabels(StimTime) = 1;
            SleepFlags = StimLabels.*SleepLabels;

        if sum(SleepFlags) == 0            
                if iscell(data) == true
                    decData{x,1} = data{a,1};
                else
                    decData(x,:) = data(a,:);
                end
                decFileIDs{x,1} = fileIDs{a,1};
                decDurations(x,1) = durations(a,1);
                decEventTimes(x,1) = eventTimes(a,1);
                x = x + 1;  
        end

    end




end

