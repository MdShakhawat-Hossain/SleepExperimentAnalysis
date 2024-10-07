function [NoSleepFileTag,CauseArousalFileTag,CauseSleepFileTag] = RemoveSleepStimData(finalOptoStimFileIDs,finalOptoStimFileEventTimes,animalID)

OptoStim_NoSleep = 1;
OptoStim_CauseArousal = 1;
OptoStim_CauseSleep = 1;

for K = 1:1:length(finalOptoStimFileIDs) % check all the files separately 
    FID = finalOptoStimFileIDs{K,1};
    FileID = [animalID '_' FID '_TrainingData.mat'];
    load(FileID)
    
    OGLabels =  trainingTable.behavState;
    SleepLabels = zeros(3120,1);
        for SL = 1:1:length(OGLabels)
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
    
        stimStartTime = round(finalOptoStimFileEventTimes(K,1)) - 5;
        stimEndTime = stimStartTime + 10;
        
        SleepBefore = 0;
        for ss = stimStartTime:1:stimStartTime+5
            if SleepLabels(ss) == 1
                SleepBefore = SleepBefore + 1;
            end
        end
    
        SleepAfter = 0;
        for ss = stimEndTime-5:1:stimEndTime
            if SleepLabels(ss) == 1
                SleepAfter = SleepAfter + 1;
            end
        end
    
        if SleepBefore == 0 % the mouse was not sleeping during optostim
            % NoSleepFileIDs{OptoStim_NoSleep,1} = finalOptoStimFileIDs{K,1};
            % NoSleepFileEventTimes{OptoStim_NoSleep,1} = finalOptoStimFileEventTimes{K,1};
            % NoSleepfinalOptoStimData{OptoStim_NoSleep,1} = finalOptoStimData{K,1};
            NoSleepFileTag(OptoStim_NoSleep) = K;
            OptoStim_NoSleep = OptoStim_NoSleep + 1;
        elseif (SleepBefore > 0) && (SleepAfter == 0)  % optostim cause arousal from sleep
            % CauseArousalFileIDs{OptoStim_CauseArousal,1} = finalOptoStimFileIDs{K,1};
            % CauseArousalFileEventTimes{OptoStim_CauseArousal,1} = finalOptoStimFileEventTimes{K,1};
            % CauseArousalfinalOptoStimData{OptoStim_CauseArousal,1} = finalOptoStimData{K,1};
            CauseArousalFileTag(OptoStim_CauseArousal) = K;
            OptoStim_CauseArousal = OptoStim_CauseArousal + 1;
        elseif (SleepBefore > 0) && (SleepAfter > 0)  % optostim does not cause arousal from sleep
            % CauseSleepFileIDs{OptoStim_CauseSleep,1} = finalOptoStimFileIDs{K,1};
            % CauseSleepFileEventTimes{OptoStim_CauseSleep,1} = finalOptoStimFileEventTimes{K,1};
            % CauseSleepfinalOptoStimData{OptoStim_CauseSleep,1} = finalOptoStimData{K,1};
            CauseSleepFileTag(OptoStim_CauseSleep) = K;
            OptoStim_CauseSleep = OptoStim_CauseSleep + 1;
        end
end

if OptoStim_NoSleep ==  1
    NoSleepFileTag = [];
end

if OptoStim_CauseArousal ==  1
    CauseArousalFileTag = [];
end

if OptoStim_CauseSleep ==  1
    CauseSleepFileTag = [];
end