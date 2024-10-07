function MicroArousal_DataSeparation(startingDirectory,trainingDirectory,baselineDirectory,MicroArousalTimeMax,MicroArousalTimeMin)

    cd(trainingDirectory)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);

    % dataTypes = {'cortical_LH','EMG','Rhodamine','GFP','Pupil'};
    dataTypes = {'Rhodamine','GFP','EMG'};
%% BLOCK PURPOSE: Create MA scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure

for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear microarousalsLogical  ArousalTransition SleepTransition ArousalTime MA_Threshold MA_Index MA_id MAStart_Index  MAEnd_Index 
    clear MAStart_Index_Data MAEnd_Index_Data indexCount BinTimes cellLength PData

    clear -regexp ^PData ^LH_ ^RH_ ^cell ^mat2Cell ^mat;
    % clear  fixedmicroarousalsIndex microarousalsCriteria microarousalsIndex
        
    microarousalsLogical = ProcData.MicroArousals.MALogical;    % Logical - ones denote potential sleep epoches (5 second bins)

    % detect microarousals that are shorter or equal to the target time
    ArousalTransition = find(diff(microarousalsLogical)==1)+1;
    SleepTransition = find(diff(microarousalsLogical)==-1)+1;
    if ArousalTransition(1) < SleepTransition(1) % arousal is detected first
        ArousalTime = SleepTransition - ArousalTransition;
    elseif ArousalTransition(1) > SleepTransition(1) % sleep is detected first
        if length(ArousalTransition) < length(SleepTransition) 
            SleepTransition(1) = []; % first value is an extra sleep. removing that
            ArousalTime = SleepTransition - ArousalTransition;
        else
            disp('code needs modification');
        end
    end

    MA_Threshold = (ArousalTime <= MicroArousalTimeMax) & (ArousalTime >= MicroArousalTimeMin) ;

    if sum(MA_Threshold >= 1) % there are some microarousals data that fits the criteria
        %% extract the microuarousals data
        MA_Index = ArousalTransition .* MA_Threshold;
        MA_id = find(MA_Index ~=0); % exclude the zero elements

        MAStart_Index = ArousalTransition(MA_id);  % the index of microarousals that is in between 2-5 seconds
        MAEnd_Index = SleepTransition(MA_id); % the index of microarousals that is in between 2-5 seconds

        % we want to extract 1 seconds before the microarousals.
        MAStart_Index_Data = MAStart_Index - 1;
        MAEnd_Index_Data = MAEnd_Index-1;
    
        for indexCount =  1:1:length(MAStart_Index_Data)
                for dn = 1:length(dataTypes)
                dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};  
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'zArea','zDiameter','eyeMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                    end
                    for sn = 1:length(subDataTypes)
                        subType  = char(subDataTypes(sn));
                        PData.(dataType).(subType){indexCount, 1} = (cat(1,ProcData.MicroArousals.parameters.(dataType).(subType){MAStart_Index_Data(indexCount):MAEnd_Index_Data(indexCount)})); %#ok<*AGROW>
                    end
                end
                % WhiskerAcceleration{indexCount, 1} = ProcData.MicroArousals.parameters.whiskerAcceleration{fixedmicroarousalsIndex(indexCount), 1};
                BinTimes{indexCount, 1} = MAEnd_Index(indexCount) - MAStart_Index(indexCount);
        end
         %% BLOCK PURPOSE: Save the data in the MicroArousalsEventData struct
        if exist('MicroArousalsData','var') == 0
%         if isfield(MicroArousalsData) == false % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(PData.GFP.Z_Ach,1)   % Loop through however many sleep epochs this file has
                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'zArea','zDiameter','eyeMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            MicroArousalsData.data.(dataType).(subType){cellLength, 1} = PData.(dataType).(subType){cellLength, 1};
                        end
                end
                
                % MicroArousalsData.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                MicroArousalsData.FileIDs{cellLength, 1} = fileID;
                MicroArousalsData.BinTimes{cellLength, 1} = BinTimes{cellLength, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(PData.GFP.Z_Ach, 1)   % Loop through however many sleep epochs this file has

                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'zArea','zDiameter','eyeMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            MicroArousalsData.data.(dataType).(subType){size(MicroArousalsData.data.(dataType).(subType), 1) + 1, 1} = PData.(dataType).(subType){cellLength, 1};

                        end
                end    
                % MicroArousalsData.data.WhiskerAcceleration{size(MicroArousalsData.data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                MicroArousalsData.FileIDs{size(MicroArousalsData.FileIDs, 1) + 1, 1} = fileID;
                MicroArousalsData.BinTimes{size(MicroArousalsData.BinTimes, 1) + 1, 1} = BinTimes{cellLength, 1};
            end
        end
    else
            disp('No microarousals data found within the criteria. Please change the criteria. Moving to next file'); disp('')
    end
        disp(['Adding microarousals epochs from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
end
%% save the data
disp(' Data added to microarousals structure.'); disp(' ')
cd(startingDirectory)
cd(baselineDirectory)
AnimalIdx = strfind(procDataFileID,'_');
animalID = procDataFileID(1:AnimalIdx(1)-1);
save([animalID '_MicroArousalsData.mat'],'MicroArousalsData')
    
end