function CreateMicroArousalsData_FP_GRABNE(startingDirectory,trainingDirectory,baselineDirectory,MicroArousalTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% Ph.D. Candidate, Department of Bioengineering
% The Pennsylvania State University
%________________________________________________________________________________________________________________________
%
%   Purpose: This function uses the sleep logicals in each ProcData file to find periods where there are 60 seconds of 
%            consecutive ones within the sleep logical (12 or more). If a ProcData file's sleep logical contains one or
%            more of these 60 second periods, each of those bins is gathered from the data and put into the SleepEventData.mat
%            struct along with the file's name. 
%________________________________________________________________________________________________________________________
%
%   Inputs: The function loops through each ProcData file within the current folder - no inputs to the function itself
%           This was done as it was easier to add to the SleepEventData struct instead of loading it and then adding to it
%           with each ProcData loop.
%
%   Outputs: SleepEventData.mat struct
%________________________________________________________________________________________________________________________

    cd(trainingDirectory)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);

%% need to make the code run through loops.
dataTypes = {'cortical_LH','EMG','Rhodamine','GFP','Pupil'};
%% BLOCK PURPOSE: Create sleep scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure
microarousalsBins = MicroArousalTime;
for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear -regexp ^PData ^LH_ ^RH_ ^cell ^mat2Cell ^mat;
clear  fixedmicroarousalsIndex microarousalsCriteria microarousalsIndex
        
    microarousalsLogical = ProcData.MicroArousals.MALogical;    % Logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1, microarousalsBins);   % Target time
    microarousalsIndex = find(conv(microarousalsLogical, targetTime) >= microarousalsBins) - (microarousalsBins - 1);   % Find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(microarousalsIndex)  % If microarousalsIndex is empty, skip this file
        % Skip file
    else
        microarousalsCriteria = (0:(microarousalsBins)-1);     % This will be used to fix the issue in microarousalsIndex
        fixedmicroarousalsIndex = unique(microarousalsIndex + microarousalsCriteria);   % MA Index now has the proper time stamps from MA logical
        for indexCount = 1:length(fixedmicroarousalsIndex)    % Loop through the length of MA Index, and pull out associated data
            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                    subDataTypes = {'Z_Ach','Z_NE'};  
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                else
                    subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                end
                for sn = 1:length(subDataTypes)
                    subType  = char(subDataTypes(sn));
                    PData.(dataType).(subType){indexCount, 1} = ProcData.MicroArousals.parameters.(dataType).(subType){fixedmicroarousalsIndex(indexCount), 1}; %#ok<*AGROW>
                end
            end
            WhiskerAcceleration{indexCount, 1} = ProcData.MicroArousals.parameters.whiskerAcceleration{fixedmicroarousalsIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 1*fixedmicroarousalsIndex(indexCount);
        end
        
        indexBreaks = find(fixedmicroarousalsIndex(2:end) - fixedmicroarousalsIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            % filtered signal bands

            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                    subDataTypes = {'Z_Ach','Z_NE'};
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                else
                    subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
            
                end
                for sn = 1:length(subDataTypes)
                    subType  = char(subDataTypes(sn));
                    % filtered signal bands
                    matData.(dataType).(subType) = cell2mat(PData.(dataType).(subType));
                    arrayData.(dataType).(subType) = reshape(matData.(dataType).(subType)', [1, size(matData.(dataType).(subType), 2)*size(matData.(dataType).(subType), 1)]);
                    cellData.(dataType).(subType) = {arrayData.(dataType).(subType)};
                end
            end

            % whisker acceleration
            for x = 1:length(WhiskerAcceleration)
                targetPoints = size(WhiskerAcceleration{1, 1}, 2);
                if size(WhiskerAcceleration{x, 1}, 2) ~= targetPoints
                    maxLength = size(WhiskerAcceleration{x, 1}, 2);
                    difference = targetPoints - size(WhiskerAcceleration{x, 1}, 2);
                    for y = 1:difference
                        WhiskerAcceleration{x, 1}(maxLength + y) = 0;
                    end
                end
            end
            
            matWhiskerAcceleration = cell2mat(WhiskerAcceleration);
            arrayWhiskerAcceleration = reshape(matWhiskerAcceleration', [1, size(matWhiskerAcceleration, 2)*size(matWhiskerAcceleration, 1)]);
            cellWhiskerAcceleration = {arrayWhiskerAcceleration};

            % bin times
            matBinTimes = cell2mat(BinTimes);
            arrayBinTimes = reshape(matBinTimes', [1, size(matBinTimes, 2)*size(matBinTimes, 1)]);
            cellBinTimes = {arrayBinTimes};
        else
            count = length(fixedmicroarousalsIndex);
            holdIndex = zeros(1, (length(indexBreaks) + 1));
            
            for indexCounter = 1:length(indexBreaks) + 1
                if indexCounter == 1
                    holdIndex(indexCounter) = indexBreaks(indexCounter);
                elseif indexCounter == length(indexBreaks) + 1
                    holdIndex(indexCounter) = count - indexBreaks(indexCounter - 1);
                else
                    holdIndex(indexCounter)= indexBreaks(indexCounter) - indexBreaks(indexCounter - 1);
                end
            end
            
            splitCounter = 1:length(PData.cortical_LH.deltaBandPower);
            convertedMat2Cell = mat2cell(splitCounter', holdIndex);
            
            for matCounter = 1:length(convertedMat2Cell)
                 for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            % filtered signal bands
                            mat2CellData.(dataType).(subType){matCounter, 1} = PData.(dataType).(subType)(convertedMat2Cell{matCounter, 1});
                        end
                 end
                mat2CellWhiskerAcceleration{matCounter, 1} = WhiskerAcceleration(convertedMat2Cell{matCounter, 1});
                mat2CellBinTimes{matCounter, 1} = BinTimes(convertedMat2Cell{matCounter, 1});
            end
            
            for cellCounter = 1:length(mat2CellData.cortical_LH.deltaBandPower)
                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            % filtered signal bands
                            matData.(dataType).(subType) = cell2mat(mat2CellData.(dataType).(subType){cellCounter, 1});
                            arrayData.(dataType).(subType) = reshape(matData.(dataType).(subType)', [1, size(matData.(dataType).(subType), 2)*size(matData.(dataType).(subType), 1)]);
                            cellData.(dataType).(subType){cellCounter, 1} = arrayData.(dataType).(subType);
                        end
                end
                
                for x = 1:size(mat2CellWhiskerAcceleration{cellCounter, 1}, 1)
                    targetPoints = size(mat2CellWhiskerAcceleration{cellCounter, 1}{1, 1}, 2);
                    if size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2) ~= targetPoints
                        maxLength = size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2);
                        difference = targetPoints - size(mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}, 2);
                        for y = 1:difference
                            mat2CellWhiskerAcceleration{cellCounter, 1}{x, 1}(maxLength + y) = 0;
                        end
                    end
                end
                
                matWhiskerAcceleration = cell2mat(mat2CellWhiskerAcceleration{cellCounter, 1});
                arrayWhiskerAcceleration = reshape(matWhiskerAcceleration', [1, size(matWhiskerAcceleration, 2)*size(matWhiskerAcceleration, 1)]);
                cellWhiskerAcceleration{cellCounter, 1} = arrayWhiskerAcceleration;
                
                matBinTimes = cell2mat(mat2CellBinTimes{cellCounter, 1});
                arrayBinTimes = reshape(matBinTimes', [1, size(matBinTimes, 2)*size(matBinTimes, 1)]);
                cellBinTimes{cellCounter, 1} = arrayBinTimes;
            end
        end
        
        %% BLOCK PURPOSE: Save the data in the MicroArousalsEventData struct
        if exist('MicroArousalsData','var') == 0
%         if isfield(MicroArousalsData) == false % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 2)   % Loop through however many sleep epochs this file has
                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            MicroArousalsData.data.(dataType).(subType){cellLength, 1} = cellData.(dataType).(subType){1, 1};
                        end
                end
                
                MicroArousalsData.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                MicroArousalsData.FileIDs{cellLength, 1} = fileID;
                MicroArousalsData.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 1)   % Loop through however many sleep epochs this file has

                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'Rhodamine') == true || strcmp(dataType,'GFP') == true
                        subDataTypes = {'Z_Ach','Z_NE'};
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            MicroArousalsData.data.(dataType).(subType){size(MicroArousalsData.data.(dataType).(subType), 1) + 1, 1} = cellData.(dataType).(subType){cellLength, 1};

                        end
                end    
                MicroArousalsData.data.WhiskerAcceleration{size(MicroArousalsData.data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                MicroArousalsData.FileIDs{size(MicroArousalsData.FileIDs, 1) + 1, 1} = fileID;
                MicroArousalsData.BinTimes{size(MicroArousalsData.BinTimes, 1) + 1, 1} = cellBinTimes{cellLength, 1};
            end
        end
    end
    
    disp(['Adding microarousals epochs from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
end

%% clear variables
clear -regexp cell^ mat^ mat2cell^ array^


disp(' Data added to microarousals structure.'); disp(' ')
cd(startingDirectory)
cd(baselineDirectory)
AnimalIdx = strfind(procDataFileID,'_');
animalID = procDataFileID(1:AnimalIdx(1)-1);
save([animalID '_MicroArousalsData.mat'],'MicroArousalsData')

end
