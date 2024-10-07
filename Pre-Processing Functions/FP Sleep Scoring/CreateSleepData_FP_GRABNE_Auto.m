function [SleepData] = CreateSleepData_FP_GRABNE_Auto(startingDirectory,trainingDirectory,baselineDirectory,NREMsleepTime,REMsleepTime,modelName,SleepData)
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

if strcmp(modelName,'Manual') == false
    cd(baselineDirectory)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
else
    cd(trainingDirectory)
    % character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
end

%% need to make the code run through loops.
% sleepType = {'NREM','REM'};
dataTypes = {'cortical_LH','cortical_RH','EMG','CBV','GFP','Pupil'};

%% BLOCK PURPOSE: Create sleep scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure
sleepBins = NREMsleepTime/5;
for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear -regexp ^PData ^LH_ ^RH_ ^cell ^mat2Cell ^mat;
        
    nremLogical = ProcData.sleep.logicals.(modelName).nremLogical;    % Logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1, sleepBins);   % Target time
    sleepIndex = find(conv(nremLogical, targetTime) >= sleepBins) - (sleepBins - 1);   % Find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(sleepIndex)  % If sleepIndex is empty, skip this file
        % Skip file
    else
        sleepCriteria = (0:(sleepBins - 1));     % This will be used to fix the issue in sleepIndex
        fixedSleepIndex = unique(sleepIndex + sleepCriteria);   % sleep Index now has the proper time stamps from sleep logical
        for indexCount = 1:length(fixedSleepIndex)    % Loop through the length of sleep Index, and pull out associated data
            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                    subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'};  
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes = {'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                else
                    subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                end
                for sn = 1:length(subDataTypes)
                    subType  = char(subDataTypes(sn));
                    PData.(dataType).(subType){indexCount, 1} = ProcData.sleep.parameters.(dataType).(subType){fixedSleepIndex(indexCount), 1}; %#ok<*AGROW>
                end
            end
            WhiskerAcceleration{indexCount, 1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 5*fixedSleepIndex(indexCount);
        end
        
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            % filtered signal bands

            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                    subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
            count = length(fixedSleepIndex);
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
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'}; % {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
        %% BLOCK PURPOSE: Save the data in the SleepEventData struct
        if isfield(SleepData,(modelName)) == false  % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 2)   % Loop through however many sleep epochs this file has
                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            SleepData.(modelName).NREM.data.(dataType).(subType){cellLength, 1} = cellData.(dataType).(subType){1, 1};
                        end
                end
                
                SleepData.(modelName).NREM.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.(modelName).NREM.FileIDs{cellLength, 1} = fileID;
                SleepData.(modelName).NREM.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 1)   % Loop through however many sleep epochs this file has

                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            SleepData.(modelName).NREM.data.(dataType).(subType){size(SleepData.(modelName).NREM.data.(dataType).(subType), 1) + 1, 1} = cellData.(dataType).(subType){cellLength, 1};

                        end
                end    
                SleepData.(modelName).NREM.data.WhiskerAcceleration{size(SleepData.(modelName).NREM.data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                SleepData.(modelName).NREM.FileIDs{size(SleepData.(modelName).NREM.FileIDs, 1) + 1, 1} = fileID;
                SleepData.(modelName).NREM.BinTimes{size(SleepData.(modelName).NREM.BinTimes, 1) + 1, 1} = cellBinTimes{cellLength, 1};
            end
        end
    end
    
    disp(['Adding NREM sleeping epochs from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
end

%% clear variables
clear -regexp cell^ mat^ mat2cell^ array^
%% BLOCK PURPOSE: Create sleep scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure
sleepBins = REMsleepTime/5;
for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear -regexp ^PData ^LH_ ^RH_ ^cell ^mat2Cell ^mat;
        
    remLogical = ProcData.sleep.logicals.(modelName).remLogical;    % Logical - ones denote potential sleep epoches (5 second bins)
    targetTime = ones(1, sleepBins);   % Target time
    sleepIndex = find(conv(remLogical, targetTime) >= sleepBins) - (sleepBins - 1);   % Find the periods of time where there are at least 11 more
    % 5 second epochs following. This is not the full list.
    if isempty(sleepIndex)  % If sleepIndex is empty, skip this file
        % Skip file
    else
        sleepCriteria = (0:(sleepBins - 1));     % This will be used to fix the issue in sleepIndex
        fixedSleepIndex = unique(sleepIndex + sleepCriteria);   % sleep Index now has the proper time stamps from sleep logical
        for indexCount = 1:length(fixedSleepIndex)    % Loop through the length of sleep Index, and pull out associated data
            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                    subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                else
                    subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
            
                end
                for sn = 1:length(subDataTypes)
                    subType  = char(subDataTypes(sn));
                    % filtered signal bands
                    PData.(dataType).(subType){indexCount, 1} = ProcData.sleep.parameters.(dataType).(subType){fixedSleepIndex(indexCount), 1}; %#ok<*AGROW>
                end
            end
            WhiskerAcceleration{indexCount, 1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 5*fixedSleepIndex(indexCount);
        end
        
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            % filtered signal bands

            for dn = 1:length(dataTypes)
            dataType = char(dataTypes(dn));
                if strcmp(dataType,'EMG') == true
                    subDataTypes = {'emg','emgSignal'};
                elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                    subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                elseif strcmp(dataType,'Pupil') == true 
                    subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
            count = length(fixedSleepIndex);
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
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
        
        %% BLOCK PURPOSE: Save the data in the SleepEventData struct
        if isfield(SleepData.(modelName),'REM') == false  % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 2)   % Loop through however many sleep epochs this file has
                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            SleepData.(modelName).REM.data.(dataType).(subType){cellLength, 1} = cellData.(dataType).(subType){1, 1};
                        end
                end
                
                SleepData.(modelName).REM.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.(modelName).REM.FileIDs{cellLength, 1} = fileID;
                SleepData.(modelName).REM.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellData.cortical_LH.deltaBandPower, 1)   % Loop through however many sleep epochs this file has

                for dn = 1:length(dataTypes)
                    dataType = char(dataTypes(dn));
                    if strcmp(dataType,'EMG') == true
                        subDataTypes = {'emg','emgSignal'};
                    elseif strcmp(dataType,'CBV') == true || strcmp(dataType,'GFP') == true || strcmp(dataType,'Isos') == true
                        subDataTypes = {'Z_ACh','Z_NE','P_ACh','P_NE'}; 
                    elseif strcmp(dataType,'Pupil') == true 
                        subDataTypes ={'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
                    else
                        subDataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower','corticalPower'};
                
                    end
                        for sn = 1:length(subDataTypes)
                            subType  = char(subDataTypes(sn));
                            SleepData.(modelName).REM.data.(dataType).(subType){size(SleepData.(modelName).REM.data.(dataType).(subType), 1) + 1, 1} = cellData.(dataType).(subType){cellLength, 1};

                        end
                end    
                SleepData.(modelName).REM.data.WhiskerAcceleration{size(SleepData.(modelName).REM.data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                SleepData.(modelName).REM.FileIDs{size(SleepData.(modelName).REM.FileIDs, 1) + 1, 1} = fileID;
                SleepData.(modelName).REM.BinTimes{size(SleepData.(modelName).REM.BinTimes, 1) + 1, 1} = cellBinTimes{cellLength, 1};
            end
        end
    end
    
    disp(['Adding REM sleeping epochs from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
end

disp([modelName ' model data added to SleepData structure.']); disp(' ')
cd(startingDirectory)

end
