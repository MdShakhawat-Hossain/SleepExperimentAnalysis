function [SleepData] = CreateSleepData_FP_EEG(startingDirectory,trainingDirectory,baselineDirectory,NREMsleepTime,REMsleepTime,modelName,SleepData)
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

%% BLOCK PURPOSE: Create sleep scored data structure.
% Identify sleep epochs and place in SleepEventData.mat structure-regexp
sleepBins = NREMsleepTime/5;
for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear -regexp ^LH_ ^RH_ ^cell ^mat2Cell ^mat;
        
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
            % filtered signal bands
            LH_deltaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.deltaBandPower{fixedSleepIndex(indexCount), 1}; %#ok<*AGROW>
            RH_deltaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.deltaBandPower{fixedSleepIndex(indexCount), 1};
            LH_thetaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.thetaBandPower{fixedSleepIndex(indexCount), 1};
            RH_thetaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.thetaBandPower{fixedSleepIndex(indexCount), 1};
            LH_alphaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.alphaBandPower{fixedSleepIndex(indexCount), 1};
            RH_alphaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.alphaBandPower{fixedSleepIndex(indexCount), 1};
            LH_betaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.betaBandPower{fixedSleepIndex(indexCount), 1};
            RH_betaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.betaBandPower{fixedSleepIndex(indexCount), 1};
            LH_gammaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.gammaBandPower{fixedSleepIndex(indexCount), 1};
            RH_gammaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.gammaBandPower{fixedSleepIndex(indexCount), 1};
            LH_EEGPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.EEGPower{fixedSleepIndex(indexCount), 1};
            RH_EEGPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.EEGPower{fixedSleepIndex(indexCount), 1};
            % Rhodamine
%             LH_Rhodamine{indexCount, 1} = ProcData.sleep.parameters.Rhodamine.LH{fixedSleepIndex(indexCount), 1};
            RH_Rhodamine{indexCount, 1} = ProcData.sleep.parameters.Rhodamine.RH{fixedSleepIndex(indexCount), 1};
%             LH_GFP{indexCount, 1} = ProcData.sleep.parameters.GFP.LH{fixedSleepIndex(indexCount), 1};
            RH_GFP{indexCount, 1} = ProcData.sleep.parameters.GFP.RH{fixedSleepIndex(indexCount), 1};   
            WhiskerAcceleration{indexCount, 1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 5*fixedSleepIndex(indexCount);
        end
        
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            % filtered signal bands
            matLH_DeltaPower = cell2mat(LH_deltaPower);
            arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
            cellLH_DeltaPower = {arrayLH_DeltaPower};
            
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            
            matLH_ThetaPower = cell2mat(LH_thetaPower);
            arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
            cellLH_ThetaPower = {arrayLH_ThetaPower};
            
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            
            matLH_AlphaPower = cell2mat(LH_alphaPower);
            arrayLH_AlphaPower = reshape(matLH_AlphaPower', [1, size(matLH_AlphaPower, 2)*size(matLH_AlphaPower, 1)]);
            cellLH_AlphaPower = {arrayLH_AlphaPower};
            
            matRH_AlphaPower = cell2mat(RH_alphaPower);
            arrayRH_AlphaPower = reshape(matRH_AlphaPower', [1, size(matRH_AlphaPower, 2)*size(matRH_AlphaPower, 1)]);
            cellRH_AlphaPower = {arrayRH_AlphaPower};
            
            matLH_BetaPower = cell2mat(LH_betaPower);
            arrayLH_BetaPower = reshape(matLH_BetaPower', [1, size(matLH_BetaPower, 2)*size(matLH_BetaPower, 1)]);
            cellLH_BetaPower = {arrayLH_BetaPower};
            
            matRH_BetaPower = cell2mat(RH_betaPower);
            arrayRH_BetaPower = reshape(matRH_BetaPower', [1, size(matRH_BetaPower, 2)*size(matRH_BetaPower, 1)]);
            cellRH_BetaPower = {arrayRH_BetaPower};
            
            matLH_GammaPower = cell2mat(LH_gammaPower);
            arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
            cellLH_GammaPower = {arrayLH_GammaPower};
            
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            
            matLH_MUAPower = cell2mat(LH_EEGPower);
            arrayLH_MUAPower = reshape(matLH_MUAPower', [1, size(matLH_MUAPower, 2)*size(matLH_MUAPower, 1)]);
            cellLH_MUAPower = {arrayLH_MUAPower};
            
            matRH_MUAPower = cell2mat(RH_EEGPower);
            arrayRH_MUAPower = reshape(matRH_MUAPower', [1, size(matRH_MUAPower, 2)*size(matRH_MUAPower, 1)]);
            cellRH_MUAPower = {arrayRH_MUAPower};

           
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
        
            % Rhodamine
%             matLH_Rhodamine = cell2mat(LH_Rhodamine);
%             arrayLH_Rhodamine = reshape(matLH_Rhodamine', [1, size(matLH_Rhodamine, 2)*size(matLH_Rhodamine, 1)]);
%             cellLH_Rhodamine = {arrayLH_Rhodamine};
            
            matRH_Rhodamine = cell2mat(RH_Rhodamine);
            arrayRH_Rhodamine = reshape(matRH_Rhodamine', [1, size(matRH_Rhodamine, 2)*size(matRH_Rhodamine, 1)]);
            cellRH_Rhodamine = {arrayRH_Rhodamine};
            
            % GFP
%             matLH_GFP = cell2mat(LH_GFP);
%             arrayLH_GFP = reshape(matLH_GFP', [1, size(matLH_GFP, 2)*size(matLH_GFP, 1)]);
%             cellLH_GFP = {arrayLH_GFP};
            
            matRH_GFP = cell2mat(RH_GFP);
            arrayRH_GFP = reshape(matRH_GFP', [1, size(matRH_GFP, 2)*size(matRH_GFP, 1)]);
            cellRH_GFP = {arrayRH_GFP};        

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
            
            splitCounter = 1:length(LH_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter', holdIndex);
            
            for matCounter = 1:length(convertedMat2Cell)
                mat2CellLH_DeltaPower{matCounter, 1} = LH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_DeltaPower{matCounter, 1} = RH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_ThetaPower{matCounter, 1} = LH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_ThetaPower{matCounter, 1} = RH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_AlphaPower{matCounter, 1} = LH_alphaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_AlphaPower{matCounter, 1} = RH_alphaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_BetaPower{matCounter, 1} = LH_betaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_BetaPower{matCounter, 1} = RH_betaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_GammaPower{matCounter, 1} = LH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_GammaPower{matCounter, 1} = RH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_MUAPower{matCounter, 1} = LH_EEGPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_MUAPower{matCounter, 1} = RH_EEGPower(convertedMat2Cell{matCounter, 1});
                
%                 mat2CellLH_Rhodamine{matCounter, 1} = LH_Rhodamine(convertedMat2Cell{matCounter, 1});
                mat2CellRH_Rhodamine{matCounter, 1} = RH_Rhodamine(convertedMat2Cell{matCounter, 1});
                
%                 mat2CellLH_GFP{matCounter, 1} = LH_GFP(convertedMat2Cell{matCounter, 1});
                mat2CellRH_GFP{matCounter, 1} = RH_GFP(convertedMat2Cell{matCounter, 1});
                
                mat2CellWhiskerAcceleration{matCounter, 1} = WhiskerAcceleration(convertedMat2Cell{matCounter, 1});
                mat2CellBinTimes{matCounter, 1} = BinTimes(convertedMat2Cell{matCounter, 1});
            end
            
            for cellCounter = 1:length(mat2CellLH_DeltaPower)
                matLH_DeltaPower = cell2mat(mat2CellLH_DeltaPower{cellCounter, 1});
                arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
                cellLH_DeltaPower{cellCounter, 1} = arrayLH_DeltaPower;
                
                matRH_DeltaPower = cell2mat(mat2CellRH_DeltaPower{cellCounter, 1});
                arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
                cellRH_DeltaPower{cellCounter, 1} = arrayRH_DeltaPower;
                
                matLH_ThetaPower = cell2mat(mat2CellLH_ThetaPower{cellCounter, 1});
                arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
                cellLH_ThetaPower{cellCounter, 1} = arrayLH_ThetaPower;
                
                matRH_ThetaPower = cell2mat(mat2CellRH_ThetaPower{cellCounter, 1});
                arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
                cellRH_ThetaPower{cellCounter, 1} = arrayRH_ThetaPower;
                
                matLH_AlphaPower = cell2mat(mat2CellLH_AlphaPower{cellCounter, 1});
                arrayLH_AlphaPower = reshape(matLH_AlphaPower', [1, size(matLH_AlphaPower, 2)*size(matLH_AlphaPower, 1)]);
                cellLH_AlphaPower{cellCounter, 1} = arrayLH_AlphaPower;
                
                matRH_AlphaPower = cell2mat(mat2CellRH_AlphaPower{cellCounter, 1});
                arrayRH_AlphaPower = reshape(matRH_AlphaPower', [1, size(matRH_AlphaPower, 2)*size(matRH_AlphaPower, 1)]);
                cellRH_AlphaPower{cellCounter, 1} = arrayRH_AlphaPower;
                
                matLH_BetaPower = cell2mat(mat2CellLH_BetaPower{cellCounter, 1});
                arrayLH_BetaPower = reshape(matLH_BetaPower', [1, size(matLH_BetaPower, 2)*size(matLH_BetaPower, 1)]);
                cellLH_BetaPower{cellCounter, 1} = arrayLH_BetaPower;
                
                matRH_BetaPower = cell2mat(mat2CellRH_BetaPower{cellCounter, 1});
                arrayRH_BetaPower = reshape(matRH_BetaPower', [1, size(matRH_BetaPower, 2)*size(matRH_BetaPower, 1)]);
                cellRH_BetaPower{cellCounter, 1} = arrayRH_BetaPower;
                
                matLH_GammaPower = cell2mat(mat2CellLH_GammaPower{cellCounter, 1});
                arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
                cellLH_GammaPower{cellCounter, 1} = arrayLH_GammaPower;
                
                matRH_GammaPower = cell2mat(mat2CellRH_GammaPower{cellCounter, 1});
                arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
                cellRH_GammaPower{cellCounter, 1} = arrayRH_GammaPower;
                
                matLH_MUAPower = cell2mat(mat2CellLH_MUAPower{cellCounter, 1});
                arrayLH_MUAPower = reshape(matLH_MUAPower', [1, size(matLH_MUAPower, 2)*size(matLH_MUAPower, 1)]);
                cellLH_MUAPower{cellCounter, 1} = arrayLH_MUAPower;
                
                matRH_MUAPower = cell2mat(mat2CellRH_MUAPower{cellCounter, 1});
                arrayRH_MUAPower = reshape(matRH_MUAPower', [1, size(matRH_MUAPower, 2)*size(matRH_MUAPower, 1)]);
                cellRH_MUAPower{cellCounter, 1} = arrayRH_MUAPower;
      
%                 matLH_Rhodamine = cell2mat(mat2CellLH_Rhodamine{cellCounter, 1});
%                 arrayLH_Rhodamine = reshape(matLH_Rhodamine', [1, size(matLH_Rhodamine, 2)*size(matLH_Rhodamine, 1)]);
%                 cellLH_Rhodamine{cellCounter, 1} = arrayLH_Rhodamine;
                
                matRH_Rhodamine = cell2mat(mat2CellRH_Rhodamine{cellCounter, 1});
                arrayRH_Rhodamine = reshape(matRH_Rhodamine', [1, size(matRH_Rhodamine, 2)*size(matRH_Rhodamine, 1)]);
                cellRH_Rhodamine{cellCounter, 1} = arrayRH_Rhodamine;
                
%                 matLH_GFP = cell2mat(mat2CellLH_GFP{cellCounter, 1});
%                 arrayLH_GFP = reshape(matLH_GFP', [1, size(matLH_GFP, 2)*size(matLH_GFP, 1)]);
%                 cellLH_GFP{cellCounter, 1} = arrayLH_GFP;
                
                matRH_GFP = cell2mat(mat2CellRH_GFP{cellCounter, 1});
                arrayRH_GFP = reshape(matRH_GFP', [1, size(matRH_GFP, 2)*size(matRH_GFP, 1)]);
                cellRH_GFP{cellCounter, 1} = arrayRH_GFP;
                
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
            for cellLength = 1:size(cellLH_DeltaPower, 2)   % Loop through however many sleep epochs this file has
                SleepData.(modelName).NREM.data.EEG_LH.deltaBandPower{cellLength, 1} = cellLH_DeltaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.deltaBandPower{cellLength, 1} = cellRH_DeltaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_LH.thetaBandPower{cellLength, 1} = cellLH_ThetaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.thetaBandPower{cellLength, 1} = cellRH_ThetaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_LH.alphaBandPower{cellLength, 1} = cellLH_AlphaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.alphaBandPower{cellLength, 1} = cellRH_AlphaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_LH.betaBandPower{cellLength, 1} = cellLH_BetaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.betaBandPower{cellLength, 1} = cellRH_BetaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_LH.gammaBandPower{cellLength, 1} = cellLH_GammaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.gammaBandPower{cellLength, 1} = cellRH_GammaPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_LH.EEGPower{cellLength, 1} = cellLH_MUAPower{1, 1};
                SleepData.(modelName).NREM.data.EEG_RH.EEGPower{cellLength, 1} = cellRH_MUAPower{1, 1};

%                 SleepData.(modelName).NREM.data.Rhodamine.LH{cellLength, 1} = cellLH_Rhodamine{1, 1};
                SleepData.(modelName).NREM.data.Rhodamine.RH{cellLength, 1} = cellRH_Rhodamine{1, 1};
                
%                 SleepData.(modelName).NREM.data.GFP.LH{cellLength, 1} = cellLH_GFP{1, 1};
                SleepData.(modelName).NREM.data.GFP.RH{cellLength, 1} = cellRH_GFP{1, 1};
                
                SleepData.(modelName).NREM.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.(modelName).NREM.FileIDs{cellLength, 1} = fileID;
                SleepData.(modelName).NREM.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellLH_DeltaPower, 1)   % Loop through however many sleep epochs this file has
                SleepData.(modelName).NREM.data.EEG_LH.deltaBandPower{size(SleepData.(modelName).NREM.data.EEG_LH.deltaBandPower, 1) + 1, 1} = cellLH_DeltaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.deltaBandPower{size(SleepData.(modelName).NREM.data.EEG_RH.deltaBandPower, 1) + 1, 1} = cellRH_DeltaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_LH.thetaBandPower{size(SleepData.(modelName).NREM.data.EEG_LH.thetaBandPower, 1) + 1, 1} = cellLH_ThetaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.thetaBandPower{size(SleepData.(modelName).NREM.data.EEG_RH.thetaBandPower, 1) + 1, 1} = cellRH_ThetaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_LH.alphaBandPower{size(SleepData.(modelName).NREM.data.EEG_LH.alphaBandPower, 1) + 1, 1} = cellLH_AlphaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.alphaBandPower{size(SleepData.(modelName).NREM.data.EEG_RH.alphaBandPower, 1) + 1, 1} = cellRH_AlphaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_LH.betaBandPower{size(SleepData.(modelName).NREM.data.EEG_LH.betaBandPower, 1) + 1, 1} = cellLH_BetaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.betaBandPower{size(SleepData.(modelName).NREM.data.EEG_RH.betaBandPower, 1) + 1, 1} = cellRH_BetaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_LH.gammaBandPower{size(SleepData.(modelName).NREM.data.EEG_LH.gammaBandPower, 1) + 1, 1} = cellLH_GammaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.gammaBandPower{size(SleepData.(modelName).NREM.data.EEG_RH.gammaBandPower, 1) + 1, 1} = cellRH_GammaPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_LH.EEGPower{size(SleepData.(modelName).NREM.data.EEG_LH.EEGPower, 1) + 1, 1} = cellLH_MUAPower{cellLength, 1};
                SleepData.(modelName).NREM.data.EEG_RH.EEGPower{size(SleepData.(modelName).NREM.data.EEG_RH.EEGPower, 1) + 1, 1} = cellRH_MUAPower{cellLength, 1};
               
%                 SleepData.(modelName).NREM.data.Rhodamine.LH{size(SleepData.(modelName).NREM.data.Rhodamine.LH, 1) + 1, 1} = cellLH_Rhodamine{cellLength, 1};
                SleepData.(modelName).NREM.data.Rhodamine.RH{size(SleepData.(modelName).NREM.data.Rhodamine.RH, 1) + 1, 1} = cellRH_Rhodamine{cellLength, 1};
                
%                 SleepData.(modelName).NREM.data.GFP.LH{size(SleepData.(modelName).NREM.data.GFP.LH, 1) + 1, 1} = cellLH_GFP{cellLength, 1};
                SleepData.(modelName).NREM.data.GFP.RH{size(SleepData.(modelName).NREM.data.GFP.RH, 1) + 1, 1} = cellRH_GFP{cellLength, 1};
 
                SleepData.(modelName).NREM.data.WhiskerAcceleration{size(SleepData.(modelName).NREM.data.WhiskerAcceleration, 1) + 1, 1} = cellWhiskerAcceleration{cellLength, 1};
                SleepData.(modelName).NREM.FileIDs{size(SleepData.(modelName).NREM.FileIDs, 1) + 1, 1} = fileID;
                SleepData.(modelName).NREM.BinTimes{size(SleepData.(modelName).NREM.BinTimes, 1) + 1, 1} = cellBinTimes{cellLength, 1};
            end
        end
    end
    
    disp(['Adding NREM sleeping epochs from ProcData file ' num2str(a) ' of ' num2str(size(procDataFileIDs, 1)) '...']); disp(' ')
end

%% REM
sleepBins = REMsleepTime/5;
for a = 1:size(procDataFileIDs, 1)           % Loop through the list of ProcData files
    procDataFileID = procDataFileIDs(a, :);    % Pull character string associated with the current file
    load(procDataFileID);                             % Load in procDataFile associated with character string
    [~,~,fileID] = GetFileInfo_FP(procDataFileID);     % Gather file info
    
    clear -regexp ^LH_ ^RH_ ^cell ^mat2Cell ^mat; 

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
            % filtered signal bands
            LH_deltaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.deltaBandPower{fixedSleepIndex(indexCount), 1};
            RH_deltaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.deltaBandPower{fixedSleepIndex(indexCount), 1};
            LH_thetaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.thetaBandPower{fixedSleepIndex(indexCount), 1};
            RH_thetaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.thetaBandPower{fixedSleepIndex(indexCount), 1};
            LH_alphaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.alphaBandPower{fixedSleepIndex(indexCount), 1};
            RH_alphaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.alphaBandPower{fixedSleepIndex(indexCount), 1};
            LH_betaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.betaBandPower{fixedSleepIndex(indexCount), 1};
            RH_betaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.betaBandPower{fixedSleepIndex(indexCount), 1};
            LH_gammaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.gammaBandPower{fixedSleepIndex(indexCount), 1};
            RH_gammaPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.gammaBandPower{fixedSleepIndex(indexCount), 1};
            LH_EEGPower{indexCount, 1} = ProcData.sleep.parameters.EEG_LH.EEGPower{fixedSleepIndex(indexCount), 1};
            RH_EEGPower{indexCount, 1} = ProcData.sleep.parameters.EEG_RH.EEGPower{fixedSleepIndex(indexCount), 1};
         
            % Rhodamine
%             LH_Rhodamine{indexCount, 1} = ProcData.sleep.parameters.Rhodamine.LH{fixedSleepIndex(indexCount), 1};
            RH_Rhodamine{indexCount, 1} = ProcData.sleep.parameters.Rhodamine.RH{fixedSleepIndex(indexCount), 1};
            
%             LH_GFP{indexCount, 1} = ProcData.sleep.parameters.GFP.LH{fixedSleepIndex(indexCount), 1};
            RH_GFP{indexCount, 1} = ProcData.sleep.parameters.GFP.RH{fixedSleepIndex(indexCount), 1};
 
            WhiskerAcceleration{indexCount, 1} = ProcData.sleep.parameters.whiskerAcceleration{fixedSleepIndex(indexCount), 1};
            BinTimes{indexCount, 1} = 5*fixedSleepIndex(indexCount);
        end
        
        indexBreaks = find(fixedSleepIndex(2:end) - fixedSleepIndex(1:end - 1) > 1);    % Find if there are numerous sleep periods
        
        if isempty(indexBreaks)   % If there is only one period of sleep in this file and not multiple
            % filtered signal bands
            matLH_DeltaPower = cell2mat(LH_deltaPower);
            arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
            cellLH_DeltaPower = {arrayLH_DeltaPower};
            
            matRH_DeltaPower = cell2mat(RH_deltaPower);
            arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
            cellRH_DeltaPower = {arrayRH_DeltaPower};
            
            matLH_ThetaPower = cell2mat(LH_thetaPower);
            arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
            cellLH_ThetaPower = {arrayLH_ThetaPower};
            
            matRH_ThetaPower = cell2mat(RH_thetaPower);
            arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
            cellRH_ThetaPower = {arrayRH_ThetaPower};
            
            matLH_AlphaPower = cell2mat(LH_alphaPower);
            arrayLH_AlphaPower = reshape(matLH_AlphaPower', [1, size(matLH_AlphaPower, 2)*size(matLH_AlphaPower, 1)]);
            cellLH_AlphaPower = {arrayLH_AlphaPower};
            
            matRH_AlphaPower = cell2mat(RH_alphaPower);
            arrayRH_AlphaPower = reshape(matRH_AlphaPower', [1, size(matRH_AlphaPower, 2)*size(matRH_AlphaPower, 1)]);
            cellRH_AlphaPower = {arrayRH_AlphaPower};
            
            matLH_BetaPower = cell2mat(LH_betaPower);
            arrayLH_BetaPower = reshape(matLH_BetaPower', [1, size(matLH_BetaPower, 2)*size(matLH_BetaPower, 1)]);
            cellLH_BetaPower = {arrayLH_BetaPower};
            
            matRH_BetaPower = cell2mat(RH_betaPower);
            arrayRH_BetaPower = reshape(matRH_BetaPower', [1, size(matRH_BetaPower, 2)*size(matRH_BetaPower, 1)]);
            cellRH_BetaPower = {arrayRH_BetaPower};
            
            matLH_GammaPower = cell2mat(LH_gammaPower);
            arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
            cellLH_GammaPower = {arrayLH_GammaPower};
            
            matRH_GammaPower = cell2mat(RH_gammaPower);
            arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
            cellRH_GammaPower = {arrayRH_GammaPower};
            
            matLH_MUAPower = cell2mat(LH_EEGPower);
            arrayLH_MUAPower = reshape(matLH_MUAPower', [1, size(matLH_MUAPower, 2)*size(matLH_MUAPower, 1)]);
            cellLH_MUAPower = {arrayLH_MUAPower};
            
            matRH_MUAPower = cell2mat(RH_EEGPower);
            arrayRH_MUAPower = reshape(matRH_MUAPower', [1, size(matRH_MUAPower, 2)*size(matRH_MUAPower, 1)]);
            cellRH_MUAPower = {arrayRH_MUAPower};
           
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
            
            % Rhodamine
%             matLH_Rhodamine = cell2mat(LH_Rhodamine);
%             arrayLH_Rhodamine = reshape(matLH_Rhodamine', [1, size(matLH_Rhodamine, 2)*size(matLH_Rhodamine, 1)]);
%             cellLH_Rhodamine = {arrayLH_Rhodamine};
            
            matRH_Rhodamine = cell2mat(RH_Rhodamine);
            arrayRH_Rhodamine = reshape(matRH_Rhodamine', [1, size(matRH_Rhodamine, 2)*size(matRH_Rhodamine, 1)]);
            cellRH_Rhodamine = {arrayRH_Rhodamine};
            
            % GFP
%             matLH_GFP = cell2mat(LH_GFP);
%             arrayLH_GFP = reshape(matLH_GFP', [1, size(matLH_GFP, 2)*size(matLH_GFP, 1)]);
%             cellLH_GFP = {arrayLH_GFP};
            
            matRH_GFP = cell2mat(RH_GFP);
            arrayRH_GFP = reshape(matRH_GFP', [1, size(matRH_GFP, 2)*size(matRH_GFP, 1)]);
            cellRH_GFP = {arrayRH_GFP};

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
            
            splitCounter = 1:length(LH_deltaPower);
            convertedMat2Cell = mat2cell(splitCounter', holdIndex);
            
            for matCounter = 1:length(convertedMat2Cell)
                mat2CellLH_DeltaPower{matCounter, 1} = LH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_DeltaPower{matCounter, 1} = RH_deltaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_ThetaPower{matCounter, 1} = LH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_ThetaPower{matCounter, 1} = RH_thetaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_AlphaPower{matCounter, 1} = LH_alphaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_AlphaPower{matCounter, 1} = RH_alphaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_BetaPower{matCounter, 1} = LH_betaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_BetaPower{matCounter, 1} = RH_betaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_GammaPower{matCounter, 1} = LH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_GammaPower{matCounter, 1} = RH_gammaPower(convertedMat2Cell{matCounter, 1});
                mat2CellLH_MUAPower{matCounter, 1} = LH_EEGPower(convertedMat2Cell{matCounter, 1});
                mat2CellRH_MUAPower{matCounter, 1} = RH_EEGPower(convertedMat2Cell{matCounter, 1});
                
%                 mat2CellLH_Rhodamine{matCounter, 1} = LH_Rhodamine(convertedMat2Cell{matCounter, 1});
                mat2CellRH_Rhodamine{matCounter, 1} = RH_Rhodamine(convertedMat2Cell{matCounter, 1});
                
%                 mat2CellLH_GFP{matCounter, 1} = LH_GFP(convertedMat2Cell{matCounter, 1});
                mat2CellRH_GFP{matCounter, 1} = RH_GFP(convertedMat2Cell{matCounter, 1});

                mat2CellWhiskerAcceleration{matCounter, 1} = WhiskerAcceleration(convertedMat2Cell{matCounter, 1});
                mat2CellBinTimes{matCounter, 1} = BinTimes(convertedMat2Cell{matCounter, 1});
            end
            
            for cellCounter = 1:length(mat2CellLH_DeltaPower)
                matLH_DeltaPower = cell2mat(mat2CellLH_DeltaPower{cellCounter, 1});
                arrayLH_DeltaPower = reshape(matLH_DeltaPower', [1, size(matLH_DeltaPower, 2)*size(matLH_DeltaPower, 1)]);
                cellLH_DeltaPower{cellCounter, 1} = arrayLH_DeltaPower;
                
                matRH_DeltaPower = cell2mat(mat2CellRH_DeltaPower{cellCounter, 1});
                arrayRH_DeltaPower = reshape(matRH_DeltaPower', [1, size(matRH_DeltaPower, 2)*size(matRH_DeltaPower, 1)]);
                cellRH_DeltaPower{cellCounter, 1} = arrayRH_DeltaPower;
                
                matLH_ThetaPower = cell2mat(mat2CellLH_ThetaPower{cellCounter, 1});
                arrayLH_ThetaPower = reshape(matLH_ThetaPower', [1, size(matLH_ThetaPower, 2)*size(matLH_ThetaPower, 1)]);
                cellLH_ThetaPower{cellCounter, 1} = arrayLH_ThetaPower;
                
                matRH_ThetaPower = cell2mat(mat2CellRH_ThetaPower{cellCounter, 1});
                arrayRH_ThetaPower = reshape(matRH_ThetaPower', [1, size(matRH_ThetaPower, 2)*size(matRH_ThetaPower, 1)]);
                cellRH_ThetaPower{cellCounter, 1} = arrayRH_ThetaPower;
                
                matLH_AlphaPower = cell2mat(mat2CellLH_AlphaPower{cellCounter, 1});
                arrayLH_AlphaPower = reshape(matLH_AlphaPower', [1, size(matLH_AlphaPower, 2)*size(matLH_AlphaPower, 1)]);
                cellLH_AlphaPower{cellCounter, 1} = arrayLH_AlphaPower;
                
                matRH_AlphaPower = cell2mat(mat2CellRH_AlphaPower{cellCounter, 1});
                arrayRH_AlphaPower = reshape(matRH_AlphaPower', [1, size(matRH_AlphaPower, 2)*size(matRH_AlphaPower, 1)]);
                cellRH_AlphaPower{cellCounter, 1} = arrayRH_AlphaPower;
                
                matLH_BetaPower = cell2mat(mat2CellLH_BetaPower{cellCounter, 1});
                arrayLH_BetaPower = reshape(matLH_BetaPower', [1, size(matLH_BetaPower, 2)*size(matLH_BetaPower, 1)]);
                cellLH_BetaPower{cellCounter, 1} = arrayLH_BetaPower;
                
                matRH_BetaPower = cell2mat(mat2CellRH_BetaPower{cellCounter, 1});
                arrayRH_BetaPower = reshape(matRH_BetaPower', [1, size(matRH_BetaPower, 2)*size(matRH_BetaPower, 1)]);
                cellRH_BetaPower{cellCounter, 1} = arrayRH_BetaPower;
                
                matLH_GammaPower = cell2mat(mat2CellLH_GammaPower{cellCounter, 1});
                arrayLH_GammaPower = reshape(matLH_GammaPower', [1, size(matLH_GammaPower, 2)*size(matLH_GammaPower, 1)]);
                cellLH_GammaPower{cellCounter, 1} = arrayLH_GammaPower;
                
                matRH_GammaPower = cell2mat(mat2CellRH_GammaPower{cellCounter, 1});
                arrayRH_GammaPower = reshape(matRH_GammaPower', [1, size(matRH_GammaPower, 2)*size(matRH_GammaPower, 1)]);
                cellRH_GammaPower{cellCounter, 1} = arrayRH_GammaPower;
                
                matLH_MUAPower = cell2mat(mat2CellLH_MUAPower{cellCounter, 1});
                arrayLH_MUAPower = reshape(matLH_MUAPower', [1, size(matLH_MUAPower, 2)*size(matLH_MUAPower, 1)]);
                cellLH_MUAPower{cellCounter, 1} = arrayLH_MUAPower;
                
                matRH_MUAPower = cell2mat(mat2CellRH_MUAPower{cellCounter, 1});
                arrayRH_MUAPower = reshape(matRH_MUAPower', [1, size(matRH_MUAPower, 2)*size(matRH_MUAPower, 1)]);
                cellRH_MUAPower{cellCounter, 1} = arrayRH_MUAPower;

%                 matLH_Rhodamine = cell2mat(mat2CellLH_Rhodamine{cellCounter, 1});
%                 arrayLH_Rhodamine = reshape(matLH_Rhodamine', [1, size(matLH_Rhodamine, 2)*size(matLH_Rhodamine, 1)]);
%                 cellLH_Rhodamine{cellCounter, 1} = arrayLH_Rhodamine;
                
                matRH_Rhodamine = cell2mat(mat2CellRH_Rhodamine{cellCounter, 1});
                arrayRH_Rhodamine = reshape(matRH_Rhodamine', [1, size(matRH_Rhodamine, 2)*size(matRH_Rhodamine, 1)]);
                cellRH_Rhodamine{cellCounter, 1} = arrayRH_Rhodamine;
                
%                 matLH_GFP = cell2mat(mat2CellLH_GFP{cellCounter, 1});
%                 arrayLH_GFP = reshape(matLH_GFP', [1, size(matLH_GFP, 2)*size(matLH_GFP, 1)]);
%                 cellLH_GFP{cellCounter, 1} = arrayLH_GFP;
                
                matRH_GFP = cell2mat(mat2CellRH_GFP{cellCounter, 1});
                arrayRH_GFP = reshape(matRH_GFP', [1, size(matRH_GFP, 2)*size(matRH_GFP, 1)]);
                cellRH_GFP{cellCounter, 1} = arrayRH_GFP;
                
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
        if isfield(SleepData.(modelName),'REM') == false % If the structure is empty, we need a special case to format the struct properly
            for cellLength = 1:size(cellLH_DeltaPower, 2)   % Loop through however many sleep epochs this file has
                SleepData.(modelName).REM.data.EEG_LH.deltaBandPower{cellLength, 1} = cellLH_DeltaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.deltaBandPower{cellLength, 1} = cellRH_DeltaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_LH.thetaBandPower{cellLength, 1} = cellLH_ThetaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.thetaBandPower{cellLength, 1} = cellRH_ThetaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_LH.alphaBandPower{cellLength, 1} = cellLH_AlphaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.alphaBandPower{cellLength, 1} = cellRH_AlphaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_LH.betaBandPower{cellLength, 1} = cellLH_BetaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.betaBandPower{cellLength, 1} = cellRH_BetaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_LH.gammaBandPower{cellLength, 1} = cellLH_GammaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.gammaBandPower{cellLength, 1} = cellRH_GammaPower{1, 1};
                SleepData.(modelName).REM.data.EEG_LH.EEGPower{cellLength, 1} = cellLH_MUAPower{1, 1};
                SleepData.(modelName).REM.data.EEG_RH.EEGPower{cellLength, 1} = cellRH_MUAPower{1, 1};

%                 SleepData.(modelName).REM.data.Rhodamine.LH{cellLength, 1} = cellLH_Rhodamine{1, 1};
                SleepData.(modelName).REM.data.Rhodamine.RH{cellLength, 1} = cellRH_Rhodamine{1, 1};
                
%                 SleepData.(modelName).REM.data.GFP.LH{cellLength, 1} = cellLH_GFP{1, 1};
                SleepData.(modelName).REM.data.GFP.RH{cellLength, 1} = cellRH_GFP{1, 1};

                SleepData.(modelName).REM.data.WhiskerAcceleration{cellLength, 1} = cellWhiskerAcceleration{1, 1};
                SleepData.(modelName).REM.FileIDs{cellLength, 1} = fileID;
                SleepData.(modelName).REM.BinTimes{cellLength, 1} = cellBinTimes{1, 1};
            end
        else    % If the struct is not empty, add each new iteration after previous data
            for cellLength = 1:size(cellLH_DeltaPower, 1)   % Loop through however many sleep epochs this file has
                SleepData.(modelName).REM.data.EEG_LH.deltaBandPower{size(SleepData.(modelName).REM.data.EEG_LH.deltaBandPower, 1) + 1, 1} = cellLH_DeltaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.deltaBandPower{size(SleepData.(modelName).REM.data.EEG_RH.deltaBandPower, 1) + 1, 1} = cellRH_DeltaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_LH.thetaBandPower{size(SleepData.(modelName).REM.data.EEG_LH.thetaBandPower, 1) + 1, 1} = cellLH_ThetaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.thetaBandPower{size(SleepData.(modelName).REM.data.EEG_RH.thetaBandPower, 1) + 1, 1} = cellRH_ThetaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_LH.alphaBandPower{size(SleepData.(modelName).REM.data.EEG_LH.alphaBandPower, 1) + 1, 1} = cellLH_AlphaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.alphaBandPower{size(SleepData.(modelName).REM.data.EEG_RH.alphaBandPower, 1) + 1, 1} = cellRH_AlphaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_LH.betaBandPower{size(SleepData.(modelName).REM.data.EEG_LH.betaBandPower, 1) + 1, 1} = cellLH_BetaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.betaBandPower{size(SleepData.(modelName).REM.data.EEG_RH.betaBandPower, 1) + 1, 1} = cellRH_BetaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_LH.gammaBandPower{size(SleepData.(modelName).REM.data.EEG_LH.gammaBandPower, 1) + 1, 1} = cellLH_GammaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.gammaBandPower{size(SleepData.(modelName).REM.data.EEG_RH.gammaBandPower, 1) + 1, 1} = cellRH_GammaPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_LH.EEGPower{size(SleepData.(modelName).REM.data.EEG_LH.EEGPower, 1) + 1, 1} = cellLH_MUAPower{cellLength, 1};
                SleepData.(modelName).REM.data.EEG_RH.EEGPower{size(SleepData.(modelName).REM.data.EEG_RH.EEGPower, 1) + 1, 1} = cellRH_MUAPower{cellLength, 1};
               
%                 SleepData.(modelName).REM.data.Rhodamine.LH{size(SleepData.(modelName).REM.data.Rhodamine.LH, 1) + 1, 1} = cellLH_Rhodamine{cellLength, 1};
                SleepData.(modelName).REM.data.Rhodamine.RH{size(SleepData.(modelName).REM.data.Rhodamine.RH, 1) + 1, 1} = cellRH_Rhodamine{cellLength, 1};
                
%                 SleepData.(modelName).REM.data.GFP.LH{size(SleepData.(modelName).REM.data.GFP.LH, 1) + 1, 1} = cellLH_GFP{cellLength, 1};
                SleepData.(modelName).REM.data.GFP.RH{size(SleepData.(modelName).REM.data.GFP.RH, 1) + 1, 1} = cellRH_GFP{cellLength, 1};
 
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
