function [] = AddMicroArousalsParameters_FP_GRABNE_SingleFiber(procDataFileIDs,RestingBaselines,baselineType,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Organize data into appropriate bins for microarousal scoring characterization
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding microarousal scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~, fileDate, ~] = GetFileInfo_FP(procDataFileID);
    strDay = ConvertDate_FP(fileDate);
    load(procDataFileID)
%     if isfield(ProcData,"Microarousals") == true
%         ProcData = rmfield(ProcData,"Microarousals");
%     end
    specDataFileID = [procDataFileID(1:end-12) 'SpecDataB.mat'];
    load(specDataFileID)
    %% BLOCK PURPOSE: Create folder for the Neural data of each electrode

    % cortical delta
    LH_Delta = ProcData.data.cortical_LH.deltaBandPower;
    LH_baselineDelta = RestingBaselines.(baselineType).cortical_LH.deltaBandPower.(strDay).mean;
    LH_DeltaNeuro = (LH_Delta-LH_baselineDelta)/LH_baselineDelta;  

    % cortical theta
    LH_Theta = ProcData.data.cortical_LH.thetaBandPower;
    LH_baselineTheta = RestingBaselines.(baselineType).cortical_LH.thetaBandPower.(strDay).mean;
    LH_ThetaNeuro = (LH_Theta-LH_baselineTheta)/LH_baselineTheta;
  
    % cortical alpha
    LH_Alpha = ProcData.data.cortical_LH.alphaBandPower;
    LH_baselineAlpha = RestingBaselines.(baselineType).cortical_LH.alphaBandPower.(strDay).mean;
    LH_AlphaNeuro = (LH_Alpha-LH_baselineAlpha)/LH_baselineAlpha;   
 
    % cortical beta
    LH_Beta = ProcData.data.cortical_LH.betaBandPower;
    LH_baselineBeta = RestingBaselines.(baselineType).cortical_LH.betaBandPower.(strDay).mean;
    LH_BetaNeuro = (LH_Beta-LH_baselineBeta)/LH_baselineBeta;

    % cortical gamma
    LH_Gamma = ProcData.data.cortical_LH.gammaBandPower;
    LH_baselineGamma = RestingBaselines.(baselineType).cortical_LH.gammaBandPower.(strDay).mean;
    LH_GammaNeuro = (LH_Gamma-LH_baselineGamma)/LH_baselineGamma;
    % cortical Cortical
    LH_Cortical = ProcData.data.cortical_LH.corticalPower;
    LH_baselineCortical = RestingBaselines.(baselineType).cortical_LH.corticalPower.(strDay).mean;
    LH_CorticalNeuro = (LH_Cortical-LH_baselineCortical)/LH_baselineCortical;
    % Divide the neural signals into five second bins and put them in a cell array

    LH_tempDeltaStruct = cell(NBins,1);
    LH_tempThetaStruct = cell(NBins,1);
    LH_tempAlphaStruct = cell(NBins,1);
    LH_tempBetaStruct = cell(NBins,1);
    LH_tempGammaStruct = cell(NBins,1);
    LH_tempCorticalStruct = cell(NBins,1);
    % loop through all samples across the full duration in 5 second bins (NBins total)
    for b = 1:NBins
        if b == 1

            % cortical
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro(b:30)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro(b:30)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro(b:30)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro(b:30)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro(b:30)};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro(b:30)};
        elseif b == NBins

            % cortical
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((30*(b - 1)) + 1)):end)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((30*(b - 1)) + 1)):end)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((30*(b - 1)) + 1)):end)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((30*(b - 1)) + 1)):end)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((30*(b - 1)) + 1)):end)};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro((((30*(b - 1)) + 1)):end)};
        else

            % cortical
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((30*(b - 1)) + 1)):(30*b))};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((30*(b - 1)) + 1)):(30*b))};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((30*(b - 1)) + 1)):(30*b))};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((30*(b - 1)) + 1)):(30*b))};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((30*(b - 1)) + 1)):(30*b))};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro((((30*(b - 1)) + 1)):(30*b))};
        end
    end

    % save cortical data under ProcData file
    ProcData.MicroArousals.parameters.cortical_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.MicroArousals.parameters.cortical_LH.thetaBandPower = LH_tempThetaStruct;
    ProcData.MicroArousals.parameters.cortical_LH.alphaBandPower = LH_tempAlphaStruct;
    ProcData.MicroArousals.parameters.cortical_LH.betaBandPower = LH_tempBetaStruct;
    ProcData.MicroArousals.parameters.cortical_LH.gammaBandPower = LH_tempGammaStruct;
    ProcData.MicroArousals.parameters.cortical_LH.corticalPower = LH_tempCorticalStruct;
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of each electrode
    trialDuration_sec = ProcData.notes.trialDuration_sec;   % sec
    offset = 0.10;   % sec
    binWidth = 1;   % sec
    T = round(SpecData.cortical_LH.T,1);
    F = SpecData.cortical_LH.F;
%     specLH = SpecData.cortical_LH.normS;
    specLH = SpecData.cortical_LH.normS;
    freqFloor = floor(F);
    % delta
    deltaLow = freqFloor == 1;%1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow,1,'first');
    deltaLowEnd = find(deltaHigh,1,'last');
    deltaSpecLH = specLH(deltaLowStart:deltaLowEnd,:);
    meanDeltaSpecLH = mean(deltaSpecLH,1);
    % theta
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 9;%10;
    thetaLowStart = find(thetaLow,1,'first');
    thetaLowEnd = find(thetaHigh,1,'last');
    thetaSpecLH = specLH(thetaLowStart:thetaLowEnd,:);
    meanThetaSpecLH = mean(thetaSpecLH,1);
    % alpha
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow,1,'first');
    alphaLowEnd = find(alphaHigh,1,'last');
    alphaSpecLH = specLH(alphaLowStart:alphaLowEnd,:);
    meanAlphaSpecLH = mean(alphaSpecLH,1);
    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow,1,'first');
    betaLowEnd = find(betaHigh,1,'last');
    betaSpecLH = specLH(betaLowStart:betaLowEnd,:);
    meanBetaSpecLH = mean(betaSpecLH,1);
    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow,1,'first');
    gammaLowEnd = find(gammaHigh,1,'last');
    gammaSpecLH = specLH(gammaLowStart:gammaLowEnd,:);
    meanGammaSpecLH = mean(gammaSpecLH,1);
    % Divide the neural signals into five second bins and put them in a cell array

    LH_tempDeltaSpecStruct = cell(NBins,1);
    LH_tempThetaSpecStruct = cell(NBins,1);
    LH_tempAlphaSpecStruct = cell(NBins,1);
    LH_tempBetaSpecStruct = cell(NBins,1);
    LH_tempGammaSpecStruct = cell(NBins,1);
    % loop through all samples across the 15 minutes in 5 second bins (NBins total)
    for c = 1:NBins
        if c == 1
            startTime = offset;
            startTime_index = find(T == startTime);
            endTime = 1;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
        elseif c == NBins
            startTime = trialDuration_sec - 1;
            [~,startTime_index] = min(abs(T - startTime));
            endTime = trialDuration_sec - offset;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
        else
            startTime = binWidth*(c - 1);
            [~,startTime_index] = min(abs(T - startTime));
            endTime = binWidth*c;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index + 1:endTime_index + 1)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index + 1:endTime_index + 1)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index + 1:endTime_index + 1)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index + 1:endTime_index + 1)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index + 1:endTime_index + 1)};
        end
    end

    % save cortical data under ProcData file
    ProcData.MicroArousals.parameters.cortical_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.MicroArousals.parameters.cortical_LH.specThetaBandPower = LH_tempThetaSpecStruct;
    ProcData.MicroArousals.parameters.cortical_LH.specAlphaBandPower = LH_tempAlphaSpecStruct;
    ProcData.MicroArousals.parameters.cortical_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.MicroArousals.parameters.cortical_LH.specGammaBandPower = LH_tempGammaSpecStruct;
    
    %% BLOCK PURPOSE: Create folder for binarized whisking and binarized force sensor
    binWhiskerAngle = ProcData.data.binWhiskerAngle;
    binForceSensor = ProcData.data.binForceSensor;
    ForceSensor = ProcData.data.forceSensor;
    whiskerAngle = ProcData.data.whiskerAngle;
    whiskerAcceleration = diff(whiskerAngle,2);
    % Find the number of whiskerBins due to frame drops.
    whiskerBinNumber = NBins;
    % Divide the signal into five second bins and put them in a cell array
    tempWhiskerStruct = cell(whiskerBinNumber,1);
    tempWhiskerAccelStruct = cell(whiskerBinNumber,1);
    tempBinWhiskerStruct = cell(whiskerBinNumber,1);
    tempForceStruct = cell(whiskerBinNumber,1);
    tempForceRStruct = cell(whiskerBinNumber,1);
    for whiskerBins = 1:whiskerBinNumber
        if whiskerBins == 1
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle(whiskerBins:30)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration(whiskerBins:30)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:30)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:30)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor(whiskerBins:30)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((30*(whiskerBins-1)) + 1)):end)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((30*(whiskerBins-1)) + 1)):end)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((30*(whiskerBins-1)) + 1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((30*(whiskerBins-1)) + 1)):end)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((30*(whiskerBins-1)) + 1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((30*(whiskerBins-1)) + 1)):(30*whiskerBins))};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((30*(whiskerBins-1)) + 1)):(30*whiskerBins))};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((30*(whiskerBins-1)) + 1)):(30*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((30*(whiskerBins-1)) + 1)):(30*whiskerBins))};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((30*(whiskerBins-1)) + 1)):(30*whiskerBins))};
        end
    end
    % save whisker and force sensor data under ProcData file
    ProcData.MicroArousals.parameters.whiskerAngle = tempWhiskerStruct;
    ProcData.MicroArousals.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
    ProcData.MicroArousals.parameters.binWhiskerAngle = tempBinWhiskerStruct;
    ProcData.MicroArousals.parameters.binForceSensor = tempForceStruct;
    ProcData.MicroArousals.parameters.ForceSensor = tempForceRStruct;
    %% add pupil parameters
    ProcData.data.Pupil.mmPerPixel = 0.018; % used to convert pixels to mm
        if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
        % create fields for the data
        dataTypes = {'pupilArea','Diameter','mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
            for aa = 1:length(dataTypes)
                dataType = dataTypes{1,aa};
                samplingRate = ProcData.notes.dsFs;
                [z,p,k] = butter(4,1/(samplingRate/2),'low');
                [sos,g] = zp2sos(z,p,k);
                if strcmp(dataType,'eyeMotion') == true
                    [z,p,k] = butter(4,10/(samplingRate/2),'low');
                    [sos,g] = zp2sos(z,p,k);
                    try
                        centroidX = filtfilt(sos,g,ProcData.data.Pupil.CentroidX - RestingBaselines.setDuration.Pupil.CentroidX.(strDay).mean);
                        centroidY = filtfilt(sos,g,ProcData.data.Pupil.CentroidY - RestingBaselines.setDuration.Pupil.CentroidY.(strDay).mean);
                    catch
                        centroidX = ProcData.data.Pupil.CentroidX - RestingBaselines.setDuration.Pupil.CentroidX.(strDay).mean;
                        centroidY = ProcData.data.Pupil.CentroidY - RestingBaselines.setDuration.Pupil.CentroidY.(strDay).mean;
                    end

                    % calculate the distance travelled by pupil
                    for dd = 1:length(centroidX)
                        if dd == 1
                            distanceTraveled(1,dd) = 0;
                        else
                            xPt = [centroidX(1,dd - 1),centroidY(1,dd - 1)];
                            yPt = [centroidX(1,dd),centroidY(1,dd)];
                            distanceTraveled(1,dd) = pdist([xPt;yPt],'euclidean');
                        end
                    end

                    data.(dataType).data = distanceTraveled*ProcData.data.Pupil.mmPerPixel;
                    ProcData.data.Pupil.distanceTraveled = distanceTraveled*ProcData.data.Pupil.mmPerPixel;
                elseif strcmp(dataType,'CentroidX') == true
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean);
                    catch
                        data.(dataType).data = ProcData.data.Pupil.CentroidX - RestingBaselines.manualSelection.Pupil.CentroidX.(strDay).mean;
                    end
                elseif strcmp(dataType,'CentroidY') == true
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean);
                    catch
                        data.(dataType).data = ProcData.data.Pupil.CentroidY - RestingBaselines.manualSelection.Pupil.CentroidY.(strDay).mean;
                    end
                elseif strcmp(dataType,'whiskerMotion')
                    data.(dataType).data = ProcData.data.whiskerAngle.^2;
                else
                    try
                        data.(dataType).data = filtfilt(sos,g,ProcData.data.Pupil.(dataType));
                    catch
                        data.(dataType).data = ProcData.data.Pupil.(dataType);
                    end
                end
                data.(dataType).struct = cell(NBins,1);
                % loop through all samples across the 15 minutes in 5 second bins (NBins total)
                for b = 1:NBins
                    if b == 1
                        data.(dataType).struct(b,1) = {data.(dataType).data(b:30)};
                    elseif b == NBins
                        data.(dataType).struct(b,1) = {data.(dataType).data((((30*(b - 1)) + 1)):(30*b))};%end)};
                    else
                        data.(dataType).struct(b,1) = {data.(dataType).data((((30*(b - 1)) + 1)):(30*b))};
                    end
                end
                ProcData.MicroArousals.parameters.Pupil.(dataType) = data.(dataType).struct;
            end
        end
    %% Create folder for the EMG
    EMG = ProcData.data.EMG.emg';
    normEMG = EMG - RestingBaselines.(baselineType).EMG.emg.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:30)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((30*(EMGBins-1)) + 1)):(30*EMGBins))};
        end
    end
    ProcData.MicroArousals.parameters.EMG.emg = tempEMGStruct; 
    %% raw EMG Signal
    EMG = ProcData.data.EMG.emgSignal';
    normEMG = EMG - RestingBaselines.(baselineType).EMG.emgSignal.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:30)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((30*(EMGBins-1)) + 1)):(30*EMGBins))};
        end
    end
    % save EMG data under ProcData file
%     ProcData.MicroArousals.parameters = rmfield(ProcData.MicroArousals.parameters,'EMG');
    ProcData.MicroArousals.parameters.EMG.emgSignal = tempEMGStruct;  
    %% BLOCK PURPOSE: Create folder for the left and right Rhodamine data
%     NE_GFP = ProcData.data.GFP.NE; 
%     NE_Rhodamine = ProcData.data.Rhodamine.NE;
% 
%     NE_NormRhodamine = (NE_Rhodamine - RestingBaselines.(baselineType).Rhodamine.NE.(strDay).mean)/RestingBaselines.(baselineType).Rhodamine.NE.(strDay).std;
%     NE_NormGFP = (NE_GFP - RestingBaselines.(baselineType).GFP.NE.(strDay).mean)/RestingBaselines.(baselineType).GFP.NE.(strDay).std;
% 
%     NE_tempGFPStruct = cell(NBins,1);  
%     RhodamineNE_tempRhodamineStruct = cell(NBins,1);
% 
%     for RhodamineBins = 1:NBins
%         if RhodamineBins == 1
%             NE_tempGFPStruct(RhodamineBins,1) = {NE_NormGFP(RhodamineBins:30)};
% 
%             RhodamineNE_tempRhodamineStruct(RhodamineBins,1) = {NE_NormRhodamine(RhodamineBins:30)};
%         else
%             NE_tempGFPStruct(RhodamineBins,1) = {NE_NormGFP((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};
%             RhodamineNE_tempRhodamineStruct(RhodamineBins,1) = {NE_NormRhodamine((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};
% 
%         end
%     end
%     % save hemodynamic data under ProcData file
%     ProcData.MicroArousals.parameters.GFP.NE = NE_tempGFPStruct;
%     ProcData.MicroArousals.parameters.Rhodamine.NE = RhodamineNE_tempRhodamineStruct;
    %% BLOCK PURPOSE: Create folder for the left and right fiber data (z Scored)
    Z_NE_GFP = ProcData.data.GFP.Z_NE; 
    Z_NE_Rhodamine = ProcData.data.Rhodamine.Z_NE;

    Z_NE_NormRhodamine = (Z_NE_Rhodamine - RestingBaselines.(baselineType).Rhodamine.Z_NE.(strDay).mean)/RestingBaselines.(baselineType).Rhodamine.Z_NE.(strDay).std;
    Z_NE_NormGFP = (Z_NE_GFP - RestingBaselines.(baselineType).GFP.Z_NE.(strDay).mean)/RestingBaselines.(baselineType).GFP.Z_NE.(strDay).std;

    Z_NE_tempGFPStruct = cell(NBins,1);  
    RhodamineZ_NE_tempRhodamineStruct = cell(NBins,1);

    for RhodamineBins = 1:NBins
        if RhodamineBins == 1
            Z_NE_tempGFPStruct(RhodamineBins,1) = {Z_NE_NormGFP(RhodamineBins:30)};

            RhodamineZ_NE_tempRhodamineStruct(RhodamineBins,1) = {Z_NE_NormRhodamine(RhodamineBins:30)};
        else
            Z_NE_tempGFPStruct(RhodamineBins,1) = {Z_NE_NormGFP((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};
            RhodamineZ_NE_tempRhodamineStruct(RhodamineBins,1) = {Z_NE_NormRhodamine((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};

        end
    end
    % save hemodynamic data under ProcData file
    ProcData.MicroArousals.parameters.GFP.Z_NE = Z_NE_tempGFPStruct;
    ProcData.MicroArousals.parameters.Rhodamine.Z_NE = RhodamineZ_NE_tempRhodamineStruct;
        %% BLOCK PURPOSE: Create folder for the left and right fiber data (z Scored)
        if  isfield(ProcData.data.GFP,'Z_Ach') == true
            Z_Ach_GFP = ProcData.data.GFP.Z_Ach; 
            Z_Ach_Rhodamine = ProcData.data.Rhodamine.Z_Ach;
        
            Z_Ach_NormRhodamine = (Z_Ach_Rhodamine - RestingBaselines.(baselineType).Rhodamine.Z_Ach.(strDay).mean)/RestingBaselines.(baselineType).Rhodamine.Z_Ach.(strDay).std;
            Z_Ach_NormGFP = (Z_Ach_GFP - RestingBaselines.(baselineType).GFP.Z_Ach.(strDay).mean)/RestingBaselines.(baselineType).GFP.Z_Ach.(strDay).std;
        
            Z_Ach_tempGFPStruct = cell(NBins,1);  
            RhodamineZ_Ach_tempRhodamineStruct = cell(NBins,1);
        
            for RhodamineBins = 1:NBins
                if RhodamineBins == 1
                    Z_Ach_tempGFPStruct(RhodamineBins,1) = {Z_Ach_NormGFP(RhodamineBins:30)};
        
                    RhodamineZ_Ach_tempRhodamineStruct(RhodamineBins,1) = {Z_Ach_NormRhodamine(RhodamineBins:30)};
                else
                    Z_Ach_tempGFPStruct(RhodamineBins,1) = {Z_Ach_NormGFP((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};
                    RhodamineZ_Ach_tempRhodamineStruct(RhodamineBins,1) = {Z_Ach_NormRhodamine((((30*(RhodamineBins-1)) + 1)):(30*RhodamineBins))};
        
                end
            end
            % save hemodynamic data under ProcData file
            ProcData.MicroArousals.parameters.GFP.Z_Ach = Z_Ach_tempGFPStruct;
            ProcData.MicroArousals.parameters.Rhodamine.Z_Ach = RhodamineZ_Ach_tempRhodamineStruct;
        end
    %%
    save(procDataFileID,'ProcData');
end

end
