function [] = AddSleepParameters_FP_GRABNE_SingleFiber(procDataFileIDs,RestingBaselines,baselineType,NBins)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Organize data into appropriate bins for sleep scoring characterization
%________________________________________________________________________________________________________________________

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    disp(['Adding sleep scoring parameters to ' procDataFileID '... (' num2str(a) '/' num2str(size(procDataFileIDs,1)) ')']); disp(' ')
    [~, fileDate, ~] = GetFileInfo_FP(procDataFileID);
    strDay = ConvertDate_FP(fileDate);
    load(procDataFileID)
    specDataFileID = [procDataFileID(1:end-12) 'SpecDataA.mat'];
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
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro(b:150)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro(b:150)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro(b:150)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro(b:150)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro(b:150)};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro(b:150)};
        elseif b == NBins

            % cortical
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((150*(b - 1)) + 1)):end)};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro((((150*(b - 1)) + 1)):end)};
        else

            % cortical
            LH_tempDeltaStruct(b,1) = {LH_DeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempThetaStruct(b,1) = {LH_ThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempAlphaStruct(b,1) = {LH_AlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempBetaStruct(b,1) = {LH_BetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempGammaStruct(b,1) = {LH_GammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            LH_tempCorticalStruct(b,1) = {LH_CorticalNeuro((((150*(b - 1)) + 1)):(150*b))};
        end
    end

    % save cortical data under ProcData file
    ProcData.sleep.parameters.cortical_LH.deltaBandPower = LH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_LH.thetaBandPower = LH_tempThetaStruct;
    ProcData.sleep.parameters.cortical_LH.alphaBandPower = LH_tempAlphaStruct;
    ProcData.sleep.parameters.cortical_LH.betaBandPower = LH_tempBetaStruct;
    ProcData.sleep.parameters.cortical_LH.gammaBandPower = LH_tempGammaStruct;
    ProcData.sleep.parameters.cortical_LH.corticalPower = LH_tempCorticalStruct;
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of each electrode
    trialDuration_sec = ProcData.notes.trialDuration_sec;   % sec
    offset = 2.5;   % sec
    binWidth = 5;   % sec
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
            endTime = 5;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            LH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecLH(startTime_index:endTime_index)};
            LH_tempThetaSpecStruct{c,1} = {meanThetaSpecLH(startTime_index:endTime_index)};
            LH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecLH(startTime_index:endTime_index)};
            LH_tempBetaSpecStruct{c,1} = {meanBetaSpecLH(startTime_index:endTime_index)};
            LH_tempGammaSpecStruct{c,1} = {meanGammaSpecLH(startTime_index:endTime_index)};
        elseif c == NBins
            startTime = trialDuration_sec - 5;
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
    ProcData.sleep.parameters.cortical_LH.specDeltaBandPower = LH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specThetaBandPower = LH_tempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specAlphaBandPower = LH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specBetaBandPower = LH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_LH.specGammaBandPower = LH_tempGammaSpecStruct;   
        %% BLOCK PURPOSE: Create folder for the Neural data of RH electrode

    % cortical delta
    RH_Delta = ProcData.data.cortical_RH.deltaBandPower;
    RH_baselineDelta = RestingBaselines.(baselineType).cortical_RH.deltaBandPower.(strDay).mean;
    RH_DeltaNeuro = (RH_Delta-RH_baselineDelta)/RH_baselineDelta;  

    % cortical theta
    RH_Theta = ProcData.data.cortical_RH.thetaBandPower;
    RH_baselineTheta = RestingBaselines.(baselineType).cortical_RH.thetaBandPower.(strDay).mean;
    RH_ThetaNeuro = (RH_Theta-RH_baselineTheta)/RH_baselineTheta;
  
    % cortical alpha
    RH_Alpha = ProcData.data.cortical_RH.alphaBandPower;
    RH_baselineAlpha = RestingBaselines.(baselineType).cortical_RH.alphaBandPower.(strDay).mean;
    RH_AlphaNeuro = (RH_Alpha-RH_baselineAlpha)/RH_baselineAlpha;   
 
    % cortical beta
    RH_Beta = ProcData.data.cortical_RH.betaBandPower;
    RH_baselineBeta = RestingBaselines.(baselineType).cortical_RH.betaBandPower.(strDay).mean;
    RH_BetaNeuro = (RH_Beta-RH_baselineBeta)/RH_baselineBeta;

    % cortical gamma
    RH_Gamma = ProcData.data.cortical_RH.gammaBandPower;
    RH_baselineGamma = RestingBaselines.(baselineType).cortical_RH.gammaBandPower.(strDay).mean;
    RH_GammaNeuro = (RH_Gamma-RH_baselineGamma)/RH_baselineGamma;
    % cortical Cortical
    RH_Cortical = ProcData.data.cortical_RH.corticalPower;
    RH_baselineCortical = RestingBaselines.(baselineType).cortical_RH.corticalPower.(strDay).mean;
    RH_CorticalNeuro = (RH_Cortical-RH_baselineCortical)/RH_baselineCortical;
    % Divide the neural signals into five second bins and put them in a cell array

    RH_tempDeltaStruct = cell(NBins,1);
    RH_tempThetaStruct = cell(NBins,1);
    RH_tempAlphaStruct = cell(NBins,1);
    RH_tempBetaStruct = cell(NBins,1);
    RH_tempGammaStruct = cell(NBins,1);
    RH_tempCorticalStruct = cell(NBins,1);
    % loop through all samples across the full duration in 5 second bins (NBins total)
    for b = 1:NBins
        if b == 1

            % cortical
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro(b:150)};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro(b:150)};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro(b:150)};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro(b:150)};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro(b:150)};
            RH_tempCorticalStruct(b,1) = {RH_CorticalNeuro(b:150)};
        elseif b == NBins

            % cortical
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro((((150*(b - 1)) + 1)):end)};
            RH_tempCorticalStruct(b,1) = {RH_CorticalNeuro((((150*(b - 1)) + 1)):end)};
        else

            % cortical
            RH_tempDeltaStruct(b,1) = {RH_DeltaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempThetaStruct(b,1) = {RH_ThetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempAlphaStruct(b,1) = {RH_AlphaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempBetaStruct(b,1) = {RH_BetaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempGammaStruct(b,1) = {RH_GammaNeuro((((150*(b - 1)) + 1)):(150*b))};
            RH_tempCorticalStruct(b,1) = {RH_CorticalNeuro((((150*(b - 1)) + 1)):(150*b))};
        end
    end

    % save cortical data under ProcData file
    ProcData.sleep.parameters.cortical_RH.deltaBandPower = RH_tempDeltaStruct;
    ProcData.sleep.parameters.cortical_RH.thetaBandPower = RH_tempThetaStruct;
    ProcData.sleep.parameters.cortical_RH.alphaBandPower = RH_tempAlphaStruct;
    ProcData.sleep.parameters.cortical_RH.betaBandPower = RH_tempBetaStruct;
    ProcData.sleep.parameters.cortical_RH.gammaBandPower = RH_tempGammaStruct;
    ProcData.sleep.parameters.cortical_RH.corticalPower = RH_tempCorticalStruct;
    %% BLOCK PURPOSE: Create folder for the Neural spectrogram data of RH electrode
    trialDuration_sec = ProcData.notes.trialDuration_sec;   % sec
    offset = 2.5;   % sec
    binWidth = 5;   % sec
    T = round(SpecData.cortical_RH.T,1);
    F = SpecData.cortical_RH.F;
    specRH = SpecData.cortical_RH.normS;
    freqFloor = floor(F);
    % delta
    deltaLow = freqFloor == 1;
    deltaHigh = freqFloor == 4;
    deltaLowStart = find(deltaLow,1,'first');
    deltaLowEnd = find(deltaHigh,1,'last');
    deltaSpecRH = specRH(deltaLowStart:deltaLowEnd,:);
    meanDeltaSpecRH = mean(deltaSpecRH,1);
    % theta
    thetaLow = freqFloor == 4;
    thetaHigh = freqFloor == 9;
    thetaLowStart = find(thetaLow,1,'first');
    thetaLowEnd = find(thetaHigh,1,'last');
    thetaSpecRH = specRH(thetaLowStart:thetaLowEnd,:);
    meanThetaSpecRH = mean(thetaSpecRH,1);
    % alpha
    alphaLow = freqFloor == 10;
    alphaHigh = freqFloor == 13;
    alphaLowStart = find(alphaLow,1,'first');
    alphaLowEnd = find(alphaHigh,1,'last');
    alphaSpecRH = specRH(alphaLowStart:alphaLowEnd,:);
    meanAlphaSpecRH = mean(alphaSpecRH,1);
    % beta
    betaLow = freqFloor == 13;
    betaHigh = freqFloor == 30;
    betaLowStart = find(betaLow,1,'first');
    betaLowEnd = find(betaHigh,1,'last');
    betaSpecRH = specRH(betaLowStart:betaLowEnd,:);
    meanBetaSpecRH = mean(betaSpecRH,1);
    % gamma
    gammaLow = freqFloor == 30;
    gammaHigh = freqFloor == 99;
    gammaLowStart = find(gammaLow,1,'first');
    gammaLowEnd = find(gammaHigh,1,'last');
    gammaSpecRH = specRH(gammaLowStart:gammaLowEnd,:);
    meanGammaSpecRH = mean(gammaSpecRH,1);
    % Divide the neural signals into five second bins and put them in a cell array

    RH_tempDeltaSpecStruct = cell(NBins,1);
    RH_tempThetaSpecStruct = cell(NBins,1);
    RH_tempAlphaSpecStruct = cell(NBins,1);
    RH_tempBetaSpecStruct = cell(NBins,1);
    RH_tempGammaSpecStruct = cell(NBins,1);
    % loop through all samples across the 15 minutes in 5 second bins (NBins total)
    for c = 1:NBins
        if c == 1
            startTime = offset;
            startTime_index = find(T == startTime);
            endTime = 5;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index:endTime_index)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index:endTime_index)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index:endTime_index)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index:endTime_index)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index:endTime_index)};
        elseif c == NBins
            startTime = trialDuration_sec - 5;
            [~,startTime_index] = min(abs(T - startTime));
            endTime = trialDuration_sec - offset;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index:endTime_index)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index:endTime_index)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index:endTime_index)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index:endTime_index)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index:endTime_index)};
        else
            startTime = binWidth*(c - 1);
            [~,startTime_index] = min(abs(T - startTime));
            endTime = binWidth*c;
            [~,endTime_index] = min(abs(T - endTime));

            % cortical
            RH_tempDeltaSpecStruct{c,1} = {meanDeltaSpecRH(startTime_index + 1:endTime_index + 1)};
            RH_tempThetaSpecStruct{c,1} = {meanThetaSpecRH(startTime_index + 1:endTime_index + 1)};
            RH_tempAlphaSpecStruct{c,1} = {meanAlphaSpecRH(startTime_index + 1:endTime_index + 1)};
            RH_tempBetaSpecStruct{c,1} = {meanBetaSpecRH(startTime_index + 1:endTime_index + 1)};
            RH_tempGammaSpecStruct{c,1} = {meanGammaSpecRH(startTime_index + 1:endTime_index + 1)};
        end
    end

    % save cortical data under ProcData file
    ProcData.sleep.parameters.cortical_RH.specDeltaBandPower = RH_tempDeltaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specThetaBandPower = RH_tempThetaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specAlphaBandPower = RH_tempAlphaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specBetaBandPower = RH_tempBetaSpecStruct;
    ProcData.sleep.parameters.cortical_RH.specGammaBandPower = RH_tempGammaSpecStruct;
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
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle(whiskerBins:150)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration(whiskerBins:150)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle(whiskerBins:150)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor(whiskerBins:150)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor(whiskerBins:150)};
        elseif whiskerBins == whiskerBinNumber
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):end)};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):end)};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):end)};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((150*(whiskerBins-1)) + 1)):end)};
        else
            tempWhiskerStruct(whiskerBins, 1) = {whiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempWhiskerAccelStruct(whiskerBins, 1) = {whiskerAcceleration((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempBinWhiskerStruct(whiskerBins, 1) = {binWhiskerAngle((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempForceStruct(whiskerBins, 1) = {binForceSensor((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
            tempForceRStruct(whiskerBins, 1) = {ForceSensor((((150*(whiskerBins-1)) + 1)):(150*whiskerBins))};
        end
    end
    % save whisker and force sensor data under ProcData file
    ProcData.sleep.parameters.whiskerAngle = tempWhiskerStruct;
    ProcData.sleep.parameters.whiskerAcceleration = tempWhiskerAccelStruct;
    ProcData.sleep.parameters.binWhiskerAngle = tempBinWhiskerStruct;
    ProcData.sleep.parameters.binForceSensor = tempForceStruct;
    ProcData.sleep.parameters.ForceSensor = tempForceRStruct;
    %% add pupil parameters
    ProcData.data.Pupil.mmPerPixel = 0.018; % used to convert pixels to mm
   
        % if strcmp(ProcData.data.Pupil.diameterCheck,'y') == true
        % create fields for the data
        dataTypes = {'zDiameter'};%{'mmarea','mmDiameter','zArea','zDiameter','eyeMotion','CentroidX','CentroidY','whiskerMotion'};
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
                % loop through all samples across the 52 minutes in 5 second bins (NBins total)
                for b = 1:NBins
                    if b == 1
                        data.(dataType).struct(b,1) = {data.(dataType).data(b:150)};
                    elseif b == NBins
                        if length(data.(dataType).data(:)) < (150*b)
                            data.(dataType).struct(b,1) = {data.(dataType).data((((150*(b - 1)))):end)};
                        elseif length(data.(dataType).data(:)) == (150*b)
                            data.(dataType).struct(b,1) = {data.(dataType).data((((150*(b - 1)) + 1)):end)};
                        end
                    else
                        data.(dataType).struct(b,1) = {data.(dataType).data((((150*(b - 1)) + 1)):(150*b))};
                    end
                end
                ProcData.sleep.parameters.Pupil.(dataType) = data.(dataType).struct;
            end
        % end
        %}
    %% Create folder for the EMG
    EMG = ProcData.data.EMG.emg;
    normEMG = EMG - RestingBaselines.(baselineType).EMG.emg.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1)) + 1)):(150*EMGBins))};
        end
    end
    ProcData.sleep.parameters.EMG.emg = tempEMGStruct; 
    %% raw EMG Signal
    EMG = ProcData.data.EMG.emgSignal;
    normEMG = EMG - RestingBaselines.(baselineType).EMG.emgSignal.(strDay).mean;  
    tempEMGStruct = cell(NBins,1);
    for EMGBins = 1:NBins
        if EMGBins == 1
            tempEMGStruct(EMGBins,1) = {normEMG(EMGBins:150)};
        else
            tempEMGStruct(EMGBins,1) = {normEMG((((150*(EMGBins-1)) + 1)):(150*EMGBins))};
        end
    end
    ProcData.sleep.parameters.EMG.emgSignal = tempEMGStruct;  
    %% BLOCK PURPOSE: Create folder for right hemisphere fiber data (z Scored)
    Z_NE_GFP = ProcData.data.GFP.Z_NE; 
    Z_NE_CBV = ProcData.data.CBV.Z_NE;

    Z_NE_NormCBV = (Z_NE_CBV - RestingBaselines.(baselineType).CBV.Z_NE.(strDay).mean);%/RestingBaselines.(baselineType).CBV.Z_NE.(strDay).std;
    Z_NE_NormGFP = (Z_NE_GFP - RestingBaselines.(baselineType).GFP.Z_NE.(strDay).mean);%/RestingBaselines.(baselineType).GFP.Z_NE.(strDay).std;

    Z_NE_tempGFPStruct = cell(NBins,1);  
    CBVZ_NE_tempCBVStruct = cell(NBins,1);

    for CBVBins = 1:NBins
        if CBVBins == 1
            Z_NE_tempGFPStruct(CBVBins,1) = {Z_NE_NormGFP(CBVBins:150)};
            CBVZ_NE_tempCBVStruct(CBVBins,1) = {Z_NE_NormCBV(CBVBins:150)};
        else
            Z_NE_tempGFPStruct(CBVBins,1) = {Z_NE_NormGFP((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
            CBVZ_NE_tempCBVStruct(CBVBins,1) = {Z_NE_NormCBV((((150*(CBVBins-1)) + 1)):(150*CBVBins))};

        end
    end
    % save hemodynamic data under ProcData file
    ProcData.sleep.parameters.GFP.Z_NE = Z_NE_tempGFPStruct;
    ProcData.sleep.parameters.CBV.Z_NE = CBVZ_NE_tempCBVStruct;
    %% BLOCK PURPOSE: Create folder for the right hemisphere fiber data (Percentage)
    P_NE_GFP = ProcData.data.GFP.P_NE; 
    P_NE_CBV = ProcData.data.CBV.P_NE;

    P_NE_NormCBV = (P_NE_CBV - RestingBaselines.(baselineType).CBV.P_NE.(strDay).mean);%/RestingBaselines.(baselineType).CBV.P_NE.(strDay).std;
    P_NE_NormGFP = (P_NE_GFP - RestingBaselines.(baselineType).GFP.P_NE.(strDay).mean);%/RestingBaselines.(baselineType).GFP.P_NE.(strDay).std;

    P_NE_tempGFPStruct = cell(NBins,1);  
    CBVP_NE_tempCBVStruct = cell(NBins,1);

    for CBVBins = 1:NBins
        if CBVBins == 1
            P_NE_tempGFPStruct(CBVBins,1) = {P_NE_NormGFP(CBVBins:150)};
            CBVP_NE_tempCBVStruct(CBVBins,1) = {P_NE_NormCBV(CBVBins:150)};
        else
            P_NE_tempGFPStruct(CBVBins,1) = {P_NE_NormGFP((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
            CBVP_NE_tempCBVStruct(CBVBins,1) = {P_NE_NormCBV((((150*(CBVBins-1)) + 1)):(150*CBVBins))};

        end
    end
    % save hemodynamic data under ProcData file
    ProcData.sleep.parameters.GFP.P_NE = P_NE_tempGFPStruct;
    ProcData.sleep.parameters.CBV.P_NE = CBVP_NE_tempCBVStruct;
        %% BLOCK PURPOSE: Create folder for the left and right fiber data (z Scored)
        if  isfield(ProcData.data.GFP,'Z_ACh') == true
            Z_ACh_GFP = ProcData.data.GFP.Z_ACh; 
            Z_ACh_CBV = ProcData.data.CBV.Z_ACh;
        
            Z_ACh_NormCBV = (Z_ACh_CBV - RestingBaselines.(baselineType).CBV.Z_ACh.(strDay).mean)/RestingBaselines.(baselineType).CBV.Z_ACh.(strDay).std;
            Z_ACh_NormGFP = (Z_ACh_GFP - RestingBaselines.(baselineType).GFP.Z_ACh.(strDay).mean)/RestingBaselines.(baselineType).GFP.Z_ACh.(strDay).std;
        
            Z_ACh_tempGFPStruct = cell(NBins,1);  
            CBVZ_ACh_tempCBVStruct = cell(NBins,1);
        
            for CBVBins = 1:NBins
                if CBVBins == 1
                    Z_ACh_tempGFPStruct(CBVBins,1) = {Z_ACh_NormGFP(CBVBins:150)};       
                    CBVZ_ACh_tempCBVStruct(CBVBins,1) = {Z_ACh_NormCBV(CBVBins:150)};
                else
                    Z_ACh_tempGFPStruct(CBVBins,1) = {Z_ACh_NormGFP((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
                    CBVZ_ACh_tempCBVStruct(CBVBins,1) = {Z_ACh_NormCBV((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
        
                end
            end
            % save hemodynamic data under ProcData file
            ProcData.sleep.parameters.GFP.Z_ACh = Z_ACh_tempGFPStruct;
            ProcData.sleep.parameters.CBV.Z_ACh = CBVZ_ACh_tempCBVStruct;
        end
        %% BLOCK PURPOSE: Create folder for the left hemisphere fiber data (z Scored)
        if  isfield(ProcData.data.GFP,'P_ACh') == true
            P_ACh_GFP = ProcData.data.GFP.P_ACh; 
            P_ACh_CBV = ProcData.data.CBV.P_ACh;
        
            P_ACh_NormCBV = (P_ACh_CBV - RestingBaselines.(baselineType).CBV.P_ACh.(strDay).mean)/RestingBaselines.(baselineType).CBV.P_ACh.(strDay).std;
            P_ACh_NormGFP = (P_ACh_GFP - RestingBaselines.(baselineType).GFP.P_ACh.(strDay).mean)/RestingBaselines.(baselineType).GFP.P_ACh.(strDay).std;
        
            P_ACh_tempGFPStruct = cell(NBins,1);  
            CBVP_ACh_tempCBVStruct = cell(NBins,1);
        
            for CBVBins = 1:NBins
                if CBVBins == 1
                    P_ACh_tempGFPStruct(CBVBins,1) = {P_ACh_NormGFP(CBVBins:150)};      
                    CBVP_ACh_tempCBVStruct(CBVBins,1) = {P_ACh_NormCBV(CBVBins:150)};
                else
                    P_ACh_tempGFPStruct(CBVBins,1) = {P_ACh_NormGFP((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
                    CBVP_ACh_tempCBVStruct(CBVBins,1) = {P_ACh_NormCBV((((150*(CBVBins-1)) + 1)):(150*CBVBins))};
        
                end
            end
            % save hemodynamic data under ProcData file
            ProcData.sleep.parameters.GFP.P_ACh = P_ACh_tempGFPStruct;
            ProcData.sleep.parameters.CBV.P_ACh = CBVP_ACh_tempCBVStruct;
        end
    %%
    save(procDataFileID,'ProcData');
end

end
