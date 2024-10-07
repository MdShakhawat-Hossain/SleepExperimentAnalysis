% load('IRFDataCell.mat')
close all
for n = 3:size(dataCell,1)
    ballData = load([dataCell{n,1}(1:31) '.txt']);
    % ballData = load('221208_006.txt');
    emgDataOnly = [ballData(:,1) ballData(:,4)];
%     emgDataOnly = dataCell{n,2};
%     movementData.targetPosition = dataCell{n,3};
%     movementData.secondsPerFrame = dataCell{n,4};
    load(dataCell{n,1});
    % load('221208_006_processe_2layerBrainInSkullDataFinal.mat');
    if movementData.frames(1) > 2 || movementData.frames(2) < 24000
        disp('movement data frames may not be entire run, resample brain motion data instead')
        continue
    end
    % Fs = 19.775;
    sampleRate = 10000;
    Fs = 39;
    emgFilt = filterEMGDataSub(emgDataOnly,sampleRate,Fs);
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(emgFilt(:,1),emgFilt(:,2))
    title('select EMG threshold')
    [~,EMGThresh] = ginput(1);
% EMGThresh = 2.2;
    close all
    EMGBinary = emgFilt(:,2)>EMGThresh;
    % moveDataX = interp1((1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame,movementData.targetPosition(:,1),(1:length(movementData.targetPosition(:,1))*.9883)*0.0256);
    % moveDataY = interp1((1:length(movementData.targetPosition(:,2)))*movementData.secondsPerFrame,movementData.targetPosition(:,2),(1:length(movementData.targetPosition(:,2))*.9883)*0.0256);
    moveDataX = interp1((1:length(movementData.targetPosition(:,1)))*movementData.secondsPerFrame,movementData.targetPosition(:,1),(1:length(movementData.targetPosition(:,1))*.984)*0.0257);
    moveDataY = interp1((1:length(movementData.targetPosition(:,2)))*movementData.secondsPerFrame,movementData.targetPosition(:,2),(1:length(movementData.targetPosition(:,2))*.984)*0.0257);
    moveDataX(isnan(moveDataX)) = [];
    moveDataY(isnan(moveDataY)) = [];
    
    if length(moveDataX) > length(EMGBinary)
        disp(['move data longer than emg by ' num2str(length(moveDataX) - length(EMGBinary))])
        moveDataX = moveDataX(1:length(EMGBinary));
        moveDataY = moveDataY(1:length(EMGBinary));
    elseif length(moveDataX) < length(EMGBinary)
        disp(['emg data longer than move by ' num2str(length(EMGBinary) - length(moveDataX))])
        %     EMGBinary = EMGBinary((length(EMGBinary) - length(moveDataX))+1:end);
        EMGBinary = EMGBinary(1:length(moveDataX));
    end
    % [~, ~, speedaveBinary] = speedProcess_2PLSM(convertBallVoltToMPS(ballData(:,2)), 10000, Fs, size(movementData.targetPosition,1));
    % [~, ~, speedaveBinary] = speedProcess_2PLSM(2*pi*0.06*(ballData(:,2))/10, 10000, Fs, size(movementData.targetPosition,1));
    
    [zeroa, poleb, gain] = butter(5,(1/(Fs/2)),'low'); % 10 Hz low pass filter design, 5th order Butterworth
    [sos,g] = zp2sos(zeroa,poleb,gain);
    moveDataX = filtfilt(sos,g,moveDataX);
    moveDataY = filtfilt(sos,g,moveDataY);
    
    if length(moveDataX) > length(EMGBinary)
        disp(['move data longer than emg by ' num2str(length(moveDataX) - length(EMGBinary))])
        moveDataX = moveDataX(1:length(EMGBinary));
        moveDataY = moveDataY(1:length(EMGBinary));
    elseif length(moveDataX) < length(EMGBinary)
        disp(['emg data longer than move by ' num2str(length(EMGBinary) - length(moveDataX))])
        EMGBinary = EMGBinary((length(EMGBinary) - length(moveDataX))+1:end);
    end
    % make baseline of emg power 0, not 1
    HRFx = OXY_HRF(double(EMGBinary),moveDataX,Fs);
    HRFx.HRF = HRFx.HRF(HRFx.HRF_time >= -1);
    HRFx.HRF_time = HRFx.HRF_time(HRFx.HRF_time >= -1);
    [paramsx, HRFx_gamma, R2x] = HRF_motion_gammaFit(HRFx.HRF, Fs, 'unconstrained');
    
    HRFy = OXY_HRF(double(EMGBinary),moveDataY,Fs);
    HRFy.HRF = HRFy.HRF(HRFy.HRF_time >= -1);
    HRFy.HRF_time = HRFy.HRF_time(HRFy.HRF_time >= -1);
    [paramsy, HRFy_gamma, R2y] = HRF_motion_gammaFit(HRFy.HRF, Fs, 'unconstrained');
    
    
%     figure(1)
%     plot(HRFx.HRF_time, HRFx.HRF,'k'); hold on;
%     plot(HRFx.HRF_time, HRFx_gamma,'r');
%     title('Medial/Lateral Impulse Response')
%     xlabel('Time (s)')
%     ylabel('Medial/Lateral Motion (\mum)')
%     
%     figure(2)
%     plot(HRFy.HRF_time, HRFy.HRF,'k'); hold on;
%     plot(HRFy.HRF_time, HRFy_gamma,'r');
%     title('Rostral/Caudal Impulse Response')
%     xlabel('Time (s)')
%     ylabel('Rostral/Caudal Motion (\mum)')
    
    dataCell{n,6} = EMGThresh;
    dataCell{n,7} = HRFx;
    dataCell{n,8} = HRFx_gamma;
    dataCell{n,9} = HRFy;
    dataCell{n,10} = HRFy_gamma;
end
save('IRFDataCell.mat','dataCell')
% 
% 
% % Gamma function fit
% 
% out.HRF = out.HRF(out.HRF_time >= 0)
% 
% % Fs = 30.2767;
% Fs = 31;
% [~, ~, speedaveBinary] = speedProcess_2PLSM(convertBallVoltToMPS(ballData(:,2)), 10000, Fs, size(movementData.emgData,1));
% 
% [zeroa, poleb, gain] = butter(5,(10/(Fs/2)),'low'); % 10 Hz low pass filter design, 5th order Butterworth
% [sos,g] = zp2sos(zeroa,poleb,gain);
% emgFilt = filtfilt(sos,g,movementData.emgData(:,1));
% 
% out = OXY_HRF(speedaveBinary,emgFilt(:,2),Fs);
% figure(3)
% plot(out.HRF_time,out.HRF)

function procData = filterEMGDataSub(emg,sampleRate,dsFs)

close all

rawEMG = emg(:,2);
t = emg(:,1);
rawEMG=filloutliers(rawEMG, 'center','mean', 'ThresholdFactor', 4);%this takes out any large artifacts
% process EMG data
% fpass = [300,(sampleRate/2)-1];   % Hz
fpass = [300,3000];   % Hz
trialDuration_sec = t(end);   % read this variable in from your data in seconds
analogSamplingRate = sampleRate;   % Hz - change if yours is different
% dsFs = 40;   % Hz - downsampled frequency
analogExpectedLength = trialDuration_sec*analogSamplingRate;
trimmedEMG = rawEMG(1:min(analogExpectedLength,length(rawEMG)));
[z,p,k] = butter(5,fpass./(analogSamplingRate/2));
[sos,g] = zp2sos(z,p,k);
filtEMG = filtfilt(sos,g,trimmedEMG - mean(trimmedEMG));
% kernelWidth = 0.5;
kernelWidth = 0.005;
smoothingKernel = gausswin(kernelWidth*analogSamplingRate)/sum(gausswin(kernelWidth*analogSamplingRate));
EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));
% EMGPwr = conv(filtEMG.^2,smoothingKernel,'same');
% EMGPwr = filtEMG.^2;
resampEMG = resample(EMGPwr,dsFs,analogSamplingRate);
procEMG = resampEMG;   % save this as your final array
procT = 0:1/dsFs:t(end);
% procT = procT(1:end-1);
if length(procT) > length(procEMG)
    procT = procT(1:length(procEMG));
elseif length(procT) < length(procEMG)
    procEMG = procEMG(1:length(procT));
end
figure('units','normalized','outerposition',[0 0 1 1])
plot(procT,procEMG)
title('Select start and stop time (t) to use to determine mean baseline')
tVals = ginput(2);
close
i = tVals(1,1) <= procT & procT <= tVals(2,1);
baseline = mean(procEMG(i));
procEMG = procEMG - (baseline - 1);
for n = length(procEMG)-100:length(procEMG)
    if procEMG(n) < 0
        procEMG(n) = procEMG(n-1);
    end
end
procData = [procT',procEMG];
% plot(procData(:,1),procData(:,2))

% subplot(2,1,1)
% plot(emg(:,1),emg(:,2))
% subplot(2,1,2)
% plot(procData(:,1),procData(:,2))
end