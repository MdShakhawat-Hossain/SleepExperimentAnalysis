clear all;
close all;
clc;
% 
load('NEACh008_230630_10_17_31_RawData.mat');

window  = 512;
fs_analog = RawData.notes.analogSamplingRate;
signal_fs = 1017;
duration = 3120;

eeg = (resample(RawData.data.cortical_LH,window,fs_analog));
emg = (resample(RawData.data.EMG,window,fs_analog));
ne = RawData.data.NE.dFF0_z.F465;
%% Process the analog signals
% eeg
    fpass = [1 100];
    [z1,p1,k1] = butter(3,fpass/(window/2));
    [sos1,g1] = zp2sos(z1,p1,k1);
    filtNeuro = filtfilt(sos1,g1,eeg - mean(eeg));
    % [z2,p2,k2] = butter(3,10/(window/2),'low');
    % [sos2,g2] = zp2sos(z2,p2,k2);
    % smoothPower = filtfilt(sos2,g2,filtNeuro.^2);
    % procNeuroPower = max(resample(smoothPower,neuroFs,analogFs),0);

% emg
    % fpass = [300,3000];
    % 
    % [z,p,k] = butter(3,fpass/(window/2));
    % [sos,g] = zp2sos(z,p,k);
    % filtEMG = filtfilt(sos,g,emg - mean(emg));
    % kernelWidth = 0.5;
    % smoothingKernel = gausswin(kernelWidth*window)/sum(gausswin(kernelWidth*window));
    % EMGPwr = log10(conv(filtEMG.^2,smoothingKernel,'same'));

    fpass = [10,100];
    [z,p,k] = butter(3,fpass/(window/2));
    [sos,g] = zp2sos(z,p,k);
    filtEMG = filtfilt(sos,g,emg - mean(emg));
%%
trial_eeg = zeros(duration,window);
trial_emg = zeros(duration,window);
trial_ne = zeros(duration,signal_fs);
for i = 1:1:duration
    trial_eeg(i,:) =  filtNeuro(1 + ((i-1)*window) : (i*window));
    trial_emg(i,:) =  filtEMG(1 + ((i-1)*window) : (i*window));
    trial_ne(i,:) =  ne(1 + ((i-1)*signal_fs) : (i*signal_fs));
end

%% training labels
load('NEACh008_230630_10_17_31_TrainingData.mat')


TableSize = height(trainingTable);
DataSize = TableSize*5;
label_t = zeros(1,DataSize);

OGLabels = trainingTable.behavState;
for SL = 1:1:TableSize
    if OGLabels(SL) == "Not Sleep"
        for ML = 1:1:5
            label_t(((SL-1)*5)+ML) = 1;
        end
    elseif OGLabels(SL) == "NREM Sleep"
        for ML = 1:1:5
            label_t(((SL-1)*5)+ML) = 2;
        end
    elseif OGLabels(SL) == "REM Sleep"
        for ML = 1:1:5
            label_t(((SL-1)*5)+ML) = 3;
        end
    end
end

clearvars -EXCEPT label_t signal_fs trial_eeg trial_emg trial_ne window
% figure; 
% subplot(311);plot(eeg);
% subplot(312);plot(emg);
% subplot(313);plot(ne);