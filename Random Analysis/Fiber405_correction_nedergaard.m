
signal_fs = FiberData.params.DataFs; 
signal2 = FiberData.NE.rawData.F465;
signal3 = FiberData.NE.rawData.F465;
control3 = FiberData.NE.rawData.F405;

%% 3) Normalize and plot FP traces

MeanFilterOrder = 1000; % for smoothing
MeanFilter = ones(MeanFilterOrder,1)/MeanFilterOrder;

fs_signal = 1:length(signal2);
sec_signal = fs_signal/signal_fs;

med_2 = median(signal2);

% deltaF/F
delta2 = ((signal2 - med_2)/med_2)*100;
figure
plot(sec_signal,delta2);
%normalization of LC-GCaMP to 405nm channel
reg = polyfit(control3, signal3, 1); 
a = reg(1);
b = reg(2);
controlFit = a.*control3 + b;
controlFit =  filtfilt(MeanFilter,1,double(controlFit));
normDat = (signal3 - controlFit)./controlFit;
delta3 = normDat * 100;

% check fit
figure
a = subplot(4,1,1);
    plot(sec_signal(1000:end), control3(1000:end));
    title('raw control');
b = subplot(4,1,2);
    plot(sec_signal(1000:end), signal3(1000:end));
    title('raw signal');
c = subplot(4,1,3);
    plot(sec_signal(1000:end), signal3(1000:end));
    hold on
    plot(sec_signal(1000:end), controlFit(1000:end));
    title('fitted control');
d = subplot(4,1,4);
    plot(sec_signal(1000:end), delta3(1000:end));
    title('normalized signal');
linkaxes([a,b,c,d],'x');

delta2_filt = filtfilt(MeanFilter,1,double(delta2));
delta3_filt = filtfilt(MeanFilter,1,double(delta3));

% downsampling traces for plotting
ds_delta2_filt = downsample(delta2_filt, 100);
ds_delta3_filt = downsample(delta3_filt, 100);

fs_signal = 1:1:length(delta2_filt);
sec_signal = fs_signal/signal_fs;
ds_sec_signal = downsample(sec_signal, 100); % to remove noise?

% the index 1000:end removes the first second of the recoding for nicer plotting
figure
a = subplot(2,1,1);
    plot(ds_sec_signal(1000:end), ds_delta2_filt(1000:end))
    title('Syn-NE2m');
b = subplot(2,1,2);
    plot(ds_sec_signal(1000:end), ds_delta3_filt(1000:end))
    title('LC-GCaMP');
linkaxes([a,b],'x');