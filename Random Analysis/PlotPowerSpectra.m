close all
clc
params.tapers = [2,3];
params.Fs = 1017;
params.fpass = [4,16];
movingwin = [3.33,1];
%%
% [S,f]=mtspectrumc(detrend(FiberData.RawData.Raw_RH(:,2),'constant'),params);
% 
% semilogx(f,S);xlabel('Frequency');ylabel('Power')
% 
% hold on
% xlim([5 15]);
% % ylim([0 1])
% yscale('log')
% % legend('Slide','Auto Fluorescence','GRAB NE','Mutant NE')
% 


%% prepare the data for spectogram
data_LC = detrend(FiberData.RawData.Raw_RH(10000:end,2));
figTime = (1:length(data_LC))/FiberData.params.DataFs;

[z2,p2,k2] = butter(4,[5 15]/(1017/2),'bandpass');
[sos2,g2] = zp2sos(z2,p2,k2);
filtData = filtfilt(sos2,g2,data_LC);

[z3,p3,k3] = butter(4,10/(1017/2),'low');
[sos3,g3] = zp2sos(z3,p3,k3);
filtData_LC = filtfilt(sos3,g3,FiberData.pressureSensor_filt(10000:end));

[Spec_LC,Tn_LC,Fn_LC] = mtspecgramc(filtData,movingwin,params);

[~,ridx] = max(Spec_LC,[],2);
HR_LC = Fn_LC(ridx);   % heart rate, in Hz
%%
data_BF = detrend(FiberData.RawData.Raw_LH(10000:end,2));
figTime = (1:length(data_BF))/FiberData.params.DataFs;

[z2,p2,k2] = butter(4,[5 15]/(1017/2),'bandpass');
[sos2,g2] = zp2sos(z2,p2,k2);
filtData = filtfilt(sos2,g2,data_BF);

[z3,p3,k3] = butter(4,10/(1017/2),'low');
[sos3,g3] = zp2sos(z3,p3,k3);
filtData_BF = filtfilt(sos3,g3,FiberData.pressureSensor_filt(10000:end));

[Spec_BF,Tn_BF,Fn_BF] = mtspecgramc(filtData,movingwin,params);

[~,ridx] = max(Spec_BF,[],2);
HR_BF = Fn_BF(ridx);   % heart rate, in Hz
%%
HR = (HR_LC+HR_BF)./2;
% patch the missing data at the beginning and end of the signal
patchedHR = horzcat(HR(1),HR,HR(end),HR(end));

% Smooth the signal with a 10 Hz low pass third-order butterworth filter
[z4,p4,k4] = butter(4,15/(1017/2),'low');
[sos4,g4] = zp2sos(z4,p4,k4);
heartRate = filtfilt(sos3,g3,patchedHR);
% heartRate = patchedHR;
%%
figure;
ax1 = subplot(412);
plot(figTime,data_LC,'color',colors('vegas gold'));
xlim([0 figTime(end)])
legend(['GRAB NE'])
xlabel('Time(s)')
ylabel('\DeltaF')

ax2 = subplot(411);
plot(figTime,filtData_LC,'b')
xlim([0 figTime(end)])
legend(['Movement'])
xlabel('Time(s)')
ylabel('force sensor')

ax3 = subplot(413);
Semilog_ImageSC(Tn_LC,Fn_LC,Spec_LC','y');
clim([0.02 0.1])
xlim([0 figTime(end)])
xlabel('Time(s)')
ylabel('frequency')
yticks([5 6 7 8 9])

ax4 = subplot(414);
figTime_2 = 1:length( heartRate);
plot(figTime_2, heartRate)
xlim([0 figTime_2(end)])
legend(['Heart Rate'])
xlabel('Time(s)')
ylabel('Frequency (Hz)')

linkaxes([ax1,ax2,ax3,ax4],'x')


