close all
load('IRFDataCell.mat')
rowNum = 4;
emgTime = dataCell{rowNum,2}(:,1);
emg = dataCell{rowNum,2}(:,2);
moveY = dataCell{rowNum,3}(:,2);
secondsPerFrame = dataCell{rowNum,4};
moveYTime = secondsPerFrame*(1:length(moveY));
hrfY = dataCell{rowNum,9};
gammaX = dataCell{rowNum,8};
gammaY = dataCell{rowNum,10};
convY = conv(emg,gammaY);
convYTrimmed = convY(1:end-(length(gammaY)-1));
figure(1)
subplot(4,1,1)
plot(emgTime,convYTrimmed)
xlim([0 610])
subplot(4,1,2)
plot(moveYTime,moveY)
xlim([0 610])
subplot(4,1,3)
plot(emgTime,emg)
xlim([0 610])
subplot(4,1,4)
plot(hrf.HRF_time,gammaY)
xlim([-10 10])

emgDataInterp=zeros(size(movementData.targetPosition));
emgDataInterp(:,1)=moveYTime;
emgDataInterp(:,2)=interp1(emgTime,emg,moveYTime,'linear');
best_lag=round(1/secondsPerFrame);%about 1 second

figure(2)
emg_bins=.5:.1:3.5;
pixel_bins=-3:.5:7;
subplot(221)
scatter(emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1))
xlabel('emg')
ylabel('shift, pixels')
subplot(222)
scatter(emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2))
hold on
plot(emg,convYTrimmed-4,'ko')
hold off
xlabel('emg')
ylabel('shift, pixels')
subplot(223)
hh1=histogram2(emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,1), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh1.Values'+1))
axis xy
xlabel('emg')
ylabel('shift, pixels')
subplot(224)
hh2=histogram2(emgDataInterp(1:(end-best_lag),2),movementData.targetPosition(best_lag:end-1,2), emg_bins, pixel_bins,...
    'DisplayStyle','tile','ShowEmptyBins','on');
imagesc(emg_bins, pixel_bins,log(hh2.Values'+1))
axis xy
colormap hot
xlabel('emg')
ylabel('shift, pixels')