resample_rate = 10;
[z2,p2,k2] = butter(4,1/(resample_rate/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);

NE_GFP = resample(FiberData.NE.expCorrected.Percentage.F465,resample_rate,FiberData.params.DataFs);
NE_GFP1 = medfilt1(NE_GFP,20);
filtNE_GFP = filtfilt(sos2,g2,NE_GFP1);
filtNE_GFP = sgolayfilt(filtNE_GFP,3,15);


figTime = (1:length(filtNE_GFP))/(resample_rate*60);
figure; 
ax5 = subplot(9,1,[3:5]);
plot(figTime,filtNE_GFP);
legend('Fluorescent Slide');
xlabel('Time(s)');
ylabel('F/F(%)')
xlim([4.5 (figTime(end)-4.5)])
ylim([-5 5])