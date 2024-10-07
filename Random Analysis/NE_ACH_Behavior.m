

ACh_Whisk_mean = mean(EventData.GFP.Z_Ach.whisk.data,1);
NE_Whisk_mean = mean(EventData.GFP.Z_NE.whisk.data,1);
LH_Whisk_mean = mean(EventData.Rhodamine.Z_Ach.whisk.data,1);
RH_Whisk_mean = mean(EventData.Rhodamine.Z_NE.whisk.data,1);

NBeans = -1:0.005:1;
figure;
NE_ACh = NE_Whisk_mean - ACh_Whisk_mean;
subplot(3,3,1);
histogram2((LH_Whisk_mean + RH_Whisk_mean)/2,NE_ACh,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE-ACh (Z)'); xlabel('mScarlet avg (Z)');axis equal
ACh_NE = ACh_Whisk_mean - NE_Whisk_mean;
subplot(3,3,2);
histogram2(RH_Whisk_mean,NE_Whisk_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('mScarlet RH (Z)'); axis equal
subplot(3,3,3);
histogram2(LH_Whisk_mean,ACh_Whisk_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('ACh (Z)'); xlabel('mScarlet LH (Z)'); axis equal
subplot(3,3,5);
histogram2(NE_Whisk_mean,ACh_Whisk_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('ACh (Z)'); axis equal
subplot(3,3,6);
histogram2(RH_Whisk_mean,LH_Whisk_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('mScarlet RH (Z)'); xlabel('mScarlet LH (Z)'); axis equal
