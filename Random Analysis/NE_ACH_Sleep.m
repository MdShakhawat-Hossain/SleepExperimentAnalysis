REMSize = size(SleepData.Manual.REM.data.GFP.Z_Ach,1);
for n = 1:REMSize

    ACh_REM(n,:) = SleepData.Manual.REM.data.GFP.Z_Ach{n,1}(end-900+1:end);
    NE_REM(n,:) = SleepData.Manual.REM.data.GFP.Z_NE{n,1}(end-900+1:end);
    LH_REM(n,:) = SleepData.Manual.REM.data.Rhodamine.Z_Ach{n,1}(end-900+1:end);
    RH_REM(n,:) = SleepData.Manual.REM.data.Rhodamine.Z_NE{n,1}(end-900+1:end);

end

ACh_REM_mean = mean(ACh_REM,1);
NE_REM_mean = mean(NE_REM,1);
LH_REM_mean = mean(LH_REM,1);
RH_REM_mean = mean(RH_REM,1);

NBeans = -2:0.05:2;
figure;
NE_ACh = NE_REM_mean - ACh_REM_mean;
subplot(3,3,1);
histogram2((LH_REM_mean + RH_REM_mean)/2,NE_ACh,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE-ACh (Z)'); xlabel('mScarlet avg (Z)');axis equal
ACh_NE = ACh_REM_mean - NE_REM_mean;
subplot(3,3,2);
histogram2(RH_REM_mean,NE_REM_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('mScarlet RH (Z)'); axis equal
subplot(3,3,3);
histogram2(LH_REM_mean,ACh_REM_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('ACh (Z)'); xlabel('mScarlet LH (Z)'); axis equal
subplot(3,3,5);
histogram2(NE_REM_mean,ACh_REM_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('ACh (Z)'); axis equal
subplot(3,3,6);
histogram2(RH_REM_mean,LH_REM_mean,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('mScarlet RH (Z)'); xlabel('mScarlet LH (Z)'); axis equal