clear all; clc
% load('NEACh001_230612_14_00_08_ProcData.mat')
close all
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

    ACh_GRAB = [];
    NE_GRAB = [];
    PupilDiameter = [];
    ACh_mScarlet = [];
    NE_mScarlet = [];

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    load(procDataFileID);

    ACh_GRAB =  [ACh_GRAB; ProcData.data.GFP.Z_Ach];
    NE_GRAB = [NE_GRAB; ProcData.data.GFP.Z_NE];
    PupilDiameter = [PupilDiameter; ProcData.data.Pupil.zDiameter'];
    ACh_mScarlet = [ACh_mScarlet; ProcData.data.Rhodamine.Z_Ach];
    NE_mScarlet = [NE_mScarlet; ProcData.data.Rhodamine.Z_NE];

end



% ACh_GRAB = ProcData.data.GFP.Z_Ach;
% NE_GRAB = ProcData.data.GFP.Z_NE;
% PupilDiameter = ProcData.data.Pupil.zDiameter';
% ACh_mScarlet = ProcData.data.Rhodamine.Z_Ach;
% NE_mScarlet = ProcData.data.Rhodamine.Z_NE;


mean_mScarlet = (NE_mScarlet+ACh_mScarlet)./2;
NBeans = -4:0.1:4;

figure;
NE_ACh = NE_GRAB - ACh_GRAB;
subplot(3,3,1);
histogram2(mean_mScarlet,NE_ACh,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE-ACh (Z)'); xlabel('mScarlet avg (Z)');
axis equal
ACh_NE = ACh_GRAB - NE_GRAB;

subplot(3,3,2);
histogram2(NE_mScarlet,NE_GRAB,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('mScarlet RH (Z)'); axis equal
subplot(3,3,3);
histogram2(ACh_mScarlet,ACh_GRAB,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('ACh (Z)'); xlabel('mScarlet LH (Z)'); axis equal
subplot(3,3,5);
histogram2(NE_GRAB,ACh_GRAB,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('ACh (Z)'); axis equal
subplot(3,3,6);
histogram2(NE_mScarlet,ACh_mScarlet,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('mScarlet RH (Z)'); xlabel('mScarlet LH (Z)'); axis equal


subplot(3,3,7);
histogram2(NE_GRAB,PupilDiameter,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('PupilDiameter (Z)'); axis equal
subplot(3,3,8);
histogram2(ACh_GRAB,PupilDiameter,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('ACh (Z)'); xlabel('PupilDiamete (Z)'); axis equal
subplot(3,3,9);
histogram2(mean_mScarlet,PupilDiameter,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('mScarlet avg (Z)'); xlabel('PupilDiamete (Z)'); axis equal

%%
%{
NE_NREM = AnalysisResults.NEACh001.Transitions.NREMtoAWAKE.GRAB_ACh(1:900); 
ACh_NREM = AnalysisResults.NEACh001.Transitions.NREMtoAWAKE.GRAB_NE(1:900);     
RH_NREM = AnalysisResults.NEACh001.Transitions.NREMtoAWAKE.NE_Rhodamine(1:900);   
LH_NREM = AnalysisResults.NEACh001.Transitions.NREMtoAWAKE.Ach_Rhodamine(1:900);  
NBeans = -4:0.05:4;
figure;
subplot(2,2,1);
histogram2(NE_NREM,ACh_NREM,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('ACh (Z)');
subplot(2,2,2);
histogram2(NE_NREM-ACh_NREM,(RH_NREM+LH_NREM)/2,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE-ACh (Z)'); xlabel('mScarlet Avg (Z)');

NE_REM = AnalysisResults.NEACh001.Transitions.REMtoAWAKE.GRAB_ACh(1:900); 
ACh_REM = AnalysisResults.NEACh001.Transitions.REMtoAWAKE.GRAB_NE(1:900);     
RH_REM = AnalysisResults.NEACh001.Transitions.REMtoAWAKE.NE_Rhodamine(1:900);   
LH_REM = AnalysisResults.NEACh001.Transitions.REMtoAWAKE.Ach_Rhodamine(1:900);
subplot(2,2,3);
histogram2(NE_REM,ACh_REM,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE (Z)'); xlabel('ACh (Z)');
subplot(2,2,4);
histogram2(NE_REM-ACh_REM,(RH_REM+LH_REM)/2,NBeans,NBeans,'DisplayStyle','tile','ShowEmptyBins','on'); ylabel('NE-ACh (Z)'); xlabel('mScarlet Avg (Z)');
%}