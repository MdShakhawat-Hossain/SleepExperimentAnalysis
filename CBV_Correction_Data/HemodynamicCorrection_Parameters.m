clear all; close all; clc
CorrectionDataFileStruct = dir('*CBVCorrectionParameters.mat');
CorrectionDataFiles = {CorrectionDataFileStruct.name}';
CorrectionDataFileIDs = char(CorrectionDataFiles);

for CP = 1:size(CorrectionDataFileIDs,1)
    FileID = CorrectionDataFileIDs(CP,:);
    load(FileID);
    Parameters.Percent.LH.Slope(CP) = CBVCorrection.Percent.LH.Slope;
    Parameters.Percent.LH.Inter(CP) = CBVCorrection.Percent.LH.Inter;

    Parameters.Percent.RH.Slope(CP) = CBVCorrection.Percent.RH.Slope;
    Parameters.Percent.RH.Inter(CP) = CBVCorrection.Percent.RH.Inter;

    Parameters.ZScored.LH.Slope(CP) = CBVCorrection.ZScored.LH.Slope;
    Parameters.ZScored.LH.Inter(CP) = CBVCorrection.ZScored.LH.Inter;

    Parameters.ZScored.RH.Slope(CP) = CBVCorrection.ZScored.RH.Slope;
    Parameters.ZScored.RH.Inter(CP) = CBVCorrection.ZScored.RH.Inter;
    clear CBVCorrection;
end

%% Plot the parameters
figure; 
subplot(211)
plot(Parameters.Percent.LH.Slope);hold on; plot(Parameters.Percent.RH.Slope); plot(Parameters.ZScored.LH.Slope); plot(Parameters.ZScored.RH.Slope);
title('CBV Correction Slope');xlabel('sessions');ylabel('slope');legend('LH_P','RH_P','LH_Z','RH_Z')
subplot(212)
plot(Parameters.Percent.LH.Inter);hold on; plot(Parameters.Percent.RH.Inter); plot(Parameters.ZScored.LH.Inter); plot(Parameters.ZScored.RH.Inter);
title('CBV Correction Intersect');xlabel('sessions');ylabel('Intersect');legend('LH_P','RH_P','LH_Z','RH_Z')

saveas(gcf,['HemodynamicCorrectionParameters_plot.fig'],'fig')
saveas(gcf,['HemodynamicCorrectionParameters_plot.tiff'],'tiff')
close all;

%% calculate the mean value for future use
CBVCorrection.ZScored.LH.Slope = nanmean(Parameters.ZScored.LH.Slope);
CBVCorrection.ZScored.LH.Inter = nanmean(Parameters.ZScored.LH.Inter);
CBVCorrection.ZScored.RH.Slope = nanmean(Parameters.ZScored.RH.Slope);
CBVCorrection.ZScored.RH.Inter = nanmean(Parameters.ZScored.RH.Inter);

CBVCorrection.Percent.LH.Slope = nanmean(Parameters.Percent.LH.Slope);
CBVCorrection.Percent.LH.Inter = nanmean(Parameters.Percent.LH.Inter);
CBVCorrection.Percent.RH.Slope = nanmean(Parameters.Percent.RH.Slope);
CBVCorrection.Percent.RH.Inter = nanmean(Parameters.Percent.RH.Inter);

save('CBVMle_Hemodynamic_Correction_Parameters.mat',"CBVCorrection");
'C:\Users\mfh5734\OneDrive - The Pennsylvania State University\Science_Research\SleepExperiment\CBV_Correction_Data\CBVMle_Hemodynamic_Correction_Parameters.mat'
