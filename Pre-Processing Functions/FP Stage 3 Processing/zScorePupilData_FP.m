function zScorePupilData_FP(procDataFileID,RestingBaselines)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Convert pupil diameter to z-units
%________________________________________________________________________________________________________________________

load(procDataFileID);
[~,fileDate,~] = GetFileInfo_JNeurosci2022(procDataFileID);
strDay = ConvertDate_JNeurosci2022(fileDate);
% try
    % areaMean = RestingBaselines.manualSelection.Pupil.mmarea.(strDay).mean;
    % areaStd = RestingBaselines.manualSelection.Pupil.mmarea.(strDay).std;
    diameterMean = RestingBaselines.manualSelection.Pupil.Diameter.(strDay).mean;
    diameterStd = RestingBaselines.manualSelection.Pupil.Diameter.(strDay).std;
    % pupilArea = ProcData.data.Pupil.mmarea;
    diameter = ProcData.data.Pupil.Diameter;
    % ProcData.data.Pupil.zArea = (pupilArea - areaMean)./areaStd;
    ProcData.data.Pupil.zDiameter = (diameter - diameterMean)./diameterStd;
% catch
%     ProcData.data.Pupil.diameterCheck = 'n';
% end
save(procDataFileID,'ProcData','-v7.3')

end