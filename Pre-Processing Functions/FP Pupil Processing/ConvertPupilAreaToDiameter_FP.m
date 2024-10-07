function ConvertPupilAreaToDiameter_FP(rawDataFileIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Convert pupil area to diameter and determine mm pixel conversion
%________________________________________________________________________________________________________________________
    for qq = 1:size(rawDataFileIDs,1)
        rawDataFileID = rawDataFileIDs(qq,:);
        disp(['Converting pupil data of file ' num2str(qq) '/' num2str(size(rawDataFileIDs,1))]); disp(' ')
        load(rawDataFileID)
        if isfield(RawData.data.Pupil,'mmDiameter') == false % check if diameter is already available
            pupilArea = RawData.data.Pupil.pupilArea;
            diameter = sqrt(pupilArea./pi)*2;
            RawData.data.Pupil.diameter = diameter;
            RawData.data.Pupil.mmPerPixel = 0.018;
            RawData.data.Pupil.mmDiameter = diameter.*RawData.data.Pupil.mmPerPixel;
            RawData.data.Pupil.mmArea = pupilArea.*(RawData.data.Pupil.mmPerPixel^2);
            save(rawDataFileID,'RawData','-v7.3')
        end
    end
end