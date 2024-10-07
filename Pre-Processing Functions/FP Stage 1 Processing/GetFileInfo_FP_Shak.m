function [animalID,fileDate,fileID] = GetFileInfo_FP_Shak(fileName)
%________________________________________________________________________________________________________________________
% Written by Md Shakhawat Hossain
% The Pennsylvania State University, Dept. of Biomedical Engineering
%
% Adapted from code written by Dr. Aaron T. Winder: https://github.com/awinde
%
% Purpose: Identify important aspects of a file name and output each individually.
%________________________________________________________________________________________________________________________

% identify the underscores
fileBreaks = strfind(fileName(1,:),'_');
% identify the extension
if length (strfind(fileName(1,:), '.mat')) == 1
    sExt = 2;
elseif length (strfind(fileName(1,:), '.bin')) == 1
    sExt = 1;
end
%
switch sExt
    case 1
        animalID = [];
        fileDate = fileName(:,1:fileBreaks(1) - 1);
        fileID = fileName(:,1:fileBreaks(4) - 1);
    case 2
        % use the known format to parse
        animalID = fileName(:,1:fileBreaks(1) - 1);
        if numel(fileBreaks) > 3
            fileDate = fileName(:,fileBreaks(1) + 1:fileBreaks(2) - 1);
            fileID = fileName(:,fileBreaks(1) + 1:fileBreaks(5) - 1);
        else
            fileDate = [];
            fileID = [];
        end
end

end
