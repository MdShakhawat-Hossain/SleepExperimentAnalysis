function [decData,decFileIDs,decBinTimes] = NRemoveStimSleepData_IOS(animalID,data,fileIDs,binTimes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Remove resting events from the various fields that aren't in the manual selection
%________________________________________________________________________________________________________________________

bb = 1;
for aa = 1:size(fileIDs,1)
    fileID = fileIDs{aa,1};
    procDataFileID = [animalID '_' fileID '_ProcData.mat'];
    load(procDataFileID,'-mat')
    % check that the event/file doesn't have stimulation
    try
            decData{bb,1} = data{aa,1};
            decFileIDs{bb,1} = fileIDs{aa,1};
            decBinTimes{bb,1} = binTimes{aa,1};
            bb = bb + 1;
    catch
            decData{bb,1} = data{aa,1};
            decFileIDs{bb,1} = fileIDs{aa,1};
            decBinTimes{bb,1} = binTimes{aa,1};
            bb = bb + 1;
    end
end

end

