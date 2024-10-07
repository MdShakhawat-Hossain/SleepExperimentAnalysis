    clc; close all; clear all
    %% character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    %%
    for i = 1:1:length(procDataFiles)
        procDataFileID = char(procDataFiles(i,:));
        load(procDataFileID);    
        %% filter the pupil movement
        if ~isfield(ProcData.data.Pupil,'filtpupilMovement')
            pupilMovement = ProcData.data.Pupil.Movement;
            pupilMovement(isnan(pupilMovement) == 1) = 0;
            filtThreshold = 10;
            filtOrder = 2;
            [z,p,k] = butter(filtOrder,filtThreshold/(ProcData.notes.dsFs/2),'low');
            [sos,g] = zp2sos(z,p,k);
            ProcData.data.Pupil.filtpupilMovement = filtfilt(sos,g,pupilMovement);
        end
        %% Binarize the Pupil Movement.
        if isfield(ProcData.data.Pupil,'PupilMovementThreshold')
            [ProcData.data.Pupil.PupilMovementThreshold] = CreateForceSensorThreshold_FP_Up(ProcData.data.Pupil.filtpupilMovement);
            ProcData.data.Pupil.binpupilMovement = BinarizeForceSensor_FP(ProcData.data.Pupil.filtpupilMovement,ProcData.data.Pupil.PupilMovementThreshold);
        
            save(procDataFileID,'ProcData','-v7.3')
        end      
        %% plot the figure
        % figure;
        % subplot(411);
        % plot(ProcData.data.Pupil.zDiameter);
        % 
        % subplot(412);
        % plot(ProcData.data.Pupil.filtpupilMovement);
        % 
        % subplot(413);
        % plot(ProcData.data.Pupil.binpupilMovement);
    end


