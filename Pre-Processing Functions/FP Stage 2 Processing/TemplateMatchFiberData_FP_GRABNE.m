function [] = TemplateMatchFiberData_FP_GRABNE(FiberDataFileIDs,rawDataFileIDs,trimTime)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Doric triggers the LabVIEW acquisition program to start recording, but there is a slight (~ 1 second) lag
%            associated with the Doric data. The force sensor is duplicated, and this function serves to correct that
%            offset by finding the peak in the cross correlation, and shifting the LabVIEW signals based on the number of
%            lags. The beginning/end of all signals are then snipped appropriately after shifting.
%________________________________________________________________________________________________________________________

for aa = 1:size(rawDataFileIDs,1)
    %% find offset between the two force sensor signals using the cross correlation
    rawDataFileID = rawDataFileIDs(aa,:);
    fiberDataFileID = FiberDataFileIDs(aa,:);
    load(fiberDataFileID);
    load(rawDataFileID);
    disp(['Correcting offset in file number ' num2str(aa) ' of ' num2str(size(fiberDataFileID, 1)) '...']); disp(' ');
    [animalID,~,fileID] = GetFileInfo_FP(rawDataFileID);
    analogSamplingRate = RawData.notes.analogSamplingRate;
    whiskCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    pupilCamSamplingRate = RawData.notes.pupilCamSamplingRate;
    tdtSamplingRate = FiberData.params.DataFs;
    whiskerCamSamplingRate = RawData.notes.whiskCamSamplingRate;
    trialDuration = RawData.notes.trialDuration_long;
    tdtTrialDuration = trialDuration - 1;
    dsFs = 100; % Hz
    %% check if there was any weird dynamics in the data that was removed. Adjust the stimulation system
    if isfield(FiberData.params,'WData')==1
        startSample = (analogSamplingRate*(FiberData.params.WData.StartDataRemove + trimTime)) + 1;
        endSample = (analogSamplingRate*(FiberData.params.WData.StartDataRemove + trimTime));
        RawData.data.stimulations_long(startSample:endSample) = 0;
    end
    %% 

    labviewSyncSignal = detrend(resample(RawData.data.forceSensor_long,dsFs,analogSamplingRate),'constant');
    labviewSyncSignal = labviewSyncSignal(dsFs*trialDuration/4:dsFs*trialDuration/2);
    fiberSyncSignal = detrend(resample(FiberData.pressureSensor,dsFs,tdtSamplingRate),'constant');
    fiberSyncSignal = fiberSyncSignal(round(dsFs*tdtTrialDuration/4):round(dsFs*tdtTrialDuration/2));
    sampleDiff = length(labviewSyncSignal) - length(fiberSyncSignal);
    fiberSyncSignal = vertcat(fiberSyncSignal,zeros(sampleDiff,1));
    % dsFs xcorr
    maxLag = 10*dsFs;
    [r,lags] = xcorr(labviewSyncSignal,fiberSyncSignal,maxLag);
    [~,index] = max(r(1:(length(r) - 1)/2) - dsFs/10);
    offset = lags(index);
    offsetTime = abs(offset/dsFs);
    disp(['LabVIEW trailed TDT by ' num2str(offsetTime) ' seconds.']); disp(' ')

    % offset
    dsOffset = round(dsFs*(abs(offset)/dsFs));
    dsFs_pad = zeros(1,abs(dsOffset));
    labviewSyncShift = horzcat(dsFs_pad,labviewSyncSignal);

    corrOffset = figure;
    % set(corrOffset,'WindowStyle','docked');
    corrOffset.WindowState = 'minimized';
    sgtitle([rawDataFileID ' cross ' fiberDataFileID])

    ax1 = subplot(3,1,1);
    plot((1:length(fiberSyncSignal))/dsFs,fiberSyncSignal,'b')
    hold on;
    plot((1:length(labviewSyncSignal))/dsFs,labviewSyncSignal,'r')
    title({[animalID ' ' fileID ' sync channel data'],'offset correction between TDT and LabVIEW DAQ'})
    legend('Original TDT','Original LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    ax2 = subplot(3,1,2); %#ok<NASGU>
    plot(lags/dsFs,r,'k')
    title('Cross Correlation between the two signals')
    ylabel('Correlation (A.U.)')
    xlabel('Lag (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    ax3 = subplot(3,1,3);
    plot((1:length(fiberSyncSignal))/dsFs,fiberSyncSignal,'k')
    hold on;
    plot((1:length(labviewSyncShift))/dsFs,labviewSyncShift,'b')
    title({'Shifted correction between TDT and LabVIEW DAQ',['Offset value: ' num2str(offset) ' samples or ~' num2str(offset/dsFs) ' seconds']})
    legend('Original TDT','Shifted LabVIEW')
    ylabel('A.U.')
    xlabel('Time (sec)')
    set(gca,'Ticklength',[0,0])
    axis tight
    linkaxes([ax1,ax3],'x')

    %% apply TDT correction to the data and trim excess time
    tdtSampleDiff = tdtSamplingRate*(trialDuration) - length(FiberData.pressureSensor);
    tdtCut = floor(trimTime*tdtSamplingRate - tdtSampleDiff);
    % remove some previous data
    if isfield(RawData.data,'NE') == 1
    RawData.data = rmfield(RawData.data,'NE');
    end
    if isfield(RawData.data,'Ach') == 1
    RawData.data = rmfield(RawData.data,'Ach');
    end
    if isfield(RawData.data,'ACh') == 1
    RawData.data = rmfield(RawData.data,'ACh');
    end

    
    tdtFields = {'NE','ACh'};
    
    tdtsubfields = {'dFF0_p','dFF0_z'};
    for cc = 1:length(tdtFields)
%         subfields = fieldnames(FiberData.(tdtFields{1,cc}))';
%         subfields = subfields(4:end);
        for dd = 1:length(tdtsubfields)
            subsubfields = fieldnames(FiberData.(tdtFields{1,cc}).(tdtsubfields{1,dd}))';
            for nn = 1:length( subsubfields)
                RawData.data.(char(tdtFields(1,cc))).(char(tdtsubfields(1,dd))).(char(subsubfields(1,nn))) = FiberData.(char(tdtFields(1,cc))).(char(tdtsubfields(1,dd))).(char(subsubfields(1,nn)))(floor(trimTime*tdtSamplingRate):end - (tdtCut + 1));
            end
        end
    end
    RawData.data.pressureSensor = FiberData.pressureSensor(floor(trimTime*tdtSamplingRate):end - (tdtCut + 1));
    %% apply labview correction to the data and trim excess time
    labviewFields = {'cortical_LH','cortical_RH','EMG','forceSensor','stimulations'};
    for cc = 1:length(labviewFields)
        labviewAnalogPad = zeros(1,round(offsetTime*analogSamplingRate));
        labviewAnalogShift = horzcat(labviewAnalogPad,RawData.data.([labviewFields{1,cc} '_long']));
        labviewAnalogSampleDiff = analogSamplingRate*trialDuration - length(labviewAnalogShift);
        labviewAnalogCut = trimTime*analogSamplingRate - labviewAnalogSampleDiff;
        RawData.data.(labviewFields{1,cc}) = labviewAnalogShift(floor(trimTime*analogSamplingRate):end - (labviewAnalogCut + 1));
    end
    
    labviewWhiskerPad = zeros(1,round(offsetTime*whiskCamSamplingRate));
    labviewWhiskerShift = horzcat(labviewWhiskerPad,RawData.data.whiskerAngle_long);
    labviewWhiskerSampleDiff = whiskCamSamplingRate*trialDuration - length(labviewWhiskerShift);
    labviewWhiskerCut = trimTime*whiskerCamSamplingRate - labviewWhiskerSampleDiff;
    RawData.data.whiskerAngle = labviewWhiskerShift(floor(trimTime*whiskCamSamplingRate):end - (labviewWhiskerCut + 1));
    
    %% pupil Data
    if isfield(RawData.data,'Pupil')==1

    labviewPupilPad = zeros(1,round(offsetTime*pupilCamSamplingRate));
    labviewPupilDiameterShift = horzcat(labviewPupilPad,RawData.data.Pupil.diameter);

    % labviewPupilAreaShift = horzcat(labviewPupilPad,RawData.data.Pupil.pupilArea);
    % labviewPupilmmAreaShift = horzcat(labviewPupilPad,RawData.data.Pupil.mmArea);
    % labviewPupilmmDiameterShift = horzcat(labviewPupilPad,RawData.data.Pupil.mmDiameter);

    % labviewPupilMajorShift = horzcat(labviewPupilPad,RawData.data.Pupil.pupilMajor);
    % labviewPupilpatchMajorShift = horzcat(labviewPupilPad,RawData.data.Pupil.patchMajorAxis);
    % labviewPupilMinorShift = horzcat(labviewPupilPad,RawData.data.Pupil.pupilMinor);
    % labviewPupilpatchMinorShift = horzcat(labviewPupilPad,RawData.data.Pupil.patchMinorAxis);
    % labviewPupilpatchCentroidX = horzcat(labviewPupilPad,RawData.data.Pupil.patchCentroidX);
    % labviewPupilpatchCentroidY = horzcat(labviewPupilPad,RawData.data.Pupil.patchCentroidY);

    % labviewPupilCPad = zeros(2,round(offsetTime*pupilCamSamplingRate));
    % labviewPupilCentroidShift = horzcat(labviewPupilCPad,RawData.data.Pupil.pupilCentroid');

    labviewPupilSampleDiff = pupilCamSamplingRate*trialDuration - length(labviewPupilDiameterShift);
    labviewPupilCut = trimTime*pupilCamSamplingRate - labviewPupilSampleDiff;
    RawData.data.pupilDiameter = labviewPupilDiameterShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));

        if isfield(RawData.data.Pupil,'movement')==1
                labviewPupilMovementShift = horzcat(labviewPupilPad,RawData.data.Pupil.movement');
                labviewPupilMovementSampleDiff = pupilCamSamplingRate*trialDuration - length(labviewPupilMovementShift);
                labviewPupilMovementCut = trimTime*pupilCamSamplingRate - labviewPupilMovementSampleDiff;
                RawData.data.pupilMovement = labviewPupilMovementShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilMovementCut + 1));
        end
    % RawData.data.pupilArea = labviewPupilAreaShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilmmArea = labviewPupilmmAreaShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilmmDiameter = labviewPupilmmDiameterShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    
    % RawData.data.pupilMajor = labviewPupilMajorShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilpatchMajor = labviewPupilpatchMajorShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilMinor = labviewPupilMinorShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilpatchMinor = labviewPupilpatchMinorShift(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilpatchCentroidX = labviewPupilpatchCentroidX(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    % RawData.data.pupilpatchCentroidY = labviewPupilpatchCentroidY(floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));

    % RawData.data.pupilCentroid = labviewPupilCentroidShift(:,floor(trimTime*pupilCamSamplingRate):end - (labviewPupilCut + 1));
    end
    RawData.notes.offsetCorrect = true;
    RawData.notes.trimTime = trimTime;
    RawData.notes.trialDuration_sec = trialDuration - 2*trimTime;
    %% check shift
    %{
    fiberSyncSignal_2 = resample(detrend(RawData.data.pressureSensor,'constant'),dsFs,tdtSamplingRate);
    labviewSyncSignal_2 = resample(detrend(RawData.data.forceSensor,'constant'),dsFs,analogSamplingRate);
    checkShift = figure;
    p1  =  plot((1:length(fiberSyncSignal_2))/dsFs,fiberSyncSignal_2 ,'b');
    hold on;
    p2 = plot((1:length(labviewSyncSignal_2))/dsFs,labviewSyncSignal_2 ,'r');
    legend([p1,p2],'TDT','LabVIEW')
    title([rawDataFileID ' cross ' fiberDataFileID])
    axis tight
    %}
    %% save files
    RawData.notes.TDT = FiberData.params;
    save(rawDataFileID,'RawData','-v7.3')
    %% Save the file to directory.
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/XCorr Shift/'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(corrOffset,[dirpath animalID '_' fileID '_XCorrShift']);
    close(corrOffset)
    % dirpath = [pathstr '/Figures/Shift Check/'];
    % if ~exist(dirpath,'dir')
    %     mkdir(dirpath);
    % end
    % savefig(checkShift,[dirpath animalID '_' fileID '_ShiftCheck']);
    % close(checkShift)
end

end
