function [Speed, time] = loadTDMSfile_Doric(filen_tdms,Fs)
%% load *.tdms file (if any)
if exist(filen_tdms,'file')==2
    TDMSdata1 = convertTDMS(0,filen_tdms);
    % two channels in TDMS file:
    %               channel 1: ball velocity signal
    %               channel 2: trigger signal
    TDMSdata = [TDMSdata1.Data.MeasuredData.Data];
    
    if size(TDMSdata,2) == 2
        Speed = TDMSdata(:,2)*2*pi*0.06;
        FrameOnset = TDMSdata(:,1); % camera trigger
    elseif size(TDMSdata,2) == 3
        
        Speed = TDMSdata(:,3)*2*pi*0.06;
        FrameOnset = TDMSdata(:,2); % camera trigger
        %     LH_560 = TDMSdata(:,1);
    end
    
    time = ((1:length(Speed))-1)/Fs;

    
%     
%     figure;
%     h(1) = subplot(311);
%     plot(time, FrameOnset);
%      h(2) = subplot(312);
%     plot(time, Speed);
%      h(3) = subplot(313);
%     plot(time, LH_560);
%     linkaxes(h,'x')
    
    
    %% we need to trim the speed data to make sure ONLY data during image acquistion will be shown
    % find start time of the data acquistion, i.e., first frame onset
    locs.collect_onset = find(FrameOnset>=4.5, 1, 'first'); % align everything using frame onset
    time = time - locs.collect_onset/Fs; % set the collection onset as time zero
    
    figure;
    h(1) = subplot(211);
    plot(time, Speed);
    h(2) = subplot(212);
    plot(time, FrameOnset);
    linkaxes(h,'x')
    
    
    
    % find positive and negative peaks to determine the on and off time for
    % each frame
    % THIS IS THE EXACT TIME FOR EACH FRAME
    FrameOnset_diff = [0; diff(FrameOnset);]; % take first derivative
    [~, idx.frame_on] = findpeaks(FrameOnset_diff, 'MinPeakDistance',200, 'MinPeakHeight',3);
    idx.frame_on = idx.frame_on-1;
    pks.frame_on = FrameOnset(idx.frame_on);
    locs.frame_on = time(idx.frame_on);
    [~, idx.frame_off] = findpeaks(-FrameOnset_diff, 'MinPeakDistance',200, 'MinPeakHeight',3);
    pks.frame_off = FrameOnset(idx.frame_off); % frame_off peaks
    locs.frame_off = time(idx.frame_off); % frame_off time stamp

    % here is the actual data collection segment for the 2PLSM image
    seg = idx.frame_on(1):idx.frame_on(end);

    
    % tructuate data to match the image
    Speed = Speed(seg);
%     LH_560 = LH_560(seg);
    time = time(seg);
    
else
    Speed = [];
%     LH_560 = [];
    time = [];
end
end