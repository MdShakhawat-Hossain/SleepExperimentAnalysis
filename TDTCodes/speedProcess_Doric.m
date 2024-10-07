function speedaccBinary = speedProcess_Doric(speed,Fs)
    % INPUTS:
    %       Fr: frame rate
    %       rawSpeed: un-processed walking speed data
    %       rawImage: Dalsa Image
    %
    % OUTPUTS:
    %       speedata: downsampled data used to plot on the front pannel
    %       speedaccBinary: Binarized accleration
    %       speedaveBinary: averaged binarized acceleration during one frame
    

        % 10 Hz low pass filter
        [zeroa, poleb, gain] = butter(2,20/(0.5*Fs), 'low');
        [sos,g] = zp2sos(zeroa,poleb, gain);
        speedfilt = filtfilt(sos,g,speed);

        % --- Continue working on this filtered speed
        % compute accleration, in m/s/point
        speedacc = [0; diff(speedfilt)]; % in Bing's program, this step reduced speed data by 1 sample point
        
        % Binarize the acceleration
        accThreshold = 0.005;%2*1e-6*30000/Fs; % 1e-6 m/s/point, i,e, 0.03 m/s^2
        speedaccBinary = (abs(speedacc)>accThreshold);
        
end