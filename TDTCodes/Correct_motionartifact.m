function [Data] =  Correct_motionartifact(Data,Params)    

    % as there may be spikes, which will totally change the polynominal fit
    % first conduct a median filter to remove potential spikes
    filt_window = 5;
    Medfit.F405 = medfilt1(Data(:,1), filt_window*Params.DataFs);
    Medfit.F465 = medfilt1(Data(:,2), filt_window*Params.DataFs);
    Medfit.F560 = medfilt1(Data(:,3), filt_window*Params.DataFs);

    Residual.F405 = Data(:,1)-Medfit.F405;
    Residual.F465 = Data(:,2)-Medfit.F465;
    Residual.F560 = Data(:,3)-Medfit.F560;

    spike_thresh = 3;

    removeSpike.F405 = Residual.F405 >= spike_thresh*(max(Medfit.F405)-min(Medfit.F405));
    removeSpike.F465 = Residual.F465 >= spike_thresh*(max(Medfit.F465)-min(Medfit.F465));
    removeSpike.F560 = Residual.F560 >= spike_thresh*(max(Medfit.F560)-min(Medfit.F560));

    if any(removeSpike.F405)
        disp('Spikes in 405 channel');
        Data(removeSpike.F405,1) = Medfit.F405(removeSpike.F405);
    end
    if any(removeSpike.F465)
        disp('Spikes in 465 channel');
        Data(removeSpike.F465,2) = Medfit.F465(removeSpike.F465);
    end
    if any(removeSpike.F560)
        disp('Spikes in 560 channel');
        Data(removeSpike.F560,3) = Medfit.F560(removeSpike.F560);
    end
    