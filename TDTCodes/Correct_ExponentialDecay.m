function [expCorrected, predicted] = Correct_ExponentialDecay(RawData,forceSensor,Params,ChannelName)

        %% remove some data to perform a better fit
            timeN = 210; % remove first 210 seconds;
            sampleN = floor(timeN*Params.DataFs);
           
            % adding dummy data to match index
            RawData = RawData(sampleN+1:end-sampleN,:);
            forceSensor = forceSensor(sampleN+1:end-sampleN,:);
%% perform exponential correction
        FitData=filtfilt(Params.sos_Fit,Params.g_Fit,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        FiltData=filtfilt(Params.sos_Low,Params.g_Low,RawData); %Low pass filter data below 0.05Hz before fitting to remove metabolic clearance/photobleaching trends
        speedaccBinary = speedProcess_Doric(forceSensor, Params.DataFs); % for detrend purpose only
        Params.ExcludeVals = find(speedaccBinary==1); % exclude running data for trend fitting
        %
        Spacing = 1:length(FiltData(:,1)); % sample index
        figTime=(1:length(FiltData(:,1)))/(Params.DataFs*60);
    %% Correct TRITC blood volume Exp
        [fitVals]=fit(Spacing',FiltData(:,3),'exp2','Exclude',Params.ExcludeVals);

%         [fitVals]=fit(Spacing',FiltData(:,3),'poly1','Exclude',Params.ExcludeVals);
%           [fitVals,gof,fitinfo]=fit(Spacing',FiltData(:,3),'poly3','Exclude',Params.ExcludeVals);


        coeffVals=coeffvalues(fitVals);

        predicted.F560=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
        
%         predicted.F560=(coeffVals(1).*Spacing)+(coeffVals(2));
        
%         predicted.F560=(coeffVals(1).*(Spacing.^2))+(coeffVals(2).*Spacing)+(coeffVals(3));
        % Plot the exponential fit
        figure; 
        h(1) = subplot(211);
        plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(predicted.F560,'constant')); plot(figTime,detrend(FitData(:,3),'constant')); %plot(figTime,WheelData);
        title([ 'Detreding of TRITC metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw TRITC brightness','Fit Line','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');

        expCorrected.F560=FiltData(:,3)-predicted.F560';

        h(2)=subplot(212);
        plot(figTime,detrend(RawData(:,3),'constant')); hold on; plot(figTime,detrend(expCorrected.F560,'constant')); 
        title(['Detreding of TRITC metabolism ' ChannelName]); xlabel('Time (min)'); legend({'Raw TRITC brightness','expCorrected TRITC'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');
        linkaxes(h,'x');
        if ~isfolder('../Figures/Corrections/')
            mkdir('../Figures/Corrections/')
        end
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_560_' ChannelName '.fig'],'fig')
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_560_' ChannelName '.tiff'],'tiff')
        close 
    %% Correct GRAB-465 Exp
        [fitVals]=fit(Spacing',FiltData(:,2),'exp2','Exclude',Params.ExcludeVals);
        coeffVals=coeffvalues(fitVals);
        predicted.F465=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));

%         [fitVals,gof,fitinfo]=fit(Spacing',FiltData(:,2),'poly2','Exclude',Params.ExcludeVals);
%         coeffVals=coeffvalues(fitVals);
%         predicted.F465=(coeffVals(1).*(Spacing.^2))+(coeffVals(2).*Spacing)+(coeffVals(3));

        % Plot the exponential fit
        figure;
        h(1)=subplot(211);
        plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(predicted.F465,'constant')); plot(figTime,detrend(FitData(:,2),'constant')); 
        title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','Exponential Fit','Low pass filtered data fit'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');

        expCorrected.F465=FiltData(:,2)-predicted.F465';

        h(2)=subplot(212);
        plot(figTime,detrend(RawData(:,2),'constant')); hold on; plot(figTime,detrend(expCorrected.F465,'constant')); 
        title(['Exponential fit of 465 ' ChannelName]); xlabel('Time (min)'); legend({'Raw 465 brightness','expCorrected 465'}); xlim([0 figTime(end)]);
        ylabel('Raw Signal Amplitude');
        linkaxes(h,'x');
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_465_' ChannelName '.fig'],'fig')
        saveas(gcf,['../Figures/Corrections/' Params.savepath 'expcorrected_465_' ChannelName '.tiff'],'tiff')
        close
    %% Correct 405 Exp
%         [fitVals]=fit(Spacing',FiltData(:,1),'exp2','Exclude',Params.ExcludeVals);
%         coeffVals=coeffvalues(fitVals);
%         predicted.F405=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));

        [fitVals]=fit(Spacing',FiltData(:,1),'exp2','Exclude',Params.ExcludeVals);
        coeffVals=coeffvalues(fitVals);
        predicted.F405=(coeffVals(1)*exp((coeffVals(2).*Spacing)))+(coeffVals(3)*exp((coeffVals(4).*Spacing)));
  
    
        expCorrected.F405=FiltData(:,1)-predicted.F405';
       
