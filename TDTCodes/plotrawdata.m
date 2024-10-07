function plotrawdata(fiberData)

            %{
            %% resample the data for plotting
            Params.LowFs = 100; % resample the data to 100 samples per second
            rawData_ch1_resampled = resample(RawData_ch1,Params.LowFs,Params.DataFs);
            rawData_ch2_resampled = resample(RawData_ch2,Params.LowFs,Params.DataFs);
            forcesensor_resampled = resample(ForceSensor,Params.LowFs,Params.DataFs); 
            % exclude data
            excludeSample = Params.DataFs*30; % remove 30 seconds
             
            figTime=(1:length(rawData_ch1_resampled(excludeSample:end-excludeSample,1)))/(Params.LowFs*60);
            figure;
            ay(1)=subplot(411);
            plot(figTime,rawData_ch1_resampled(excludeSample:end-excludeSample,1))
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('405')
        
            ay(2)=subplot(412);
             plot(figTime,rawData_ch1_resampled(excludeSample:end-excludeSample,2))
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('465')
        
            ay(3)=subplot(413);
            plot(figTime,rawData_ch1_resampled(excludeSample:end-excludeSample,3)); 
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('560')

            ay(4)=subplot(414);
            plot(figTime,forcesensor_resampled(excludeSample:end-excludeSample))
            xlabel('Time (min)');ylabel('Activity');
            legend('Ball Motion')
            linkaxes(ay,'x');
        
            figure;
            az(1)=subplot(411);
            plot(figTime,rawData_ch2_resampled(excludeSample:end-excludeSample,1))
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('405')
        
            az(2)=subplot(412);
            plot(figTime,rawData_ch2_resampled(excludeSample:end-excludeSample,2))
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('465')
        
            az(3)=subplot(413);
            plot(figTime,rawData_ch2_resampled(excludeSample:end-excludeSample,3)); 
            xlabel('Time (min)');ylabel('Raw Signal Amplitude');
            title('Raw data before corrections');
            legend('560')
        
            az(4)=subplot(414);
            plot(figTime,forcesensor_resampled(excludeSample:end-excludeSample))
            xlabel('Time (min)');ylabel('Activity');
            legend('Ball Motion')
            linkaxes(az,'x');

            %}