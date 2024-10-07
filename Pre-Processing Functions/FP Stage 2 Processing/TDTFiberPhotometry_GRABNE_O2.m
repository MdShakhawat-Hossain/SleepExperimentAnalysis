% close all;
% clear all; 
% clc
function TDTFiberPhotometry_GRABNE_O2(filepath,rawDataFilespath)
Data = TDTbin2mat(filepath);
%%
    Params.DataFs=round(Data.streams.x560B.fs);
    Params.RawFs = Data.streams.Fitr.fs;
    Params.low_Freq = 1;
    Params.Fit_Freq=0.05;
    Params.plot_Freq=0.1;
    Params.ball_Freq=10;
    Params.decay_Freq=0.001; % 0.001
    %% Extract the onsets 
    Params.Gas.Onset = Data.epocs.Note.onset; % get the time when gas was either turned off or on
    Params.Gas.notes = Data.epocs.Note.notes; % the note saying whether it was on of off
    Params.DataOnset = Data.epocs.Top_.onset; % get the onset from the top epoc
    Params.DataOffset = Data.epocs.Bot_.offset;% get the offset from the bot epoc
%% Filter Parameters
    [z,p,k]=butter(3,Params.low_Freq/(0.5*Params.DataFs),'low'); %Low pass for optical data to physiologically relevant range
    [Params.sos_Low,Params.g_Low]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.ball_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for locomotion data
    [Params.sos_ball,Params.g_ball]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.Fit_Freq/(0.5*Params.DataFs),'low'); %design lowpass filter for hemodynamic correction
    [Params.sos_Fit,Params.g_Fit]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.plot_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for plot
    [Params.sos_plot,Params.g_plot]=zp2sos(z,p,k);

    [z,p,k]=butter(3,Params.decay_Freq/(0.5*Params.DataFs),'low'); %Low pass filter for long term decay
    [Params.sos_decay,Params.g_decay]=zp2sos(z,p,k);

    for trialNum = 1:1:length(Params.DataOnset)
        %% generate file save path for each trials
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            Params.savepath =  [saveExt(1:strIdx-1) 'FiberData'];
        %% Fiber signals
            Raw_LH = [double(Data.streams.x405A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x465A.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x560B.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
            Raw_RH = [double(Data.streams.x405C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x465C.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));
                double(Data.streams.x560D.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))));]';
            
            ForceSensor_raw = double(Data.streams.Forc.data(round(Params.DataOnset(trialNum)*Params.DataFs)+1:round(Params.DataFs*Params.DataOffset(trialNum))))';            
            %% Collected signal 
            Mod_Raw_LH_long = [double(Data.streams.Fitr.data(:,round(Params.DataOnset(trialNum)*Data.streams.Fitr.fs)+1:round(Data.streams.Fitr.fs*Params.DataOffset(trialNum))))]';
            Mod_Raw_RH_long = [double(Data.streams.Fipr.data(:,round(Params.DataOnset(trialNum)*Data.streams.Fitr.fs)+1:round(Data.streams.Fitr.fs*Params.DataOffset(trialNum))))]';
            resample_freq = 100;
            Mod_Raw_RH = detrend(resample(Mod_Raw_RH_long(round(180*Data.streams.Fitr.fs):end-round(180*Data.streams.Fitr.fs),:),resample_freq, round(Data.streams.Fitr.fs)),'constant');
            Mod_Raw_LH = detrend(resample(Mod_Raw_LH_long(round(180*Data.streams.Fitr.fs):end-round(180*Data.streams.Fitr.fs),:),resample_freq, round(Data.streams.Fitr.fs)),'constant');
            ForceSensor_raw_N = detrend(resample(ForceSensor_raw(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');
%% 
            if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
            end
            MplotTime = (1:length(Mod_Raw_LH(:,1)))/(resample_freq*60);
            NPlotTime = (1:length(ForceSensor_raw_N))/(resample_freq*60);
            figure;
            h(1)=subplot(511);
            plot(NPlotTime,ForceSensor_raw_N)
            legend('Force Sensor')
            h(2)=subplot(512);
            plot(MplotTime,Mod_Raw_LH(:,1))
            legend('LH Raw 465')
            title('Raw signals recorded from the fibers')
            h(3)= subplot(513);
            plot(MplotTime,Mod_Raw_LH(:,2))
            legend('LH Raw 560')
            h(4)=subplot(514);
            plot(MplotTime,Mod_Raw_RH(:,1))
            legend('RH Raw 465')
            h(5)=subplot(515);
            plot(MplotTime,Mod_Raw_RH(:,2))
            legend('RH Raw 560')
            ylabel('Raw signals')
            xlabel('Time (min)')
            linkaxes(h,'x');
             xlim([1 54])
%             xlim([1 60])
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawSignal.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawSignal.tiff'],'tiff')
            close 
            
            figure;

            resample_freq = 100;
            Raw_LH_N = detrend(resample(Raw_LH(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');
            Raw_RH_N = detrend(resample(Raw_RH(round(180*Params.DataFs):end-round(180*Params.DataFs),:),resample_freq, round(Params.DataFs)),'constant');
            NPlotTime = (1:length(ForceSensor_raw_N))/(resample_freq*60);
            h(1)=subplot(711);
            plot(NPlotTime,ForceSensor_raw_N)
            legend('Force Sensor')
            title('Raw signals recorded from the fibers')
            h(2)=subplot(712);
            plot(NPlotTime,Raw_LH_N(:,1))
            legend('LH Raw 405')
            h(3)=subplot(713);
            plot(NPlotTime,Raw_LH_N(:,2))
            legend('LH Raw 465')
            h(4)= subplot(714);
            plot(NPlotTime,Raw_LH_N(:,3))
            legend('LH Raw 560')

            h(5)=subplot(715);
            plot(NPlotTime,Raw_RH_N(:,1))
            legend('RH Raw 405')
            h(6)=subplot(716);
            plot(NPlotTime,Raw_RH_N(:,2))
            legend('RH Raw 465')
            h(7)=subplot(717);
            plot(NPlotTime,Raw_RH_N(:,3))
            legend('RH Raw 560')
            ylabel('Raw signals')
            xlabel('Time (min)')
            linkaxes(h,'x');
            xlim([1 55])
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawData.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath '_rawData.tiff'],'tiff')
            close 
            %% remove data from the front to fix the exponential decay
            timeN = 1; % do not remove any data;
            sampleN = floor(timeN*Params.DataFs);
            RawData_LH = Raw_LH;
            RawData_RH = Raw_RH;
            ForceSensor = ForceSensor_raw;
            % adding dummy data to match index
            RawData_LH(1:sampleN,:) = Raw_LH(sampleN+1:sampleN+sampleN,:);
            RawData_RH(1:sampleN,:) = Raw_RH(sampleN+1:sampleN+sampleN,:);
            ForceSensor(1:sampleN,:) = ForceSensor_raw(sampleN+1:sampleN+sampleN,:);

            RawData_LH(sampleN+1:end,:) = Raw_LH(sampleN+1:end,:);
            RawData_RH(sampleN+1:end,:) = Raw_RH(sampleN+1:end,:);
            ForceSensor(sampleN+1:end,:) = ForceSensor_raw(sampleN+1:end,:);
            % add the data back to first 
            ForceSensor_filt = filtfilt(Params.sos_ball,Params.g_ball,ForceSensor);
 %% detect any sudden spikes in the data
             time_N = (1:length(RawData_LH(:,2)))/Params.DataFs;
             %	Outliers are defined as elements more than three 
             % local scaled MAD from the local median over a window length specified by window. This method is also known as a Hampel filter.
             [OR_RawData_LH,TF_Outliers_LH] = filloutliers(RawData_LH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N); % slide using a 0.25 second window
             time_N = (1:length(RawData_RH(:,2)))/Params.DataFs;
             [OR_RawData_RH,TF_Outliers_RH] = filloutliers(RawData_RH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N);
            %% Correct exponential decay
            removalTime = 30; % the duration of the data you want to remove from analysis
            [expCorrected_LH, ~] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_LH,ForceSensor,Params,'Ach',removalTime);
            [expCorrected_RH, ~] = Correct_ExponentialDecay_HighPass_Update(OR_RawData_RH,ForceSensor,Params,'NE',removalTime);
            %% Low pass the data 
            lowPassData_LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F405);
            lowPassData_LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F465);
            lowPassData_LH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F560);
            lowPassData_RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F405);
            lowPassData_RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F465);
            lowPassData_RH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F560);
                %% zScore the signal
            % determine the resting state baseline to calculate ZScore
            if isfield(Params,'baselineStartTime') == false 
                %{
                FrcSensor = ForceSensor((210*Params.DataFs)+1:end-(210*Params.DataFs));
                ImageCheck = figure;
                ImageCheck.WindowState = 'maximized';
                drawnow;
                L(1) = subplot(3,1,1);
                plot((1:length(FrcSensor))/Params.DataFs,FrcSensor);
                hold on 
                yyaxis right
                plot((1:length(lowPassData_RH.F405))/Params.DataFs,lowPassData_RH.F405);
                legend('Force Sensor', 'RH 405')
    
                xlim([0 length(FrcSensor)/Params.DataFs]);
                ylabel('forceSensor');
                L(2) = subplot(3,1,2);
                plot((1:length(lowPassData_LH.F560))/Params.DataFs,lowPassData_LH.F560);
                ylabel('dF LH 560');            
                xlim([0 length(lowPassData_RH.F560)/Params.DataFs]);
    
                L(3) = subplot(3,1,3);
                plot((1:length(lowPassData_RH.F560))/Params.DataFs,lowPassData_RH.F560);
                ylabel('dF RH 560');
                xlim([0 length(lowPassData_RH.F560)/Params.DataFs]);
                linkaxes(L,'x');
                commandwindow;
                Params.baselineStartTime = input('Input the start time for resting data: '); disp(' ')
                Params.baselineEndTime = input('Input the end time for resting data: '); disp(' ')
                close(ImageCheck) 
                %}
                Params.baselineStartTime = 500;
                Params.baselineEndTime = 3000; 
            end

            [zScored_LH] =  ZScoreFiberData(lowPassData_LH,Params,'LH');
            [zScored_RH] =  ZScoreFiberData(lowPassData_RH,Params,'RH');            
            %% remove motion related signal
            [motionCorrect_LH] =  Correct_motionartifact(zScored_LH,Params);  
            [motionCorrect_RH] =  Correct_motionartifact(zScored_RH,Params); 
            %% remove isosbestic changes related signal
            [isosCorrect_LH] = control405Correction_Update(motionCorrect_LH,Params,'LH');
            [isosCorrect_RH] = control405Correction_Update(motionCorrect_RH,Params,'RH');
            
            Fiberpath = [Params.savepath '.mat'];
            [Params.CorrectSlope_LH , Corrected_F465_LH ] = HemodynamicsCalculation(isosCorrect_LH,Params,Fiberpath,'LH');
            [Params.CorrectSlope_RH , Corrected_F465_RH ] = HemodynamicsCalculation(isosCorrect_RH,Params,Fiberpath,'RH');

            CBVCorrection_LH = isosCorrect_LH;
            CBVCorrection_RH = isosCorrect_RH;

            CBVCorrection_LH.F465 = Corrected_F465_LH;
            CBVCorrection_RH.F465 = Corrected_F465_RH;

   
            figTime=(1:length(CBVCorrection_LH.F465))/(Params.DataFs*60);

            ImageCheck = figure;
            ImageCheck.WindowState = 'minimized';

            k(1) = subplot(4,1,1);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_RH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_RH.F560))
            legend('RH 465','RH 560')
            title('Pre-Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\delta F rescaled');

            k(2) = subplot(4,1,3);
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_LH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,isosCorrect_LH.F560))
            legend('LH 465','LH 560')
            title('Before Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\delta F rescaled');

            k(3)=subplot(4,1,2); 
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F560))
            legend('RH 465','RH 560')
            title('After Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F rescaled');

            k(4)=subplot(4,1,4); 
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F465))
            hold on;
            plot(figTime,filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F560))
            legend('LH 465','LH 560')
            title('After Hemodynamic Correction'); xlabel('Time (min)'); xlim([0 figTime(end)]);
            ylabel('\Delta F rescaled');

            linkaxes(k,'x');

            saveas(gcf,['../Figures/Corrections/' Params.savepath 'HemoCorrection.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'HemoCorrection.tiff'],'tiff')
            close(ImageCheck)
            %% house keeping
            timeN = removalTime; % remove first few seconds;
            sampleN = floor(timeN*Params.DataFs);
          
           % add some dummy values to match the index
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                     expCorrected_LH.(fieldNames{corI}) = [expCorrected_LH.(fieldNames{corI})(1:sampleN); expCorrected_LH.(fieldNames{corI})(1:end); expCorrected_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    expCorrected_RH.(fieldNames{corI}) = [expCorrected_RH.(fieldNames{corI})(1:sampleN); expCorrected_RH.(fieldNames{corI})(1:end); expCorrected_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    lowPassData_LH.(fieldNames{corI}) = [lowPassData_LH.(fieldNames{corI})(1:sampleN); lowPassData_LH.(fieldNames{corI})(1:end); lowPassData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    lowPassData_RH.(fieldNames{corI}) = [lowPassData_RH.(fieldNames{corI})(1:sampleN); lowPassData_RH.(fieldNames{corI})(1:end); lowPassData_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    zScored_LH.(fieldNames{corI}) = [zScored_LH.(fieldNames{corI})(1:sampleN); zScored_LH.(fieldNames{corI})(1:end); zScored_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    zScored_RH.(fieldNames{corI}) = [zScored_RH.(fieldNames{corI})(1:sampleN); zScored_RH.(fieldNames{corI})(1:end); zScored_RH.(fieldNames{corI})(end-sampleN-1:end)];
            end
            
            fieldNames = {'F465','F560'};
            for corI=1:1:length(fieldNames)
                    motionCorrect_LH.(fieldNames{corI}) = [motionCorrect_LH.(fieldNames{corI})(1:sampleN); motionCorrect_LH.(fieldNames{corI})(1:end); motionCorrect_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    motionCorrect_RH.(fieldNames{corI}) = [motionCorrect_RH.(fieldNames{corI})(1:sampleN); motionCorrect_RH.(fieldNames{corI})(1:end); motionCorrect_RH.(fieldNames{corI})(end-sampleN-1:end)];
                    
                    CBVCorrection_LH.(fieldNames{corI}) = [CBVCorrection_LH.(fieldNames{corI})(1:sampleN); CBVCorrection_LH.(fieldNames{corI})(1:end); CBVCorrection_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    CBVCorrection_RH.(fieldNames{corI}) = [CBVCorrection_RH.(fieldNames{corI})(1:sampleN); CBVCorrection_RH.(fieldNames{corI})(1:end); CBVCorrection_RH.(fieldNames{corI})(end-sampleN-1:end)];
            end
            %% save the data
            FiberData.Ach.ModrawData.F465 = Mod_Raw_LH_long(:,1);
            FiberData.NE.ModrawData.F465 = Mod_Raw_RH_long(:,1);
            FiberData.Ach.ModrawData.F560 = Mod_Raw_LH_long(:,2);
            FiberData.NE.ModrawData.F560 = Mod_Raw_RH_long(:,2);


            FiberData.Ach.OrawData.F405 = RawData_LH(:,1);
            FiberData.Ach.OrawData.F465 = RawData_LH(:,2);
            FiberData.Ach.OrawData.F560 = RawData_LH(:,3);
            
            FiberData.NE.OrawData.F405 = RawData_RH(:,1);
            FiberData.NE.OrawData.F465 = RawData_RH(:,2);
            FiberData.NE.OrawData.F560 = RawData_RH(:,3);

            FiberData.Ach.rawData.F405 = OR_RawData_LH(:,1);
            FiberData.Ach.rawData.F465 = OR_RawData_LH(:,2);
            FiberData.Ach.rawData.F560 = OR_RawData_LH(:,3);
            FiberData.NE.rawData.F405 = OR_RawData_RH(:,1);
            FiberData.NE.rawData.F465 = OR_RawData_RH(:,2);
            FiberData.NE.rawData.F560 = OR_RawData_RH(:,3);
            
            FiberData.pressureSensor  = ForceSensor;
            FiberData.pressureSensor_filt  = ForceSensor_filt;


            FiberData.Ach.expCorrected = expCorrected_LH;
            FiberData.NE.expCorrected = expCorrected_RH;

            FiberData.Ach.lowPassData = lowPassData_LH;
            FiberData.NE.lowPassData = lowPassData_RH;

            FiberData.Ach.zScored= zScored_LH;
            FiberData.NE.zScored = zScored_RH;

            FiberData.Ach.motionCorrect= motionCorrect_LH;
            FiberData.NE.motionCorrect = motionCorrect_RH;

            FiberData.Ach.RhodamineCorrection = CBVCorrection_LH;
            FiberData.NE.RhodamineCorrection = CBVCorrection_RH;

            FiberData.Ach.dFF0_z = CBVCorrection_LH;
            FiberData.NE.dFF0_z = CBVCorrection_RH;

            Params.outliersPos.Ach = TF_Outliers_LH;
            Params.outliersPos.NE = TF_Outliers_RH;

            FiberData.params = Params;
            
            clearvars -except -regexp Data$ path$ Params$ trialNum$
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            savepathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            save(savepathmat,"FiberData",'-v7.3')
    end