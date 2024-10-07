% close all;
% clear all; 
% clc
function TDTFiberPhotometry_bilateralGCaMP7s(filepath,rawDataFilespath)
Data = TDTbin2mat(filepath);
%%
    Params.DataFs=round(Data.streams.x560B.fs);
    Params.RawFs = Data.streams.Fitr.fs;
    Params.low_Freq = 1;
    Params.Fit_Freq=0.05;
    Params.plot_Freq=0.1;
    Params.ball_Freq=10;
    Params.decay_Freq=0.001; % 0.0014
    %% Extract the onsets 
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
            timeN = 1; % remove first 15 seconds;
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
            % add 15s back to first 
            ForceSensor_filt = filtfilt(Params.sos_ball,Params.g_ball,ForceSensor);

 %% detect any sudden spikes in the data
             time_N = (1:length(RawData_LH(:,2)))/Params.DataFs;
             %	Outliers are defined as elements more than three 
             % local scaled MAD from the local median over a window length specified by window. This method is also known as a Hampel filter.
             [OR_RawData_LH,TF_Outliers_LH] = filloutliers(RawData_LH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N); % slide using a 0.25 second window
             time_N = (1:length(RawData_RH(:,2)))/Params.DataFs;
             [OR_RawData_RH,TF_Outliers_RH] = filloutliers(RawData_RH,"linear","movmedian",round(1*Params.DataFs),"SamplePoints",time_N);
            %% Correct exponential decay
            [expCorrected_LH, predicted_LH] = Correct_ExponentialDecay_HighPass(OR_RawData_LH,ForceSensor,Params,'LH');
            [expCorrected_RH, predicted_RH] = Correct_ExponentialDecay_HighPass(OR_RawData_RH,ForceSensor,Params,'RH');
            %% Low pass the data 
            lowPassData_LH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F405);
            lowPassData_LH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F465);
            lowPassData_LH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_LH.F560);
            lowPassData_RH.F405=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F405);
            lowPassData_RH.F465=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F465);
            lowPassData_RH.F560=filtfilt(Params.sos_Low,Params.g_Low,expCorrected_RH.F560);
            %% rescale the data from 0 to 1
            for qN=1:size(lowPassData_LH,2)
            rescaleData_LH.F405(:,qN)=rescale(lowPassData_LH.F405(:,qN),0,1); %rescale all data between 0 to 1
            rescaleData_LH.F465(:,qN)=rescale(lowPassData_LH.F465(:,qN),0,1); %rescale all data between 0 to 1
            rescaleData_LH.F560(:,qN)=rescale(lowPassData_LH.F560(:,qN),0,1); %rescale all data between 0 to 1
            end
            for pN=1:size(lowPassData_RH,2)
            rescaleData_RH.F405(:,pN)=rescale(lowPassData_RH.F405(:,pN),0,1); %rescale all data between 0 to 1
            rescaleData_RH.F465(:,pN)=rescale(lowPassData_RH.F465(:,pN),0,1); %rescale all data between 0 to 1
            rescaleData_RH.F560(:,pN)=rescale(lowPassData_RH.F560(:,pN),0,1); %rescale all data between 0 to 1
            end
            %% remove motion related signal
            [motionCorrect_LH] =  Correct_motionartifact(rescaleData_LH,Params);  
            [motionCorrect_RH] =  Correct_motionartifact(rescaleData_RH,Params); 
            %% hemodynamic Corrections
            Params.correctionConstant = -0.23; % used from Kyles Data -0.23
            CBVCorrection_LH = motionCorrect_LH;
            CBVCorrection_RH = motionCorrect_RH;
            CBVCorrection_LH.F465 = CBVCorrection_LH.F465-(Params.correctionConstant*CBVCorrection_LH.F560);
            CBVCorrection_RH.F465 = CBVCorrection_RH.F465-(Params.correctionConstant*CBVCorrection_RH.F560);
            %
            figTime=(1:length(CBVCorrection_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_LH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals LH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,CBVCorrection_RH.F465),'constant')); 
            title('Hemodynamic Corrected GCaMP7s Signals RH'); xlabel('Time (min)'); legend({'GCaMP7s ExpCorrect ','GCaMP7s TRITCCorrect'}); xlim([0 figTime(end)]);
            ylabel('\Delta F');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'CBVCorrectionCompare.tiff'],'tiff')
            close 
%             CBV_attenuation_correction_KLT_001
            %% Correct 465 and 560 using 405
%             [dF_LH, dFF0_LH,dFF0_z_LH] =  control405Correction(CBVCorrection_LH,predicted_LH,Params,'LH');
%             [dF_RH, dFF0_RH,dFF0_z_RH] =  control405Correction(CBVCorrection_RH,predicted_RH,Params,'RH');
            [dF_LH] =  control405Correction(CBVCorrection_LH,Params,'LH');
            [dF_RH] =  control405Correction(CBVCorrection_RH,Params,'RH');
            %% plot Data comparison
            figTime=(1:length(expCorrected_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F465),'constant')); 
            title('Exp Corrected GCaMP7s Signals'); xlabel('Time (min)'); legend({'LH GCaMP7s','RH GCaMP7s'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_LH.F560),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,expCorrected_RH.F560),'constant')); 
            title('Exp Corrected TRITC Signals'); xlabel('Time (min)'); legend({'LH TRITC','RH TRITC'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'ExpCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'ExpCorrectionCompare.tiff'],'tiff')
            close 
            
            %
            figTime=(1:length(dF_LH.F465))/(Params.DataFs*60);
            figure; 
            h(1) = subplot(211);
             plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_LH.F465),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_RH.F465),'constant')); 
            title('Isosbestic Correction GCaMP7s Signals'); xlabel('Time (min)'); legend({'LH GCaMP7s','RH GCaMP7s'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');

            h(2)=subplot(212);
            plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_LH.F560),'constant')); hold on; plot(figTime,detrend(filtfilt(Params.sos_plot,Params.g_plot,dF_RH.F560),'constant')); 
            title('Isosbestic Correction TRITC Signals'); xlabel('Time (min)'); legend({'LH TRITC','RH TRITC'}); xlim([0 figTime(end)]);
            ylabel('\DeltaF');
            linkaxes(h,'x');
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrectionCompare.fig'],'fig')
            saveas(gcf,['../Figures/Corrections/' Params.savepath 'IsosbesticCorrectionCompare.tiff'],'tiff')
            close 
            %% house keeping
            timeN = 210; % remove first 15 seconds;
            sampleN = floor(timeN*Params.DataFs);
          
           % add some dummy values to match the index
            fieldNames = {'F405','F465','F560'};
            for corI=1:1:length(fieldNames)
                    predicted_LH.(fieldNames{corI}) = [predicted_LH.(fieldNames{corI})(1:sampleN); predicted_LH.(fieldNames{corI});predicted_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    predicted_RH.(fieldNames{corI}) = [predicted_RH.(fieldNames{corI})(1:sampleN); predicted_RH.(fieldNames{corI});predicted_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    expCorrected_LH.(fieldNames{corI}) = [expCorrected_LH.(fieldNames{corI})(1:sampleN); expCorrected_LH.(fieldNames{corI})(1:end); expCorrected_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    expCorrected_RH.(fieldNames{corI}) = [expCorrected_RH.(fieldNames{corI})(1:sampleN); expCorrected_RH.(fieldNames{corI})(1:end); expCorrected_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    lowPassData_LH.(fieldNames{corI}) = [lowPassData_LH.(fieldNames{corI})(1:sampleN); lowPassData_LH.(fieldNames{corI})(1:end); lowPassData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    lowPassData_RH.(fieldNames{corI}) = [lowPassData_RH.(fieldNames{corI})(1:sampleN); lowPassData_RH.(fieldNames{corI})(1:end); lowPassData_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    rescaleData_LH.(fieldNames{corI}) = [rescaleData_LH.(fieldNames{corI})(1:sampleN); rescaleData_LH.(fieldNames{corI})(1:end); rescaleData_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    rescaleData_RH.(fieldNames{corI}) = [rescaleData_RH.(fieldNames{corI})(1:sampleN); rescaleData_RH.(fieldNames{corI})(1:end); rescaleData_RH.(fieldNames{corI})(end-sampleN-1:end)];
            end
            
            fieldNames = {'F465','F560'};
            for corI=1:1:length(fieldNames)
                    motionCorrect_LH.(fieldNames{corI}) = [motionCorrect_LH.(fieldNames{corI})(1:sampleN); motionCorrect_LH.(fieldNames{corI})(1:end); motionCorrect_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    motionCorrect_RH.(fieldNames{corI}) = [motionCorrect_RH.(fieldNames{corI})(1:sampleN); motionCorrect_RH.(fieldNames{corI})(1:end); motionCorrect_RH.(fieldNames{corI})(end-sampleN-1:end)];
                    
                    CBVCorrection_LH.(fieldNames{corI}) = [CBVCorrection_LH.(fieldNames{corI})(1:sampleN); CBVCorrection_LH.(fieldNames{corI})(1:end); CBVCorrection_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    CBVCorrection_RH.(fieldNames{corI}) = [CBVCorrection_RH.(fieldNames{corI})(1:sampleN); CBVCorrection_RH.(fieldNames{corI})(1:end); CBVCorrection_RH.(fieldNames{corI})(end-sampleN-1:end)];

                    dF_LH.(fieldNames{corI}) = [dF_LH.(fieldNames{corI})(1:sampleN); dF_LH.(fieldNames{corI})(1:end); dF_LH.(fieldNames{corI})(end-sampleN-1:end)];
                    dF_RH.(fieldNames{corI}) = [dF_RH.(fieldNames{corI})(1:sampleN); dF_RH.(fieldNames{corI})(1:end); dF_RH.(fieldNames{corI})(end-sampleN-1:end)];

            end
            %% save the data
            
            FiberData.LH.ModrawData.F465 = Mod_Raw_LH_long(:,1);
            FiberData.RH.ModrawData.F465 = Mod_Raw_RH_long(:,1);
            FiberData.LH.ModrawData.F560 = Mod_Raw_LH_long(:,2);
            FiberData.RH.ModrawData.F560 = Mod_Raw_RH_long(:,2);


            FiberData.LH.OrawData.F405 = RawData_LH(:,1);
            FiberData.LH.OrawData.F465 = RawData_LH(:,2);
            FiberData.LH.OrawData.F560 = RawData_LH(:,3);
            
            FiberData.RH.OrawData.F405 = RawData_RH(:,1);
            FiberData.RH.OrawData.F465 = RawData_RH(:,2);
            FiberData.RH.OrawData.F560 = RawData_RH(:,3);

            FiberData.LH.rawData.F405 = OR_RawData_LH(:,1);
            FiberData.LH.rawData.F465 = OR_RawData_LH(:,2);
            FiberData.LH.rawData.F560 = OR_RawData_LH(:,3);
            FiberData.RH.rawData.F405 = OR_RawData_RH(:,1);
            FiberData.RH.rawData.F465 = OR_RawData_RH(:,2);
            FiberData.RH.rawData.F560 = OR_RawData_RH(:,3);
            
            FiberData.pressureSensor  = ForceSensor;
            FiberData.pressureSensor_filt  = ForceSensor_filt;
            FiberData.LH.expCorrected = expCorrected_LH;
            FiberData.RH.expCorrected = expCorrected_RH;
            FiberData.LH.predicted = predicted_LH;
            FiberData.RH.predicted = predicted_RH;
            FiberData.LH.lowPassData = lowPassData_LH;
            FiberData.RH.lowPassData = lowPassData_RH;
            FiberData.LH.rescaleData= rescaleData_LH;
            FiberData.RH.rescaleData = rescaleData_RH;
            FiberData.LH.motionCorrect= motionCorrect_LH;
            FiberData.RH.motionCorrect = motionCorrect_RH;
            FiberData.LH.TRITCCorrection = CBVCorrection_LH;
            FiberData.RH.TRITCCorrection = CBVCorrection_RH;

            FiberData.LH.dF = dF_LH;
            FiberData.RH.dF = dF_RH;
            Params.outliersPos.LH = TF_Outliers_LH;
            Params.outliersPos.RH = TF_Outliers_RH;
            FiberData.params = Params;
            clearvars -except -regexp Data$ path$ Params$ trialNum$
            saveExt = [rawDataFilespath{trialNum}];
            strIdx = strfind(saveExt,'RawData');
            savepathmat =  [saveExt(1:strIdx-1) 'FiberData.mat'];
            save(savepathmat,"FiberData",'-v7.3')
    end