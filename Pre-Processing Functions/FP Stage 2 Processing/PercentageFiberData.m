function [dFF0_percentage] =  PercentageFiberData(Corrected,Params,ChannelName,RawData,timeN)

    % rmStartIndex = Params.DataFs*Params.baselineStartTime;
    % rmEndIndex = Params.DataFs*Params.baselineEndTime;
     sampleN = floor(timeN*Params.DataFs);
           
    % adding dummy data to match index
    RawData = RawData(sampleN+1:end-sampleN,:);
    %%
    rF.F405 = RawData(:,1);
    rF.F465 = RawData(:,2);
    rF.F560 = RawData(:,3);

    dF.F405 = Corrected.F405;   
    dF.F465 = Corrected.F465;
    dF.F560 = Corrected.F560;

    DataLength = length(dF.F405);
    DataSize =  Params.DataFs*60*5; % 5 minutes bin
    for Dn = 1:1:ceil(DataLength/DataSize)
        if Dn == 1
            DataHolder = 0;
        end
        DataStart = 1 + DataHolder ;
        DataEnd = DataStart + DataSize;
        if DataEnd > DataLength
            DataEnd = DataLength;
        end
        DataHolder = DataEnd;

        % F405
        dFF0_percentage.F405(DataStart:DataEnd) = 100 + ((100.*(dF.F405(DataStart:DataEnd)-mean(rF.F405(DataStart:DataEnd))))./mean(rF.F405(DataStart:DataEnd)));
        % F465
        % temp_F465 = dF.F465(DataStart:DataEnd)./dF.F405(DataStart:DataEnd);

        % dFF0_percentage.F465(DataStart:DataEnd) = (100.*(temp_F465-mean(temp_F465)))./mean(temp_F465);
        dFF0_percentage.F465(DataStart:DataEnd) = 100 + ((100.*(dF.F465(DataStart:DataEnd)-mean(rF.F465(DataStart:DataEnd))))./mean(rF.F465(DataStart:DataEnd)));

        % F560    
        dFF0_percentage.F560(DataStart:DataEnd) = 100 + ((100.*(dF.F560(DataStart:DataEnd)-mean(rF.F560(DataStart:DataEnd))))./mean(rF.F560(DataStart:DataEnd)));

    end
    % %% F405
    % dFF0_percentage.F405 = (100.*(dF.F405-mean(dF.F405)))./mean(dF.F405);
    % %% F465
    % temp_F465 = dF.F465./dF.F405;
    % dFF0_percentage.F465 = (100.*(temp_F465-mean(temp_F465)))./mean(temp_F465);
    % %% F560    
    % dFF0_percentage.F560 = (100.*(dF.F560-mean(dF.F560)))./mean(dF.F560);

   
    
    
    %% plot the data
    figTime=(1:length(dFF0_percentage.F465))/(Params.DataFs*60);
    % figure; 
    h(1) = subplot(311);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F560)); 
    title('Exp Corrected percentagescored CBV'); xlabel('Time (min)'); xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF');

    h(2) = subplot(312);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F465)); 
    title('Exp Corrected percentagescored Green FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF');   

    h(3) = subplot(313);
    plot(figTime,filtfilt(Params.sos_Low,Params.g_Low,dFF0_percentage.F405)); 
    title('Exp Corrected percentagescored 405 FP'); xlabel('Time (min)');  xlim([0 figTime(end)]);
    ylabel('percentageScored \DeltaF'); 

    linkaxes(h);
    saveas(gcf,['../Figures/Corrections/' Params.savepath 'Percentage_' ChannelName '.fig'],'fig')
    close 