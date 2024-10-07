% function DLCArtery_Processing()    
clear all; clc; close all;
% Auhtor: Md Shakhawat Hossain
% Updated: March 2024
% purpose: Process the tracked locations from CSV files


    % if ~isfolder('../Figures/Corrections/')
    %     mkdir('../Figures/Corrections/')
    % end
    %% list eh CSV files in the folder
    CSVList = dir('*.csv');
    % FileLocation = CSVList(1).folder;
    % LocationFind = strfind(FileLocation,'\');
    % AnimalID = FileLocation(LocationFind(end-2)+1:LocationFind(end-1)-1);
    %%
    for CN =  1:1:length(CSVList)
        ArteryFile = CSVList(CN).name;
        NameFind = strfind(ArteryFile,'DLC');
        % RawFileName = [AnimalID '_' ArteryFile(1:NameFind(1)-1) 'RawData.mat'];
        FigFileName = [ArteryFile(1:NameFind(1)-1)];
        %%
        DLCData = readtable(ArteryFile);
        mid = table2array(DLCData(:,2:3));
        top = table2array(DLCData(:,5:6));
        right = table2array(DLCData(:,8:9));
        bot = table2array(DLCData(:,11:12));
        left = table2array(DLCData(:,14:15));
    
        for i = 1:1:length(top)
        diameter_tb(i) = (norm(top(i,:) - bot(i,:)));
        diameter_rl(i) = (norm(right(i,:) - left(i,:)));   
        diameter_mt(i) = 2*(norm(mid(i,:) - top(i,:)));
        diameter_mr(i) = 2*(norm(mid(i,:) - right(i,:)));
        diameter_mb(i) = 2*(norm(mid(i,:) - bot(i,:)));
        diameter_ml(i) = 2*(norm(mid(i,:) - left(i,:)));
        end
        %% load the rawData
        % load(RawFileName); % load the rawdata to get frame rate and append the Artery tracking information
        FrameRate = 10;%RawData.notes.ArteryCamSamplingRate;
        Session_time = 930;%RawData.notes.trialDuration_long;
        %% plot the data and check
        ImageCheck = figure;
        plot_time =  (1:length(diameter_tb))./FrameRate; % seconds
        ax1 = subplot(311);
        plot(plot_time,diameter_tb);
        hold on
        plot(plot_time,diameter_rl);
        plot(plot_time,diameter_mt);
        plot(plot_time,diameter_mr);
        plot(plot_time,diameter_mb);
        plot(plot_time,diameter_ml);
        xlim([0 Session_time]);
        
        legend('tb','rl','mt','mr','mb','ml');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        title([FigFileName ' DLC Artery Tracking'])
        %% average Artery diameter
        ax2 = subplot(312);
        avg_diameter = (diameter_tb+diameter_rl+diameter_mt+diameter_mr+diameter_mb+diameter_ml)/6;
        plot(plot_time,avg_diameter);
        xlim([0 Session_time]);
        
        legend('Avg Diameter');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        %% filter the Artery diameter
        Diameter = medfilt1(avg_diameter, 5);
        ax3 = subplot(313);
        plot(plot_time,Diameter);
        xlim([0 Session_time]);
        
        legend('Artery Diameter');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        linkaxes([ax1,ax2,ax3],'x');
        %%  
        saveas(gcf,[FigFileName  '_ArteryDiameter.fig'],'fig')
        saveas(gcf,[FigFileName  '_ArteryDiameter.tiff'],'tiff')
        close(ImageCheck) 

        load('VesselDiameterSummary.mat');

        ImageCompare = figure;
        ax1 = subplot(211);
        percentageDiameter = (100*(Diameter - mean(Diameter)))./mean(Diameter);
        plot(plot_time,percentageDiameter);
        xlim([0 Session_time]);
        hold on
        VesselSummary.mscanDataFiles
        FigFileName
        Filenumber = input('What is the radon file number: '); disp(' ')
        radondia = medfilt1(VesselSummary.RadonDiameter{Filenumber, 1},5);
        radonDiameter = (100*(radondia - mean(radondia)))./mean(radondia);
        plot(plot_time,radonDiameter);
        xlim([0 Session_time]);
        
        legend('DLC Diameter','Radon Diameter');
        xlabel('Time (s)');
        ylabel('Diameter (percentage)')

        saveas(gcf,[FigFileName  '_ArteryDiameter_Compare.fig'],'fig')
        saveas(gcf,[FigFileName  '_ArteryDiameter_Compare.tiff'],'tiff')
        close(ImageCompare) 

        %%
        % figure;
        % ax1 = subplot(311);
        % plot(plot_time,Diameter);
        % xlim([0 Session_time]);
        % 
        % xlabel('Time (s)');
        % ylabel('Diameter (pixels)')
        % 
        % ax2 = subplot(312);
        % plot(plot_time,RawData.data.Artery.ArteryArea);
        % xlim([0 Session_time]);
        % 
        % xlabel('Time (s)');
        % ylabel('Diameter (pixels)')
        % 
        % ax3 = subplot(313);
        % plot((1:length(FiberData.NE.dFF0_p.F465))./FiberData.params.DataFs,FiberData.NE.dFF0_p.F465);
        % xlim([0 Session_time]);
        % xlabel('Time (s)');
        % ylabel('NE')
        % linkaxes([ax1,ax2],'x');
        %% add the new Artery Diameter to RawData
        ArteryDiameter.Diameter = Diameter;
        ArteryDiameter.PercentDiameter = percentageDiameter;
        save([FigFileName '_Diameter'],'','-v7.3')
    
        clearvars -except CSVList CN
    end


% end