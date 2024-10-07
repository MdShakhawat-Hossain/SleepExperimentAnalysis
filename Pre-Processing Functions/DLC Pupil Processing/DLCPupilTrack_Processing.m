function DLCPupilTrack_Processing(AnimalID)    
% clear all; clc; close all;
% Auhtor: Md Shakhawat Hossain
% Updated: March 2024
% purpose: Process the tracked locations from CSV files


    if ~isfolder('../Figures/Corrections/')
        mkdir('../Figures/Corrections/')
    end
    %% list eh CSV files in the folder
    CSVList = dir('*.csv');
    % FileLocation = CSVList(1).folder;
    % LocationFind = strfind(FileLocation,'\');
    % AnimalID = FileLocation(LocationFind(end-2)+1:LocationFind(end-1)-1);
    %%
    for CN =  1:1:length(CSVList)
        PupilFile = CSVList(CN).name;
        NameFind = strfind(PupilFile,'PupilCam');
        RawFileName = [AnimalID '_' PupilFile(1:NameFind(1)-1) 'RawData.mat'];
        FigFileName = [PupilFile(1:NameFind(1)-1)];
        %%
        DLCData = readtable(PupilFile);
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
        load(RawFileName); % load the rawdata to get frame rate and append the pupil tracking information
        FrameRate = RawData.notes.pupilCamSamplingRate;
        Session_time = RawData.notes.trialDuration_long;
        %% plot the data and check
        ImageCheck = figure;
        plot_time =  (1:length(diameter_tb))./30; % seconds
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
        title([FigFileName ' DLC Pupil Tracking'])
        %% average pupil diameter
        ax2 = subplot(312);
        avg_diameter = (diameter_tb+diameter_rl+diameter_mt+diameter_mr+diameter_mb+diameter_ml)/6;
        plot(plot_time,avg_diameter);
        xlim([0 Session_time]);
        
        legend('Avg Diameter');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        %% filter the pupil diameter
        Diameter = medfilt1(avg_diameter, 5);
        ax3 = subplot(313);
        plot(plot_time,Diameter);
        xlim([0 Session_time]);
        
        legend('Pupil Diameter');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        linkaxes([ax1,ax2,ax3],'x');
        %%  
        saveas(gcf,['../Figures/Corrections/' FigFileName  '_PupilDiameter.fig'],'fig')
        saveas(gcf,['../Figures/Corrections/' FigFileName  '_PupilDiameter.tiff'],'tiff')
        close(ImageCheck) 
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
        % plot(plot_time,RawData.data.Pupil.pupilArea);
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
        %% add the new Pupil Diameter to RawData
        RawData.data.Pupil.diameter = Diameter;
        save(RawFileName,'RawData','-v7.3')
    
        clearvars -except CSVList CN AnimalID
    end


end