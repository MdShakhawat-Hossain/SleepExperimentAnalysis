function DLCPupilTrack_Processing_8point(AnimalID)    
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
        TR = table2array(DLCData(:,8:9));
        right = table2array(DLCData(:,11:12));
        RB = table2array(DLCData(:,14:15));
        bot = table2array(DLCData(:,17:18));
        BL = table2array(DLCData(:,20:21));
        left = table2array(DLCData(:,22:23));
        LT = table2array(DLCData(:,14:15));
    
        for i = 1:1:length(top)
        diameter_top(i) = (norm(mid(i,:) - top(i,:)));
        diameter_TR(i) = (norm(mid(i,:) - TR(i,:)));  
        diameter_right(i) = (norm(mid(i,:) - right(i,:)));
        diameter_RB(i) = (norm(mid(i,:) - RB(i,:)));  
        diameter_bot(i) = (norm(mid(i,:) - bot(i,:)));
        diameter_BL(i) = (norm(mid(i,:) - BL(i,:))); 
        diameter_left(i) = (norm(mid(i,:) - left(i,:)));
        diameter_LT(i) = (norm(mid(i,:) - LT(i,:)));  
        end
        %% load the rawData
        load(RawFileName); % load the rawdata to get frame rate and append the pupil tracking information
        FrameRate = RawData.notes.pupilCamSamplingRate;
        Session_time = RawData.notes.trialDuration_long;
        %% plot the data and check
        ImageCheck = figure;
        % set(ImageCheck,'WindowStyle','docked');
        ImageCheck.WindowState = 'minimized';
        plot_time =  (1:length(diameter_top))./30; % seconds
        ax1 = subplot(411);
        plot(plot_time,diameter_top);
        hold on
        plot(plot_time,diameter_TR);
        plot(plot_time,diameter_right);
        plot(plot_time,diameter_RB);
        plot(plot_time,diameter_bot);
        plot(plot_time,diameter_BL);
        plot(plot_time,diameter_left);
        plot(plot_time,diameter_LT);
        xlim([0 Session_time]);
        
        legend('top','TR','right','RB','bot','BL','left','LT');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)')
        title([FigFileName ' DLC Pupil Tracking'])
        %% average pupil diameter
        ax2 = subplot(412);
        diameter = [diameter_top; diameter_TR; diameter_right; diameter_RB; diameter_bot; diameter_BL; diameter_left; diameter_LT];
        
        for il = 1:1:length(diameter)
            diamterNew(:,il) = filloutliers(diameter(:,il),"linear","median"); 
        end

        plot(plot_time,diamterNew);
        xlim([0 Session_time]);

        legend('Diameter (outlier removed)');
        xlabel('Time (s)');
        ylabel('Diameter (pixels)');
        
        avg_diameter = median(diamterNew);    %average diameter    
        %% filter the pupil diameter
        Diameter = medfilt1(avg_diameter, 5);
        ax3 = subplot(413);
        plot(plot_time,Diameter);

        hold on
        Diameter =  filloutliers(Diameter,"linear","movmean",10*FrameRate);
        plot(plot_time,Diameter);

        xlim([0 Session_time]);
        xlabel('Time (s)');
        ylabel('Diameter (pixels)') 
        legend('Pupil Diameter','Outliers Removed');
        %% movement calculation
        Movement = zeros(length(mid),1);
        for i = 2:1:length(mid)
            Movement(i) = (norm(mid(i,:) - mid(i-1,:)));
        end
        Movement = medfilt1(Movement, 15);
        
        ax4 = subplot(414);
        plot(plot_time,Movement);
        
        hold on
        Movement = filloutliers(Movement,"linear","movmean",10*FrameRate);

        plot(plot_time,Movement);
        xlim([0 Session_time]);
        ylim([min(Movement) max(Movement)])
        
        legend('Pupil Displacement','Outlier removed');
        xlabel('Time (s)');
        ylabel('Displacement (pixels)')
        linkaxes([ax1,ax2,ax3,ax4],'x');
        %% save the figure
        saveas(gcf,['../Figures/Corrections/' FigFileName  '_PupilDiameter.fig'],'fig')
        saveas(gcf,['../Figures/Corrections/' FigFileName  '_PupilDiameter.tiff'],'tiff')
        close(ImageCheck) 
        %% add the new Pupil Diameter to RawData
        RawData.data.Pupil.diameter = Diameter;
        RawData.data.Pupil.movement = Movement;
        
        save(RawFileName,'RawData','-v7.3')
    
        clearvars -except CSVList CN AnimalID
    end


end