% function Analyze_microarousal_triggeredResponse()

clear all; clc; close all
load('NEACh001_MicroArousalsData.mat');

DataSize = size(MicroArousalsData.data.GFP.Z_NE,1);

ExpectedLength = 2*30; % 5 seconds Data
for MA = 1:1:DataSize
    MAData_NE(:,MA) = (MicroArousalsData.data.GFP.Z_NE{MA,1}(1:ExpectedLength) - mean(MicroArousalsData.data.GFP.Z_NE{MA,1}(1:15)));%/std(MicroArousalsData.data.GFP.Z_NE{MA,1}(1:15));
    MAData_ACh(:,MA) = (MicroArousalsData.data.GFP.Z_Ach{MA,1}(1:ExpectedLength) - mean(MicroArousalsData.data.GFP.Z_Ach{MA,1}(1:15)));%/std(MicroArousalsData.data.GFP.Z_Ach{MA,1}(1:15));
    MAData_CBV(:,MA) = (MicroArousalsData.data.Rhodamine.Z_Ach{MA,1}(1:ExpectedLength)  - mean(MicroArousalsData.data.Rhodamine.Z_Ach{MA,1}(1:15)));%/std(MicroArousalsData.data.Rhodamine.Z_Ach{MA,1}(1:15));
    MAData_EMG(:,MA) = (MicroArousalsData.data.EMG.emgSignal{MA,1}(1:ExpectedLength) - mean(MicroArousalsData.data.EMG.emgSignal{MA,1}(1:15)));%/std(MicroArousalsData.data.EMG.emgSignal{MA,1}(1:15));
end

figure;
pTime = (1:ExpectedLength)/30;
subplot(411)
for MA = 1:1:DataSize
    hold on
    plot(pTime,MAData_NE(:,MA),'color',[0.4660 0.6740 0.1880]);
    % plot(MAData_ACh(:,MA),'color',[0 0.4470 0.7410]);
    % plot(MAData_CBV(:,MA),'color',[0.6350 0.0780 0.1840]);
    xline(1);
    xlabel('time(s)');
    ylabel('NE');
end
subplot(412)
for MA = 1:1:DataSize
    hold on
    % plot(MAData_NE(:,MA),'color',[0.4660 0.6740 0.1880]);
    plot(pTime,MAData_ACh(:,MA),'color',[0 0.4470 0.7410]);
    % plot(MAData_CBV(:,MA),'color',[0.6350 0.0780 0.1840]);
    xline(1);
    ylabel('ACh');
end
subplot(413)
for MA = 1:1:DataSize
    hold on
    % plot(MAData_NE(:,MA),'color',[0.4660 0.6740 0.1880]);
    % plot(MAData_ACh(:,MA),'color',[0 0.4470 0.7410]);
    plot(pTime,MAData_CBV(:,MA),'color',[0.6350 0.0780 0.1840]);
    xline(1);
    ylabel('CBV');
end
subplot(414)
for MA = 5:1:DataSize
    hold on
    % plot(MAData_NE(:,MA),'color',[0.4660 0.6740 0.1880]);
    % plot(MAData_ACh(:,MA),'color',[0 0.4470 0.7410]);
    plot(pTime,MAData_EMG(:,MA),'color','k');
    xline(1);
    ylabel('EMG');
end

% figure;
% 
% 
%     plot(mean(MAData_NE),'color',[0.4660 0.6740 0.1880]);
%     hold on
%     plot(mean(MAData_ACh),'color',[0 0.4470 0.7410]);
%     plot(mean(MAData_CBV),'color',[0.6350 0.0780 0.1840]);
%     plot(mean(MAData_EMG),'color','k');
%     xline(30);
