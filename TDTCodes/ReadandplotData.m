close all; clear
Data = TDTbin2mat('N:\TDT_Tanks\2022-06-24_GRAB-NE_Mouse_004_waterbath\GRAB_NE_004-220624-114320');
AnalogWheel = double(Data.streams.Ball.data);
plotTime =  (1:length(Data.streams.Ball.data))/Data.streams.Ball.fs;
ax(1) = subplot(411);
plot(plotTime,AnalogWheel); 
xlabel('Time(S)')
legend('Ball')


ax(4)=subplot(414);

Markersss = ones(1,length(Data.epocs.Note.onset));

 plot(Data.epocs.Note.onset,Markersss*Data.streams.Ball.fs,'*')
 hold on;
 plot(Data.epocs.Note.offset,Markersss*Data.streams.Ball.fs,'|')
 legend('Gas On', 'Gas Off')
 


ax(2) = subplot(412);

 plot(plotTime,Data.streams.x465A.data)
xlabel('Time(S)')
legend('465 EX')
title('S1 Cortex Fiberphotometry')

ax(3) = subplot(413);
 plot(plotTime,Data.streams.x560B.data)
xlabel('Time(S)')
legend('560EX')
title('S1 Cortex Fiberphotometry')

linkaxes(ax,'x');