figure;

subplot(412)
plotTime = (1:length(FiberData.RawData.Raw_RH(:,2)))./(FiberData.params.DataFs*60);
plot(plotTime,FiberData.RawData.Raw_RH(:,2));
xlabel('Time (min)');
ylabel('Autofluorescence');
legend('Autofluorescence');
xlim([1 60])

subplot(414)
plotTime = (1:length(FiberData.RawData.Raw_RH(:,2)))./(FiberData.params.DataFs*60);
plot(plotTime,FiberData.RawData.Raw_RH(:,2));
xlabel('Time (min)');
ylabel('Raw Fluorescence');
legend('GRAB NE');
xlim([1 60])

subplot(413)
plotTime = (1:length(FiberData.RawData.Raw_RH(:,2)))./(FiberData.params.DataFs*60);
plot(plotTime,FiberData.RawData.Raw_RH(:,2));
xlabel('Time (min)');
ylabel('Raw Fluorescence');
legend('GFP');
xlim([1 60])


subplot(411)
plotTime = (1:length(FiberData.RawData.Raw_RH(:,2)))./(FiberData.params.DataFs*60);
plot(plotTime,FiberData.RawData.Raw_RH(:,2));
xlabel('Time (min)');
ylabel('Raw Fluorescence');
legend('Slides');
xlim([1 60])