function [CorrectSlope , Corrected_F465 ] = HemodynamicsCalculation(Data,Params,FiberPath,HemiSphere)
% clc; close all; clear all

    saveNameIx = strfind(FiberPath,'.');
    saveName = FiberPath(1:saveNameIx(1)-1);
    FigName = [saveName '_' HemiSphere '_hemodynamicCorrection.fig'];
    if ~isfolder('../Figures/Corrections/')
                mkdir('../Figures/Corrections/')
    end 
%     AniamIx = strfind(FiberPath,'_');
%     AnimalID = FiberPath(1:AniamIx(1)-1);
    %%
    figTime=(1:length(Data.F560))/(Params.DataFs*60);

        
    p = polyfit(Data.F560, Data.F465,1);
    predicted.F465 = 10*polyval(p,Data.F560);

    figure; plot(predicted.F465);
        
    % [fitVals]=fit(IsosCorrectData.F560,IsosCorrectData.F465,'exp2');

    % coeffVals=coeffvalues(fitVals);

    % predicted.F465=(coeffVals(1)*exp((coeffVals(2).*IsosCorrectData.F560)))+(coeffVals(3)*exp((coeffVals(4).*IsosCorrectData.F560)));
        

    Corrected_F465=Data.F465-predicted.F465;
    ImageCheck = figure;
    ImageCheck.WindowState = 'minimized'; 

    h(1) = subplot(311);plot(figTime,Data.F560); hold on;plot(figTime,Data.F465); 
    ylabel('Time (minutes)');
    xlabel('\Delta F/F (Z)');
    legend('F560','F465');
    title('Pre-correction')
    xlim([0 figTime(end)])

    h(2) = subplot(312);plot(figTime,Data.F560); hold on;plot(figTime,Corrected_F465); 
   xlabel('Time (minutes)');
   ylabel('\Delta F/F (Z)');
    legend('F560','F465');
    title('Post-correction')
    xlim([0 figTime(end)])

    h(3) = subplot(313);plot(figTime,Data.F465); hold on;plot(figTime,Corrected_F465); 
    ylabel('Time (minutes)');
    xlabel('\Delta F/F (Z)');
    legend('Pre-F465','Post-F465');
    title('465 comparison')
    xlim([0 figTime(end)])

    % eq_str = sprintf('y = %.2f*exp(%.2f*x) + %.2f*exp(%.2f*x)', coeffVals);
    % text(5, 4, eq_str, ... 
    % 'FontSize', 10)

    figure; 

    linkaxes(h,'x')

    saveas(gcf,['../Figures/Corrections/' FigName],'fig')
    saveas(gcf,['../Figures/Corrections/' FigName],'tiff')
    
    close all
    CorrectSlope = coeffVals;
% end
