%% Model hemodynamic response attenuation of RFP signal
figure(99);scatter(RescaleData(:,2),RescaleData(:,1));
figure(99);scatter(LowPassData(:,2),LowPassData(:,1));
figure(303);scatter(ZscoredFiberData(:,2),ZscoredFiberData(:,1))

%% Hemodynamic response correction code using Mle
BinEdges=(0:0.001:1);
figure(100);
DataHist=histogram2(RescaleData(:,2),RescaleData(:,1),BinEdges,BinEdges);

BinCounts=DataHist.BinCounts;
% ColCounts=sum(BinCounts,2);
allCounts=sum(BinCounts,'all');
% NormMat=repmat(ColCounts,1,size(BinCounts,2));
% NormHist=(BinCounts./NormMat)*100;
allHist=(BinCounts/allCounts)*100;
XCounts=sum(BinCounts,2);
mleFit.allHist=allHist;
for binNum=2:size(BinCounts,1)
    if binNum==2
        percData(1)=(XCounts(1)/sum(XCounts))*100;
    end
    percData(binNum)=(sum(XCounts(1:binNum))/sum(XCounts))*100;
end
startInd=find(percData<=2,1,'last')+1; %exclude lower 2% of data points
endInd=find(percData>=98,1,'first')-1; %exclude upper 2% of data points
mleFit.startInd=startInd;
mleFit.endInd=endInd;

startInt=BinEdges(startInd);
belowThresh=find(RescaleData(:,2)<startInt);
endInt=BinEdges(endInd);
aboveThresh=find(RescaleData(:,2)>endInt);

blankVals=union(belowThresh,aboveThresh);
lm_Data=RescaleData;
lm_Data(blankVals,:)=[];
for binNum=startInd:endInd
    index=(binNum-startInd)+1;
    leftEdge=BinEdges(binNum);
    rightEdge=BinEdges(binNum+1);
    
        dataInds=RescaleData(:,2)>=leftEdge...
            & RescaleData(:,2)<=rightEdge;
        theData=RescaleData(dataInds,1);
    

    
    [phat,pci]=mle(theData,'distrbution','normal');
    figure(604); normHist=histogram(theData,BinEdges);
    normCounts=normHist.BinCounts./sum(normHist.BinCounts);
    theFit=pdf('normal',normHist.BinEdges(1:(end-1)),phat(1),phat(2));
    theFit=theFit./sum(theFit);
    
        %     figure(index);plot(normHist.BinEdges(1:(end-1)),normCounts);
        %     hold on; plot(normHist.BinEdges(1:(end-1)),theFit);
    rsqr=1-(sum((normCounts-theFit).^2)/sum((normCounts-mean(normCounts)).^2));
    
    mleFit.fitParams.phat(index,:)=phat;
    mleFit.fitParams.pci(:,:,index)=pci;
    mleFit.goodness.rsqr(index)=rsqr;
    mleFit.fitData.eGFPData{index}=theData;
    mleFit.fitData.BinEdges(index,:)=[leftEdge,rightEdge];
    mleFit.fitData.normHist(index,:)=normCounts;
    mleFit.fitData.HistFit(index,:)=theFit;
    mleFit.fitData.fitHistEdges(index,:)=normHist.BinEdges;
end

fitPeaks(1:length(BinEdges))=NaN;
fitposStan(1:length(BinEdges))=NaN;
fitnegStan(1:length(BinEdges))=NaN;
fitPeaks(startInd:endInd)=mleFit.fitParams.phat(:,1);
fitposStan(startInd:endInd)=mleFit.fitParams.phat(:,1)+mleFit.fitParams.phat(:,2);
fitnegStan(startInd:endInd)=mleFit.fitParams.phat(:,1)-mleFit.fitParams.phat(:,2);

[linFit_2]=fitlm(BinEdges(startInd:endInd)',mleFit.fitParams.phat(:,1),'RobustOpts','on');
[linFit]=fitlm(lm_Data(:,2),lm_Data(:,1),'RobustOpts','on');
    % linVals=coeffvalues(linFit);
linSlope=table2array(linFit.Coefficients(2,1));
linInt=table2array(linFit.Coefficients(1,1));
slope_pVal=table2array(linFit.Coefficients(2,4));
linPlot=linSlope*BinEdges+linInt;
mleFit.linFit=linFit;

linSlope_2=table2array(linFit_2.Coefficients(2,1));
linInt_2=table2array(linFit_2.Coefficients(1,1));
slope_pVal_2=table2array(linFit_2.Coefficients(2,4));
linPlot_2=linSlope_2*BinEdges+linInt_2;
mleFit.linFit_2=linFit_2;
 

figure(204);subplot(3,3,1); hold on;
imagesc(BinEdges,BinEdges,allHist');
caxis([0 0.002]);
h=colorbar('eastoutside');
axis xy;
ylabel('\Delta RFP (%)');
xlabel('\Delta FITC (%)');

title(' \Delta FITC vs \Delta RFP');
plot(BinEdges,fitPeaks,'r','LineWidth',2);
plot(BinEdges,fitposStan,'--r','LineWidth',1);
plot(BinEdges,fitnegStan,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot,'w','LineWidth',2);
xlim([0 1]);
ylim([0 1]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
slopeTxt=num2str(linSlope);
intTxt=num2str(linInt);
pvalTxt=num2str(slope_pVal);
rsqrTxt=num2str(linFit.Rsquared.Adjusted);
all_an=annotation('textbox',[0.13, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt 'x+' intTxt ' rsqr=' rsqrTxt],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);


figure(204);subplot(3,3,2); hold on;
imagesc(BinEdges,BinEdges,allHist');
caxis([0 0.002]);
h=colorbar('eastoutside');
axis xy;
ylabel('\Delta RFP (%)');
xlabel('\Delta FITC (%)');

title(' \Delta FITC vs \Delta RFP');
plot(BinEdges,fitPeaks,'r','LineWidth',2);
plot(BinEdges,fitposStan,'--r','LineWidth',1);
plot(BinEdges,fitnegStan,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot_2,'w','LineWidth',2);
xlim([0 1]);
ylim([0 1]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
slopeTxt_2=num2str(linSlope_2);
intTxt_2=num2str(linInt_2);
pvalTxt_2=num2str(slope_pVal_2);
rsqrTxt_2=num2str(linFit_2.Rsquared.Adjusted);
MLE_an=annotation('textbox',[0.43, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt_2 'x+' intTxt_2 ' rsqr=' rsqrTxt_2],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);
%%

coeffVals=table2array(mleFit.linFit_2.Coefficients);

y_int=table2array(mleFit.linFit_2.Coefficients(1,1));
the_slope=table2array(mleFit.linFit_2.Coefficients(2,1));
FinalRFP=RescaleData(:,1)-(the_slope*RescaleData(:,2)+y_int);
FinalRFP_Z=(FinalRFP-mean(FinalRFP))/std(FinalRFP);
figure(50);plot(ZscoredFiberData(:,1));hold on; plot(FinalRFP_Z);plot(ZscoredFiberData(:,2));
legend({'Raw tdTomato','HbT Corrected tdTomato','FITC blood plasma'});