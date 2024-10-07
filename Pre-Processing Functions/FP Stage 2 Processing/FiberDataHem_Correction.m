close all; clc
clearvars -Except FiberData
[z,p,k]=butter(3,1/(0.5*FiberData.params.DataFs),'low');
[sos,g]=zp2sos(z,p,k);
BinEdges=(-2:0.001:2);
figure;DataHist_FC=histogram2(FiberData.NE.rescaleData.F560,FiberData.NE.rescaleData.F465,BinEdges,BinEdges,'DisplayStyle','tile','LineStyle',':');colorbar
BinCounts_FC=DataHist_FC.BinCounts;
allCounts=sum(BinCounts_FC,'all');
allHist_FC=(BinCounts_FC/allCounts)*100;
XCounts=sum(BinCounts_FC,2);
mleFit.FC.allHist_FC=allHist_FC;
for binNum=2:size(BinCounts_FC,1)
    if binNum==2
        percData(1)=(XCounts(1)/sum(XCounts))*100;
    end
    percData(binNum)=(sum(XCounts(1:binNum))/sum(XCounts))*100;
end
startInd_FC=find(percData<=5,1,'last')+1; %exclude lower 2% of data points
endInd_FC=find(percData>=95,1,'first')-1; %exclude upper 2% of data points
mleFit.FC.startInd=startInd_FC;
mleFit.FC.endInd=endInd_FC;

startInt=BinEdges(startInd_FC);
belowThresh=find(FiberData.NE.expCorrected.F560<startInt);
endInt=BinEdges(endInd_FC);
aboveThresh=find(FiberData.NE.expCorrected.F560>endInt);

blankVals=union(belowThresh,aboveThresh);
lm_Data.F560=FiberData.NE.expCorrected.F560;
lm_Data.F465=FiberData.NE.expCorrected.F465;
lm_Data.F560(blankVals)=[];
lm_Data.F465(blankVals)=[]; 
for binNum=startInd_FC:endInd_FC
    index=(binNum-startInd_FC)+1;
    leftEdge=BinEdges(binNum);
    rightEdge=BinEdges(binNum+1);
    
        dataInds=FiberData.NE.expCorrected.F560>=leftEdge...
            & FiberData.NE.expCorrected.F560<=rightEdge;
        theData=FiberData.NE.expCorrected.F465(dataInds);
    
    [phat,pci]=mle(theData,'Distribution','Normal');
%     figure(10); normHist_FC=histogram(theData,BinEdges);
    [normHist_FC.BinCounts, normHist_FC.BinEdges] = histcounts(theData,BinEdges);
    normCounts_FC=normHist_FC.BinCounts./sum(normHist_FC.BinCounts);
    theFit_FC=pdf('normal',normHist_FC.BinEdges(1:(end-1)),phat(1),phat(2));
    theFit_FC=theFit_FC./sum(theFit_FC);
    
    rsqr=1-(sum((normCounts_FC-theFit_FC).^2)/sum((normCounts_FC-mean(normCounts_FC)).^2));
    
    mleFit.FC.fitparams.phat(index,:)=phat;
    mleFit.FC.fitparams.pci(:,:,index)=pci;
    mleFit.FC.goodness.rsqr(index)=rsqr;
    mleFit.FC.fitData.eGFPData{index}=theData;
    mleFit.FC.fitData.BinEdges(index,:)=[leftEdge,rightEdge];
    mleFit.FC.fitData.normHist(index,:)=normCounts_FC;
    mleFit.FC.fitData.HistFit(index,:)=theFit_FC;
    mleFit.FC.fitData.fitHistEdges(index,:)=normHist_FC.BinEdges;
end

fitPeaks(1:length(BinEdges))=NaN;
fitposStan(1:length(BinEdges))=NaN;
fitnegStan(1:length(BinEdges))=NaN;
fitPeaks(startInd_FC:endInd_FC)=mleFit.FC.fitparams.phat(:,1);
fitposStan(startInd_FC:endInd_FC)=mleFit.FC.fitparams.phat(:,1)+mleFit.FC.fitparams.phat(:,2);
fitnegStan(startInd_FC:endInd_FC)=mleFit.FC.fitparams.phat(:,1)-mleFit.FC.fitparams.phat(:,2);

[linFit]=fitlm(lm_Data.F560,lm_Data.F465,'RobustOpts','on');
% linVals=coeffvalues(linFit);
linSlope=table2array(linFit.Coefficients(2,1));
linInt=table2array(linFit.Coefficients(1,1));
slope_pVal=table2array(linFit.Coefficients(2,4));
linPlot=linSlope*BinEdges+linInt;
mleFit.FC.linFit=linFit;
% mleFit.FC.goodness.linGOF=gof;

figure;
hold on;
imagesc(BinEdges,BinEdges,allHist_FC');
caxis([0 0.002]);
h=colorbar('eastoutside');
axis xy;
ylabel('\Delta GFP (%)');
xlabel('\Delta Rhodamine (%)');

title('\Delta Rhodamine vs \Delta GFP');
plot(BinEdges,fitPeaks,'r','LineWidth',2);
plot(BinEdges,fitposStan,'--r','LineWidth',1);
plot(BinEdges,fitnegStan,'--r','LineWidth',1,'HandleVisibility','off');
plot(BinEdges,linPlot,'w','LineWidth',2);
xlim([-2 2]);
ylim([-2 2]);
legend({'MLE Peaks','+/-1Std','linear fit of all data'},'TextColor','white','Location','northeast');
legend('boxoff');
slopeTxt=num2str(linSlope);
intTxt=num2str(linInt);
pvalTxt=num2str(slope_pVal);
rsqrTxt=num2str(linFit.Rsquared.Adjusted);
FC500an=annotation('textbox',[0.13, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt 'x+' intTxt ' rsqr=' rsqrTxt],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);