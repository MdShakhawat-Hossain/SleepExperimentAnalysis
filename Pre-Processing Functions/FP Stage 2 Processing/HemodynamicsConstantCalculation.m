function CorrectSlope = HemodynamicsConstantCalculation(IsosCorrectData,FiberPath,HemiSphere)
% clc; close all; clear all
% Params = FiberData.params;
% IsosCorrectData = control405Correction_Update(FiberData.NE.motionCorrect,Params,'RH');

    saveNameIx = strfind(FiberPath,'.');
    saveName = FiberPath(1:saveNameIx(1)-1);
    FigName = [saveName '_' HemiSphere '_hemodynamicCorrection.fig'];
%     AniamIx = strfind(FiberPath,'_');
%     AnimalID = FiberPath(1:AniamIx(1)-1);
    %%
%     MaxData = max(max(IsosCorrectData.F465),max(IsosCorrectData.F465));
%     MinData = min(min(IsosCorrectData.F465),min(IsosCorrectData.F465));
%     BinEdges=(MinData:0.01:MaxData);
    BinEdges=(-6:0.1:6);
%     figure;DataHist_FC=histogram2(IsosCorrectData.F560,IsosCorrectData.F465,BinEdges,BinEdges,'DisplayStyle','tile','ShowEmptyBins','on');
% 
%     figure;DataHist_FC=histogram2(IsosCorrectData.F465,IsosCorrectData.F560,BinEdges,BinEdges,'DisplayStyle','tile','ShowEmptyBins','on');
%     colorbar
%     xlabel('560');ylabel('465');
 
    [DataHist_FC.BinCounts,~,~]=histcounts2(IsosCorrectData.F560,IsosCorrectData.F465,BinEdges,BinEdges);
    
    BinCounts_FC=DataHist_FC.BinCounts;
    allCounts=sum(BinCounts_FC,'all');
    
    allHist_FC=(BinCounts_FC/allCounts)*100;
    XCounts=sum(BinCounts_FC,2);
    mleFit.FC.allHist_FC=allHist_FC;
    %%
    for binNum=2:size(BinCounts_FC,1)
        if binNum==2
            percData(1)=(XCounts(1)/sum(XCounts))*100;
        end
        percData(binNum)=(sum(XCounts(1:binNum))/sum(XCounts))*100;
    end
    %%
    startInd_FC=find(percData<=2,1,'last')+1; %exclude lower 2% of data points
    endInd_FC=find(percData>=98,1,'first')-1; %exclude upper 2% of data points
    mleFit.FC.startInd=startInd_FC;
    mleFit.FC.endInd=endInd_FC;
    
    startInt=BinEdges(startInd_FC);
    belowThresh=find(IsosCorrectData.F560<startInt);
    endInt=BinEdges(endInd_FC);
    aboveThresh=find(IsosCorrectData.F560>endInt);
    
    blankVals=union(belowThresh,aboveThresh);
    lm_data.F560 = IsosCorrectData.F560;
    lm_data.F465 = IsosCorrectData.F465;
    
    % [B,TF] = rmoutliers(lm_data_F560,'percentiles',[2 98]);
    % removeVals = double(TF);
    % lm_Data_F465=IsosCorrectData.F465(~removeVals);
    
    lm_data.F560(blankVals)=[];
    lm_data.F465(blankVals)=[];
    
    for binNum=startInd_FC:endInd_FC
        index=(binNum-startInd_FC)+1;
        leftEdge=BinEdges(binNum);
        rightEdge=BinEdges(binNum+1);
        
%             dataInds=IsosCorrectData.F560>=leftEdge...
%                 & IsosCorrectData.F560<=rightEdge;

            theData= lm_data.F465; %IsosCorrectData.F465;
        
        
        [phat,pci]=mle(theData); %,'Distrbution','normal'

        [normHist_FC.BinCounts,normHist_FC.BinEdges]=histcounts(theData,BinEdges);
        normCounts_FC=normHist_FC.BinCounts./sum(normHist_FC.BinCounts);
        theFit_FC=pdf('normal',normHist_FC.BinEdges(1:(end-1)),phat(1),phat(2));
        theFit_FC=theFit_FC./sum(theFit_FC);
        
        rsqr=1-(sum((normCounts_FC-theFit_FC).^2)/sum((normCounts_FC-mean(normCounts_FC)).^2));
        
        mleFit.FC.fitParams.phat(index,:)=phat;
        mleFit.FC.fitParams.pci(:,:,index)=pci;
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
    fitPeaks(startInd_FC:endInd_FC)=mleFit.FC.fitParams.phat(:,1);
    fitposStan(startInd_FC:endInd_FC)=mleFit.FC.fitParams.phat(:,1)+mleFit.FC.fitParams.phat(:,2);
    fitnegStan(startInd_FC:endInd_FC)=mleFit.FC.fitParams.phat(:,1)-mleFit.FC.fitParams.phat(:,2);
    
    [linFit]=fitlm(lm_data.F560,lm_data.F465,'RobustOpts','on');
    linSlope=table2array(linFit.Coefficients(2,1));
    linInt=table2array(linFit.Coefficients(1,1));
    slope_pVal=table2array(linFit.Coefficients(2,4));
    linPlot=linSlope*BinEdges+linInt;
    mleFit.FC.linFit=linFit;
    %%
    ImageCheck = figure;
    ImageCheck.WindowState = 'minimized'; 
    
    hold on;

    histogram2(IsosCorrectData.F560,IsosCorrectData.F465,BinEdges,BinEdges,'DisplayStyle','tile','ShowEmptyBins','on');
%     h=colorbar('eastoutside');
%     caxis([MinData MaxData]);
    
    axis xy;
    ylabel('\Delta GFP (%)');
    xlabel('\Delta Rhodamine (%)');
    
    title(' \Delta Rhodamine vs \Delta GFP');

    plot(BinEdges,fitPeaks,'r','LineWidth',2);
    plot(BinEdges,fitposStan,'--r','LineWidth',1);
    plot(BinEdges,fitnegStan,'--r','LineWidth',1,'HandleVisibility','off');
    plot(BinEdges,linPlot,'k','LineWidth',2);
%     xlim([0 1]);
%     ylim([0 1]);
    legend({'MLE Peaks','linear fit of all data','+/-1Std'},'Location','best');
    slopeTxt=num2str(linSlope);
    intTxt=num2str(linInt);
    pvalTxt=num2str(slope_pVal);
    rsqrTxt=num2str(linFit.Rsquared.Adjusted);
    annotation('textbox',[0.13, 0.64, 0.1, 0.1],'String',['Linear fit y=' slopeTxt 'x+' intTxt ' rsqr=' rsqrTxt],'FitBoxToText','on','LineStyle','none','Color',[1 1 1]);
    %%
    savefig(FigName);
    
    close all
    CorrectSlope = linSlope;
end
