% ACh =  AnalysisResults.(animalID).Transitions.AWAKEtoNREM.GRAB_ACh;
% NE =  AnalysisResults.(animalID).Transitions.AWAKEtoNREM.GRAB_NE;
% EMG = AnalysisResults.(animalID).Transitions.AWAKEtoNREM.EMG;
% Force = AnalysisResults.(animalID).Transitions.AWAKEtoNREM.force;
% CBV = AnalysisResults.(animalID).Transitions.AWAKEtoNREM.NE_Rhodamine;
% Whisk = AnalysisResults.(animalID).Transitions.AWAKEtoNREM.whisk  ;
% 
% tableSize = cat(1,ACh,NE,EMG,Force,CBV,Whisk);
% GLM_Table = table('Size',[size(tableSize,2),6],'VariableTypes',{'double','double','double','double','double','double'},'VariableNames',{'ACh','NE','EMG','Force','CBV','Whisk'});
% 
% GLM_Table.ACh = ACh';
% GLM_Table.NE = NE';
% GLM_Table.EMG = EMG';
% GLM_Table.Force = Force';
% GLM_Table.CBV = CBV';
% GLM_Table.Whisk = Whisk';
% 
% modelspec = 'CBV ~ ACh+NE';
% mdl = fitglm(GLM_Table,modelspec,'Distribution','normal')

clear all; clc; close all
load('GLMDataMatrix.mat')
AnimalID = {'NEACh007'};
animalID = AnimalID{1,1};
%% REM
CBV_LH =  mean(GLMDataMatrix.(animalID).GLMModelData.REM.Ach_Rhodamine,2)';
CBV_RH =  mean(GLMDataMatrix.(animalID).GLMModelData.REM.NE_Rhodamine,2)';


NE =  mean(GLMDataMatrix.(animalID).GLMModelData.REM.NE_GFP,2)';
ACh = mean(GLMDataMatrix.(animalID).GLMModelData.REM.Ach_GFP,2)';

figure;plot(CBV_LH,'r'); hold on; plot(ACh,'c');hold on; plot(NE,'g');
tableSize = cat(1,ACh,NE,CBV_LH,CBV_RH);
GLM_Table3_REM = table('Size',[size(tableSize,2),4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'ACh','NE','CBVRH','CBVLH'});

GLM_Table3_REM.ACh = ACh';
GLM_Table3_REM.NE = NE';
GLM_Table3_REM.CBVLH = CBV_LH';
GLM_Table3_REM.CBVRH = CBV_RH';

modelspec = 'CBVLH ~ ACh+NE';
mdl_REM = fitglm(GLM_Table3_REM,modelspec,'Distribution','normal')

%% NREM
CBV_LH =  mean(GLMDataMatrix.(animalID).GLMModelData.NREM.Ach_Rhodamine,2)';
CBV_RH =  mean(GLMDataMatrix.(animalID).GLMModelData.NREM.NE_Rhodamine,2)';


NE =  mean(GLMDataMatrix.(animalID).GLMModelData.NREM.NE_GFP,2)';
ACh = mean(GLMDataMatrix.(animalID).GLMModelData.NREM.Ach_GFP,2)';

figure;plot(CBV_LH,'r'); hold on; plot(ACh,'c');hold on; plot(NE,'g');title('NREM');
tableSize = cat(1,ACh,NE,CBV_LH,CBV_RH);
GLM_Table3_NREM = table('Size',[size(tableSize,2),4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'ACh','NE','CBVRH','CBVLH'});

GLM_Table3_NREM.ACh = ACh';
GLM_Table3_NREM.NE = NE';
GLM_Table3_NREM.CBVLH = CBV_LH';
GLM_Table3_NREM.CBVRH = CBV_RH';

modelspec = 'CBVLH ~ ACh+NE';
mdl_NREM = fitglm(GLM_Table3_NREM,modelspec,'Distribution','normal')

%% Rest
CBV_LH =  mean(GLMDataMatrix.(animalID).GLMModelData.Rest.Ach_Rhodamine,2)';
CBV_RH =  mean(GLMDataMatrix.(animalID).GLMModelData.Rest.NE_Rhodamine,2)';


NE =  mean(GLMDataMatrix.(animalID).GLMModelData.Rest.NE_GFP,2)';
ACh = mean(GLMDataMatrix.(animalID).GLMModelData.Rest.Ach_GFP,2)';

figure;plot(CBV_LH,'r'); hold on; plot(ACh,'c');hold on; plot(NE,'g');
tableSize = cat(1,ACh,NE,CBV_LH,CBV_RH);
GLM_Table3_Rest = table('Size',[size(tableSize,2),4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'ACh','NE','CBVRH','CBVLH'});

GLM_Table3_Rest.ACh = ACh';
GLM_Table3_Rest.NE = NE';
GLM_Table3_Rest.CBVLH = CBV_LH';
GLM_Table3_Rest.CBVRH = CBV_RH';

modelspec = 'CBVLH ~ ACh+NE';
mdl_Rest = fitglm(GLM_Table3_Rest,modelspec,'Distribution','normal')

%% All

CBV_LH =  mean(GLMDataMatrix.(animalID).GLMModelData.All.Ach_Rhodamine,2)';
CBV_RH =  mean(GLMDataMatrix.(animalID).GLMModelData.All.NE_Rhodamine,2)';


NE =  mean(GLMDataMatrix.(animalID).GLMModelData.All.NE_GFP,2)';
ACh = mean(GLMDataMatrix.(animalID).GLMModelData.All.Ach_GFP,2)';

figure;plot(CBV_LH,'r'); hold on; plot(ACh,'c');hold on; plot(NE,'g');
tableSize = cat(1,ACh,NE,CBV_LH,CBV_RH);
GLM_Table3_All = table('Size',[size(tableSize,2),4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'ACh','NE','CBVRH','CBVLH'});

GLM_Table3_All.ACh = ACh';
GLM_Table3_All.NE = NE';
GLM_Table3_All.CBVLH = CBV_LH';
GLM_Table3_All.CBVRH = CBV_RH';

modelspec = 'CBVLH ~ ACh+NE';
mdl_All = fitglm(GLM_Table3_All,modelspec,'Distribution','normal')