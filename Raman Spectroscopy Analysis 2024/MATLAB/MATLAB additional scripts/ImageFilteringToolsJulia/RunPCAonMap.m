function RunPCAonMap

global SpectraArray;
global PCARun
%      SpectraArray.Spec{ii}=Spec; SpectraArray.mapname{ii}=filename{ii}; SpectraArray.pathname{ii}=pathname;
%      SpectraArray.mask{ii}=mask;SpectraArray.Results{ii}=Results;SpectraArray.ReshapeSize{ii}=ReshapeSize;
%      SpectraArray.xj{ii}=XList;
close all force
addpath('C:\Users\ramandemo\Desktop\ImageFilteringToolsJulia\inputdlg')
addpath('E:\Renishaw_C_Data\Martin\Matlab\GlobalWDFTools\MappingTools')
clc

InputforPCA

PCARun.xj=cell2mat(SpectraArray(1).xj(1));
PCARun.dataforPCA=cell2mat(SpectraArray(1).Spec(1));
if Answer.LowWN~=0
xjnew=Answer.LowWN:Answer.HighWN;
PCARun.dataforPCA=interp1(PCARun.xj, PCARun.dataforPCA',xjnew,'linear','extrap')'; 
PCARun.xj=xjnew;
end

PCARun.masktouse=cell2mat(SpectraArray(1).mask(1))';
PCARun.goodspectra=ones(size(PCARun.dataforPCA,1),1);PCARun.goodspectra(PCARun.masktouse)=0;%PCARun.goodspectra=PCARun.goodspectra';
PCARun.ReshapeSize=cell2mat(SpectraArray(1).ReshapeSize(1));

if strmatch(Answer.BaseLineType,'Polynomial')
disp(['PolyNomial Order to use: ',num2str(Answer.PolynomialOrder)])

looksgood=0;
CheckBaselinedSpectra
if looksgood==1
baselinesubtractfunc   
end

disp('Polynomial')
elseif strmatch(Answer.BaseLineType,'EMSC')
disp(['PolyNomial Order to use: ',num2str(Answer.PolynomialOrder)])    
performEMSCfunc    
disp('EMSC')
end
% figure;plot(PCARun.xj,PCARun.dataforPCA(PCARun.goodspectra==1,:))

if Answer.NormaliseMode==1
PCARun.dataforPCA(PCARun.goodspectra==1,:)=bsxfun(@rdivide,PCARun.dataforPCA(PCARun.goodspectra==1,:),sqrt(sum(PCARun.dataforPCA(PCARun.goodspectra==1,:).^2,2)));
disp('Normalise data')
end

if Answer.MeanCentreMode==1
PCARun.dataforPCA(PCARun.goodspectra==1,:)=bsxfun(@minus,PCARun.dataforPCA(PCARun.goodspectra==1,:),mean(PCARun.dataforPCA(PCARun.goodspectra==1,:),1));
disp('Mean Centre data')
end

PerformPCA_waitbar = waitbar(0.05,'Initializing and performing PCA....Please wait.....');

disp(['Number of PCs: ',num2str(Answer.PCs)])
Answer.PCs;

if strmatch(Answer.PCAType,'NIPALS')
disp('Running NIPALS')
[tempT,PCARun.PCA.P,PCARun.PCA.e]=PCA(PCARun.dataforPCA(PCARun.goodspectra==1,:),Answer.PCs); %CHECK
PCARun.PCA.T=zeros(size(PCARun.dataforPCA,1),Answer.PCs);
n=0;
for jj=1:size(PCARun.dataforPCA,1)
    if PCARun.goodspectra(jj)==1
        n=n+1;
        PCARun.PCA.T(jj,:)=tempT(n,:);               %CHECK
    end
end
PCARun.PCA.T=PCARun.dataforPCA*PCARun.PCA.P';               %CHECK
elseif strmatch(Answer.PCAType,'SVD')
disp('Running SVD')
[U,S,V]=svd(PCARun.dataforPCA(PCARun.goodspectra==1,:)'*PCARun.dataforPCA(PCARun.goodspectra==1,:),'econ');
PCARun.PCA.P=V(:,1:Answer.PCs)';
PCARun.PCA.T=zeros(size(PCARun.dataforPCA,1),Answer.PCs);
PCARun.PCA.T=PCARun.dataforPCA*V(:,1:Answer.PCs);
T2=sum(PCARun.PCA.T.*PCARun.PCA.T);
X2=sum(sum(PCARun.dataforPCA.*PCARun.dataforPCA));
PCARun.PCA.e=T2./X2;
end

waitbar(0.5,PerformPCA_waitbar,'PCA complete. Generating PCA score maps.....Please wait.');

PCARun.PCAScores=reshape(PCARun.PCA.T,PCARun.ReshapeSize(1),PCARun.ReshapeSize(2),Answer.PCs);
PCAfig=figure('name','PCA Analysis');
PCARun.PCAfig=PCAfig;
subplot(2,1,1)
imagesc(PCARun.PCAScores(:,:,1));
title(['PC Score: ',num2str(1)])
axis image
colorbar
colormap winter

PCARun.PCALoadings=PCARun.PCA.P;%
subplot(2,1,2)
plot(PCARun.xj,PCARun.PCALoadings(1,:),'tag','PCloading');
hold on
plot(PCARun.xj,zeros(size(PCARun.xj,2)),'k')
hold off
title('PCloading')

PCAnoselectiontext = uicontrol('Parent',PCAfig,...
              'Style','text',...
              'String','PC Score: ',...
              'units','normalized',...
              'Position',[0.01 0.899 0.08 0.1]);
PCAnoselection = uicontrol('Parent',PCAfig,...
              'Tag', 'PCAnoselection',...
              'Style','popupmenu',...
              'String','PCA',...
              'units','normalized',...
              'Position',[0.1 0.9 0.08 0.1]);
          peaklabelpush = uicontrol('Parent',PCAfig,...
              'Tag', 'peaklabelpush',...
              'Style','pushbutton',...
              'String','Peak Label Loading',...
              'units','normalized',...
              'callback','PeakLabelLoading',...
              'Position',[0.8 0.95 0.1 0.02]);
abc=1:Answer.PCs;abc=abc';
set(findobj(PCAfig,'Tag','PCAnoselection'),'String', abc);
set(findobj(PCAfig,'Tag','PCAnoselection'),'callback','changePCAplot2018');
% set(findobj(PCAfig,'Tag','PCAnoselection'),'Value', abc);

close(PerformPCA_waitbar)
end

function baselinesubtractfunc

global PCARun


%RaMP.UnProcessedMapdata_baselinesub=zeros(size(RaMP.UnProcessedMapdata));
Baseline=zeros(size(PCARun.dataforPCA));

pp = waitbar(0,'Please wait.....');
for i=1:size(PCARun.dataforPCA,1)
clc
Statusmessage=(['Performing baseline subtraction... ',num2str(round((i/size(PCARun.dataforPCA,1)*100))),'%']);drawnow;disp(Statusmessage);
% OPT: 'Poly order' VAL: 1 TYPE: double
% OPT: 'Force +ve' VAL: 0 TYPE: double
% OPT: 'Iterations' VAL: 100 TYPE: double
% OPT: 'Binning' VAL: 1 TYPE: double
% OPT: 'Known' VAL: 'none' TYPE: getfile
if PCARun.goodspectra(i)==1
%     [PCARun.dataforPCA(i,:),Baseline(i,:)]= removeBaseline(PCARun.xj,PCARun.dataforPCA(i,:),PCARun.PolynomialOrder,0,50,1,'none');
%       waitbar(i/size(PCARun.dataforPCA,1),pp);
    [PCARun.dataforPCA(i,:),Baseline(i,:)]=f_baseline_corr(PCARun.xj,...
        PCARun.dataforPCA(i,:),5,600,5);
end
end
close(pp);
end

function performEMSCfunc
global PCARun
EMSC_waitbar = waitbar(0,'Please wait..Performing EMSC correction...');
waitbar(0.5,EMSC_waitbar,'Please wait..Performing EMSC correction...');
M=mean(PCARun.dataforPCA(PCARun.goodspectra==1,:));
[~,~,~,~,~,PCARun.dataforPCA(PCARun.goodspectra==1,:),~,~,~,~]=EMSC_martin(PCARun.dataforPCA(PCARun.goodspectra==1,:),PCARun.xj,M,[],PCARun.PolynomialOrder);  
waitbar(1,EMSC_waitbar,'EMSC correction complete.');close(EMSC_waitbar);
end

