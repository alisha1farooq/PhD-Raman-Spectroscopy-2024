close all
close all hidden
close all force

clc
addpath('E:\OneDrive - University of Leeds\University of Leeds\Matlab scripts\Renishaw Mfiles for loading'); %Folder where your WDFreader scripts and requried subfunctions are located (usually something like: C:\Users\mi141418\Documents)
%addpath('E:\Renishaw_C_Data\Martin\Collaborations\James Livermore\Classification Tool\DataClassificationTool_160218');
Defaultscriptpath='E:\OneDrive - University of Leeds\University of Leeds\Data\Raman\20180326_Candelas_TomatoCallose';
addpath(Defaultscriptpath)
global rememberpath
rememberpath=Defaultscriptpath;
if isempty(findstr(rememberpath(end),'\'))
rememberpath=[rememberpath,'\'];
end
    
ButtonName = questdlg('Retrieve XList from file or input manually?', ...
                         'XList Question', ...
                         'From File', 'Input Manually','From File');
   switch ButtonName,
     case 'From File',
[filename,pathfolder] = uigetfile([rememberpath,'*.wdf'],'Select WDF file for XList');
wdf = WdfReader([pathfolder,'\',filename]);  %Load wdf
rememberpath=pathfolder;
xList = wdf.GetXList;
wdf.Close();
% ReshapeSize  = wdf.GetWMapblock().numPoints;
% ReshapeSize(3) = [];
% ReshapeSize  = ReshapeSize';
       case 'Input Manually',
   prompt={'Input low wavenumber:','Input high wavenumber:'};
   name='Input for low and high wavenumber';
   numlines=1;
   defaultanswer={'450','1800'};
   options.Resize='on';
   options.WindowStyle='normal';
   options.Interpreter='tex';
   answer=inputdlg(prompt,name,numlines,defaultanswer,options);   
   xList = str2num(answer{1}):str2num(answer{2});           
   end % switch

% PreprocessSpectraFromFiles(XList,MedianFilterSize,S2NLimit,SDfilter,savetoRaMP);
prompt={'Median filter Size (0=no filter):','S2N Limit (high number = strong filter;0=no filter):','SD Filter threshold Size(low number = strong filter;0=no filter):','Wavenumber for Chemical maps (0=Total intensity)','Save to RaMP (1=yes;2=no):','Plot Figures (1=yes;2=no):'};
   name='Input for Peaks function';
   numlines=1;
   defaultanswer={'0','4','5','1002','0','0'};
   answer=inputdlg(prompt,name,numlines,defaultanswer);
   
global chemicalwavenumber
chemicalwavenumber=str2num(answer{4});
global plotfigures
plotfigures=str2num(answer{6});   
   
[~,mapoutput,~] = PreprocessSpectraFromFiles(xList,str2num(answer{1}),str2num(answer{2}),str2num(answer{3}),str2num(answer{5}));
disp('Filtering and saving complete')
% cd(Defaultscriptpath)


%%
   ButtonName = questdlg('Do you want to run PCA on map?', ...
                         'PCA Question', ...
                         'Yes', 'No', 'No');
   switch ButtonName,
     case 'Yes',
      RunPCAonMap
%       disp('Your favorite color is Yes');
     case 'No',
       disp('PCA on map not selected')
   end % switch