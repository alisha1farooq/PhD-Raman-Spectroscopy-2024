function [Spectra,mapoutput,mask] = PreprocessSpectraFromFiles(XList,MedianFilterSize,S2NLimit,SDfilter,savetoRaMP);
global rememberpath;
%close all; clc;
Spectra = []; 
strOut = []; 
cancel = false; 

if MedianFilterSize < 0 | MedianFilterSize > 2 
    errordlg('Error: MedianFilterSize parameter must be 0, 1 or 2','Error')
    return
end

%Get filenames to process
if length(rememberpath)> 0
[filename,pathname] = uigetfile([rememberpath,'*.wdf'],'Select WDF files to preprocess',...
                                'MultiSelect','on');
                            rememberpath=pathname;
else 
    
    [filename,pathnamne] = uigetfile('*.wdf','Select WDF files to preprocess',...
                                'MultiSelect','on');
   rememberpath=pathname;
end
%Allow user to cancel wdf file loading operation, or convert a single
%string into a cell
if ~iscell(filename)    %ie only one file selected
    if filename == 0;   %Cancel selected
        return
    else 
        filename = cellstr(filename); %Turn strings into cells
    end
end

% try

h = waitbar(0,['Loading and preprocessing spectra dataset 1/' num2str(length(filename))],...
            'Name','Performing Preprocessing',...
            'CreateCancelBtn','setappdata(gcbf,''canceling'',1)','pointer','watch');
setappdata(h,'canceling',0)
set(h,'CloseRequestFcn','');

%Go through the wdfFiles. Load them in and perform required preprocessing 
%operations on them.
global SpectraArray;
SpectraArray=cell(0);
for ii = 1:size(filename,2)
    fullFilename = [pathname filename{ii}];

   clc;disp(['1. Loading from WDF. File ' num2str(ii) '/' num2str(size(filename,2))]); 
   waitbar(1/7,h,['Loading from WDF. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar

    wdf = WdfReader(fullFilename);  %Load wdf
    xList = wdf.GetXList;

%     a=findstr(Filename,'\');a=max(a);filenametxt=Filename(a+1:end);
    global plotfigures
    if plotfigures==1
    clc;disp(['2. Plotting Raw data' num2str(ii) '/' num2str(size(filename,2))]);
    waitbar(2/7,h,['Plotting Raw data' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar    
        
        tempSpectra = wdf.GetSpectra(1,wdf.Count);

        figure; plot (xList,tempSpectra);title([filename{ii},' (Uncorrected)']);
        clear tempSpectra xList
    end
    
    waitbar(3/7,h,['Generating Raw Chemical map. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar
    clc;disp(['3. Generating Raw Chemical map. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar
    global chemicalwavenumber
    if chemicalwavenumber==0
        [ChemicalMap2D,ChemicalMap3D] = GetInteAreaMap(wdf);
        titletype='Raw Total Intensity Map';
        % TotalIntMap1D = trapz(RawSpectra,2);%sum(z2);%cumtrapz(z,2);
        % TotalIntMap2D = reshape(TotalIntMap1D,ReshapeSize(1),ReshapeSize(2),1);
        % figure;imagesc(TotalIntMap2D);title('Raw Total Intensity Map')
    else
        [~,ChemicalMap2D,ChemicalMap3D] = GetChemicalMap(wdf,chemicalwavenumber);
        titletype=['Raw Chemical Map at ',num2str(chemicalwavenumber)];
    end

    wdf.Close();    
    imagetoshow=ChemicalMap3D;
    figure;imagesc(imagetoshow);title(titletype)
        try
            qlo=quantile(ChemicalMap2D,0.05);
            qhi=quantile(ChemicalMap2D,0.95);    
        catch
            X=ChemicalMap2D;    
            P = 0.05;       %# Your probability
            S = sort(X);    %# Sort the columns of your data X
            N = size(X,1);  %# The number of rows of X
            Y = interp1q([0 (0.5:(N-0.5))./N 1]',S([1 1:N N],:),P);  %'# Get the quantiles
            qlo=Y;
            P = 0.95;       %# Your probability
            S = sort(X);    %# Sort the columns of your data X
            N = size(X,1);  %# The number of rows of X
            Y = interp1q([0 (0.5:(N-0.5))./N 1]',S([1 1:N N],:),P);  %'# Get the quantiles
            qhi=Y;
        end  
        set(gca,'CLim',[qlo,qhi]);

    
    %Perform initial preprocessing - spatial filtering,
    %removal of saturated data and O2 peak, interpolation. 
    clc;disp(['4. Filtering map file. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar
    waitbar(4/7,h,['Filtering map file. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar

    [Spec,mask,Results] = GetPreprocessedSpectra(fullFilename,XList,MedianFilterSize,S2NLimit,SDfilter);  

    %Check ptMaps
    [Spec,Map,ReshapeSize,mask]=CheckFixPtMap(Spec,fullFilename,mask);

    if chemicalwavenumber==0
        ChemicalMap2D = trapz(Spec,2);
        titletype='Filtered Total Intensity Map';
    else
        a=find(XList>chemicalwavenumber-2&XList<chemicalwavenumber+2);a=a(1);
        ChemicalMap2D=Spec(:,a);
        titletype=['Filtered Chemical Map at ',num2str(chemicalwavenumber)];
    end
    ChemicalMap2D(mask',:)=0;
    imagetoshow = reshape(ChemicalMap2D,ReshapeSize(1),ReshapeSize(2),1);
    
    disp(['5. Generating Filtered Chemical map. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar
    waitbar(5/7,h,['Generating Filtered Chemical map. File ' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar

    figure;imagesc(imagetoshow);title({titletype,Results});
        try
            qlo=quantile(ChemicalMap2D,0.05);
            qhi=quantile(ChemicalMap2D,0.95);    
        catch
            X=ChemicalMap2D;    
            P = 0.05;       %# Your probability
            S = sort(X);    %# Sort the columns of your data X
            N = size(X,1);  %# The number of rows of X
            Y = interp1q([0 (0.5:(N-0.5))./N 1]',S([1 1:N N],:),P);  %'# Get the quantiles
            qlo=Y;
            P = 0.95;       %# Your probability
            S = sort(X);    %# Sort the columns of your data X
            N = size(X,1);  %# The number of rows of X
            Y = interp1q([0 (0.5:(N-0.5))./N 1]',S([1 1:N N],:),P);  %'# Get the quantiles
            qhi=Y;
        end  
        set(gca,'CLim',[qlo,qhi]);

    
    if plotfigures==1
    disp(['6. Plotting Cleaned data' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar        
    waitbar(6/7,h,['Plotting Cleaned data' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar    

    FilSpectra=Spec;
    FilSpectra(mask',:)=[];
    figure; plot (XList,FilSpectra);
    %title([filename{ii},Results]);
    title({filename{ii},Results})
    end
    
    %Add to Spectra array.
    Spectra = [Spectra;Spec];
     SpectraArray.Spec{ii}=Spec; SpectraArray.mapname{ii}=filename{ii}; SpectraArray.pathname{ii}=pathname;
     SpectraArray.mask{ii}=mask;SpectraArray.Results{ii}=Results;SpectraArray.ReshapeSize{ii}=ReshapeSize;
     SpectraArray.xj{ii}=XList;
     mapoutput=SpectraArray;
    if getappdata(h,'cancelling')
        cancel = true;
        break
    end
    
        
    if savetoRaMP==1
    disp(['7. Saving Cleaned data' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar    
    waitbar(7/7,h,['Saving Cleaned data' num2str(ii) '/' num2str(size(filename,2))]);       %Update the waitbar    
    %fullFilename = [pathname filename{ii}];
    RaMP = WriteDatatoRaMP(XList,Map,ReshapeSize,mask);
    newfilename=strrep([filename{ii}],'%','');
    newfilename=strrep(newfilename,' ','_');
    newfilename=strrep(newfilename,'-','_');
    newfilename=strrep(newfilename,'.wdf','.mat');
%     cd(pathname)
    save([pathname,newfilename],'RaMP');
    cd(pathname)
    end
    
end

%==========================================================================
% Error handling and clean-up
% catch ME
%     delete(h)
%     errordlg(['Error occurred: preprocessing spectra from files failed.', char(10) ME.message])
% end

delete(h)  %Remove waitbar

%If preprocessing was cancelled, return nothing and keep the handles
%structure unchanged.
if cancel
    Spectra = [];
    return
end

end


