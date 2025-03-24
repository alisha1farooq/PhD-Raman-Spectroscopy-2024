function [Spectra,mask,Results] = GetPreprocessedSpectra(Filename,MasterXList,MedianFilterSize,S2NLimit,SDfilter)
%Preprocess spectra from file specified by Filename.  Returns preprocessed 
% Spectra.
% 1 Replaces any saturated spectra with neighbours (stops propagation 
%   through median filtering for heavily saturated maps).
% 2 Performs median filtering of radius MedianFilterSize. 
% 3 Removes the edge spectra to a depth of MedianFilterSize and any
%   saturated or below S2NLimit.
% 4 Interpolates sectra onto MasterXist.
    %Open wdf file whose filename was supplied and get data and various parameters from file
    disp('4.1 Filtering map file......Reading Raw Spectra')
    wdf = WdfReader(Filename);  %Load wdf
    xList = wdf.GetXList;
    Spectra = wdf.GetSpectra(1,wdf.Count);
    [Xi,Yi] = wdf.GetOriginCoords(1,wdf.Count);
    rows = sum(Yi==min(Yi));              %Number of rows
    cols = sum(Xi==min(Xi));              %Number of Columns
    saturated = wdf.GetOriginFlags(1,wdf.Count);
    clear Xi Yi;
    
Found = [];    
if SDfilter > 0
disp('4.2 Filtering map file......STD filtering')
% Size=2;
Found = zeros(wdf.Count(),1);
h=waitbar(0,'Flagging cosmic rays...');
TotalReadInSoFar = 1;
wdf.StartChunkwiseReading();
while (wdf.AreMoreChunks())
    z = wdf.GetNextDataChunk();
    Diffy = diff(z');
    MeanDiffy = mean(Diffy);
    stdDiffy = std(Diffy);
    NumberReadIn = size(z,1);
    for i = 1:NumberReadIn
        if(any(Diffy(:,i) > SDfilter*stdDiffy(i) + MeanDiffy(i)))
            Found(i+TotalReadInSoFar-1) = true;
        end
    end
    TotalReadInSoFar=TotalReadInSoFar+NumberReadIn;
    waitbar(TotalReadInSoFar/wdf.Count());
end
close(h);
Found=find(Found==1);
end
%temp=Found(Found==1); 
%Found=temp;
%else 
 %   Found=[];
%end
   mu = MapUtils(rows,cols);

    %Perform median filtering. This will only work if all data is in from 
    %2D Raman maps.  Code requires the image processing toolbox.
    if MedianFilterSize>0
disp('4.3 Filtering map file......Median filtering')
        
          %Replace any saturated spectra with a non-saturated neighbour.  Ensures 
    %saturated artefacts do not propage through the median filtering step.
            if not(isempty(saturated))
                Spectra = mu.ReplaceSpectraWithNeighbours(Spectra,saturated);
            end;    
        
    Spectra = GN_SpatialFilters(Spectra,rows,cols,false,true,0,MedianFilterSize);
    
    %Get the list of spectra coresponding to the edge 0width to be removed.
    %N.B. list will be empty if EdgeWidth = 0
    edgeMask = mu.EdgeMask(MedianFilterSize);
    else
        edgeMask =[];
    end
    
    %Get the list of below signal-to-noise threshold spectra.  Requires the
    %signal processing toolbox rms function
disp('4.4 Filtering map file......SNR filtering')
    s2NConfig.Active = true;
    s2NConfig.NoiseRegion1 = 1653; s2NConfig.NoiseRegion2 = 1692;
    s2NConfig.SignalRegion1 = 1574; s2NConfig.SignalRegion2 =1643;
    s2NConfig.Threshold = S2NLimit;
    belowS2N = mu.UnderS2nThreshold(xList,Spectra,s2NConfig);
    clear mu;

    %Build combined mask and remove from the training set. This can't be 
    %done before spatial filtering.
   % mask = unique([edgeMask;saturated;belowS2N]);
    mask = unique([edgeMask;saturated;belowS2N;Found]);

    outStr = sprintf(['Spectra from file %s: included = %d, edge = %d, saturated = %d, below S2N(<',num2str(S2NLimit),')= %d, foundSTD(>',num2str(SDfilter),')= %d'],...
        Filename, size(Spectra,1), size(edgeMask,1), size(saturated,1), size(belowS2N,1),size(Found,1));
    disp(outStr)
    Results=sprintf(['included = %d, edge = %d, saturated = %d, below S2N(<',num2str(S2NLimit),')= %d, foundSTD(>',num2str(SDfilter),')= %d'],...
        size(Spectra,1), size(edgeMask,1), size(saturated,1), size(belowS2N,1),size(Found,1));
    
    if size(mask,1)==size(Spectra,1)
        warndlg(sprintf('Warning: data filtering has excluded all data from file %s', Filename),'Warning');
    elseif size(mask,1)==0
        warndlg(sprintf('Warning: No Spectra removed from file %s', Filename),'Warning');
    else
        
    Spectra =  (interp1(xList,Spectra',MasterXList,'spline'))';    %Cubic spline interpolation
    disp('4.5 Replacing flagged Spectra with zero')
    for j=mask'
    Spectra(j,:) = zeros(size(MasterXList));
    end
    disp('4.6 Map filtering complete')
%      if isempty(Spectra) 
%         warndlg(sprintf('Warning: data filtering has excluded all data from file %s', Filename),'Warning');
%     else
%         Spectra =  (interp1(xList,Spectra',MasterXList,'spline'))';    %Cubic spline interpolation
%     end
    end
    %Spectra(mask,:) = [];
%     outStr = sprintf('Spectra from file %s: included = %d, edge = %d, saturated = %d, below S2N = %d',...
%         Filename, size(Spectra,1), size(edgeMask,1), size(saturated,1), size(belowS2N,1));
%     disp(outStr)
   
end

