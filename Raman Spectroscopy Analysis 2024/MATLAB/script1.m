%% Preparation of the datasets to load onto the other script

% Specify locations of .wdf data files, Matlab scripts, and output path. List all of the folders that have data to be processed. Just to for this
% script; can add new ones following the same format. Where there is an
% edit required I have added an X.

% directory where the data is 
dirPath={...
    'X';...  
};

% Matlab scripts required for this analysis. Should be in a folder labelled
% MATLAB additional scripts
addpath( "X" );

% The path used to save the output from this script. 
savePath = 'X';



% Specify each of the data sets.

% For each of repeats, need to provide ('i' refers to the data set;
% starting from 1, and ad one for each new data set):
% - start and end of all filenames as Name{i,1} and Name{i,2}, such that the
%   actual filename is: "<Name{i,1}><repeat number><Name{i,2}>.wdf
% - the range of indices for the repeats as Index{i}=first:last;
% - start and end of all filenames for the background data, as BName{i,1}
%   and BName{i,2}, used similar to Name{} above.
% - the range of start and end of the background data runs, sim. to above.
% - the prefix for the output filename as "NameSaving{i}".
% - a "calibration peak" as CalibrationPeak{i}".

%%save each repeat with the same beginning (e.g. cpa / cpb), 
%%add an underscore _  and then each repeat/new spectrum is just given a 
%%number in ascending order. All spectrum names must start the same 
%%so just the number at the end changes. 
 
%Set 1:
Name{1,1}='X';
Name{1,2}='';      % No text after the repeat index for this set of runs.
Index{1}=[X:X];       % i.e. repeat indices 0, 1, 2, 3, 4 and 5.

selectedDataSets = [1];       % Starts from 1. Can do ranges as '=[start,end];'

%% Load each of the selected data sets
for j=selectedDataSets
    
    %
    % Load all repeats for this data set.
    %
    fprintf( "Loading data from data set '%s' (index %i)\n", dirPath{j}, j );

    % Move to the directory for this dataset, and keep a record of where we
    % were.
    oldPath = cd( dirPath{j} );

    % At the end of the loops below, m will be the number of data sets,
    % and b will be the number of background repeats. They don't seem to
    % be used anywhere, though - maybe this was orignally meant as a check?
    m=1;
    b=1;

    % Loop over all repeats for this data set.
    for i=Index{j}

        % Get the filename for a single repeat.
        FileName{j,m}=strcat(Name{j,1},num2str(i),Name{j,2},'.wdf');
        fprintf( "- loading repeat: %s\n", FileName{j,m} ); % Check.

        % Import the ".wdf" file for this repeat.
        wdf=WdfReader(FileName{j,m});
        RamanSpectra{m}(:,1)=wdf.GetXList.';
        RamanSpectra{m}(:,2)=wdf.GetSpectra(1,1).';
        wdf.Close()

        m=m+1;      % Count the number of data repeats.
    end
    
    % Get the calibration peak for this dataset.
    CalibrationPeakTemp = CalibrationPeak{j};
    
    % Save the data. Move back to the original path first of all, i.e.
    % prior to cd'ing to the data directory.
    cd( oldPath );
    oldPath = cd( savePath );
    save([NameSaving{j} 'datasetraw'],...
       'RamanSpectra','CalibrationPeakTemp') 

    % Move back to the original directory before finishing.
    cd( oldPath );
%% X alignment of the spectra respect to the Calibration peak position in Si

% Every day of measuring a calibrating spectra is taken and it is aligned
% as good as it can to the 520.5 cm^-1 peak of Silicon. However, sometimes
% the peak seems to be around 520.3 or 520.8. This peak will be used to
% align the shift of the spectra and reduce the variability between days.

% As the PCA only considers the Y-values for analysis,
% a spline function should be fit and the Y-values should be guaranteed to
% have exactly the same X.
fprintf('Aligning spectra and fitting to the selected X values using a spline \n')

tic   % Start timing from here.

Xint_st = round(RamanSpectra{1}(1,1))
Xintvalues = Xint_st:-0.5:round(RamanSpectra{1}(end,1))-4;
% Xintvalues = Xint_st + Xint_en;
% Xintvalues = round(RamanSpectra{1}(1,1)) + 4:0.5:round(RamanSpectra{1}(end,1))-4;
    % Array of X values for the interpolation. RamanSpectra{1}(1,1) is the
    % start point, and 5:0.5:round(...)-4 gives values separated by 0.5
    % (the rounding ensures all points are an integer or half-integer,
    % and evenly spaced),

% Spline interpolation of the actual data.
for k=1:length(RamanSpectra)
   RamanSpectraInt{k}=spline(RamanSpectra{k}(:,1)+(520.5-CalibrationPeak{j}),RamanSpectra{k}(:,2),Xintvalues);
end

toc   % 'tic', 'toc' used for timing how long operations take; toc displays the result.

%% Subtract the baseline
fprintf('Subtracting the baseline \n')

tic   % Start timing from here.

RamanSpectraBaseline = cell( size(RamanSpectraInt) );

% Two sets of parameters here for the two baseline corrections; see below.
smoothwidth1=40;
smoothwidth2=60;

frame1=500;
frame2=800;

iterations1=15;
iterations2=3;

for i=1:size(RamanSpectra,2)   
    disp(['Correcting spectrum ' num2str(i) ' of ' num2str(length(RamanSpectraBaseline))])

    % 1. Remove the maximum possiible of the baseline; parameters ending '1'
    [RamanSpectraB{i}, varargout{i}] = f_baseline_corr(Xintvalues,...
        RamanSpectraInt{i},smoothwidth1,frame1,iterations1);
        % By adding a final argument 'plot2D' as suggested by the comment,
        % will generate plots showing the baseline corrections.

    % 2. Help remove the zero baseline ; parameters ending '2'.
    [RamanSpectraB2{i}, varargout{i}] = f_baseline_corr(Xintvalues,...
        RamanSpectraB{i},smoothwidth2,frame2,iterations2);

    RamanSpectraBaseline{i}=RamanSpectraB2{i}.';
end

toc   % Display time from the previous 'tic' command.


%% Smoothing of the spectra
fprintf('Smoothing the spectra \n')

tic

%%can change the iterations and polynomial here if need be (from 2,17)
for i=1:size(RamanSpectraBaseline,2)
    RamanSpectraSmooth{i}=sgolayfilt(RamanSpectraBaseline{i},2,17);
end

toc

%% Truncating the spectra
%Could truncate down to fingerprint region here. May be easier than zapping
%in the next section (zapping gets rid of section from 1800 - 2800, keeping
%the high freq region - not sure if this is necessary, i find it hassle.)
fprintf('Truncating the spectra \n')

tic

lowerlim  = 500;  % units of cm^{-1}
higherlim = 3200;

% Get the array indices corresponding to the above values.
%May need to switch these round depending on if xdata is ascending or
%descending
higherlimIndex  = find(floor(Xintvalues)==lowerlim,1);
lowerlimIndex = find(floor(Xintvalues)==higherlim,1);

% Remove values outside of the above index range.
XvaluesTrunc = [Xintvalues(lowerlimIndex:higherlimIndex)];

% Same for all of the smoothed spectra.
for i=1:length(RamanSpectraSmooth)
    RamanSpectraTrunc{i}=RamanSpectraSmooth{i}(lowerlimIndex:higherlimIndex);
end

toc


%% Normalization of the spectra
fprintf('Normalizing the spectra \n')

tic

% Intensity of the peaks may vary due to different setups, so only relative
% intensity between peaks can be considered real data.
for i=1:length(RamanSpectraTrunc)
        
      % For normalization to the amide peak
     x2=find(round(XvaluesTrunc)==1600); 
      x1=find(round(XvaluesTrunc)==1700);
      
 
      
      % If we assume that the baseline subtraction is good, then we can
      % just normalize by the maximum
      RamanSpectraNorm{i}=RamanSpectraTrunc{i}/... %was RamanSpectraZap but trying out no zap
          max(sgolayfilt(RamanSpectraTrunc{i}(x1(1):x2(1)),2,31));       
end

%%%toc

%% Data Matrix
% Get the colums of intensities to analyse in a matrix:
fprintf('Creating the Data Matrix \n')

tic

DataMatrix=[];
for i=1:size(RamanSpectraNorm,2)
    DataMatrix(:,end+1)=RamanSpectraNorm{i};
    
end

toc

%% Calculate the Average and STD
Average = mean(DataMatrix,2);
Error   = std (DataMatrix,0,2);


%% Plot figures
fprintf('Plotting Figure \n')
figure = figure
subplot(2,1,1)

plot(XvaluesTrunc, Average)
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity ')

%hold on
subplot(2,1,2)
plot(XvaluesTrunc,Average)
xlim([500 1800])

xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity ')
xticks(500:100:1800)
%hold on


%% Saving Data
    disp('saving data')
    cd(savePath)
    filename = [NameSaving{j} '.tif']
    save([NameSaving{j} 'datasetanalysed'],...
    'DataMatrix','Average','Error','XvaluesTrunc','RamanSpectra','CalibrationPeakTemp')
    saveas(figure, filename)
%waitforbuttonpress
    
% Clear data for next dataset
    clear RamanSpectraB RamanBackgroundSpectra ...
        RamanBSpectraInt RamanSpectraBaseline RamanSpectraNorm ...
        RamanSpectraQuartz RamanSpectraSmooth RamanSpectraTrunc ...
        RamanBSpectraInt RamanSpectraInt RamanSpectraZero...
        BackgroundAverage BackgroundAverageZero...
        DataMatrix RamanSpectraZap

    %%
end

