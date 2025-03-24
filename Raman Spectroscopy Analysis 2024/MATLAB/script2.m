%% Preparation of the datasets to load into data analysis script

%% WDF FILE IMPORTING
clear all
close all
clc

% All lines to be edited have been labelled with an X

% add path for emsc folder & raw data
dirPath={'X'...
               };     
           
% Add MATLAB additional files          
addpath('X')    
% EMSC function path
addpath('X')
% For the background, add the path where the background data is (reference
% dataset and analyseddataset from previous matlab script is also in this folder!!):
addpath('X')


savingPath='X';

% The file names are always an identifier + a number: (I'm assuming this
% means namesaving)

%raw data path for testing EMSC function
cd('X');
wdf_files = dir('*.wdf');

Name{1,1}='X';
Name{1,2}='';      % No text after the repeat index for this set of runs.
Index{1}=[X:X];
CalibrationPeak{1}=520.5;
NameSaving{1}='X'; %the analysed file you want to run has to start with this name also and keep it in the same folder as raw data


selectedDataset=1;  


%% Optimizing the parameters
% This section gives a summary plot of what the spectra will look like
% after the different stages of dataset creation.
% This should help you choose the right parameters before you run the whole
% thing 

datasettesting=1; % Choose which dataset you are testing (if only 1 then change to 1)
cd(dirPath{datasettesting});
indextesting=1; % Choose one spectrum to test on. 
% Creating the structure file with the map inside
MapTesting=WdfReader([Name{datasettesting} '.wdf']);
% Obtaining the spectra from the structure file, where each row is a
% spectrum
RamanSpectratesting=MapTesting.GetSpectra(1,MapTesting.Count);
% Changing it to columns because that is how I usually work with them
RamanSpectratesting=RamanSpectratesting';
Wavenumbers=MapTesting.GetXList'; %added ' at the end to flip values to 1948x1
CalibrationPeaktesting=CalibrationPeak{datasettesting};
% Choose the final Xintvalues that you want to use for the spline
% interpolation.
% Make sure you choose values withing the range you have acquired, or the
% spline will go very weird. Typically go for 750-1800 for FP region
Xintvalues=[500:0.5:3200]; 


RamanSpectraInttesting=spline(Wavenumbers+(520.5-CalibrationPeak{datasettesting}),...
    RamanSpectratesting(:,indextesting),Xintvalues);



% Import the background and reference data for the correction with EMSC
% For background (b)
load('X_datasetraw.mat')
% Variables called "Raman Spectra","Wavenumbers" and "CalibrationPeakTemp"
for k=1:size(RamanSpectra,2)
   b(:,k)=spline(Wavenumbers +(520.5-CalibrationPeakTemp),RamanSpectra{k}(:,2),Xintvalues);
end

% For biological reference (r)
load('X_datasetraw')
% Variables called "Ramtest3=test2'an Spectra","Wavenumbers" and "CalibrationPeakTemp"
for k=1:size(RamanSpectra,2)
   r(:,k)=spline(Wavenumbers+(520.5-CalibrationPeakTemp),RamanSpectra{k}(:,2),Xintvalues);
end

% Choose the parameters for the EMSC correction

N=1;
r_EMSC=r.';
b_EMSC=b.';
[RamanSpectraEMSCtesting,background,c_r,c_b,B_N] = EMSC(RamanSpectraInttesting,r_EMSC,b_EMSC,N);

% Choose the parameters for the baseline subtraction
addpath('X')
smoothwidth=10;
bwidth=400;
iterations=10;
[RamanSpectraBaselinetesting, baseline]=f_baseline_corr(Xintvalues,...
        RamanSpectraEMSCtesting,smoothwidth,bwidth,iterations);  

% Plot the whole thing so you get to change the parameters :)
figure

subplot(2,2,1)
title('Raw data')
hold on
plot(Wavenumbers,RamanSpectratesting(:,indextesting))
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity / counts')
axis tight

subplot(2,2,2)
title('EMSC corrected data')
hold on
plot(Xintvalues,RamanSpectraInttesting,'k')
plot(Xintvalues,RamanSpectraEMSCtesting,'b')
plot(Xintvalues,background,'r')
legend('Average Raw Raman Spectra','EMSC-corrected Raman Spectra',...
    'Subtracted PDMS Background')
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity / counts')
axis tight

subplot(2,2,3)
title('Baseline')
hold on
plot(Xintvalues,RamanSpectraEMSCtesting)
plot(Xintvalues,baseline,'r')
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity / counts')
legend('EMSC-corrected Raman Spectra','Baseline fitted')
axis tight

subplot(2,2,4)
title('Baseline subtracted data')
hold on
plot(Xintvalues,RamanSpectraBaselinetesting)
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity / counts')
axis tight

MapTesting.Close
clear MapTesting

%% Loading the data
fprintf('Loading the Data \n')
tic

for j=selectedDataset
    cd(savingPath) 
    try load([NameSaving{j} 'datasetraw' '.mat']) %move from output to the same folder as the PDMS and biological reference data  
    catch
        warning('Raw dataset could not be found, proceeding to load raw wdf files');
    end
    fprintf(['Loading data from dataset' num2str(j) '\n']);
    cd(dirPath{j});
   
       m=1;


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
    
    % Creating the structure file with the map inside
    Map=WdfReader([Name{j} '.wdf']);   
    if exist('RamanSpectra')==0
        % Obtaining the spectra from the structure file, where each row is a
        % spectrum
        RamanSpectra=Map.GetSpectra(1,Map.Count);
        % Changing it to columns because that is how I usually work with them
        RamanSpectra=RamanSpectra';
        Wavenumbers=Map.GetXList;
        CalibrationPeakTemp=CalibrationPeak{j};
        cd(savingPath) 
        save([NameSaving{j} '_datasetraw_emsc9623'],...
            'RamanSpectra','Wavenumbers','CalibrationPeakTemp')
    end
    toc

%% X alignment of the spectra respect to the Calibration peak position in Si

% Prior to each experiment, the system should be aligned as closely as possible 
% to the 520.5 cm^-1 peak of Silicon. Sometimes the peak can be around 520.3 or 520.8. 
% This peak will be used to align the Shift of the spectra and reduce the 
% variability between days.

% As the PCA only considers the Y-values for analysis,
% a spline function should be fit and the Y-values should be guarantee to
% have exactly the same X.
fprintf('Aligning spectra and fitting to the selected X values using a spline \n')
tic
RamanSpectraInt=NaN(length(Xintvalues),Map.Count);

for k=1:length(RamanSpectra)
  RamanSpectraInt(:,k)=spline(RamanSpectra{k}(:,1)+(520.5-CalibrationPeak{j}),RamanSpectra{k}(:,2),Xintvalues);
end
toc

clear RamanSpectra Wavenumbers


%% Full correction of all spectra with EMSC

% Usign the EMSC function for correction with multiple background spectra.

disp('EMSC correcting')

tic
r_EMSC=r.';
b_EMSC=b.';

% smoothing of the reference and background files to reduce noise in the
% EMSC output
r_EMSC=[sgolayfilt(r_EMSC.',2,17)].';
b_EMSC=[sgolayfilt(b_EMSC.',2,17)].';

RamanSpectraEMSC=NaN(size(RamanSpectraInt));

for i=1:size(RamanSpectraEMSC,2)
    [RamanSpectraEMSC(:,i),background,c_r,c_b,B_N] = EMSC(RamanSpectraInt(:,i).',r_EMSC,b_EMSC,N);
end
toc

%% Subtract the baseline
fprintf('Subtracting the baseline \n')

tic
RamanSpectraBaseline=NaN(size(RamanSpectraEMSC));
for i=1:size(RamanSpectraEMSC,2)
    if mod(i,20)==0
        disp(['Baseline correcting spectrum ' num2str(i) ' of ' num2str(size(RamanSpectraEMSC,2))])
    end
    % First step to remove the maximum possible of the baseline
    if i==1
        tic
    end
    [RamanSpectraBB, ~]=f_baseline_corr(Xintvalues,...
        RamanSpectraEMSC(:,i),smoothwidth,bwidth,iterations);%,'plot2D');    
    RamanSpectraBaseline(:,i)=RamanSpectraBB.';
    if i==10
        elapsedTime = toc;
    end
    if mod(i,20)==0
        disp(['Time left ' num2str(round(elapsedTime*(Map.Count-i)/600)) 'min'])
    end
    clear RamanSpectraBB
end
toc
clear RamanSpectraB RamanSpectraEMSC


%% Data Matrix
fprintf('Creating the Data Matrix \n')
tic
DataMatrix=RamanSpectraBaseline;

%clear RamanSpectraBaseline

%% Calculate the Average and STD

Average=mean(DataMatrix,2);
Error=std(DataMatrix,0,2);
Average2=flip(Average, 2);

%% Plot figures
fprintf('Plotting Figure \n')
figure
subplot(2,1,1)

plot(Xintvalues, Average2)
xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity ')

%hold on
subplot(2,1,2)
plot(Xintvalues,Average2)
xlim([500 1800])

xlabel('Wavenumber / cm^{-1}')
ylabel('Raman Intensity ')
xticks(500:100:1800)
%hold on

%% Saving Data
    cd(savingPath)
    save([NameSaving{j} 'X'],...
        'DataMatrix','Average','Error',...
        'CalibrationPeakTemp','Xintvalues')

        
end

disp('Done! :)')