classdef WdfUtils < handle
    % Attach the wdf object to the WdfUtils object, but only allow set
    % access, to prevent accessing all other methods within the wdf object
    % from here.
    %
    % To use, set up a wdf object using 
    % wdf = WdfReader(WdfName);
    % Then set up a WdfUtils object as 
    % Utils = WdfUtils;
    % Set the wdf in the utils object to the required one as
    % Utils.wdf = wdf;
    % We can then call the methods in here using (for example)
    % MeanSpectrum = Utils.GetMeanSpectrum(false)
    % Which will generate the mean spectrum, without a waitbar.
    
    properties (SetAccess = public, GetAccess = private)
        wdf = '';
    end;
    methods (Access = public)
        
        % Calculates the mean spectrum of the data. Use the argument to
        % optionally show a progress bar.
        function MeanSpec = GetMeanSpectrum(this, showProgress)
            if(showProgress)
                h = waitbar(0,'Calculating mean...');
            end
            MeanSpec = zeros(1,this.wdf.PointsPerSpectrum());
            TotalReadInSoFar = 0;
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                Weight = (TotalReadInSoFar/(TotalReadInSoFar+NumberReadIn));
                MeanSpec = (1-Weight).*mean(z)+Weight.*MeanSpec;
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        % Calculates the mean spectrum of the data. Use the input to
        % optionally pass in the mean spectrum (if previously calculated)
        % or show a progress bar.
        function VarianceSpec = GetVarianceSpectrum(this, varargin)
            if(nargin == 1)
                MeanSpec = this.wdf.GetMeanSpectrum(false);
                showProgress = false;
            elseif(nargin == 2)
                MeanSpec = varargin{1};
                showProgress = false;
            elseif(nargin == 3)
                MeanSpec = varargin{1};
                showProgress = varargin{2};
            end
            if(showProgress)
                h = waitbar(0,'Calculating variance...');
            end
            VarianceSpec = zeros(1,this.wdf.PointsPerSpectrum());
            this.wdf.StartChunkwiseReading();
            TotalReadInSoFar = 0;
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                VarianceSpec = VarianceSpec + sum((z-MeanSpec).*(z-MeanSpec));
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        
        % Generates a map of the total intensity of each spectrum (sum of
        % the spectrum). Note this total intensity map is trivially related
        % to the mean value. 
        % Optionally show a progress bar by passing in true or false.
        function TotalIntensityMap = CalculateIntegratedIntensity(this, showProgress)
            TotalIntensityMap = zeros(this.wdf.Count(),1);
            TotalReadInSoFar = 1;
            if(showProgress)
                h = waitbar(0,'Calculating intensity map...');
            end
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                TotalIntensityMap(TotalReadInSoFar:TotalReadInSoFar+NumberReadIn-1) = sum(z,2);
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
        % Generates a map of the value of the spectrum at the nearest
        % wavenumber in the xlist to the chosen one.
        % Optionally show a progress bar by passing in true or false along
        % with the desired wavenumber.
        function TotalIntensityMap = PeakIntensityMap(this, Wavenumber, showProgress)
            xList = this.wdf.GetXList();
            Position = find(xList < Wavenumber, 1);
            TotalIntensityMap = zeros(this.wdf.Count(),1);
            TotalReadInSoFar = 1;
            if(showProgress)
                h = waitbar(0,'Calculating peak intensity map...');
            end
            this.wdf.StartChunkwiseReading();
            while (this.wdf.AreMoreChunks())
                z = this.wdf.GetNextDataChunk();
                NumberReadIn = size(z,1);
                TotalIntensityMap(TotalReadInSoFar:TotalReadInSoFar+NumberReadIn-1) = z(:,Position);
                TotalReadInSoFar = TotalReadInSoFar + NumberReadIn;
                if(showProgress)
                    waitbar(TotalReadInSoFar/this.wdf.Count());
                end
            end
            if(showProgress)
                close(h);
            end
        end
    end
end