classdef MapUtils
    %Class for exposing various indexing operations on the Map.  Assumes a
    %rectangular map with indexing running over columns first. 
    
    properties (GetAccess = public, SetAccess = protected)
        Rows = nan;
        Columns = nan;
    end;
    
    methods (Access = public)
        % Constructor; initialise the MapUtils object with # Rows and Columns
        function this = MapUtils(Rows, Columns)
            this.Rows = Rows;
            this.Columns = Columns;
        end;
        
        % Destructor: just calls the Close() method.
        function delete(this)
            this.Close();
        end;
        
        %Returns a sorted spectrum number mask consisting of 1-based 
        %indices of all spectra at the edge of a rectangular area map of 
        %dims Rows x Cols where are within NumberOfSpectra. 
        %If EdgeWidth <> 1 or 2 then an empty array is returned.  
        %ASSUMES RA802 DATA ORDERING / RASTER MAPS.
        function [Result] = EdgeMask(this, EdgeWidth, NumberOfSpectra)
            N = this.Rows*this.Columns;
            Result = [];
            if EdgeWidth > 0 && EdgeWidth < 3
                Result = [linspace(2,this.Columns-1,this.Columns-2),linspace(1,N-this.Columns+1,this.Rows)];
                Result = [Result,linspace(this.Columns,N,this.Rows),linspace(N-this.Columns+2,N-1,this.Columns-2)];
                if EdgeWidth == 2
                    Result = [Result,linspace(this.Columns+3,2*this.Columns-2,this.Columns-4),linspace(this.Columns+2,N-2*this.Columns+2,this.Rows-2)];
                    Result = [Result,linspace(2*this.Columns-1,N-this.Columns-1,this.Rows-2),linspace(N-2*this.Columns+3,N-this.Columns-2,this.Columns-4)];
                end;
                Result = sort(Result(Result <= NumberOfSpectra))';
            end;
        end;
        
        %Returns a sorted spectrum number mask consisting of 1-based
        %indices of all spectra for which all elements are zero.
        function [Result] = ZeroSpectra(this, Spectra)
            Result = sort(find(all(~Spectra')))';
        end

        %Replace the indicated spectra with most proximal neighbour that a) exists
        %and b) is not in the ReplaceIndices list.  Search extends out to map edge  
        %if necessary.  Note there is no assessment as to the similarity of the 
        %spectrum, the function just selects the closest available geographically.  
        %This fuction may be used to replace saturated spectra with unsaturated
        %neighbours.  
        function [ SpectraOut ] = ReplaceSpectraWithNeighbours(this, SpectraIn, ReplaceIndices)
            SpectraOut = SpectraIn;
            maxRadius = max(this.Rows, this.Columns);
            for i = 1:size(ReplaceIndices,1)
                for irad = 1:maxRadius
                    %get list of neighbours with Radius irad in nearness order, 
                    %taking account of their existence (i.e. lying within the 
                    %map dims and excluding any )
                    neighbours = this.ExistingRadiusNeighbours(ReplaceIndices(i), ReplaceIndices, irad);
                    %exit the loop if we have found a neighbour
                    if not(isempty(neighbours)) 
                        break;
                    end;
                end;
                %update the returned spectrum array with the first neighbour
                if not(isempty(neighbours))
                    SpectraOut(ReplaceIndices(i),:) = SpectraIn(neighbours(1),:);
                end;
            end;
        end; 
        
        %Insert the masked indexes data back into the predictions matrix 
        %with the value passed-in. 
        function [ PredictionsOut ] = InsertMaskedIndices(this, PredictionsIn, Mask, Value)
            if isempty(Mask)
                PredictionsOut = PredictionsIn;
            else
                Mask = sort(Mask);  %following assumes that the mask is sorted in ascending order.
                PredictionsOut=ones(1,length(PredictionsIn)+length(Mask));
                m=length(Mask);
                PredictionsOut(1:Mask(1)-1)=PredictionsIn(1:Mask(1)-1);
                for ii=1:m-1;
                    PredictionsOut(Mask(ii))=Value;
                    PredictionsOut(Mask(ii)+1:Mask(ii+1)-1)=PredictionsIn(Mask(ii)+1-ii:Mask(ii+1)-1-ii);
                end;
                PredictionsOut(Mask(m))=Value;
                PredictionsOut(Mask(m)+1:end)=PredictionsIn(Mask(m)+1-m:end);
            end;
        end;

        %Returns a sorted spectrum number mask consisting of all spectral  
        %indices that DO NOT reach the specified signal to noise criterion  
        function [Result] = UnderS2nThreshold(this, XList, Spectra, Config)
            if Config.Active 
                %Sort out region limits
                noiseLow = min(Config.NoiseRegion1, Config.NoiseRegion2);
                noiseHigh = max(Config.NoiseRegion1, Config.NoiseRegion2);
                signalLow = min(Config.SignalRegion1, Config.SignalRegion2);
                signalHigh = max(Config.SignalRegion1, Config.SignalRegion2);
                %Check that data covers signal region 
                nx=size(XList,2);
                if (signalLow - XList(1))*(signalLow - XList(nx)) > 0 || (signalHigh - XList(1))*(signalHigh - XList(nx)) > 0
                    error('Signal-to-noise filter cannot be applied because the data, after any truncation, does not cover the signal region.') 
                end
                %Calculate rms pp noise over noise evaluation region
                noiseSpec = Spectra(:,XList > noiseLow & XList < noiseHigh);
                %Check noise region contains at least some points
                if size(noiseSpec,2) < 2
                    error('Signal-to-noise filter cannot be applied because there are insufficient data in the noise region.') 
                end
                %Following avoids Signal Processing toolbox - verified vs. rms(diff(noiseSpec,1,2),2)
                rmsNoise = sqrt(sum(diff(noiseSpec,1,2).^2,2)/(size(noiseSpec,2)-1));
                %Get CSpline interpolated values at signalLow and signalHigh -
                %no transpose, so avoid transpose in next interp1 call
                baseSignal = interp1(XList,Spectra',[signalLow,signalHigh],'spline');  
                %Limited X list and spectra over max signal region
                signalXList = XList(XList > signalLow & XList < signalHigh);
                signalSpec = Spectra(:,XList > signalLow & XList < signalHigh);
                %Linear baselines between these pairs of points over this region 
                baselines = (interp1([signalLow,signalHigh],baseSignal,signalXList,'linear'))';
                %Remove linear baselines
                signalSpec = signalSpec - baselines;
                %Get signal max above linear baseline
                maxSignal = max(signalSpec,[],2);            
                %Calculate the signal to rms pt-pt noise ratio
                s2n = maxSignal ./ rmsNoise;
                %Get list of s2n vector indeces below s2nThreshold
                Result = find(s2n < Config.Threshold);    
            else
                Result = [];
            end;
        end;

    end; %end public methods
    
    methods (Access = protected)
        
        %Construct an array of neighbours to the point MapIndex, that are 
        %at radius = Radius.  Exclude any indices lying outside the map or in 
        %the ExcludeIndices array.
        function [ Neighbours ] = ExistingRadiusNeighbours(this, MapIndex, ExcludeIndices, Radius)
            mapSize = this.Rows * this.Columns;
            %rooks move - first priority
            Neighbours = [MapIndex-Radius*this.Columns; MapIndex-Radius; 
                MapIndex+Radius; MapIndex+Radius*this.Columns];
            %everything inbetween
            for i = 1:Radius-1
                Neighbours = [Neighbours; MapIndex-Radius*this.Columns-i; 
                    MapIndex-Radius*this.Columns+i; MapIndex-i*this.Columns-Radius;
                    MapIndex-i*this.Columns+Radius; MapIndex+i*this.Columns-Radius; 
                    MapIndex+i*this.Columns+Radius; MapIndex+Radius*this.Columns-i;
                    MapIndex+Radius*this.Columns+i];
            end;
            %corners / bishops move
            Neighbours = [Neighbours; MapIndex-Radius*this.Columns-Radius; 
                MapIndex-Radius*this.Columns+Radius; MapIndex+Radius*this.Columns-Radius
                MapIndex+Radius*this.Columns+Radius]; 
            %remove indices lying outside the map
            Neighbours((Neighbours < 1) | (Neighbours > mapSize)) = [];
            %remove indices in the ExcludeIndices map
            Neighbours = Neighbours .* ~ismember(Neighbours,ExcludeIndices);
            Neighbours(Neighbours == 0) = [];            
        end;
        
    end; %end protected methods
    
end

