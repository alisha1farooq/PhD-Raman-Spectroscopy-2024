% WdfFile  Provides access to selected data in WDF files.
%
% Renishaw's WDF file format is designed for storing Raman spectral data.
% The WdfFile class provides read-only access to a subset of the data in
% a WDF file, including:
%   * The spectral data;
%   * Data origin lists (see below for more details); and
%   * Key file metadata (such as measurement title and username).
%   * Saturated spectra and cosmic ray removed spectra
%   * Loadings and explained variance from principle components analysis.
% It allows modified spectra to be written back to the file.
%
% wdf = WdfFile(FILENAME) creates a new WdfFile object, opening the
% WDF file specified by FILENAME.  
%
% wdf.Close() closes the file.
%
%
% Accessing spectral data
% =======================
%
% SPECTRA = wdf.GetSpectra(IndexStart, IndexEnd) reads the specified
% spectra from the file.  Indices are 1-based.
%
% wdf.WriteSpectra(IndexStart, Data) writes the specified spectra back to
% the file, starting at the specified spectrum index.
%
% Spectral datasets are associated with an 'X-list' and 'Y-list'; for
% normal spectra the Y-list contains only one element and can be ignored.
% The meaning of the data in these lists can be determined by the relevant
% properties of the WdfFile object.
%
% XLIST = wdf.GetXList() retrieves the X-list data.
% YLIST = wdf.GetYList() retrieves the Y-list data.
%
% For large files that cannot fit into memory, data can be processed in
% manageable consecutive chunks using the following set of functions:
%    wdf.StartChunkwiseReading();
%    [X, indices] = wdf.GetNextChunk(numberOfSpectra);
%    trueOrFlase  = wdf.AreMoreChunks();
%    fraction     = wdf.GetChunkwiseProgress();
% (note that numberOfSpectra is optional for GetNextChunk()).
%
%
% Accessing Data Origin List values
% =================================
%
% WDF files may contain one or more Data Origin Lists.  Each list stores a
% unique type of data (such as time, or position in the X-axis) for each
% spectrum / dataset in the file.  Each Data Origin List is identified by
% the WiREDataType of the values it stores.  Lists are either 'primary'
% (important) or 'alternate' (less important).  The data in each list can
% either be stored as floating-point or integer values, and the caller must
% use the correct method to access the data (otherwise meaningless values
% will be returned); this can typically be inferred from the WiREDataType
% value associated with the list.
% 
% INFO = wdf.GetOriginListInfo() returns an Nx4 cell-array where each row
% describes a Data Origin List in the file.  The columns contain:
%   (1) A LOGICAL value indicating if the list is a 'primary' origin list;
%   (2) The list data type (as a WiREDataType value);
%   (3) The units of measurement (as a WiREDataUnit value); and
%   (4) The list name.
%
% V = wdf.GetOriginListValues(WiREDataType, IndexStart, IndexEnd) reads the
% specified values from the Data Origin List with the specified data type,
% treating the binary data in the origin list as double-precision floating
% point values.  Indices are 1-based.
%
% V = wdf.GetOriginListValuesInt(WiREDataType, IndexStart, IndexEnd) is
% equivalent to GetOriginListValues(), but treats the binary data in the
% origin list as 64-bit integer values.
%
%
% WDF file metadata
% =================
%
% WdfFile exposes the following properties:
%    Title               The measurement title.
%    Username            The name of the user who acquired the data.
%    MeasurementType     The measurement type of the data in the file.
%    ScanType            The scan type used to acquire the data.
%    LaserWavenumber     The wavenumber of the laser used.
%    Count               The actual number of spectra that are stored in
%                        the file (the number of spectra collected).
%    SpectralUnits       A WiREDataUnit value indicating the units of
%                        measurement associated with the spectral data.
%    XListType           A WiREDataType value indicating the type (or
%                        meaning) of the values stored in the X-list.
%    XListUnits          A WiREDataUnit value indicating the units of
%                        measurement associated with the X-list data.
%    YListType           A WiREDataType value indicating the type (or
%                        meaning) of the values stored in the Y-list.
%    YListUnits          A WiREDataUnit value indicating the units of
%                        measurement associated with the Y-list data.
%    PointsPerSpectrum   The number of elements in each spectrum / dataset.
%                        Equal to (XListLength * YListLength).
%    DataOriginCount     The number of Data Origin Lists in the file.
%    Capacity            The maximum number of spectra that can be stored
%                        in the file.
%    ApplicationName     The name of the software used to acquire the data.
%    ApplicationVersion  The software version used to acquire the data (as
%                        a 4-element vector: [major, minor, patch, build]).
%    XListLength         The number of elements in the X-list.
%    YListLength         The number of elements in the Y-list.
%    AccumulationCount   The number of accumulations co-added for each
%                        spectrum.
%
% Accessing Spectra with Saturation and Cosmic Ray Removal
% The WMAP block for formatting the orientation
% ========================================================
%
% [Saturated CosmicRay] = wdf.GetOriginFlags(IndexStart,IndexEnd) will return
% two lists. The first is a list of all spectra that have been flagged as 
% suffering from saturation. The second is a list of all spectra that have
% had cosmic rays removed. INDEXSTART and INDEXEND define the range of
% spectra to consider, where INDEXSTART is 1 or greater (ie corresponds to 
% spectrum 0).
%
% SATURATED and COSMICRAY return vectors of numbers that correspond to the
% Nth spectrum. For instance an appearance of the number 59 means that the
% 59th spectrum has been affected by the relevant flag.
% 
% Accessing the Loadings Data from Principle Component Analysis
% =============================================================
%
% [Loadings VarianceExplained] = wdf.getLoadingsData(PCANumber) will
% return both the loadings for all principle components computed in the
% WiRE software and the corresponding percentage of the overall variance
% that each principle component represents.
% 
% %LOADINGS is a 2D array in which the first column corresponds to the
% Raman shift in wavenumbers, and each subsequent column the corresponding
% loadings for each principle component at each Raman shift. The second
% column is PC1, the third PC2 etc, listed in order of decreasing
% percentage of variance explained. 
%
% VARIANCEXPLANINED is a vector listing each of the percentage variance 
% explained values. THe first entry corresponds to PC1, the second to PC2
% etc.
%
% PCANUMBER corresponds to the specific PCA loadings data that are required.
% If one has performed PCA three times, there will be three sets of PCA data
% stored. If the user wanted to extract the loadings from the first
% analysis, PCANUMBER would be equal to 1.
%
% Unsupported methods
% ===================
%
% The following methods are exposed by WdfFile but are intended for use
% only by Renishaw:
%    GetBlockID(),  GetFullBlockID(),  GetNextAvailableBlockUID(),
%    GetPrimaryOriginLists().

% Copywrite Notice:
% (c) 2012 - 2014 Renishaw plc. All rights reserved.

classdef WdfFile < handle
    %% Private fields
    properties (Access = protected)
        m_NextSpectrumIndex = 1;
        Handle = -1;
        XListOffset = int64(-1);
        YListOffset = int64(-1);
        OriginsOffset = int64(-1);
        DataOffset = int64(-1);
        OriginListInfo = cell(0, 4);
        PCAOffset = int64(-1);
        MapOffset = int64(-1);
    end;
    
    %% Public properties relating to WDF metadata
    properties (GetAccess = public, SetAccess = protected)
        Title = '';
        Username = '';
        MeasurementType = WiREMeasurementType.Unspecified;
        ScanType = WiREScanType(0);
        LaserWavenumber = nan;
        Count = nan;
        SpectralUnits = WiREDataUnit.Arbitrary;
        XListType = WiREDataType.Arbitrary;
        XListUnits = WiREDataUnit.Arbitrary;
        YListType = WiREDataType.Arbitrary;
        YListUnits = WiREDataUnit.Arbitrary;
        PointsPerSpectrum = nan;
        DataOriginCount = nan;
        Capacity = nan;
        ApplicationName = '';
        ApplicationVersion = nan(1, 4);
        XListLength = nan;
        YListLength = nan;
        AccumulationCount = nan;
        OPut = cell(1,2);
        WMAPFlag = struct;
    end;
    
    %% Constructor, destructor and Close
    methods
        % Constructor; opens an existing WDF file with read-only access
        % using the protected constructor.
        function this = WdfFile(fileName)
            this.ProtectedConstructor(fileName, 'r+b');
        end;
        
        % Closes the file.
        function Close(this)
            if (this.Handle ~= -1)
                fclose(this.Handle);
                this.Handle = -1;
            end;
        end;
        
        % Destructor: just calls the Close() method.
        function delete(this)
            this.Close();
        end;
    end;
    
    %% Protected helper methods
    methods (Access = protected)
        % Opens a WDF file with the specified access mode, reads the WDF
        % header, and validates the DATA, XLST, YLST and ORGN blocks.
        function ProtectedConstructor(this, fileName, fileAccess)
            % Open the file
            this.Handle = fopen(fileName, fileAccess);
            if (this.Handle == -1)
                throw(WdfError('Cannot open file "%s".', fileName));
            end;
            
            % Parse the header, locate the data block and check its size is
            % consistent with expected value.
            this.ReadWdfHeader();
            [this.DataOffset, ~, dataSize] = LocateBlock(this, 'DATA');
            if (this.DataOffset == -1)
                throw(WdfError('Cannot locate Spectral Data block.'));
            end;
            if (dataSize < (16 + (4 * this.PointsPerSpectrum * this.Capacity)))
                throw(WdfError('Spectral Data block size inconsistent with WDF header.'));
            end;
            
            % Locate the X-list and Y-list blocks, and validate their sizes.
            [this.XListOffset, ~, dataSize] = LocateBlock(this, 'XLST');
            if (this.XListOffset == -1)
                throw(WdfError('Cannot locate X-list block.'));
            end;
            if (dataSize < (24 + (4 * this.XListLength)))
                throw(WdfError('X-list block size inconsistent with WDF header.'));
            end;
            
            [this.YListOffset, ~, dataSize] = LocateBlock(this, 'YLST');
            if (this.YListOffset == -1)
                throw(WdfError('Cannot locate Y-list block.'));
            end;
            if (dataSize < (24 + (4 * this.YListLength)))
                throw(WdfError('Y-list block size inconsistent with WDF header.'));
            end;
            
            % Read the X and Y list types / units
            this.XListType = GetXListType(this);
            this.YListType = GetYListType(this);
            this.XListUnits = GetXListUnits(this);
            this.YListUnits = GetYListUnits(this);
            
            % If the Data Origin List count is non-zero, attempt to locate
            % then size-check the ORGN block.  Finally, read in the data
            % origin list info.
            if (this.DataOriginCount ~= 0)
                [this.OriginsOffset, ~, dataSize] = LocateBlock(this, 'ORGN');
                if (this.OriginsOffset == -1)
                    throw(WdfError('Cannot locate Data Origin List block.'));
                end;
                if (dataSize < (20 + (this.DataOriginCount * (24 + (this.Capacity * 8)))))
                    throw(WdfError('Data Origin List block size inconsistent with WDF header.'));
                end;
            end;
            this.ReadOriginListInfo();
            
        this.WMAPFlag = WMapblock(this);
        
        end;
           
        % Reads the WDF header block, extracting key header fields and
        % storing them in class properties.
        function ReadWdfHeader(this)
            % Confirm that the signature, version and size fields
            % contain the expected values.
            this.AssertSeek(0, 'bof');
            blockID = this.AssertRead(1, 'uint32');
            blockUID = this.AssertRead(1, 'uint32');
            blockLength = this.AssertRead(1, 'uint64');
            if (~isequal(blockID, this.GetBlockID('WDF1')) || ...
                    ~(isequal(blockUID, 0) || isequal(blockUID, 1)) || ...
                    ~isequal(blockLength, 512))
                throw(WdfError('File does not use a recognised WDF format / version.'));
            end;
            
            % Read key fields from the main file header
            this.AssertSeek(60, 'bof');
            this.PointsPerSpectrum = this.AssertRead(1, 'uint32');
            this.Capacity = this.AssertRead(1, 'uint64');
            this.Count = this.AssertRead(1, 'uint64');
            this.AccumulationCount = this.AssertRead(1, 'uint32');
            this.YListLength = this.AssertRead(1, 'uint32');
            this.XListLength = this.AssertRead(1, 'uint32');
            this.DataOriginCount = this.AssertRead(1, 'uint32');
            this.ApplicationName = this.ReadUtf8String(24);
            this.ApplicationVersion = this.AssertRead([1 4], 'uint16');
            this.ScanType = WiREScanType(this.AssertRead(1, '*uint32'));
            this.MeasurementType = WiREMeasurementType(this.AssertRead(1, '*uint32'));
            this.AssertSeek(152, 'bof');
            this.SpectralUnits = WiREDataUnit(this.AssertRead(1, '*uint32'));
            this.LaserWavenumber = this.AssertRead(1, 'float32');
            this.AssertSeek(208, 'bof');
            this.Username = this.ReadUtf8String(32);
            this.Title = this.ReadUtf8String(160);
        end;
        

        
        % Searches for a specific block within the WDF file, by ID and
        % optionally UID.  Returns the location (offset), UID and length
        % (in bytes) of the first matching block found, or -1 if no
        % matching block was located.
        function [Location, BlockUID, BlockSize] = LocateBlock(this, TargetID, TargetUID)
            % Set default return values if block is not found
            Location = -1;
            BlockUID = [];
            BlockSize = uint64(0);
            
            % Initialise search
            found = false;
            offset = int64(512);
            TargetID = this.GetBlockID(TargetID);
            
            % Iteratively walk over all data-blocks in the file
            while ((~found) && (feof(this.Handle) == 0))
                % Attemp to jump to next block
                if (fseek(this.Handle, offset, 'bof') ~= 0)
                    return;
                end;
                
                % Read block header, and abort search if header fields not
                % read successfully
                blockID = fread(this.Handle, 1, 'uint32');
                uid = fread(this.Handle, 1, 'uint32');
                bSize = fread(this.Handle, 1, 'uint64');
                if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                    return;
                end;
                
                % Check if we have found the requested block
                if (nargin == 3)
                    found = (isequal(TargetID, blockID) && isequal(TargetUID, uid));
                else
                    found = isequal(TargetID, blockID);
                end;
                
                % Store block info if match found, otherwise prepare to
                % look for the next block
                if (found)
                    Location = offset;
                    BlockUID = uid;
                    BlockSize = bSize;
                else
                    offset = offset + bSize;
                end;
            end;
        end;
        
        % Reads the specified number of bytes and converts them from UTF8
        % to a Matlab string, with trailing NULLs and whitespace trimmed.
        function s = ReadUtf8String(this, Length)
            s = this.AssertRead([1 Length], 'uint8=>char');
            s = deblank(char(unicode2native(s, 'UTF-8')));
        end;
        
       
        % Seeks to the requested position in the WDF file, raising an error
        % if the seek operation is unsuccessful.
        function AssertSeek(this, Position, Origin)
            if (fseek(this.Handle, Position, Origin) ~= 0)
                throwAsCaller(WdfError('Failed to seek to requested position within file.'));
            end;
        end;
        
        % Reads the requested data from the WDF file, raising an error if
        % the actual number of elements read is fewer than requested.
        function [Data] = AssertRead(this, RDims, Precision)
            [Data, readCount] = fread(this.Handle, RDims, Precision);
            if (readCount ~= prod(RDims))
                throwAsCaller(WdfError('Failed to read requested data from file.'));
            end;
        end;
        
        % Writes the specified data to the WDF file, in column order, and
        % raising an error if the actual number of elements written is
        % fewer than requested.
        function AssertWrite(this, Data, Precision)
            [writeCount] = fwrite(this.Handle, Data, Precision);
            if (writeCount ~= numel(Data))
                throwAsCaller(WdfError('Failed to write requested data to file.'));
            end;
        end;
        
        % Helper function used to read data origin list values.
        function [Data] = ReadOriginListData(this, ListType, IndexStart, IndexEnd, Precision)
            % Validate inputs.
            if (nargin ~= 5)
                throwAsCaller(WdfError('ListType, IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throwAsCaller(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throwAsCaller(WdfError('IndexEnd is out-of-range.'));
            end;
            listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                throwAsCaller(WdfError('WDF file does not contain a Data Origin List with the specified type.'));
            end;
            
            % Read and return the data.
            count = (IndexEnd - IndexStart) + 1;
            listSize = 24 + (8 * this.Capacity);
            offset = this.OriginsOffset + 20 + ((listIndex - 1) * listSize) + 24 + ((IndexStart - 1) * 8);
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(count, Precision);
        end;
        
        %Check that the list origin actually exists. Currently used to
        %check the origin list for the presence of (x,y) coordinates.
        function Check = CheckOriginListData(this, ListType, IndexStart, IndexEnd, Precision)
            % Validate inputs.
            if (nargin ~= 5)
                throwAsCaller(WdfError('ListType, IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throwAsCaller(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throwAsCaller(WdfError('IndexEnd is out-of-range.'));
            end;
            listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
            if (isempty(listIndex))
                Check = false;
                return
            end;
            Check = true;
        end;
    end;
    

   
    
    %% Non-public methods that access key WDF data
    methods (Access = public)
        % Get the X-list type
        function [Result] = GetXListType(this)
            this.AssertSeek(this.XListOffset + 16, 'bof');
            Result = WiREDataType(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the X-list units
        function [Result] = GetXListUnits(this)
            this.AssertSeek(this.XListOffset + 20, 'bof');
            Result = WiREDataUnit(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the Y-list type
        function [Result] = GetYListType(this)
            this.AssertSeek(this.YListOffset + 16, 'bof');
            Result = WiREDataType(this.AssertRead(1, '*uint32'));
        end;
        
        % Get the Y-list units
        function [Result] = GetYListUnits(this)
            this.AssertSeek(this.YListOffset + 20, 'bof');
            Result = WiREDataUnit(this.AssertRead(1, '*uint32'));
        end;
        
        % Creates and stores a cell-array containing information about the
        % Data Origin Lists in the WDF file (one row per data origin list).
        function ReadOriginListInfo(this)
            ListInfo = cell(this.DataOriginCount, 4);
            listSize = 24 + (8 * this.Capacity);
            for n = 1:this.DataOriginCount
                this.AssertSeek(this.OriginsOffset + 20 + ((n - 1) * listSize), 'bof');
                listInfo = this.AssertRead([1 2], '*uint32');
                ListInfo{n, 1} = bitget(listInfo(1), 32) ~= 0;
                ListInfo{n, 2} = WiREDataType(bitset(listInfo(1), 32, 0));
                ListInfo{n, 3} = WiREDataUnit(listInfo(2));
                ListInfo{n, 4} = this.ReadUtf8String(16);
            end;
            this.OriginListInfo = ListInfo;
        end;
   end; 
   

    
    %% Access to spectra in the file, plus X- and Y-list data
    methods
        
        function [WL,x,y]=GetWL(this)
            this.AssertSeek(0, 'bof'); % go to beginning of file
            [wl_offset, ~, dataSize] = LocateBlock(this, 'WHTL'); % find the white light location and get data size
            if (wl_offset == -1) % if no data found throw error
                throw(WdfError('Cannot locate WL image.'));
            end;
            this.AssertSeek(wl_offset+16, 'bof'); % go to the image position
            imgdata = fread(this.Handle, dataSize, '*uint8'); % read all the data as bytes
            % make use of imread function to avoid having to write a function to decompress jpgs
            fid = fopen('img_temp.jpg', 'w'); % open a temporary jpg file with write access
            fwrite(fid, imgdata,'*uint8'); % write the data as bytes
            fclose(fid); % close the temporary file
            WL=imread('img_temp.jpg'); % read the temporary file in using imread function
            info=imfinfo('img_temp.jpg');
            % create x and y axes from info
            FOV=info.UnknownTags(2).Value./info.UnknownTags(3).Value;
            x=linspace(info.UnknownTags(1).Value(1),info.UnknownTags(1).Value(1)+FOV(1),info.Width);
            y=linspace(info.UnknownTags(1).Value(2),info.UnknownTags(1).Value(2)+FOV(2),info.Height);
        end;
        
        % Read the X-list
        function [XList] = GetXList(this)
            this.AssertSeek(this.XListOffset + 24, 'bof');
            XList = double(this.AssertRead([1 this.XListLength], 'single'));
        end;
        
        % Read the Y-list
        function [YList] = GetYList(this)
            this.AssertSeek(this.YListOffset + 24, 'bof');
            YList = double(this.AssertRead([1 this.YListLength], 'single'));
        end;
        
        % Import a range of spectra, specified by start & end indices.
        function [Data] = GetSpectra(this, IndexStart, IndexEnd,type)
            % Validate inputs.
            if (nargin < 3)
                throw(WdfError('IndexStart and IndexEnd must be specified.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throw(WdfError('IndexStart is out-of-range.'));
            end;
            if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
                throw(WdfError('IndexEnd is out-of-range.'));
            end;
            
            % Read the data, and confirm the expected quantity was read
            offset = (this.DataOffset + 16) + ((IndexStart - 1) * 4 * this.PointsPerSpectrum);
            nRowsToRead = (IndexEnd - IndexStart) + 1;
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(nRowsToRead * this.PointsPerSpectrum, 'single=>single');
            if nargin < 4 || strcmp(type,'double') 
                % Re-shape the data and convert to double-type
                Data = double(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            elseif strcmp(type,'single')
                 Data = single(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            else
                 throw(WdfError('Type must be single or double'))
            end
        end;
        
        % Write spectra back to the file, starting at the specified start
        % spectrum index.
        function WriteSpectra(this, IndexStart, Data)
            if (nargin ~= 3)
                throw(WdfError('IndexStart and Data must be specified.'));
            end;
            if (size(Data, 2) ~= this.PointsPerSpectrum)
                throw(WdfError('Number of columns in Data must match PointsPerSpectrum.'));
            end;
            if ((IndexStart < 1) || (IndexStart > this.Count))
                throw(WdfError('IndexStart is out-of-range.'));
            end;
            if (IndexStart + size(Data, 1) - 1 > this.Count)
                throw(WdfError('IndexStart and number of rows in Data would result in writing beyond end of file.'));
            end;
            
            % Write the data.
            offset = (this.DataOffset + 16) + ((IndexStart - 1) * 4 * this.PointsPerSpectrum);
            this.AssertSeek(offset, 'bof');
            this.AssertWrite(Data', 'float32');
        end;
        
        
        % Initiates chunk-wise reading of the spectral data
        function StartChunkwiseReading(this)
            this.AssertSeek(this.DataOffset + 16, 'bof');
            this.m_NextSpectrumIndex = 1;
        end;
        
        % Reads the next chunk of spectra from the file; NumberOfSpectra is
        % optional.
        function [Data, Indices] = GetNextDataChunk(this, NumberOfSpectra)
            % Determine how many rows (spectra) to read into this chunk.
            if (nargin == 2)
                if (NumberOfSpectra < 1)
                    throw(WdfError('NumberOfSpectra must be greater-than-or-equal-to 1.'));
                end;
                nRowsToRead = NumberOfSpectra;
            else
                nRowsToRead = 4096;
            end;
            
            % If necessary, limit the number of rows to read to the number
            % of remaining available rows.
            if (this.m_NextSpectrumIndex + nRowsToRead > (this.Count + 1))
               nRowsToRead = (this.Count - this.m_NextSpectrumIndex) + 1;
            end;
            
            % Check there are *some* more spectra remaining.
            if (nRowsToRead < 1)
                throw(WdfError('All spectra have already been read in previous chunks.'));
            end;
            
            % Read the data.
            offset = (this.DataOffset + 16) + ((this.m_NextSpectrumIndex - 1) * this.PointsPerSpectrum * 4);
            this.AssertSeek(offset, 'bof');
            Data = this.AssertRead(nRowsToRead * this.PointsPerSpectrum, 'single');
            
            % Calculate the spectrum indices of the data we have just read.
            if (nargout > 1)
                Indices = this.m_NextSpectrumIndex + colon(0, nRowsToRead - 1);
            end;
                
            % Re-shape the data, convert to double type, and update the
            % next read-position marker.
            Data = double(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
            this.m_NextSpectrumIndex = this.m_NextSpectrumIndex + nRowsToRead;
        end;
        
        % Returns a boolean indicating if there are more chunks of spectra
        % to be read from the file.
        function [MoreChunks] = AreMoreChunks(this)
            MoreChunks = (this.m_NextSpectrumIndex <= this.Count);
        end;
        
        % Returns a value in the range [0-1] that indicates the fraction of
        % spectra already imported via chunk-wise reading.
        function [Fraction] = GetChunkwiseProgress(this)
            Fraction = (this.m_NextSpectrumIndex - 1) / this.Count;
        end;
    end;
    
    %% Access to Data Origin Lists
    methods
        % Returns a cell-array containing information about the Data Origin
        % Lists in the WDF file (one row per data origin list).
        function [ListInfo] = GetOriginListInfo(this)
            ListInfo = this.OriginListInfo;
        end;
        
        % Gets a cell-array containing a list of primary (non-alternate)
        % data origin list data types.
        function [PrimaryOriginLists] = GetPrimaryOriginLists(this)
            PrimaryOriginLists = cell(sum(cell2mat(this.OriginListInfo(:, 1))), 2);
            
            % Cycle over the data origin lists, and gather all primary
            % lists into the result.
            resultIndex = 1;
            for n = 1:this.DataOriginCount
                if (this.OriginListInfo{n, 1})
                    PrimaryOriginLists{resultIndex, 1} = sprintf('DataList%d', resultIndex - 1);
                    PrimaryOriginLists{resultIndex, 2} = int32(this.OriginListInfo{n, 2});
                    resultIndex = resultIndex + 1;
                end;
            end;
        end;
        
        % Gets a range of values from a Data Origin List, specified by list
        % data type, treating the binary data as double-precision floating
        % point values.
        function [Data] = GetOriginListValues(this, ListType, IndexStart, IndexEnd)
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
        end;
        
        % Gets a range of values from a Data Origin List, specified by list
        % data type, treating the binary data as 64-bit integer values.
        function [Data] = GetOriginListValuesInt(this, ListType, IndexStart, IndexEnd)
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
        end;
        
        %Finds out which spectra have a saturated spectrum or cosmic ray 
        %removal flag associated with them. Returns the indeces of the
        %flagged datasets corresponding to the (1-based) data array 
        %returned by the class constructor.
        function [Saturated, CosmicRay] = GetOriginFlags(this,IndexStart,IndexEnd)
            ListType = WiREDataType.Flags;
            Data = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
            S = 7*ones(length(Data),1); Data = bitand(Data,S);
            %Find any byte with bit 1 set - DataSaturated flag
            A = find(Data==1); B = find(Data==3); C = find(Data==5); D = find(Data==7); 
            Saturated = sort([(A); (B); (C); (D)]); clear A B C D
            %Find any byte with bit 4 set - Cosmic Ray Removed flag
            A = find(Data==4); B = find(Data==5); C = find(Data==6); D = find(Data==7);
            CosmicRay = sort([(A); (B); (C); (D)]); clear A B C D
        end;
        
        %Extract all x and y coordinates. Note that streamline and
        %streamline HR have transposed (x,y) coordinates relative to each
        %other.
        function [Xcoord, Ycoord] = GetOriginCoords(this,IndexStart,IndexEnd)
            %if(~isempty(strfind(this.WMapblock.flag,'ColumnMajor')))
            %    ListType = WiREDataType.SpatialY;
            %    Check = this.CheckOriginListData(ListType,IndexStart,IndexEnd,'double');
            %    if Check
            %        Xcoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
            %        ListType = WiREDataType.SpatialX;
            %        Ycoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
            %    else
            %        Xcoord = []; Ycoord = [];
            %    end
            %elseif(~isempty(strfind(this.WMapblock.flag,'RowMajor')))
                ListType = WiREDataType.SpatialX;
                Check = this.CheckOriginListData(ListType,IndexStart,IndexEnd,'double');
                if Check
                    Xcoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
                    ListType = WiREDataType.SpatialY;
                    Ycoord = this.ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
                else
                    Xcoord = []; Ycoord = [];
                end
            %end
        end
    end;
    
    %% Access to loadings data stored in PCA maps.
    
   methods
        
        %Access the loadings from PCA analysis data. This set of methods
        %makes use of the new PSet parser. 
        function [Loadings, VarEx] = getLoadingsData(this,PCAno)
            qq = 2;
            NM = this.getNumberofID('MAP ');
            for ii = 1:NM
                [this.MapOffset, ~, ~] = LocateBlock(this,'MAP ',ii);
                this.OPut = this.ReadPSetData();
                CP = this.FindParent('PCAR',PCAno);   %Only needed for PCA maps
                Disc = this.ReadCellValue('overlayXListType');
                if CP == 1 && Disc~=-1
                    V(qq-1) = this.ReadCellValue('%VarianceExplained');
                    L(:,1) = this.ReadCellValue('overlayXList');
                    L(:,qq) = this.ReadCellValue('overlaySpectrum');
                    qq = qq + 1;
                end;
            end;
            if exist('V','var')
                VarEx = V; Loadings = L;
            else 
                throw(WdfError('Invalid PCANUMBER'))
            end
        end;
        
         %Parse the data from the PCA Map psets.
         function OPt = ReadPSetData(this)
             OPt = this.PSetParser(20);
         end;
         
         %The actual parser. Set up in this manner to permit recursive
         %calls to itself. This is necessary for nested PSets, which
         %certain properties that are themselves PSets.
         function OP = PSetParser(this,nbytes)
             key = cell(1,2); OP = cell(1,2); 
             for xx = 1:2      %Parse twice: 1st for key, 2nd for output
             this.AssertSeek(this.MapOffset + nbytes,'bof')
             TotSI = this.AssertRead(1,'uint32'); TotS = TotSI; %Size of PSet Block Content
                 n = 1;
                 while TotS>0
                     [Type, F, k] = this.ReadFlagKey();
                     switch Type
                         case 99; Txt = char('uint8');  by = 1;
                         case 115; Txt = char('uint16');by = 2;
                         case 105; Txt = char('uint32');by = 4;
                         case 119; Txt = char('uint64');by = 8;
                         case 114; Txt = char('single');by = 4;
                         case 113; Txt = char('double');by = 8;
                         case 116; this.AssertRead(1,'uint32'); this.AssertRead(1,'uint32'); by = 8;
                         case 117; [A, by] = this.ReadPSetString();
                         case 107; [A, by] = this.ReadPSetString();
                                   if xx == 1
                                       key{n,1} = k; key{n,2} = A; 
                                       n=n+1;
                                   end;
                         case 112; by = this.AssertRead(1,'uint32')+4;
                                   A = this.PSetParser(TotSI-TotS+28);
                         otherwise; throw(WdfError('PSet data type not recognised'));
                     end;
                     if (Type==99)||(Type==115)||(Type==105)||(Type)==119||(Type==114)||(Type==113) 
                         if F==0
                             A = this.AssertRead(1,Txt);
                         elseif F==128  
                             A = this.FillPSetArray(Txt); by = 4 + by*length(A);
                         end;
                     end;
                     by = 4 + by; TotS = TotS-by;
                     if xx==2 && Type~=107     %Discard the key values
                         LkUp1 = cell2mat(key(:,1)); LkUp2 = key(:,2);
                         f = LkUp1==k;
                         OP{n,1} = LkUp2(f); OP{n,2} = A;
                         n=n+1;
                     end;
                 end;
             end;
         end;
         
         %Function to read the flag and key each time. These occur for 
         %every PSet data type.
         function [Tq, Fq, kq] = ReadFlagKey(this)
             Tq = this.AssertRead(1,'uint8');
             Fq = this.AssertRead(1,'uint8');
             kq = this.AssertRead(1,'uint16');
         end;
         
         %Function to enable certain data types read in strings
         function [Str, by] = ReadPSetString(this)
             Len = this.AssertRead(1,'uint32');
             Str = this.ReadUtf8String(Len); by = Len+4;
         end;
         
         %Function to read in PSet data to an array. Reads the flags to
         %determine whether data is scalar or vector.
         function A = FillPSetArray(this,Tx)
              Leng = this.AssertRead(1,'uint32'); A = zeros(1,Leng);
              for ii = 1:Leng
                  A(ii) = this.AssertRead(1,Tx);
              end 
         end;
   
         %Function to determine whether or not the current MAP being 
         %parsed contains the correct parent PSet. Compares the parent ID
         %that we want with the actual PSet parent ID
         function  CorrectParent = FindParent(this,TargetID,TargetUID)
             PIndex = 0; CorrectParent = 0;
             TargetID = dec2hex(this.GetBlockID(TargetID));
             WantPID = uint64(str2double(strcat(num2str(TargetUID),num2str(TargetID)))); %In hex
             Strings = this.OPut(:,1);
             for x = 1:length(this.OPut(:,1))
                 St = Strings{x};
                 if strcmp(St,'parent') == 1; PIndex = x; end
             end;
             if PIndex ~=0
                 ActualPID = uint64(str2double(dec2hex(this.OPut{PIndex,2})));
                 CorrectParent = isequal(WantPID,ActualPID);
             end;
         end;
         
         %Helper function to read out the value of a cell, given the 
         %corresponding PSet Property. CellValue must be a string. If the 
         %CellValue field does not exist, then the function will return a 
         %value of -1.
         function  MapT = ReadCellValue(this,CellValue)
             Strings = this.OPut(:,1); PIndex = -1;
             for x = 1:length(this.OPut(:,1))
                 St = Strings{x};
                 if strcmp(St,CellValue) == 1; PIndex = x; end
             end;
             if PIndex == -1; MapT = -1;
             else MapT = (this.OPut{PIndex,2}); end;
         end
         
         %Function to determine the number of PSets that have a given ID in
         %the code
         function NP = getNumberofID(this,ID)
             count = 0; datasize = -1;       %Initialise variables.
             [Anchor, ~, ~] = LocateBlock(this,'YLST');
             this.AssertSeek(Anchor,'bof');
             while datasize ~= 0
                 count = count + 1;
                 [~, ~, datasize] = LocateBlock(this,ID,count);                
             end;
             NP = count - 1;
         end;
    end;
    
    %% Write a discriminant analysis map to the wdf file. 
    %Currently, this is written for discriminant analysis, but could
    %probably be generalised for other kinds of PSets.
    methods
        function WriteDiscriminantMap(this,Data) 
    
        end
    end
    
    %% Get the white light image from the file (GRL 05/10/2015)
    methods
        function [WL,x,y]=GetWLImage(this)
            this.AssertSeek(0, 'bof');                              %go to beginning of file
            [wl_offset, ~, dataSize] = LocateBlock(this, 'WHTL');   %find the white light location and get data size
            if (wl_offset == -1)                                    %if no data found throw error
                throw(WdfError('Cannot locate WL image.'));
            end;
            this.AssertSeek(wl_offset+16, 'bof');                   %go to the image position
            imgdata = fread(this.Handle, dataSize, '*uint8');       %read all the data as bytes
            %make use of imread function to avoid having to write a function to decompress jpgs
            fid = fopen('img_temp.jpg', 'w');                       %open a temporary jpg file with write access
            fwrite(fid, imgdata,'*uint8');                          %write the data as bytes
            fclose(fid);                                            %close the temporary file
            WL=imread('img_temp.jpg');                              %read the temporary file in using imread function
            % read the EXIF data and create x and y axes
            info=imfinfo('img_temp.jpg');                           
            FOV=info.UnknownTags(2).Value./info.UnknownTags(3).Value;
            x=linspace(info.UnknownTags(1).Value(1),info.UnknownTags(1).Value(1)+FOV(1),info.Width);
            y=linspace(info.UnknownTags(1).Value(2),info.UnknownTags(1).Value(2)+FOV(2),info.Height);
        end;        
    end;
    %% Read WMAP block
        methods
        function Returned = WMapblock(this)
            this.AssertSeek(0, 'bof');                              %go to beginning of file
            [wmap_offset, ~, ~] = LocateBlock(this, 'WMAP');        %find the WMAPBlock location and get data size
            if (wmap_offset == -1)                                  %if no data found throw error
                throw(WdfError('Cannot locate WMap section.'));
            end;
            x = this.AssertRead(1, 'uint32');
            this.AssertSeek(wmap_offset+24, 'bof');
            %unused = this.AssertRead(1, 'uint32');
            Returned.Location = this.AssertRead(3, 'float32');
            Returned.StepSize = this.AssertRead(3, 'float32');
            Returned.numPoints = this.AssertRead(3, 'uint32');
            Returned.linefocus_size = this.AssertRead(3, 'uint32');
            flag2 = binary2vector(x,8);
            Returned.flag = ConvertFlagToString(flag2);
        end;
    end;
    
    %% Query for next available block UID
    methods
        % Get the UID of the next available block given its ID
        function [UID] = GetNextAvailableBlockUID(this, BlockID)
            % Setup search parameters
            offset = int64(512);
            highestFoundUID = 0;
            targetID = this.GetBlockID(BlockID);
            
            % Iteratively walk over all blocks in the file
            while (feof(this.Handle) == 0)
                % Attemp to jump to next block
                this.AssertSeek(offset, 'bof');
                
                % Read block header, and throw an error if block header
                % fields could not be read
                blockID = fread(this.Handle, 1, 'uint32');
                uid = fread(this.Handle, 1, 'uint32');
                bSize = fread(this.Handle, 1, 'uint64');
                if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                    if (feof(this.Handle))
                        break;
                    else
                        throw(WdfError('Error whilst searching for next available block UID.'));
                    end;
                end;
                
                % Check if this block is of the target type ID, and if so,
                % record the UID if it is larger than any previously-found
                % block UID
                if (isequal(targetID, blockID))
                    if (uid > highestFoundUID)
                        highestFoundUID = uid;
                    end;
                end;
                
                % Increment the offset
                offset = offset + bSize;
            end;
            
            % We've searched through the whole file, so now we can safely
            % conclude that we have found the largest existing block UID,
            % and return the next available value.
            UID = int32(highestFoundUID + 1);
        end;
    end;
    
    
    
    
    %% Public static helper methods
    methods (Static, Access = public)
        % Converts a 4-character block-ID to its UINT32 value
        function [ID] = GetBlockID(String)
            ID = uint32(sum(String .* (256 .^ (0:3))));
        end;
        
        % Converts a 4-character block-ID and UINT32 UID to a 8-byte unique
        % block identifier
        function [FullID] = GetFullBlockID(BlockIDString, BlockUID)
            blockID = uint64(sum(BlockIDString .* (256 .^ (0:3))));
            FullID = int64(blockID + bitshift(uint64(BlockUID), 32));
        end;
    end;
end