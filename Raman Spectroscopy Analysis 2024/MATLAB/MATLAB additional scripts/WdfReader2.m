function [data,X,Y,zAxis,Saturated,CosmicRay,Z]=WdfReader2(fileName)
this=[];
% wdf = WdfReader(fileName);
WdfReader();
data=GetSpectra(1, this.Count);
try
    X = GetOriginListValues(WiREDataType('SpatialX'), 1, this.Count);
    Y = GetOriginListValues(WiREDataType('SpatialY'), 1, this.Count);
    try
        Z = GetOriginListValues(WiREDataType('SpatialZ'), 1, this.Count);
    catch
        Z=[];
    end
    catch
    X=[];
    Y=[];
    Z=[];
%     warning('no xy coordinate detected. possibly not an image file.');
end
zAxis=GetXList();
try
[Saturated, CosmicRay] = GetOriginFlags(1,this.Count);
catch
    Saturated=[];
    CosmicRay=[];
%     warning('no cosmic ray or saturated data detected');
end
Close();


    function v=WiREScanType(bitFlags)
        
        MultitrackStitched = uint32(256);
        MultitrackDiscrete = uint32(512);
        LineFocusMapping   = uint32(1024);
        BasicType = WiREScanBasicType('Unspecified');
        IsMultitrackStitched = false;
        IsMultitrackDiscrete = false;
        IsLineFocusMapping = false;
        
        bitFlags = uint32(bitFlags);
        v.BasicType = WiREScanBasicType(bitand(uint32(255), bitFlags));
        v.IsMultitrackStitched = [];%bitand(WiREScanType('MultitrackStitched'), bitFlags) ~= 0;
        v.IsMultitrackDiscrete = [];%bitand(WiREScanType('MultitrackDiscrete'), bitFlags) ~= 0;
        v.IsLineFocusMapping = [];%bitand(WiREScanType('LineFocusMapping'), bitFlags) ~= 0;
    end;
    
    function val=WiREScanBasicType(str)
        enums={
            'Unspecified';
            'Static';
            'Continuous';
            'StepRepeat';
            'FilterScan';
            'FilterImage';
            'StreamLine';
            'StreamLineHR';
            'PointDetector';
            };
        if ischar(str)
            val=find(strcmp(str,enums));
            if ~isempty(val)
                val=enums{val};
            end
        elseif length(str)>1
            val=find(strcmp(char(str),enums));
            if ~isempty(val)
                val=enums{val};
            end
        else
            val=enums{str+1};
        end
    end;
    
    function val= WiREDataUnit(str)
        enums={
            'Arbitrary';
            'RamanShift';
            'Wavenumber';
            'Nanometre';
            'ElectronVolt';
            'Micron';
            'Counts';
            'Electrons';
            'Millimetres';
            'Metres';
            'Kelvin';
            'Pascal';
            'Seconds';
            'Milliseconds';
            'Hours';
            'Days';
            'Pixels';
            'Intensity';
            'RelativeIntensity';
            'Degrees';
            'Radians';
            'Celsius';
            'Fahrenheit';
            'KelvinPerMinute';
            'FileTime';};
        if ischar(str)
            val=find(strcmp(str,enums));
            if ~isempty(val)
                val=enums{val};
            end
        else
            val=enums{str+1};
        end
    end;
    
    
    function val=WiREDataType(str)
        enums={
            'Arbitrary';
            'Frequency';
            'Intensity';
            'SpatialX'; % X axis position
            'SpatialY'; % Y axis position
            'SpatialZ'; % Z axis (vertical) position
            'SpatialR'; % rotary stage R axis position
            'SpatialTheta'; % rotary stage theta angle
            'SpatialPhi'; % rotary stage phi angle
            'Temperature';
            'Pressure';
            'Time';
            'Derived'; % derivative type
            'Polarization';
            'FocusTrack'; % focus track Z position
            'RampRate'; % temperature ramp rate
            'Checksum';
            'Flags'; % bit flags
            'ElapsedTime'; % elapsed time interval
            };
        if ischar(str)
            val=find(strcmp(str,enums));
            if ~isempty(val)
                val=enums{val};
            end
        else
            val=enums{str+1};
        end
    end

    function val=WiREMeasurementType(str)
        enums={
            'Unspecified';
            'Single';
            'Series';
            'Map';
            };
        if ischar(str)
            val=find(strcmp(str,enums));
            if ~isempty(val)
                val=enums{val};
            end
        else
            val=enums{str+1};
        end
        
    end

%% Constructor, destructor and Close
% Constructor; opens an existing WDF file with read-only access
% using the protected constructor.
    function WdfReader()
        ProtectedConstructor('rb');
    end;
    
    % Closes the file.
    function Close()
        if (this.Handle ~= -1)
            fclose(this.Handle);
            this.Handle = -1;
        end;
    end;
    
    
    % Opens a WDF file with the specified access mode, reads the WDF
    % header, and validates the DATA, XLST, YLST and ORGN blocks.
    function ProtectedConstructor(fileAccess)
        % Open the file
        this.Handle = fopen(fileName, fileAccess);
        if (this.Handle == -1)
            throw(error('Cannot open file "%s".', fileName));
        end;
        
        % Parse the header, locate the data block and check its size is
        % consistent with expected value.
        ReadWdfHeader();
        [this.DataOffset, meh, dataSize] = LocateBlock('DATA');
        if (this.DataOffset == -1)
            throw(error('Cannot locate Spectral Data block.'));
        end;
        if (dataSize < (16 + (4 * this.PointsPerSpectrum * this.Capacity)))
            throw(error('Spectral Data block size inconsistent with WDF header.'));
        end;
        
        % Locate the X-list and Y-list blocks, and validate their sizes.
        [this.XListOffset, meh, dataSize] = LocateBlock('XLST');
        if (this.XListOffset == -1)
            throw(error('Cannot locate X-list block.'));
        end;
        if (dataSize < (24 + (4 * this.XListLength)))
            throw(error('X-list block size inconsistent with WDF header.'));
        end;
        
        [this.YListOffset, meh, dataSize] = LocateBlock('YLST');
        if (this.YListOffset == -1)
            throw(error('Cannot locate Y-list block.'));
        end;
        if (dataSize < (24 + (4 * this.YListLength)))
            throw(error('Y-list block size inconsistent with WDF header.'));
        end;
        
        % Read the X and Y list types / units
        this.XListType = GetXListType();
        this.YListType = GetYListType();
        this.XListUnits = GetXListUnits();
        this.YListUnits = GetYListUnits();
        
        % If the Data Origin List count is non-zero, attempt to locate
        % then size-check the ORGN block.  Finally, read in the data
        % origin list info.
        if (this.DataOriginCount ~= 0)
            [this.OriginsOffset, meh, dataSize] = LocateBlock('ORGN');
            if (this.OriginsOffset == -1)
                throw(error('Cannot locate Data Origin List block.'));
            end;
            if (dataSize < (20 + (this.DataOriginCount * (24 + (this.Capacity * 8)))))
                throw(error('Data Origin List block size inconsistent with WDF header.'));
            end;
        end;
        ReadOriginListInfo();
    end;
    
    % Reads the WDF header block, extracting key header fields and
    % storing them in class properties.
    function ReadWdfHeader()
        % Confirm that the signature, version and size fields
        % contain the expected values.
        AssertSeek(0, 'bof');
        blockID = AssertRead(1, 'uint32');
        blockUID = AssertRead(1, 'uint32');
        blockLength = AssertRead(1, 'uint64');
        if (~isequal(blockID, GetBlockID('WDF1')) || ...
                ~(isequal(blockUID, 0) || isequal(blockUID, 1)) || ...
                ~isequal(blockLength, 512))
            throw(error('File does not use a recognised WDF format / version.'));
        end;
        
        % Read key fields from the main file header
        AssertSeek(60, 'bof');
        this.PointsPerSpectrum = AssertRead(1, 'uint32');
        this.Capacity = AssertRead(1, 'uint64');
        this.Count = AssertRead(1, 'uint64');
        this.AccumulationCount = AssertRead(1, 'uint32');
        this.YListLength = AssertRead(1, 'uint32');
        this.XListLength = AssertRead(1, 'uint32');
        this.DataOriginCount = AssertRead(1, 'uint32');
        this.ApplicationName = ReadUtf8String(24);
        this.ApplicationVersion = AssertRead([1 4], 'uint16');
        this.ScanType = WiREScanType(AssertRead(1, '*uint32'));
        this.MeasurementType = WiREMeasurementType(AssertRead(1, '*uint32'));
        AssertSeek(152, 'bof');
        this.SpectralUnits = WiREDataUnit(AssertRead(1, '*uint32'));
        this.LaserWavenumber = AssertRead(1, 'float32');
        AssertSeek(208, 'bof');
        this.Username = ReadUtf8String(32);
        this.Title = ReadUtf8String(160);
    end;
    
    % Searches for a specific block within the WDF file, by ID and
    % optionally UID.  Returns the location (offset), UID and length
    % (in bytes) of the first matching block found, or -1 if no
    % matching block was located.
    function [Location, BlockUID, BlockSize] = LocateBlock(TargetID, TargetUID)
        % Set default return values if block is not found
        Location = -1;
        BlockUID = [];
        BlockSize = uint64(0);
        
        % Initialise search
        found = false;
        offset = double(512);
        TargetID = GetBlockID(TargetID);
        
        % Iteratively walk over all data-blocks in the file
        while ((~found) && (feof(this.Handle) == 0))
            % Attemp to jump to next block
            if (fseek(this.Handle, offset, 'bof') ~= 0)
                return;
            end;
            
            % Read block header, and abort search if header fields not
            % read successfully
            blockID = fread(this.Handle, 1, '*uint32');
            uid = fread(this.Handle, 1, '*uint32');
            bSize = fread(this.Handle, 1, 'uint64=>double');
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
    function s = ReadUtf8String(Length)
        s = AssertRead([1 Length], 'uint8=>char');
        s = deblank(char(unicode2native(s, 'UTF-8')));
    end;
    
    % Seeks to the requested position in the WDF file, raising an error
    % if the seek operation is unsuccessful.
    function AssertSeek(Position, Origin)
        if (fseek(this.Handle, Position, Origin) ~= 0)
            (error('Failed to seek to requested position within file.'));
        end;
    end;
    
    % Reads the requested data from the WDF file, raising an error if
    % the actual number of elements read is fewer than requested.
    function [Data] = AssertRead(RDims, Precision)
        [Data, readCount] = fread(this.Handle, RDims, Precision);
        if (readCount ~= prod(RDims))
            (error('Failed to read requested data from file.'));
        end;
    end;
    
    % Helper function used to read data origin list values.
    function [Data] = ReadOriginListData(ListType, IndexStart, IndexEnd, Precision)
        % Validate inputs.
        if (nargin ~= 4)
            (error('ListType, IndexStart and IndexEnd must be specified.'));
        end;
        if ((IndexStart < 1) || (IndexStart > this.Count))
            (error('IndexStart is out-of-range.'));
        end;
        if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
            (error('IndexEnd is out-of-range.'));
        end;
        listIndex = find(cellfun(@(x) isequal(x, ListType), this.OriginListInfo(:, 2)), 1);
        if (isempty(listIndex))
            (error('WDF file does not contain a Data Origin List with the specified type.'));
        end;
        
        % Read and return the data.
        count = (IndexEnd - IndexStart) + 1;
        listSize = 24 + (8 * this.Capacity);
        offset = this.OriginsOffset + 20 + ((listIndex - 1) * listSize) + 24 + ((IndexStart - 1) * 8);
        AssertSeek(offset, 'bof');
        Data = AssertRead(count, Precision);
    end;
    
    
    %% Non-public methods that access key WDF data
    % Get the X-list type
    function [Result] = GetXListType()
        AssertSeek(this.XListOffset + 16, 'bof');
        Result = WiREDataType(AssertRead(1, '*uint32'));
    end;
    
    % Get the X-list units
    function [Result] = GetXListUnits()
        AssertSeek(this.XListOffset + 20, 'bof');
        Result = WiREDataUnit(AssertRead(1, '*uint32'));
    end;
    
    % Get the Y-list type
    function [Result] = GetYListType()
        AssertSeek(this.YListOffset + 16, 'bof');
        Result = WiREDataType(AssertRead(1, '*uint32'));
    end;
    
    % Get the Y-list units
    function [Result] = GetYListUnits()
        AssertSeek(this.YListOffset + 20, 'bof');
        Result = WiREDataUnit(AssertRead(1, '*uint32'));
    end;
    
    % Creates and stores a cell-array containing information about the
    % Data Origin Lists in the WDF file (one row per data origin list).
    function ReadOriginListInfo()
        ListInfo = cell(this.DataOriginCount, 4);
        listSize = 24 + (8 * this.Capacity);
        for n = 1:this.DataOriginCount
            AssertSeek(this.OriginsOffset + 20 + ((n - 1) * listSize), 'bof');
            listInfo = AssertRead([1 2], '*uint32');
            ListInfo{n, 1} = bitget(listInfo(1), 32) ~= 0;
            ListInfo{n, 2} = WiREDataType(bitset(listInfo(1), 32, 0));
            ListInfo{n, 3} = WiREDataUnit(listInfo(2));
            ListInfo{n, 4} = ReadUtf8String(16);
        end;
        this.OriginListInfo = ListInfo;
    end;
    
    
    %% Access to spectra in the file, plus X- and Y-list data
    % Read the X-list
    function [XList] = GetXList()
        AssertSeek(this.XListOffset + 24, 'bof');
        XList = double(AssertRead([1 this.XListLength], 'single'));
    end;
    
    % Read the Y-list
    function [YList] = GetYList()
        AssertSeek(this.YListOffset + 24, 'bof');
        YList = double(AssertRead([1 this.YListLength], 'single'));
    end;
    
    % Import a range of spectra, specified by start & end indices.
    function [Data] = GetSpectra(IndexStart, IndexEnd)
        % Validate inputs.
        if (nargin ~= 2)
            throw(error('IndexStart and IndexEnd must be specified.'));
        end;
        if ((IndexStart < 1) || (IndexStart > this.Count))
            throw(error('IndexStart is out-of-range.'));
        end;
        if ((IndexEnd < IndexStart) || (IndexEnd > this.Count))
            throw(error('IndexEnd is out-of-range.'));
        end;
        
        % Read the data, and confirm the expected quantity was read
        offset = (this.DataOffset + 16) + ((IndexStart - 1) * 4 * this.PointsPerSpectrum);
        nRowsToRead = (IndexEnd - IndexStart) + 1;
        AssertSeek(offset, 'bof');
        Data = AssertRead([this.PointsPerSpectrum,nRowsToRead], 'single=>double');
        Data=Data';
        % Re-shape the data and convert to double-type
        %             Data = double(reshape(Data, this.PointsPerSpectrum, nRowsToRead)');
    end;
    
    % Initiates chunk-wise reading of the spectral data
    function StartChunkwiseReading()
        AssertSeek(this.DataOffset + 16, 'bof');
        this.m_NextSpectrumIndex = 1;
    end;
    
    % Reads the next chunk of spectra from the file; NumberOfSpectra is
    % optional.
    function [Data, Indices] = GetNextDataChunk(NumberOfSpectra)
        % Determine how many rows (spectra) to read into this chunk.
        if (nargin == 2)
            if (NumberOfSpectra < 1)
                throw(error('NumberOfSpectra must be greater-than-or-equal-to 1.'));
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
            throw(error('All spectra have already been read in previous chunks.'));
        end;
        
        % Read the data.
        offset = (this.DataOffset + 16) + ((this.m_NextSpectrumIndex - 1) * this.PointsPerSpectrum * 4);
        AssertSeek(offset, 'bof');
        Data = AssertRead(nRowsToRead * this.PointsPerSpectrum, 'single');
        
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
    function [MoreChunks] = AreMoreChunks()
        MoreChunks = (this.m_NextSpectrumIndex <= this.Count);
    end;
    
    % Returns a value in the range [0-1] that indicates the fraction of
    % spectra already imported via chunk-wise reading.
    function [Fraction] = GetChunkwiseProgress()
        Fraction = (this.m_NextSpectrumIndex - 1) / this.Count;
    end;
    
    
    %% Access to Data Origin Lists
    % Returns a cell-array containing information about the Data Origin
    % Lists in the WDF file (one row per data origin list).
    function [ListInfo] = GetOriginListInfo()
        ListInfo = this.OriginListInfo;
    end;
    
    % Gets a cell-array containing a list of primary (non-alternate)
    % data origin list data types.
    function [PrimaryOriginLists] = GetPrimaryOriginLists()
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
    function [Data] = GetOriginListValues(ListType, IndexStart, IndexEnd)
        Data = ReadOriginListData(ListType, IndexStart, IndexEnd, 'double');
    end;
    
    % Gets a range of values from a Data Origin List, specified by list
    % data type, treating the binary data as 64-bit integer values.
    function [Data] = GetOriginListValuesInt(ListType, IndexStart, IndexEnd)
        Data = ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
    end;
    
    %Finds out which spectra have a saturated spectrum or cosmic ray
    %removal flag associated with them.
    function [Saturated CosmicRay] = GetOriginFlags(IndexStart,IndexEnd)
        ListType = WiREDataType('Flags');
        Data = ReadOriginListData(ListType, IndexStart, IndexEnd, 'int64=>double');
        S = 7*ones(length(Data),1); Data = bitand(Data,S);
        A = find(Data==1); B = find(Data==3); C = find(Data==5); D = find(Data==7);
        Saturated = sort([(A);(B);(C);(D)]); clear A B C D
        A = find(Data==4); B = find(Data==5); C = find(Data==6); D = find(Data==7);
        CosmicRay = sort([(A);(B);(C);(D)]); clear A B C D
    end;
    
    %% Query for next available block UID
    
    % Get the UID of the next available block given its ID
    function [UID] = GetNextAvailableBlockUID(BlockID)
        % Setup search parameters
        offset = int64(512);
        highestFoundUID = 0;
        targetID = GetBlockID(BlockID);
        
        % Iteratively walk over all blocks in the file
        while (feof(this.Handle) == 0)
            % Attemp to jump to next block
            AssertSeek(offset, 'bof');
            
            % Read block header, and throw an error if block header
            % fields could not be read
            blockID = fread(this.Handle, 1, 'uint32');
            uid = fread(this.Handle, 1, 'uint32');
            bSize = fread(this.Handle, 1, 'uint64');
            if (isempty(blockID) || isempty(uid) || (numel(bSize) ~= 1) || (bSize < 16))
                if (feof(this.Handle))
                    break;
                else
                    throw(error('Error whilst searching for next available block UID.'));
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
    
    
    %% Public static helper methods
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
    
end