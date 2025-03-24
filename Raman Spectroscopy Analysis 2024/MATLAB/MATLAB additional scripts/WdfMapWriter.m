% WdfMapWriter  Helper class for writing WDF map-analysis results to file
%
% writer = WdfMapWriter(filename, blockID, blockUID) creates a new
% WdfMapWriter object, which is designed to simplify writing results from a
% mapping analysis of WDF data to file, especially when used with
% WdfReader.  The resulting file will contain an Analysis / history-item
% block with the specified block ID and UID, and related MAP blocks can be
% appended.  This file can then be directly appended onto an existing WDF
% file in order to add the results to the WDF file.
%
% writer.FinaliseBlock() should be called after all data has been written
% to the block, and records the total block size.  This method should be
% called before starting to write additional MAP blocks.
%
% writer.Close() closes the file.
%
% writer.WritePSet(cellArray) stores the cell-array as a WDF property-set
% in the file.  The first column contains property names, and the second
% column contains property values.  Nested cell-arrays are supported.
% Values can be of single, double, or any signed-integer type, or a string
% (single-line character array).
%
% Simple arrays can be stored in column-major ordering using various
% data-types:
%    writer.WriteSingleArray(array)
%    writer.WriteDoubleArray(array)
%    writer.WriteInt32Array(array)
%    writer.WriteInt64Array(array)
% 
% To allow chunk-wise processing of large files when producing MAPs, the
% following routines provide interleaved chunk-wise MAP block writing.  The
% "MapProperties" input should be a cell array containing one p-set per map
% to store in each MAP block.
%    writer.BeginInterleavedMapWriting(FirstMapUID, MapSize, MapProperties)
%    writer.WriteInterleavedMapChunk(IndexOfFirstRow, DataChunk)
%    writer.EndInterleavedMapWriting()

classdef (Sealed) WdfMapWriter < handle
    %% Private fields
    properties (Access = private)
        m_PSetKeyIndex = 33000;
        m_BlockOffset = 0;
        m_MapPsetOffsets = [];
        m_MapDataOffsets = [];
        m_MapBlocksEndOffset = [];
        Handle = -1;
    end;
    
    %% Constructor, destructor, Close and FinaliseBlock
    methods
        % Constructor: creates a new object, writing data to a
        % newly-created file, and write a WDF block header.
        function this = WdfMapWriter(fileName, BlockID, UID)
            % Open the file, checking it does not already exist
            if (exist(fileName, 'file') ~= 0)
                throw(WdfError('File already exists "%s"', fileName));
            end;
            this.Handle = fopen(fileName, 'w+b');
            if (this.Handle == -1)
                throw(WdfError('Cannot open file "%s"', fileName));
            end;
            
            % Write the block header (with a temporary block-size of 0)
            fwrite(this.Handle, WdfReader.GetBlockID(BlockID), 'uint32');
            fwrite(this.Handle, uint32(UID), 'uint32');
            fwrite(this.Handle, 0, 'uint64');
            if (ftell(this.Handle) ~= (this.m_BlockOffset + 16))
                throw(WdfError('Failed to write result block header'));
            end;
        end;
        
        % Finalises the block, by writing the block size to the block
        % header.
        function FinaliseBlock(this)
            if (fseek(this.Handle, 0, 'eof') ~= 0)
                this.ThrowSeekError();
            end;
            blockSize = ftell(this.Handle) - this.m_BlockOffset;
            if (fseek(this.Handle, 8, 'bof') ~= 0)
                this.ThrowSeekError();
            end;
            if (fwrite(this.Handle, blockSize, 'uint64') ~= 1)
                this.ThrowWriteFailedError();
            end;
            if (fseek(this.Handle, blockSize + this.m_BlockOffset, 'bof') ~= 0)
                this.ThrowSeekError();
            end;
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
    
    %% Property-set writing (public methods)
    methods (Access = public)
        % Exports the contents of a Matlab cell-array to a WDF property-set
        % at the current position in the WDF block.  Nested structures are
        % supported.  The returned PSetSize includes the 8 bytes required
        % for the 'PSET' header and p-set length fields.
        function [PSetSize] = WritePSet(this, Properties)
            % Validate the input cell-array.
            if (~iscell(Properties) || (size(Properties, 2) ~= 2))
                throw(GenericError('Invalid Properties format (expected a two-column cell-array)'));
            end;
            
            % Write a p-set marker, then call the private p-set writing
            % function.
            if (fwrite(this.Handle, 'PSET', 'char') ~= 4)
                throw(WdfError('Failed to write p-set marker'));
            end;
            PSetSize = this.WritePSetPrivate(Properties);
            
            % The total number of bytes occupied by the pset includes the
            % four bytes used by the "PSET" marker, and the p-set length
            % field, as well as the p-set data.
            PSetSize = PSetSize + 8;
        end;
        
        % Re-exports the contents of a Matlab cell-array of cell-arrays to
        % the Map-block p-sets regions.  This is intended to be called
        % after the map data range is known, in order to update the
        % 'dataRange' p-set field.  The p-set size before and after this
        % call must be identical.
        function UpdateMapPSets(this, Properties)
            % Validate the input cell-array and the object state.
            if (~iscell(Properties))
                throw(GenericError('Invalid Properties format (expected a cell-array).'));
            elseif (length(Properties) ~= length(this.m_MapPsetOffsets))
                throw(GenericError('Number of p-sets does not match number of map blocks previously written.'));
            else
                for n = 1:length(Properties)
                    if (~iscell(Properties{n}) || (size(Properties{n}, 2) ~= 2))
                        throw(GenericError('Invalid p-set format (expected a two-column cell-array)'));
                    end;
                end;
            end;
            
            % Process the p-set for each map in turn.
            for n = 1:length(Properties)
                if (fseek(this.Handle, this.m_MapPsetOffsets(n), 'bof') ~= 0)
                    this.ThrowSeekError();
                end;
                if (~isequal('PSET', fread(this.Handle, [1 4], '*char')))
                    throw(WdfError('PSET has moved / been overwritten since map was created.'));
                end;
                [oldSize, nRead] = fread(this.Handle, 1, 'uint32');
                if (nRead ~= 1)
                    throw(WdfError('Failed to read PSET size from Map block.'));
                end;
                if (fseek(this.Handle, -4, 'cof') ~= 0)
                    this.ThrowSeekError();
                end;
                newSize = this.WritePSetPrivate(Properties{n});
                if (oldSize ~= newSize)
                    throw(GenericError('P-set size has changed; map data may now be corrupted.'));
                end;
            end;
        end;
    end;
    
    %% Property-set writing (private helper methods)
    methods (Access = private)
        % Exports the contents of a Matlab cell-array to a WDF property-set
        % at the current position in the WDF block.  This private routine
        % does not mark the start of the p-set with the four-byte "PSET"
        % header, so can be used for nested p-sets or called from the
        % public function above.  The returned PSetSize is the length of
        % the p-set data, EXCLUDING the 4-byte length field.
        function [PSetSize] = WritePSetPrivate(this, Properties)
            % Note the offset into file of the p-set size header, but write
            % a dummy value for now.
            psetSizeOffset = ftell(this.Handle);
            if (fwrite(this.Handle, 0, 'uint32') ~= 1)
                throw(WdfError('Failed to write p-set header'));
            end;
            
            % Write each row of the input cell-array to a separate
            % property in turn.
            for n = 1:size(Properties, 1)
                % To check data is written properly, we will compare the
                % expected and actual number of bytes written.
                expectedSize = 0;
                fieldOffset = ftell(this.Handle);
                
                % Write the property name.
                name = Properties{n, 1};
                fwrite(this.Handle, 'k', 'char');
                fwrite(this.Handle, 0, 'uint8');
                fwrite(this.Handle, this.m_PSetKeyIndex, 'uint16');
                utf8str = native2unicode(name, 'UTF-8');
                fwrite(this.Handle, numel(utf8str), 'uint32');
                fwrite(this.Handle, utf8str, 'uint8');
                expectedSize = expectedSize + 8 + numel(utf8str);
                
                % Write the property value.
                expectedSize = expectedSize + this.WritePSetValue(this.m_PSetKeyIndex, name, Properties{n, 2});
                
                % Finally, check that the number of bytes actually written
                % matches the expected number of bytes written.
                if ((ftell(this.Handle) - fieldOffset) ~= expectedSize)
                    throw(WdfError('Internal error whilst writing p-set data'));
                end;
                
                % Increment the property-set custom key index number
                this.m_PSetKeyIndex = this.m_PSetKeyIndex + 1;
            end
            
            % Having successfully written all fields of the p-set, we now
            % calculate how many bytes the pset occupies, then write this
            % to the pset "header" area.  The size does NOT include the
            % size itself.
            PSetSize = ftell(this.Handle) - (psetSizeOffset + 4);
            if (fseek(this.Handle, psetSizeOffset, 'bof') ~= 0)
                this.ThrowSeekError();
            end;
            if (fwrite(this.Handle, PSetSize, 'uint32') ~= 1)
                this.Throw(WdfError('Failed to record p-set size'));
            end;
            if (fseek(this.Handle, psetSizeOffset + 4 + PSetSize, 'bof') ~= 0)
                this.ThrowSeekError();
            end;
        end;
        
        % Writes a single p-set value to the WDF block, the property name
        % having already been stored.  The returned FieldSize is the length
        % of the serialized field, including the key name.
        function FieldSize = WritePSetValue(this, Index, Name, Value)
            % Write the property value, according to the field data type .           
            if (iscell(Value))
                % For nested cell-arrays, we create a nested p-set.
                if (size(Value, 2) ~= 2)
                    throw(WdfError('Nested cell-arrays must have exactly two columns'));
                end;
                fwrite(this.Handle, 'p', 'char');
                fwrite(this.Handle, 0, 'uint8');
                fwrite(this.Handle, Index, 'uint16');
                FieldSize = 8 + this.WritePSetPrivate(Value);
            elseif (ischar(Value) && (ndims(Value) == 2) && (size(Value, 1) == 1))
                % For strings (single-line char arrays), write as a
                % UTF8-encoded string.
                fwrite(this.Handle, 'u', 'char');
                fwrite(this.Handle, 0, 'uint8');
                fwrite(this.Handle, Index, 'uint16');
                utf8str = native2unicode(Value, 'UTF-8');
                fwrite(this.Handle, numel(utf8str), 'uint32');
                fwrite(this.Handle, utf8str, 'uint8');
                FieldSize = 8 + (numel(utf8str));
            elseif (islogical(Value))
                % Boolean values.
                fwrite(this.Handle, '?', 'char');
                if (numel(Value) > 1)
                    fwrite(this.Handle, 128, 'uint8');
                    fwrite(this.Handle, Index, 'uint16');
                    fwrite(this.Handle, numel(Value), 'uint32');
                    FieldSize = 8 + numel(Value);
                else
                    fwrite(this.Handle, 0, 'uint8');
                    fwrite(this.Handle, Index, 'uint16');
                    FieldSize = 5;
                end;
                fwrite(this.Handle, int8(Value), 'int8');
            elseif (isnumeric(Value) && isreal(Value))
                % For non-complex numeric data, write using the input
                % data-type.  Unsigned integer values are currently
                % unsupported.
                className = class(Value);
                if (isequal(Name, 'Time'))
                    className = 'filetime';
                end;
                bytesPerElement = 0;
                switch (className)
                    case 'int8'
                        typeCode = 'c';
                        bytesPerElement = 1;
                    case 'int16'
                        typeCode = 's';
                        bytesPerElement = 2;
                    case 'int32'
                        typeCode = 'i';
                        bytesPerElement = 4;
                    case 'int64'
                        typeCode = 'w';
                        bytesPerElement = 8;
                    case 'single'
                        typeCode = 'r';
                        bytesPerElement = 4;
                    case 'double'
                        typeCode = 'q';
                        bytesPerElement = 8;
                    case 'filetime',
                        typeCode = 't';
                        bytesPerElement = 8;
                        className = 'uint64';
                    otherwise
                        this.ThrowUnsupportedPsetTypeError(Name);
                end;
                
                % Type and sizeof() has been determined; now write to file.
                fwrite(this.Handle, typeCode, 'char');
                if (numel(Value) > 1)
                    fwrite(this.Handle, 128, 'uint8');
                    fwrite(this.Handle, Index, 'uint16');
                    fwrite(this.Handle, numel(Value), 'uint32');
                    FieldSize = 8 + (numel(Value) * bytesPerElement);
                else
                    fwrite(this.Handle, 0, 'uint8');
                    fwrite(this.Handle, Index, 'uint16');
                    FieldSize = 4 + bytesPerElement;
                end;
                fwrite(this.Handle, Value, className);
            else
                % Anything else is unsupported.
                this.ThrowUnsupportedPsetTypeError(Name);
            end;
        end;
    end;
    
    %% Writing simple matrices
    methods
        % Writes a matrix, in column-major order, to the file using single
        % precision floating-point.
        function WriteSingleArray(this, Array)
            if (fwrite(this.Handle, Array, 'single') ~= numel(Array))
                this.ThrowWriteFailedError();
            end;
        end;
        
        % Writes a matrix, in column-major order, to the file using double
        % precision floating-point.
        function WriteDoubleArray(this, Array)
            if (fwrite(this.Handle, Array, 'double') ~= numel(Array))
                this.ThrowWriteFailedError();
            end;
        end;
        
        % Writes a matrix, in column-major order, to the file as INT32s.
        function WriteInt32Array(this, Array)
            if (fwrite(this.Handle, Array, 'int32') ~= numel(Array))
                this.ThrowWriteFailedError();
            end;
        end;
        
        % Writes a matrix, in column-major order, to the file as INT64s.
        function WriteInt64Array(this, Array)
            if (fwrite(this.Handle, Array, 'int64') ~= numel(Array))
                this.ThrowWriteFailedError();
            end;
        end;
    end;
    
    %% Chunk-wise interleaved MAP block writing
    methods
        % Prepares for chunk-wise interleaved writing of several MAP
        % blocks, by creating MAP blocks, initialising property-sets and
        % reserving space for each map.  The map data data is initialised
        % to zeros.
        function BeginInterleavedMapWriting(this, FirstMapUID, MapSize, MapProperties)
            % How many maps should we create?
            nMaps = numel(MapProperties);
            
            % Pre-allocate a vector of map block data offsets
            this.m_MapPsetOffsets = nan(1, nMaps);
            this.m_MapDataOffsets = nan(1, nMaps);
            
            % Create an empty map block for each map
            for n = 1:nMaps
                % Write the block header
                blockStartOffset = ftell(this.Handle);
                fwrite(this.Handle, WdfReader.GetBlockID('MAP '), 'uint32');
                fwrite(this.Handle, FirstMapUID + (n - 1), 'uint32');
                fwrite(this.Handle, 0, 'uint64');
                if ((ftell(this.Handle) - blockStartOffset) ~= 16)
                    this.ThrowWriteFailedError();
                end;
                
                % Write the property-set for this map block
                this.m_MapPsetOffsets(n) = ftell(this.Handle);
                this.WritePSet(MapProperties{n});
                
                % Write the number of elements in the map block as a UINT64
                if (fwrite(this.Handle, MapSize, 'uint64') ~= 1)
                    this.ThrowWriteFailedError();
                end;
                
                % Record the offset to the start of the map block data
                this.m_MapDataOffsets(n) = ftell(this.Handle);
                
                % Reserve space in the file for the map data
                nToReserve = MapSize;
                while (nToReserve > 0)
                    nElementsToWrite = min(nToReserve, 2^18);
                    if (fwrite(this.Handle, zeros(nElementsToWrite, 1), 'single') ~= nElementsToWrite)
                        this.ThrowWriteFailedError();
                    end;
                    nToReserve = nToReserve - nElementsToWrite;
                end;
            
                % At this point we have finished reserving space for this
                % map, so we now go back and write the actual map block
                % size in the map block header, before returning the file
                % pointer to the end of the map block.
                blockEndOffset = ftell(this.Handle);
                if (fseek(this.Handle, blockStartOffset + 8, 'bof') ~= 0)
                    this.ThrowSeekError();
                end;
                if (fwrite(this.Handle, blockEndOffset - blockStartOffset, 'uint64') ~= 1)
                    this.ThrowWriteFailedError();
                end;
                if (fseek(this.Handle, blockEndOffset, 'bof') ~= 0)
                    this.ThrowSeekError();
                end;
            end;
            
            % Record the offset of the end of the last map (the place at
            % which the next block could begin)
            this.m_MapBlocksEndOffset = ftell(this.Handle);
        end;
        
        % Writes a chunk of map data, corresponding to a limited range of
        % rows from a larger array, to the file, the space having been
        % previously reserved using BeginInterleavedMapWriting().
        function WriteInterleavedMapChunk(this, FirstRow, Data)
            rowOffset = 4 * (FirstRow - 1);
            [nRows, nCols] = size(Data);
            for n = 1:nCols
                if (fseek(this.Handle, this.m_MapDataOffsets(n) + rowOffset, 'bof') ~= 0)
                    this.ThrowSeekError();
                end;
                if (fwrite(this.Handle, Data(:, n), 'single') ~= nRows)
                    this.ThrowWriteFailedError();
                end;
            end;
        end;
        
        % Finishes the chunk-wise writing of a very large array, simply by
        % moving the file stream position to the end of the reserved
        % space.
        function EndInterleavedMapWriting(this)
            this.m_MapDataOffsets = [];
            if (fseek(this.Handle, this.m_MapBlocksEndOffset, 'bof') ~= 0)
                this.ThrowSeekError();
            end;
            this.m_MapBlocksEndOffset = [];
        end;
    end;
    
    %% Private static helper methods
    methods (Static, Access = private)
        % Throws a generic error for a failed FSEEK operation.
        function ThrowSeekError()
            throwAsCaller(WdfError('Failed to seek to requested position within file'));
        end;
        
        % Throws an error indicating that a p-set field used an invalid
        % datatype.
        function ThrowUnsupportedPsetTypeError(FieldName)
            throwAsCaller(WdfError('P-set field "%s" uses an unsupported data type', FieldName));
        end;
        
        % Throws a generic error for a failed FWRITE operation.
        function ThrowWriteFailedError()
            throwAsCaller(WdfError('Failed to write data to file'));
        end;
    end;
end