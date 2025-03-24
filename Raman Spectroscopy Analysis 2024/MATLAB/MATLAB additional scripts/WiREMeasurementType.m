% WiREMeasurementType  Enumeration indicating WiRE measurement type
%
% There are three types of WiRE measurement:
%   Single    The file contains a single spectrum / dataset.
%   Series    The file contains a series of spectra / datasets, but these
%             spectra are related via something other than spatial position
%             (for example: time, temperature, pressure, etc).
%   Map       The file contains multiple spectra / datasets that relate to
%             different spatial positions in the sample.

% Copyright Notice:
% (c) 2012 - 2013 Renishaw plc. All rights reserved.

classdef WiREMeasurementType < uint32
    enumeration
        Unspecified        (0),
        Single             (1),
        Series             (2),
        Map                (3)
    end;
end