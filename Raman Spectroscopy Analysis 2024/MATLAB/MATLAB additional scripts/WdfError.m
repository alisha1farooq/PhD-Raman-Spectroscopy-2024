function Exception = WdfError(Msg, varargin)

% WdfError  Initialise an MException with "WdfError" message-ID
%
% WE = WdfError(ERRMSG, V1, V2, ... VN) captures information about a
% specific error and stores it in WdfError object WE, which has the
% message-ID "Renishaw:SPD:WiRE:WdfError".
%
% See also: MException

% Copyright Notice:
% (c) 2012 - 2013 Renishaw plc. All rights reserved.

Exception = MException('Renishaw:SPD:WiRE:WdfError', Msg, varargin{:});