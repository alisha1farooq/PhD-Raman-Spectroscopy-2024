% WireLog  Send formatted messages to the WiREDiagnostics log
%
% The WireLog class allows messages to be sent from Matlab code to the
% WiREDiagnostics message log.  Messages in the log are associated with an
% Application Name and a Priority level.  The four available priority
% levels are Error, Warning, Function and Information; WireLog provides
% methods to send a message of each type to the log:
%    WireLog.Error(message, ...)
%    WireLog.Warning(message, ...)
%    WireLog.Function(message, ...)
%    WireLog.Info(message, ...)
% where <message> is a string and can use SPRINTF syntax.
%
% By default, messages sent from Matlab will have the application name
% "Matlab"; this can be changed by constructing a new WireLog object
% using either of the following syntaxes:
%    WireLog(name)
%    WireLog(name, message, ...)
% or by using the method
%    WireLog.SetAppName(name)
% NOTE: there is only ever one "name" for messages sent from Matlab; if
% messages need to be sent with different application names, it will be
% necessary to call SetAppName() before sending each message.
%
% With the exception of constructor-like syntax above, all methods are
% static, so an object-orientated syntax can optionally be used, e.g.:
%    log = WireLog(name)
%    log.Info(message, ...)
%
% WireLog requires the LogWireMsg MEX file.

classdef WireLog < handle
    % Constructor
    methods (Access = public)
        function this = WireLog(name, message, varargin)
            LogWireMsg(name);
            if (nargin < 2)
                LogWireMsg(4, sprintf('Logging started for [%s]', name));
            else
                LogWireMsg(4, sprintf(message, varargin{:}));
            end;
        end;
    end;
        
    methods (Static, Access = public)
        function SetAppName(name)
            LogWireMsg(name);
        end;
        
        function Error(message, varargin)
            if (nargin > 1)
                LogWireMsg(1, sprintf(message, varargin{:}));
            else
                LogWireMsg(1, message);
            end;
        end;
        
        function Warning(message, varargin)
            if (nargin > 1)
                LogWireMsg(2, sprintf(message, varargin{:}));
            else
                LogWireMsg(2, message);
            end;
        end;
        
        function Function(message, varargin)
            if (nargin > 1)
                LogWireMsg(4, sprintf(message, varargin{:}));
            else
                LogWireMsg(4, message);
            end;
        end;
        
        function Info(message, varargin)
            if (nargin > 1)
                LogWireMsg(8, sprintf(message, varargin{:}));
            else
                LogWireMsg(8, message);
            end;
        end;
    end
end

