Title = 'Input Dialog for PCA';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'on';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 4; % Horizontal dimension in fields
Prompt = {};
Formats = {};
DefAns = struct([]);

% Prompt(1,:) = {'Choose Currency','MoneyUnit',[]};
% Formats(1,1).type = 'list';
% Formats(1,1).format = 'text';
% Formats(1,1).style = 'radiobutton';
% Formats(1,1).items = {'U.S. Dollar' 'Euro';'Japanese Yen' ''};
% DefAns.MoneyUnit = 'Japanese Yen';%3; % yen

Prompt(1,:) = {['This Dialogue will allow the user to enter the PCA parameters ' ... 
   'and demonstrates how Formats input can be used to layout these controls.'],[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
Formats(1,1).span = [1 2]; % item is 1 field x 4 fields
n=1;

n=n+1;
Prompt(end+1,:) = {'User''s Name', 'Name',[]};
Formats(n,1).type = 'edit';
Formats(n,1).format = 'text';
Formats(n,1).size = 200; % automatically assign the height
DefAns(1).Name = 'Martin Isabelle';

n=n+1;
Prompt(end+1,:) = {'Low Wavenumber (0 = use original)', 'LowWN',[]};
Formats(n,1).type = 'edit';
Formats(n,1).format = 'integer';
Formats(n,1).limits = [0 9999]; % 9-digits (positive #)
Formats(n,1).size = 80;
Formats(n,1).unitsloc = 'bottomleft';
DefAns.LowWN = 450;

n=n+1;
Prompt(end+1,:) = {'High Wavenumber (0 = use original)', 'HighWN',[]};
Formats(n,1).type = 'edit';
Formats(n,1).format = 'integer';
Formats(n,1).limits = [0 9999]; % 9-digits (positive #)
Formats(n,1).size = 80;
Formats(n,1).unitsloc = 'bottomleft';
DefAns.HighWN = 1800;

n=n+1;
Prompt(end+1,:) = {'PCA Type','PCAType',[]};
Formats(n,2).type = 'list';
Formats(n,2).format = 'text';
Formats(n,2).style = 'radiobutton';
Formats(n,2).items = {'NIPALS' 'SVD' ''};
DefAns.PCAType = 'SVD';

n=n+1;
Prompt(end+1,:) = {'Normalise' 'NormaliseMode',[]};
Formats(n,1).type = 'check';
DefAns.NormaliseMode = false;

n=n+1;
Prompt(end+1,:) = {'Mean Centre' 'MeanCentreMode',[]};
Formats(n,1).type = 'check';
DefAns.MeanCentreMode = true;

n=n+1;
Prompt(end+1,:) = {'Number of PCs', 'PCs',[]};
Formats(n,1).type = 'edit';
Formats(n,1).format = 'integer';
Formats(n,1).limits = [0 999]; % 9-digits (positive #)
Formats(n,1).size = 80;
Formats(n,1).unitsloc = 'bottomleft';
DefAns.PCs = 20;

n=n+1;
Prompt(end+1,:) = {'Choose Baseline Type','BaseLineType',[]};
Formats(n,1).type = 'list';
Formats(n,1).format = 'text';
Formats(n,1).style = 'radiobutton';
Formats(n,1).items = {'Polynomial' 'EMSC';'None' ''};
DefAns.BaseLineType = 'None';

n=n+1;
Prompt(end+1,:) = {'EMSC/Polynomial Order', 'PolynomialOrder',[]};
Formats(n,1).type = 'edit';
Formats(n,1).format = 'integer';
Formats(n,1).limits = [0 999]; % 9-digits (positive #)
Formats(n,1).size = 80;
Formats(n,1).unitsloc = 'bottomleft';
DefAns.PolynomialOrder = 4;

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

PCARun.PolynomialOrder=Answer.PolynomialOrder;
