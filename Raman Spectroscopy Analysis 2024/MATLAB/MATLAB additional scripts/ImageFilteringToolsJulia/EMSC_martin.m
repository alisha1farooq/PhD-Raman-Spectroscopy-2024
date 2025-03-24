function [c,m,h,g,E,Zcor,P,D,V,coeff]=EMSC_martin(Zraw,v,Zref,model,pow)
% EMSC algorithm
% Zraw = c + mv + hZref + gModel + E
%
% INPUTS
% Zraw  = raw spectrum (I x J)
% v     = wavenumbers  (1 x J)
% Zref  = reference spectrum (1 x J)
% model = dataset to create a model of for correction e.g. parafin dataset,gg spectra (K x J)
%         NB if no model then use a row of 1xJ zeros;
% N     = if more than one spectrum in model, then PCA can be used with N components
%         NB if no model then use N=0;
% pow   = polynomial to model for v (if excluded = 1)
%
% OUTPUTS
% c,m,h,g,E,P and D are modelling parameters
% Zcor (I x J) contains the EMSC corrected spectra
%
if isempty(model)
    model=zeros(size(Zref));
end
P=model;

% use least squares to estimate parameters

% I = number of spectra to correct
% J = number of wavenumbers
[I,J]=size(Zraw);

% construct design matrix
V=zeros(J,pow);
for i=1:pow
    V(:,i)=v(:).^i;
end

D=[ones(J,1),V,Zref(:),P'];
% D=[ones(J,1),v(:),P'];
% least squares
coeff=pinv(D)*Zraw';

c=coeff(1,:); % offset
m=coeff(2:2+pow-1,:); % baseline
h=coeff(2+pow,:); % ref
% h=0;
g=coeff(2+pow+1:end,:); % model
% g=coeff(3:end,:);
Zcor=D*coeff;
E=Zraw'-Zcor;
G=g'*P;
M=m'*V';
Zcor=(Zraw'-(repmat(c,J,1)+M'+G'))./repmat(h,J,1); %this is the corrected matrix
% Zcor=(Zraw'-repmat(c,J,1)-G')./repmat(h,J,1);
% Zcor=(Zraw'+E)./repmat(h,J,1);
Zcor=Zcor'; %this is the corrected matrix

% Principal Component Analysis using the NIPALS algorithm
function [T,P,eig]=PCA(X,A,tol)
% by Gavin Rhys Lloyd 12/01/10
%
% Inputs:
% X   = preprocessed data matrix (I samples x J variables)
% A   = number of components to calculate
% tol = tolerance (optional, default = 0.001)
%
% Outputs:
% T = Scores matrix (IxA)
% P = Loadings matrix (AxJ)

% if tolerance not specified use default
if nargin<3
    tol=1e-6;
end

% prepare output matrices...
T=zeros(size(X,1),A); % scores
P=zeros(A,size(X,2)); % loadings
eig=zeros(1,A);
X2=sum(sum(X.^2));

% initial guess of t (column of X with max sum of squares)
[hi,ind]=max(sum(X.^2));
t_hat=X(:,ind);

% for each component...
for i=1:A
%     i
    carryon=true; % while loop terminator
    while carryon
        % project X onto t to get loadings...
        p_hat=(X'*t_hat)/(t_hat'*t_hat);      
        % normalise the loadings...
        p_hat=p_hat/sqrt(p_hat'*p_hat);
        % project X onto P to get scores...
        t=(X*p_hat)/(p_hat'*p_hat);
               
        % check for convergence...
        if sqrt((t_hat-t)'*(t_hat-t))<tol
            T(:,i)=t;
            P(i,:)=p_hat;
            X=X-(t*p_hat');
            carryon=false;
        end
        % if not converged, try again using current scores estimate...
        t_hat=t;
    end
end

eig=(sum(T.^2))./X2;




