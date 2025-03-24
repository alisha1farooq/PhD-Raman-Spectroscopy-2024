function [spectra,y_hat]=removeBaseline(x,y,order,pos,max_iter,binning,known)
% Function to correct the baseline of Raman spectra
% by Gavin Rhys Lloyd 30.04.14
% REF: Applied Spectroscopy, Volume 57, Issue 11,Pages 320A-340A and 1317-1453 (November 2003) , pp. 1363-1367(5)
tol=1e-8;
cont=true;
counter=0;
x_orig=x;
y_orig=y;
known_orig=known;

if strcmpi('none',known)
    known_orig=zeros(size(x));
    known=known_orig;
end

x=x(1:binning:length(x));
y=y(1:binning:length(y));
known=known(1:binning:length(known));

while cont
    if counter>max_iter
        cont=false;
        %                 disp('max number of iterations reached');
    end
    [na,y_hat,b]=fitpoly(x,y,order,[],known);
    
    % replace data with min of poly and data
    y=min([y_hat(:),y(:)]');
    
    % check for convergence...
    if (y-y_hat')*(y-y_hat')'<tol %| counter>max_iter
        cont=false;
    end
    counter=counter+1;
end
[na,y_hat,b]=fitpoly(x,y,order,x_orig,known,known_orig);

spectra=(y_orig-y_hat');
if pos
    spectra(spectra<0)=0;
    %     spectra=spectra-min(spectra);
end
end

function [x_hat,y_hat,b]=fitpoly(x,y,order,x_new,knowns,known_orig)
% function to fit a polynomial to data
% by Gavin Rhys Lloyd 14.01.10
if nargin<4
    x_new=x;
end
if isempty(x_new)
    x_new=x;
end
m=mean(x);
y=y(:);
x=x(:)-m;

p=1:order;
p=repmat(p,length(x),1);
D=[ones(size(x)),repmat(x,1,order).^p,knowns'];
b=pinv(D)*y;
if nargin==6
    x_new=x_new(:)-m;
    p=1:order;
    p=repmat(p,length(x_new),1);
    D=[ones(size(x_new)),repmat(x_new,1,order).^p,known_orig'];
end
y_hat=D*b;
x_hat=x;
end
