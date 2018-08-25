function [H, H0, V, W, Lvalues, P_num, Hprime, H0prime] = myridge(Y, X, varargin)
%% This function performs a ridge regression between X and Y
%   Y is a vector of responses of size N, X is a matrix of stimuli with N
% rows of observations and P columns of parameters
% [H, H0, V, W, Lvalues, P_num] = MYRIDGE(Y, X) returns the M by P matrix H
% (the M filters given by the ridge algorithm for the M values of lambda
% tested and the P parameters of the stimuli); H0 is the 1xM vector of bias
% parameter; V and W correspond to the matrices obtained by the svd X such
% that X = U*W*V'; Lvalues is the 1xM vector of Lambda values tested as
% ridge parameters; P_num is the number of eigen values retained in the svd.
%
%   [...] = MYRIDGE(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies lambda
%   parameters and their values. If nothing is specified, default values of
% Lambda are calculated as:
% CO.*mean(diag(W).^2)
% with CO=[10^-20 10^-5 5*10^-5 10^-4 5*10^-4 0.001 0.005 0.01 0.05 0.1 0.5 1 10 20 30 40 50 75 100 150 200 300 600 1e3 5e3 1e4]
% Valid parameters of lambdas are the following:
%        Parameter      Value
%         'Lambda Indices'    The indices of CO that the function
%                             should test (any entire value between 1 and 26.
%                             The default is all indices.
%         'Lambda Values'     The values of Lambda that the function
%                             should test.

FigFlag=0; % Set to 1 to see debugging figures
P_thresh = 0.99;%Threshold to choose the number of parameters or eigen values in the singular value decomposition 

Xridge = X - repmat(mean(X), size(X,1),1);
Yridge = Y - mean(Y);
[U_local,W_local,V_local]=svd(Xridge, 'econ');
P_num = min(find(cumsum(power(diag(W_local),2)./sum(power(diag(W_local),2)))>P_thresh));

if FigFlag==1
    figure(1)
     plot(1:size(W_local,1),cumsum(power(diag(W_local),2)./sum(power(diag(W_local),2))))
     hold on
     hline(0.99)
     hold off
     xlabel('# of parameters')
     ylabel('Cumulative sum of the normalized eigen values')
     title(sprintf('number of eigen values choosen=%d', P_num));
     pause(1);
end


[U, W, V]=svds(Xridge, P_num, 'L');

% Define Lambdas
CO=[10^-20 10^-5 5*10^-5 10^-4 5*10^-4 0.001 0.005 0.01 0.05 0.1 0.5 1 10 20 30 40 50 75 100 150 200 300 600 1e3 5e3 1e4];
pnames = {'Lambda Indices'  'Lambda Values'};
dflts  = {[]    CO.*mean(diag(W).^2)};
[LIndices,Lvalues] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if ~isempty(LIndices)
    Lvalues=CO(LIndices).*mean(diag(W).^2);
end


H=nan(length(Lvalues), size(X,2));
Hprime=nan(length(Lvalues), size(V,2));
H0=nan(length(Lvalues), 1);
H0prime=nan(length(Lvalues), 1);


for ll=1:length(Lvalues)
    if FigFlag==1
        fprintf(1,'Calculation of filter H: %d/%d lambda\n',ll,length(Lvalues));
    end
    WL = 1./(diag(W).^2 + Lvalues(ll));
    WL_h = WL.^(0.5);
    X_prime = X*V*diag(WL_h)';
    %Hprime(ll, :) = (diag(sqrt(WL))*V'*Xridge'*Yridge)';
    % H(ll,:) = V*diag(sqrt(WL)*Hprime(ll,:);
    H(ll,:) = V*diag(WL)*V'*Xridge'*Yridge;
    H0(ll) = mean(Y) - H(ll,:)*mean(X)';
    Hprime(ll,:)=H(ll,:)*V*(diag((diag(W).^2 + Lvalues(ll)).^(1/2)));
    H0prime(ll) = mean(Y) - Hprime(ll,:)*mean(X_prime)';
end
end