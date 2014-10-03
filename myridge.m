function [H, H0, V, W, L, P_num] = myridge(Y, X, L)
FigFlag=1;
P_thresh = 0.99;%Threshold to choose the number of parameters or eigen values in the singular value decomposition 

Xridge = X - repmat(mean(X), size(X,1),1);
Yridge = Y - mean(Y);
[U_local,W_local,V_local]=svd(Xridge, 'econ');    

if FigFlag==1
    figure()
     plot(1:size(W_local,1),cumsum(power(diag(W_local),2)./sum(power(diag(W_local),2))))
     hold on
     hline(0.99)
     hold off
     xlabel('# of parameters')
     ylabel('Cumulative sum of the normalized eigen values')
end

P_num = min(find(cumsum(power(diag(W_local),2)./sum(power(diag(W_local),2)))>P_thresh));
[U, W, V]=svds(Xridge, P_num, 'L');

% Define Lambdas if they were not given
if nargin<3
    CO=[0.001 0.01 0.1 1 10 20 30 40 50 75 100 150 200 300 600 1e3 5e3 1e4 5e4 1e5 5e5 1e6 5e7 1e8 5e8 1e9];
    L=CO.*mean(diag(W).^2);
end

H=nan(length(L), size(X,2));
H0=nan(length(L), 1);

for ll=1:length(L)
    if FigFlag==1
        fprintf(1,'Calculation of filter H: %d/%d lambda\n',ll,length(L));
    end
    WL = 1./(diag(W).^2 + L(ll));
    H(ll,:) = V*diag(WL)*V'*Xridge'*Yridge;
    H0(ll) = mean(Y) - H(ll,:)*mean(X)';
end
end