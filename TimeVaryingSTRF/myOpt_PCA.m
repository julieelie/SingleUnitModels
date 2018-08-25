function [H, H0, P_num, SSE, PC_invest] = myOpt_PCA(Y, X, varargin)
%% This function performs a PC regression between X and Y
%   Y is a vector of responses of size N, X is a matrix of stimuli with N
% rows of observations and P columns of parameters
% [H, H0, P_num] = MYRIDGE(Y, X) returns the 1 X P_num vector H
% (the filter given by the PC algorithm for the P parameters of the stimuli); 
% P_num is the number of eigen values retained in the svd.
%
%   [...] = myOpt_PCA(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies validating
% dataset. If nothing is specified, R2 and SSE are calculated on the
% training dataset
% Valid parameters of lambdas are the following:
%        Parameter      Value
%         'Validating Y'    The vector of responses of size M.
%                             The default is all indices.
%         'Validating X'     The matrix of stimuli with M
%                           rows of observations and P columns of parameters
%         'Dim'             the number of singular values along which X should
%                           be decomposed in other words the number of principal
%                           components in the PCA

FigFlag=1; % Set to 1 to see debugging figures
pnames = {'Validating Y'  'Validating X' 'Dim'};
dflts  = {Y    X []};
[Y_val,X_val,Nb_PC] = internal.stats.parseArgs(pnames,dflts,varargin{:});
if isempty(Nb_PC)
    SVD=1;% princomp (pca) is taking too much memory use svd by default:
%[U_local,W_local,V_local]=svd(Xridge, 'econ');
else
    SVD=2;
end

%% Run PCA on stimuli to find eigenvalues (W_local) and eigen vectors (V_localT)
Xcentered = X - repmat(mean(X), size(X,1),1);
Ycentered = Y - mean(Y);
if SVD==0
    [V_local,SCORE,W_local2] = pca(X);
    Nb_1=size(X,1)-1;
    W_local2 = W_local2.*Nb_1;
    %% Find out all the filters by progressively using the PCs of the PCA and
    %Calculate the SSE on the training or validating dataset for each of the filter
    nPC = length(W_local2);
    if nPC>40
        PC_invest = 5:1:40;
    else
        PC_invest = 5:1:nPC;
    end
    nPC_invest=length(PC_invest);
    H_local = nan(nPC_invest,size(X,2));
    H0_local = nan(nPC_invest,1);
    SSE = nan(nPC_invest,2);
    for pc=1:length(PC_invest)
        pp=PC_invest(pc);
        %fprintf(1,'Calculating filter H: %d/%d\n', pp, nPC);
        W_local_regress = diag([1./(W_local2(1:pp)); zeros((nPC-pp),1)]);
        H_local(pc,:) = V_local * W_local_regress * V_local' * Xcentered' * Ycentered;
        H0_local(pc) = mean(Y) - H_local(pc,:) * mean(X)';
        YY_val = Y_val- (H0_local(pc) + H_local(pc,:)*X_val')';
        YY2_val = YY_val.^2;
        SSE(pc,1) = sum(YY2_val);
        YY_tra = Y- (H0_local(pc) + H_local(pc,:)*X')';
        YY2_tra = YY_tra.^2;
        SSE(pc,2) = sum(YY2_tra);
    end
    %% Ckecks
    Score = Xcentered * V_local;
    if sum(sum(SCORE-Score))>10^(-6)
        fprintf(1,'WARNING check out the calculations!!!!');
    end
    P_num = find(SSE(:,1)==min(SSE(:,1)));
    H=H_local(P_num,:);
    H0=H0_local(P_num);
    
elseif SVD==1
    [U_localsvd,W_local,V_local]=svd(Xcentered, 'econ');
    W_localsvdd=diag(W_local).^2;
    %% Find out all the filters by progressively using the PCs of the PCA and
    %Calculate the SSE on the training or validating dataset for each of the filter
    nPC = length(W_localsvdd);
    if nPC>40
        PC_invest = 5:1:40;
    else
        PC_invest = 5:1:(nPC-1);
    end
    nPC_invest=length(PC_invest);
    H_local = nan(nPC_invest,size(X,2));
    H0_local = nan(nPC_invest,1);
    SSE = nan(nPC_invest,2);
    for pc=1:nPC_invest
        pp=PC_invest(pc);
        %fprintf(1,'Calculating filter H: %d/%d\n', pp, nPC);
        W_local_regress = diag([1./(W_localsvdd(1:pp)); zeros((nPC-pp),1)]);
        H_local(pc,:) = V_local * W_local_regress * V_local' * Xcentered' * Ycentered;
        H0_local(pc) = mean(Y) - H_local(pc,:) * mean(X)';
        YY_val = Y_val- (H0_local(pc) + H_local(pc,:)*X_val')';
        YY2_val = YY_val.^2;
        SSE(pc,1) = sum(YY2_val);
        YY_tra = Y- (H0_local(pc) + H_local(pc,:)*X')';
        YY2_tra = YY_tra.^2;
        SSE(pc,2) = sum(YY2_tra);
    end
    P_num = find(SSE(:,1)==min(SSE(:,1)));
    H=H_local(P_num,:);
    H0=H0_local(P_num);
elseif SVD==2
    [U_localsvd,W_local,V_local]=svds(Xcentered, Nb_PC,'L');
    W_localsvdd=diag(W_local).^2;
    %% Find out all the filter and Calculate the SSE on the training or validating dataset
    %fprintf(1,'Calculating filter H: %d/%d\n', pp, nPC);
    SSE=nan(1,2);
    W_local_regress = diag(1./W_localsvdd);
    H = V_local * W_local_regress * V_local' * Xcentered' * Ycentered;
    H0 = mean(Y) - H' * mean(X)';
    YY_val = Y_val- (H0 + H'*X_val')';
    YY2_val = YY_val.^2;
    SSE(1) = sum(YY2_val);
    YY_tra = Y- (H0 + H'*X')';
    YY2_tra = YY_tra.^2;
    SSE(2) = sum(YY2_tra);
    P_num=Nb_PC;
end

if FigFlag==1 && SVD<2
    figure(1)
    plot(PC_invest,SSE(:,1),PC_invest,SSE(:,2))
    legend('SSE validating set','SSE training set')
     xlabel('# of PCs')
     ylabel('SSE')
     pause(1);
end

end