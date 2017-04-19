function [Icum_Full, HY_Full, HYgivenS_Full, Icum_JK, HY_JK, HYgivenS_JK, N_MC_tot, Icum_std]=info_cumulative_model_Calculus_MCJK(P_YgivenS_Full,P_YgivenS_JK, NTrials, varargin)
tic
%% This function calculates the cumulative information of the neural response
% P_YgivenS is a cell array of length equals to the number of time points t in ...
... the neural responses. Each cell contains a c*N matrix of probabilities,
    ... the conditional probability of a given neural spike count c at the
    ... particular time point t for the stimulus N. The cumulative
    ... information is calculated/estimated for the last time point from
    ... the first time point. Icum is a unique floating point (value of ...
    ... cumulative information at time t). The algorithm uses a MonteCarlo
    ... estimation with a weight correction, the number of samples used 
    ... for the estimation is given by MCParameter.
    
%   [...] =
%   info_cumulative_model_Calculus(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   specifies various parameters of the calculation.
% Valid parameters are the following:
%  Parameter      Value
%         
%   'MaxMCParameter' is the maximum number of samples in the Monte Carlo estimation
%                    of the entropy of the response defaulted to 10^6
%   'IncrMCParameter' is the increment by which the number of samples is
%                   increased in the MonteCarlo estimation of the entropies
%   'ConvThresh'  is the targeted maximum error in bits on the calculation of
%                   cumulative information at each bin
%   'MinProbThresh'  set to 1 to discard probabilities lower than 1/life of
%                    a zebra finch in the lab, set to 0 to take into
%                    acccount all probabilities. Default 0.
%   'DebugFig'      set to 1 to see some debugging figures
%   'MinProb'       threshold of probability under which values are set to
%                   zero. Default is 1/number of 20ms bins in a life of a
%                   zebra finch in the lab (8 years)
%   'ScaleY'        Set to 1 to scale the distribution of conditional
%                   probabilities given the stimulus the same way for both
%                   entropy calculations
%   'Verbose'       set to 1 to see a lot of output, 0 for none default is
%                   none

%% Sorting input arguments
pnames = {'MaxMCParameter','IncrMCParameter','ConvThresh','MinProbThresh', 'DebugFig','MinProb','ScaleY', 'Verbose'};

% Calculating default values of input arguments
NStim_local = size(P_YgivenS_Full{end},2);


MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab

% Get input arguments
dflts  = {10^6 10^4 0.3 0 0 MinProb 1 0};
[N_MC_max, N_MC, ConvThresh, MinProbThresh, DebugFig, MinProb, ScaleY, Verbose] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%% Set some parameters
win = length(P_YgivenS_Full);
Nb_Boot = size(P_YgivenS_JK,1);

%% Construct the data set of individual p for all JK Bootstrap and also for the full dataset
P_YgivenS_all = cell(win,Nb_Boot+1);
%parfor
parfor bb=1:(Nb_Boot+1)
    if bb==1
        P_YgivenS_local = P_YgivenS_Full;
    else
        bb_local = bb-1;
        P_YgivenS_local = P_YgivenS_JK(bb_local,:);
    end
    % Gather for each window before the current window win the probabilities of
    % the stimuli that are still part of the calculations at the window win
    for ww=1:(win)
        % Don't keep the probabilities corresponding to least probable spike
        % count values for all stimuli at the time point win
        if MinProbThresh
            Badyy = find(sum((P_YgivenS_local{ww} < MinProb),2)==NStim_local);
            Goodyy=setdiff(1:size(P_YgivenS_local{1},1), Badyy);
            if ScaleY
                P_YgivenS_all{ww,bb}=P_YgivenS_local{ww}(Goodyy,:)./repmat(sum(P_YgivenS_local{ww}(Goodyy,:),1),size(P_YgivenS_local{ww},1),1);
            else
                P_YgivenS_all{ww,bb}=P_YgivenS_local{ww}(Goodyy,:);
            end
        else
            if ScaleY
                P_YgivenS_all{ww,bb}=P_YgivenS_local{ww}(:,:)./repmat(sum(P_YgivenS_local{ww}(:,:),1),size(P_YgivenS_local{ww},1),1);
            else
                P_YgivenS_all{ww,bb}=P_YgivenS_local{ww};
            end
        end
    end
end


%% Information

HY_Full = 0;
HYgivenS_Full = 0;
HYgivenS_JK = zeros(1,Nb_Boot);
HY_JK = zeros(1, Nb_Boot);
N_MC_tot = 0;
Icum_std = ConvThresh+1; % set the error on information to an arbitrary value above the threshold

while Icum_std>ConvThresh && N_MC_tot<N_MC_max
    [HY_Full_local,HYgivenS_Full_local,HY_JK_local,HYgivenS_JK_local]= cuminfo_MC(P_YgivenS_all, N_MC);

    N_MC_tot_new = N_MC + N_MC_tot;
    if N_MC_tot ~= 0
        HY_Full = (HY_Full * N_MC_tot + HY_Full_local * N_MC) / N_MC_tot_new;
        HYgivenS_Full = (HYgivenS_Full * N_MC_tot + HYgivenS_Full_local * N_MC) / N_MC_tot_new;
        HYgivenS_JK = (HYgivenS_JK .* N_MC_tot + HYgivenS_JK_local .* N_MC) ./ N_MC_tot_new;
        HY_JK = (HY_JK .* N_MC_tot + HY_JK_local.*N_MC) ./ N_MC_tot_new;
    else
        HY_Full = HY_Full_local;
        HYgivenS_Full = HYgivenS_Full_local;
        HYgivenS_JK = HYgivenS_JK_local;
        HY_JK = HY_JK_local;
    end
    N_MC_tot = N_MC_tot_new;

    Icum_Full = HY_Full-HYgivenS_Full;
    Icum_JK = HY_JK-HYgivenS_JK;

    Icum_JKcorrected = NTrials * repmat(Icum_Full,1,Nb_Boot) - (NTrials-1) * Icum_JK(1:Nb_Boot);
    Icum_std = std(Icum_JKcorrected,1);
end


if Verbose
    fprintf(1,'Cumulative information = %f bits\n',Icum_Full);
    ElapsedTime = toc;
    fprintf(1,'info_cumulative_model_Calculus_MCJK run for %d seconds with %d samples\n',ElapsedTime, N_MC_tot);
    %fprintf('info_cumulative_model_Calculus run for %d seconds or %d seconds per path for %d paths\n',ElapsedTime, ElapsedTime/NbP, NbP);
end




%% Strategy: Monte Carlo Estimate of the entropy of the response
function[HY_Full,HYgivenS_Full,HY_JK,HYgivenS_JK]= cuminfo_MC(P_YgivenS_all, N_MC)

if Verbose
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples\n',N_MC);
end

% Calculate the cdf of the marginal probabilities (probabilities of
% responses at each time point) of the dataset using all trials
P_Yt=cell(win,1);
for www=1:win
    P_Yt{www} = cumsum(mean(P_YgivenS_all{www,1},2)./sum(mean(P_YgivenS_all{www,1},2)));
end

%% Set the number of MonteCarlo samples N_MC and their values
MC_Randu = rand(N_MC, win);
Resp_MC = nan(N_MC, win);

%% First calculate for the full dataset as we need QY_MC of that dataset to
% implement in the calculation of JK datasets
PYgivenS_MC = nan(N_MC,NStim_local); %exact joint probability of MC sequences of responses for each stimulus
PY_MC = nan(N_MC,1); %exact joint probability of MC sequences of responses averaged over the stimuli
QY_MC = nan(N_MC,1); %estimation of the joint probability of MC sequences of responses using the marginals

% Loop through the number of samples
for ss=1:N_MC
    if sum(ss==(1:10)*N_MC/10) && Verbose
        fprintf('%d/%d MC samples Full dataset\n',ss,N_MC);
    end
    % Probabilistically determine a sequence of responses using the marginal
    % probabilities
    % for each time point, take a random number from the uniform
    % distribution 0-1 and find the neural response that correspond to
    % that probability in the full dataset
    P_YgivenS_local_Resp = nan(win,NStim_local);
    for www=1:win
        u=MC_Randu(ss,www);
        Resp_MC(ss, www) = sum(P_Yt{www}<u)+1;
        P_YgivenS_local_Resp(www,:) = P_YgivenS_all{www,1}(Resp_MC(ss,www),:);
    end

    % Calculate the exact joint probability of that sequence of responses
    % and store it
    PY_MC(ss) = mean(prod(P_YgivenS_local_Resp,1));% Joint proba over bins (product of elements percolumn) then average over stimuli
    QY_MC(ss) = prod(mean(P_YgivenS_local_Resp,2));% Mean over stims (kind of getting the marginal) then product of the probability of the path
    PYgivenS_MC(ss,:) = prod(P_YgivenS_local_Resp,1);
end
% Calculate the MC estimate of the conditional response entropy to the
% stimulus
HYgivenS_Full = - sum(sum(PYgivenS_MC./repmat(QY_MC,1,NStim_local).*log2(PYgivenS_MC),1))/(N_MC*NStim_local);
if Verbose
    fprintf(1, 'Conditional entropy (entropy of the neural response given the stimulus): %f\n', HYgivenS_Full);
end

% Calculate the MC estimate of the response entropy
% Weights
HY_Full = sum(-PY_MC./QY_MC.*log2(PY_MC))/N_MC;

%% Now calculate for JackKnife sets using the distribution of probabilities calculated for the full dataset to choose the MC samples
HYgivenS_JK = nan(1,Nb_Boot);
HY_JK = nan(1, Nb_Boot);

%parfor
parfor bout=1:(Nb_Boot)
    PYgivenS_MC_bb = nan(N_MC,NStim_local); %exact joint probability of MC sequences of responses for each stimulus
    PY_MC_bb = nan(N_MC,1); %exact joint probability of MC sequences of responses averaged over the stimuli
    % Loop through the number of samples
    for ss=1:N_MC
        if sum(ss==(1:10)*N_MC/10) && Verbose
            fprintf('%d/%d MC samples JK bootstrap %d/%d\n',ss,N_MC,bout-1,Nb_Boot);
        end
        % Probabilistically determine a sequence of responses using the marginal
        % probabilities
        % for each time point, take a random number from the uniform
        % distribution 0-1 and find the neural response that correspond to
        % that probability in the full dataset
        P_YgivenS_local_Resp = nan(win,NStim_local);
        for www=1:win
            P_YgivenS_local_Resp(www,:) = P_YgivenS_all{www,bout}(Resp_MC(ss, www),:);
        end
    
        % Calculate the exact joint probability of that sequence of responses
        % and store it
        PY_MC_bb(ss) = mean(prod(P_YgivenS_local_Resp,1));% Joint proba over bins (product of elements percolumn) then average over stimuli
        PYgivenS_MC_bb(ss,:) = prod(P_YgivenS_local_Resp,1);
    end


    % Calculate the MC estimate of the conditional response entropy to the
    % stimulus
    HYgivenS_JK(bout) = - sum(sum(PYgivenS_MC_bb./repmat(QY_MC,1,NStim_local).*log2(PYgivenS_MC_bb),1))/(N_MC*NStim_local);
    if Verbose
        fprintf(1, 'JK Bootstrap %d: Conditional entropy (entropy of the neural response given the stimulus): %f\n',bout, HYgivenS_JK(bout));
    end

    % Calculate the MC estimate of the response entropy
    % Weights
    HY_JK(bout) = sum(-PY_MC_bb./QY_MC.*log2(PY_MC_bb))/N_MC;

    % % Estimate the Error on MC estimates (based on the asymptotics method
    % % of Koehler et al. 2009 Am. Stat.
    % PYgivenS_MC_PY = PYgivenS_MC./repmat(PY_MC, 1,length(Stim_local));
    % Icum_supp = sum(PY_MC./QY_MC.*sum(PYgivenS_MC_PY.*log2(PYgivenS_MC_PY),2))/(length(Stim_local)*N_MC);
    % Icum_MCE = ((sum(PY_MC./QY_MC.*(sum(PYgivenS_MC_PY.*log2(PYgivenS_MC_PY),2)./length(Stim_local) - Icum(1)).^2)).^0.5)/N_MC;


    if DebugFig
        figure()
        plot(1:100, QY_MC(1:100), 1:100,PY_MC(1:100))
        legend('Product of probas','exact joint probas')
        xlabel('Monte Carlo samples')
        ylabel('probability')
    end


    if Verbose
        fprintf(1, 'JK Bootstrap %d: Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples: %f\n',bout,N_MC, HY_JK(bout));
    end

end

end
end
