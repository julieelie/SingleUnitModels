function [I,P_YgivenS_all,H_y, H_ymu]=info_model_Calculus(Mu_in)
%% Calculates the information obtain from N Poisson with different means.

% This functions is used to calculate the mutual information between n
% events (stimuli, categories of stimuli, motor output...) and
% a poisson response defined by its n means. The maximum value is log2(n).
% The minum value is zero and will be obtained when the n means are identical.
% This function assumes a Poisson noise for the neural response with a mean
% between 0 and 50.
% Write your own routine for values higher than 50 using a the Normal approximation
% for the Poisson distribution.


% Requires:
% Mu_in         vector 1xn that specifies the mean (in count/bin i.e. no
%               units here) of the n Poisson distributions.


% Returns: 
% I             the mutual information in bits between n events (stimuli,
%               categories of stimuli, motor output...) and a poisson
%               response defined by its means mu_in.

% P_YgivenS_all a matrix ym x n where each column gives the conditional
%               probability distribution of the response given the stimulus
%               i.e. the input mean mu_in(n). ym corresponds to the number 
%               of values of y (0:Ymax) that are used to estimate entropies.

% H_y           a scalar that specifies the entropy of the response 

% H_y_mu        a scalar that specifies the entropy of the response given
%               the stimulus i.e., the input means Mu_in.

%% Sanitary checks and parameters setting
% Check for max value of Mu_in to be below upper bound for a Poisson
% distribution assumption.
MAXMEANVALUE = 50;
if ~sum(Mu_in <= MAXMEANVALUE)
    ERRMSG= 'This function assumes Poisson distribution means inferior or equal to %d spikes per bin. Write your own routine for values higher than 50 using a the Normal approximation for the Poisson distribution.';
    error(ERRMSG, MAXMEANVALUE)
end


% Set a lower bound for Mu_in
% We need to attend the case of mu=0 which ends up in abberant calculations
% because log(0) = -Inf.
% The natural increment for mu is 1 so let's have log(mu)
% almost linear around 0 (around mu=1) and fix muLims=1/20.
MuLims=1/20;
% Here we know the natural scale for mu: the neuron can spike between
%0 and 1 time per ms so depending on the window size Win (ms) the range
% of values is 0 to Win spikes. Make the lower limit for log(mu) as
% - the max limit by replacing zero values in mu by muLims=1/Win.
% Other idea: 
% In glmfit and here if no scale is provided, other choice because they don't know the
% natural scale of mu. Their choice keeps mu^4 from underflowing.No upper limit.
% dataClass = superiorfloat(mu,y);
% muLims = realmin(dataClass).^.25; %Choice in glmfit


% Y_world defines the range of values that are used to calculate the entropies
% Ymax corresponds to the maximum number of spikes you can observed per bin
% given the iput mean. We can indeed reasonably limit the calculation
% of entropy from 0 up to the maximum number of spikes expected instead of
% summing to + infinity.
% Ymax is three standard deviations above the maximum variance (note that
% the mean equals the variance for Poisson distributions).
Ymax = 3.0*sqrt(max(Mu_in)) + max(Mu_in);
Y_world = 0:Ymax;

% To save time, calculations are only done for unique values of Mu_in    
[Mu_unique, ~, I_mu_unique] = unique(Mu_in); 
Nb_Unq_mu = length(Mu_unique);


%% Calculate the conditional entropy H(y/mu) assuming a Poisson distribution
% Initialize output variables
N_mu = length(Mu_in);
P_YgivenS_all = nan(length(Y_world),N_mu);
P_YgivenS_all_rescaled = nan(length(Y_world),N_mu);

% To save time calculation of conditional entropy is only done for unique values of mu.
H_ymu_unique = nan(size(Mu_unique));

% Y factorial
Fac_y = factorial(Y_world);   % Saved to save computational time.

% Now loop through unique mu and calculate the conditional entropy for each
% of them
for mm=1:Nb_Unq_mu
    
    % Local value of mean
    Mu = Mu_unique(mm);
    
    % Calculate the log of the probability for numerical reasons.
    Log_P = nan(size(Y_world));   % Initialize
    
    % Deal with cases where mu is equal to zero and/or is below muLims
    if Mu == 0
        Log_P(1) = 1;  % The value for y=0
        Mu = MuLims;
        Log_P(2:end)= Y_world(2:end).*log2(Mu) - log2(Fac_y(2:end)) - Mu/log(2);
    elseif Mu < MuLims
        Mu = MuLims;
        Log_P = Y_world.*log2(Mu) - log2(Fac_y) - Mu/log(2);
    else
        Log_P = Y_world.*log2(Mu) - log2(Fac_y) - Mu/log(2);
    end
    
    % From logP to P
    P = 2.^(Log_P);
    P_YgivenS_all(:,(I_mu_unique==mm))=repmat(P',1,sum(I_mu_unique==mm));
   
    % Check that sum(P) sum to 1 and rescale if not to have a more
    % realistic distribution.
    P = P ./sum(P);
    P_YgivenS_all_rescaled(:,(I_mu_unique==mm))=repmat(P',1,sum(I_mu_unique==mm));
    % conditional entropy value for each unique mu (each stim) across all y of the world
    H_ymu_unique(mm) = sum(-P.*log2(P + (P==0))); % Note that adding P==0 in the logarithmic expression ensures that P*log2(P) = 0 when P=0.
end

% Given that all stimuli have equal probability: H(y|s) = expectation(H(y|si)) = sum(p(si)*H(y|si)) = sum(H(y|si))/(number of stims)
H_ymu = sum(H_ymu_unique(I_mu_unique))/N_mu;


%% Calculate the entropy of the response y based on Poisson distributions
% Probability of responses (marginal)
P_y_i = sum(P_YgivenS_all_rescaled,2)./N_mu;

% Entropy of the response
H_y = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% conditional entropy value for each unique y of the world across all mu of the model

%% Calculate mutual information
I = (H_y - H_ymu);

end



