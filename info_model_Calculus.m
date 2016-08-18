function [I,P_YgivenS_all1,P_YgivenS_all2]=info_model_Calculus(mu_in,mu_max, ymax,y,Scale)
%% This function calculate how much information there is in the prediction of a model
% Note that there is a trade off in the information value calculation
% between the values that mu can get and how widely the entropy of the model
% can be estimated over the possible values of y. if ymax is fixed to 170
% then mu has to be lower than 65 for now.
% mu_in are the predictions from the model for every single stim and every
% single trial
% ymax is the maximum number of spike you can observed in your dataset. We
% can indeed reasonably uperbound the calculation of entropy on all
% possible y but the maximum nuber of spikes you can observed given the
% size of your observation window (like 20ms -> 20 spikes) instead of
% summing to + infinity.
% y is optional and is the observed data if you want to estimate the
% entropy of the neural response as a histogram (Bialeck kind of
% estimation)

% mu_max is the maximum value expected for the model prediction, anything
% above it does not give any information as the model of the response is
% obviously wrong. For those values of mu>mu_max, the probability is set to
% zero

% Here we know the natural scale for mu: the neuron can spike between
%0 and 1 time per ms so depending on the window size Win (ms) the range
% of values is 0 to Win spikes. Make the lower limit for log(mu) as
% - the max limit by replacing zero values in mu by muLims=1/Win.
% Other idea: the natural increment for mu is 1 so let's have log(mu)
% almost linear around 0 (around mu=1) and fix muLims=1/20.
% In glmfit and here if no scale is provided, other choice because they don't know the
% natural scale of mu. Their choice keeps mu^4 from underflowing.No upper limit.

SameDataScaleEntropy=1; %set this to zero if you prefer to scale data differently but that sounds weird!!

if nargin<2
    mu_max = 100*median(mu_in);
end

if any(mu_in>(10*mu_max))
    OkPredict = find(mu_in<(10*mu_max));% These are the mu that give information the others are abberant predictions that are non informative
    mu_in_checked = mu_in(OkPredict); 
    fprintf('%d/%d mu were suppressed from the entry as those are abberant\n',length(mu_in) - length(mu_in_checked),length(mu_in));
else
    mu_in_checked = mu_in;
    OkPredict = 1:length(mu_in);
end

if nargin<4
%     dataClass = superiorfloat(mu,y);
%     muLims = realmin(dataClass).^.25; %Choice in glmfit
    Scale=20;%Best value for our data
    muLims = 1/Scale;
else
    %muLims = 1/2;%make the log function linear just around 1 (log(1)=0)
    muLims = 1/Scale;%exp(-log(Win))
end

if nargin<4
    Yhist_Flag = 0;
else
    Yhist_Flag=1;
end

if nargin<3
    ymax = 20; % For julie's data spike count are within 20ms windows!!
    %ymax = ceil(max(mu_in_checked*10));
end
if ymax>170
    ymax=170; % factorial is reaching Inf if ymax>170
end


    
N_pres = length(mu_in_checked);
y_world = 0:ymax;
[mu_unique, ~, I_mu_unique] = unique(mu_in_checked); % mu most likely has repetitions du to the trial presentations 
Nb_mu = length(mu_unique);


%% Calculate the conditional entropy (H(y/mu))
H_ymu_unique = nan(size(mu_unique)); 
Fac_y = factorial(y_world);
RescaledP=0;
P_YgivenS_all1_local = nan(length(y_world),N_pres);
P_YgivenS_all1_rescaled = nan(length(y_world),N_pres);
for mm=1:Nb_mu
    mu_local = repmat(mu_unique(mm),size(y_world));
    %P = (mu_local .^ y_world) .* exp(-mu_local) ./ Fac_y;
    
    % Calculate log of probability
    % identify indices for which ylog(y)=0 when y=0
    mu_temp = mu_local;
    mu_temp(intersect((y_world-mu_local)==0, mu_local==0)) = 1;
    %apply a minimum value for mu, cannot be equal to zero
    mu_temp = max(mu_temp,muLims);
    Log_P= y_world.*log2(mu_temp) - log2(Fac_y) - mu_local/log(2);   %threshold mu so the log calculation does not blow up when mu=0
    
    P = 2.^(Log_P);
    P_YgivenS_all1_local(:,find(I_mu_unique==mm))=repmat(P',1,length(find(I_mu_unique==mm)));

    
    % Check that sum(P) sum to 1 and rescale if not to have a more
    % realistic distribution.
    P = P ./sum(P);
    P_YgivenS_all1_rescaled(:,find(I_mu_unique==mm))=repmat(P',1,length(find(I_mu_unique==mm)));
    H_ymu_unique(mm) = sum(-P.*log2(P + (P==0)));% conditional entropy value for each unique mu (each stim) across all y of the world
end
H_ymu = sum(H_ymu_unique(I_mu_unique))/N_pres;
%H_ymu = sum(sum(-P_YgivenS_all1_rescaled.*log2(P_YgivenS_all1_rescaled+(P_YgivenS_all1_rescaled==0))))/N_pres;

%% Calculate the entropy of y based on Poisson distributions
% This calculation is not using the data as scaled per stim as used for the
% calculation of the exact entropy
if ~SameDataScaleEntropy
    P_YgivenS_all2_local = nan(length(y_world),N_pres);
    for ii=1:(ymax +1)
        yy = y_world(ii);
        y_local = repmat(yy,size(mu_in_checked));
        Fac_y = repmat(factorial(yy),size(mu_in_checked));
        if any(mu_in_checked < muLims) % This is just to be consistent with the rules regarding mu that were applied higher
            mu_temp = max(mu_in_checked,muLims);
            P_YgivenS_all2_local(ii,:) = (mu_temp.^y_local) .* exp(-mu_temp) ./ Fac_y;
        else
            P_YgivenS_all2_local(ii,:) = (mu_in_checked.^y_local) .* exp(-mu_in_checked) ./ Fac_y;
        end
    end
    
    P_y_i = sum(P_YgivenS_all2_local,2)./N_pres;
    % Check that sum(P) sum to 1 and rescale if not to have a more
    % realistic distribution.
    if sum(P_y_i)~=1
        P_y_i = P_y_i ./sum(P_y_i);
        fprintf('rescaled the distribution of y\n');
    end
    
else
    % Calculation using the same data as for the conditional entropy
    P_y_i = sum(P_YgivenS_all1_rescaled,2)./N_pres;
end
% Calculate log probability
H_y = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% conditional entropy value for each unique y of the world across all mu of the model

%% Calculate the model information
I = (H_y - H_ymu);

%% OUtput matrices of conditional probabilities for all mu in, even the one that wee aberrant by setting values to 0
P_YgivenS_all1 = zeros(length(y_world),length(mu_in));
P_YgivenS_all1(:,OkPredict) = P_YgivenS_all1_local;

P_YgivenS_all2 = zeros(length(y_world),length(mu_in));
P_YgivenS_all2(:,OkPredict) = P_YgivenS_all2_local;

%% Calculate H(y) using the actual dataset
if Yhist_Flag
    P_y = nan(ymax +1,1);
    for ii = 1:(ymax +1)
        yy = y_world(ii);
        P_y(ii) = sum(y==yy)./length(y);
    end
    H_y_histo = sum(-P_y.*log(P_y + (P_y==0)));
    I_histo = H_y_histo - H_ymu;
    I = [I I_histo];
end

    


end



