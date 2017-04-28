function I=info_models_Calculus(y,mu1,mu2,Scale)
%%  Calculates the information between the prediction of two models
% y is the observed data
% mu1 are the prediction from the null model
% mu2 are the prediction from the full model
% y, mu1 and mu2 should all have the same dimensions
% Here we know the natural scale for mu: the neuron can spike between
%0 and 1 time per ms so depending on the window size Win (ms) the range
% of values is 0 to Win spikes. Make the lower limit for log(mu) as
% - the max limit by replacing zero values in mu by muLims=1/Win.
% Other idea: the natural increment for mu is 1 so let's have log(mu)
% almost linear around 0 (around mu=1) and fix muLims=1/20.
% In glmfit and here if no scale is provided, other choice because they don't know the
% natural scale of mu. Their choice keeps mu^4 from underflowing.No upper limit.
if nargin<4
%     dataClass = superiorfloat(mu,y);
%     muLims = realmin(dataClass).^.25; %Choice in glmfit
    Scale=20;%Best value for our data
    muLims = 1/Scale;
else
    %muLims = 1/2;%make the log function linear just around 1 (log(1)=0)
    muLims = 1/Scale;%exp(-log(Win))
end


% Calculate poisson probability for each model
Fac_y = factorial(y);
P1 = (mu1 .^ y) .* exp(-mu1) ./ Fac_y;
P2 = (mu2 .^ y) .* exp(-mu2) ./ Fac_y;

% Calculate log of probability
if sum((y-mu1)==0)
    Log_P1= y.*log2(mu1+(mu1==0)) - log2(Fac_y) - mu1/log(2); %ylog(y)=0 when y=0
elseif any(mu1 < muLims)
    mu_temp = max(mu1,muLims);
    Log_P1= y.*log2(mu_temp) - log2(Fac_y) - mu1/log(2);   %threshold mu so the log calculation does not blow up when mu=0
else
    Log_P1= y.*log2(mu1) -log2(Fac_y) - mu1/log(2);
end

if sum((y-mu2)==0)
    Log_P2= y.*log2(mu2+(mu2==0)) - log2(Fac_y) - mu2/log(2); %ylog(y)=0 when y=0
elseif any(mu2 < muLims)
    mu_temp = max(mu2,muLims);
    Log_P2= y.*log2(mu_temp) - log2(Fac_y) - mu2/log(2);   %threshold mu so the log calculation does not blow up when mu=0
else
    Log_P2= y.*log2(mu2) -log2(Fac_y) - mu2/log(2);
end

% Calculate the entropies
H1 = - P1.* Log_P1;
H2 = - P2.* Log_P2;

% Calculate the mutual information
I = sum(H1 - H2);

end

