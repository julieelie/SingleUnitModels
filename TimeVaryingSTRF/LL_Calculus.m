function LL=LL_Calculus(y,mu,Scale,SUM)
if nargin<2
    mu=y;
end
% Here we know the natural scale for mu: the neuron can spike between
%0 and 1 time per ms so depending on the window size Win (ms) the range
% of values is 0 to Win spikes. Make the lower limit for log(mu) as
% - the max limit by replacing zero values in mu by muLims=1/Win.
% Other idea: the natural increment for mu is 1 so let's have log(mu)
% almost linear around 0 (around mu=1) and fix muLims=1/4.
% In glmfit and here if no scale is provided, other choice because they don't know the
% natural scale of mu. Their choice keeps mu^4 from underflowing.No upper limit.
if nargin<3
%     dataClass = superiorfloat(mu,y);
%     muLims = realmin(dataClass).^.25; %Choice in glmfit
    Scale=20;%Best value for our data according to graph that can be obtained in the commented section below
    muLims = 1/Scale;
else
    %muLims = 1/2;%make the log function linear just around 1 (log(1)=0)
    muLims = 1/Scale;%exp(-log(Win))
end

if nargin<4
    SUM=1;%By default, sum the likelihood values
end

if sum((y-mu)==0)
    LL=y .* log(mu+(mu==0))-log(factorial(y))-mu; %ylog(y)=0 when y=0
elseif any(mu < muLims(1))
    mu2 = max(mu,muLims(1));
    LL= y.*log(mu2) - log(factorial(y)) - mu;   %threshold mu so the log calculation does not blow up when mu=0
else
    LL= y.*log(mu) -log(factorial(y)) - mu;
end
if SUM
    LL=sum(LL);
end

%% CODE TO FIND muLims
% %% Investigating how to set a floor value for mu in the LL Calculation
% M=0:5;%Hypothetical spike count value for mu
% R=0:2:10;% Hypothetical spike count for y
% SCALE=2:2:20;%Scales
% L=cell(length(SCALE),1);
% for ss=1:length(SCALE)
%     L{ss}=sprintf('1/%s',num2str(SCALE(ss)));
% end
% figure()
% for ii=1:length(R)
%     Res=repmat(R(ii),1,length(M));
%     subplot(2,3,ii)
%     hold on
%     for ss=1:length(SCALE)
%         LL=LL_Calculus(Res,M,SCALE(ss),0);
%         plot(M,LL,'.-')
%         pause(1)
%     end
%     xlabel('Mu spike count')
%     ylabel('LogLikelihood')
%     title(sprintf('R=%d spikes\ncolors=floor value for mu',R(ii)))
%     legend(L)
%     hold off
%     
% end

end

