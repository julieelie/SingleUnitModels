load('DC_128ResultsAlphaLambdaDev.mat')
% Find Max min values of Deviance in the dataset
BT=25;
MaxDev=nan(length(Deviance_All)*BT,1);
MinDev=nan(length(Deviance_All)*BT,1);
LambdaMinDev = nan(length(Deviance_All)*BT,1);
dd=0;
Lambda_values=nan(length(Deviance_All)*BT*100,1);
Deviance_values = Lambda_values;
Alpha_values = Lambda_values;
vv=0;
for AA=1:length(Deviance_All)
    Deviance_local=Deviance_All{AA};
    Lambda_local=Lambda_All{AA};
    for BB=1:length(Deviance_local)
        dd=dd+1;
        MaxDev(dd)=max(Deviance_local{BB});
        MinDev(dd)=min(Deviance_local{BB});
        LambdaMinDev(dd)=log(Lambda_local{BB}(Deviance_local{BB}==min(Deviance_local{BB})));
        for ll=1:length(Deviance_local{BB})
            vv=vv+1;
            Deviance_values(vv)=Deviance_local{BB}(ll);
            Lambda_values(vv)=log10(Lambda_local{BB}(ll));
            Alpha_values(vv)=Alphas(AA);
        end
    end
end
MinDev=MinDev(1:dd);
LambdaMinDev=LambdaMinDev(1:dd);
MAXDev = max(MaxDev);
MINDev = min(MinDev);
Alpha_values=Alpha_values(1:vv);
Deviance_values=Deviance_values(1:vv);
Lambda_values=Lambda_values(1:vv);
%cubehelix_niceplot(Alpha_values,Lambda_values,Deviance_values)
figure()
plot3(Alpha_values,Lambda_values,Deviance_values,'.')
xlabel('Alpha')
ylabel('Log lambda')
zlabel('Deviance')
title(sprintf('DC_128, %d bootstraps 100 lambda values',BT))
% Noise in lambda minimum values and deviance minimumu values
N_alpha = length(Deviance_All);
Ind=1:BT:N_alpha*BT;
Mean_Dev=nan(N_alpha,1);
Dev_perAlpha=nan(BT,N_alpha);
Std_Dev=nan(N_alpha,1);
Mean_Lambda = nan(N_alpha,1);
Std_Lambda = nan(N_alpha,1);
Lambda_perAlpha=nan(BT,N_alpha);
for ii=1:5
    nn=Ind(ii);
    Mean_Dev=mean(MinDev(nn:nn+BT-1));
    Std_Dev=std(MinDev(nn:nn+BT-1));
    Dev_perAlpha(:,ii)=MinDev(nn:nn+BT-1);
    Mean_Lambda=mean(LambdaMinDev(nn:nn+BT-1));
    Std_Lambda=std(LambdaMinDev(nn:nn+BT-1));
    Lambda_perAlpha(:,ii)=LambdaMinDev(nn:nn+BT-1);
end
figure()
ss=subplot(1,2,1);
boxplot(Dev_perAlpha)
xlabel('Alpha')
ylabel('Minimum Deviance')
set(ss,'XTickLabel',[0.1:0.2:0.8 1])
ss=subplot(1,2,2);
boxplot(Lambda_perAlpha)
xlabel('Alpha')
ylabel('Best log(Lambda)')
set(ss,'XTickLabel',[0.1:0.2:0.8 1])
