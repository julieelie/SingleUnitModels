function [R2A, SSres, SSexp, SStot, ModelPredict, LL, NEC, PvalLRatio, HLRatio, NeuroRes, VOC, NbOptPC, Pvalue, Wins, NeuralResponse, STRF_time, STRF_to, STRF_fo, ModSem] = GrowingModelsRidge(Spectro, VocType, PSTH, MinWin, MaxWin, Increment, ResDelay, NeuroRes)
FIG=1; % set to 1 for debugging figures
Check=1;%set to 1 to compare ridge results with Linear model results
if nargin<8
    NeuroRes = 'mean';
end
if nargin<7
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<6
    Increment = 5; %increase the size of the spectro window with a 5ms pace
end
if nargin<5
    MaxWin = 600; %maximum values the window of analysis can reach
end
if nargin<4
    MinWin = 40; %minimum size of the window of analysis from the begining and also size of analysis of spike rate
end
Flow = 8000;%spectrograms are low passed at 8Khz for the calculations

%define the increasing size of the window of the spectrogram
Wins = MinWin:Increment:MaxWin;

% # of models to run on the data
modNum = length(Wins);

% Number of stims in the data set
NbStim = length(VocType);

% define the range of lambda values for the ridge regression to investigate
% WARNING!! WHAT RANGE OF LAMBDAS SHOULD WE TEST????
%Lambdas = 0:1e-5:1e-3;
%Lambdas = [1 10 20 50 100 10^3 10^4 10^5 10^6 10^7 10^8];


% Initialize a bunch of output variables
% R2A.Acoustic = nan(modNum,1);% adjsuted R squared
% R2A.Semantic = nan(modNum,1);
% R2A.AcSem = nan(modNum,1);
RidgeLambdas = nan(modNum,1);
RidgeB = cell(modNum,1);
% LL = R2A;%Loglikelihood
% Pvalue = R2A;%pvalue of the anova on the model
% NEC = R2A;%Number of estimated coeeficients in the model
% ModelPredict.Acoustic = cell(modNum,1);
% ModelPredict.Semantic = cell(modNum,1);
% ModelPredict.AcSem = cell(modNum,1);
% NeuralResponse = cell(modNum,1);
% STRF_time = cell(modNum,1);
% STRF_to = cell(modNum,1);
% STRF_fo = cell(modNum,1);
% ModSem = cell(modNum,1);
% NbOptPC = nan(1,modNum);
% PvalLRatio.AcAcSem = nan(modNum,1);
% PvalLRatio.SemAcSem = nan(modNum,1);
% HLRatio.AcAcSem = nan(modNum,1);
% HLRatio.SemAcSem = nan(modNum,1);
VOC = cell(modNum,1);

% SSres.Acoustic = nan(modNum,1);
% SSres.Semantic = nan(modNum,1);
% SSres.AcSem = nan(modNum,1);
% SSexp = SSres;
% SStot = nan(modNum,1);



%% Now loop through window sizes and calculate models
for mm = 1:modNum
    fprintf(1,'%d/%d models\n', mm, modNum)
    Win = Wins(mm);
    % define new dataset depending on the size of the window of the model
    % loop through the stims and only keep the Win first ms of them when
    % they are longer than Win ms or disgard
    
    duration = nan(NbStim,1);
    for ss = 1:NbStim
        duration(ss)=Spectro.to{ss}(end)*1000; %converting s in ms here
    end
    Stim_local = find(duration >= (Win+ResDelay));% here we add ResDelay because we need to get sounds with corresponding psth that go ResDelay beyond the spectrogram of size Win
    NbStim_local = length(Stim_local);
    if NbStim_local<20
        sprintf('Only %d stims long enough to run the model: no model is run with window size %dms\n', NbStim_local, Win);
        break
    end
    %NBPC = [1:9 10:5:(NbStim_local*0.8)]; % This is not a good solution since the R2A profile of most cells show stairs-like structure with increasing number of PC we need to apply a ridge regression after the PCA.
    Dt = sum((1000.*Spectro.to{Stim_local(1)})<= Win);
    Df = sum(Spectro.fo{Stim_local(1)}<= Flow);
    x= nan(NbStim_local,Df*Dt);%this matrix will contain the vectors of spectrograms for all the stims for that window size
    y = nan(NbStim_local,1);%this matrix will contain the average spike rate in spikes/ms at that precise position and for all the stims choosen for that run
    VOC{mm} = VocType(Stim_local);
    for ss = 1:NbStim_local
        dd=Stim_local(ss);
        %new spectro
        MatSpec = reshape(Spectro.spec{dd}, length(Spectro.fo{dd}), length(Spectro.to{dd}));
        FreqBelowFlow = find(Spectro.fo{dd}<=Flow);
        EndFreq = FreqBelowFlow(end);
        NFreq=length(FreqBelowFlow);
        if NFreq~=Df
            sprintf('WARNING!! Trouble with the size of the spectros for stim %d\n', dd);
        end
        TimeBelowWin = find((1000.*Spectro.to{dd})<= Win);
        EndTime = TimeBelowWin(end);
        NTime = length(TimeBelowWin);
        if NTime~=Dt
            sprintf('WARNING!! Trouble with the size of the spectros for stim %d\n', dd);
        end
        Newspectro=MatSpec(1:EndFreq,1:EndTime);
        x(ss,:)=reshape(Newspectro, 1, NFreq*NTime);
        
        % Values of max spike rate and mean spike rate within the window
        if strcmp(NeuroRes, 'max')
            y(ss) = max(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));
        elseif strcmp(NeuroRes, 'mean')
            y(ss) = mean(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', NeuroRes);
    
        end
    end
    
    
    % Take the log of the spectro and ground the output to supress -Inf
    % values
    x = 20*log10(abs(x));
    MAXI = max(max(x));
    x(find(x<(MAXI-80)))=MAXI-80;
    
    % Run Ridge and find the best parameter (Lambdas) to
    % optimize the Costfunction of the ridge. WARNING YOU DON'T WANT TO RUN
    % ON ALL DATASET!!!
    fprintf(1,'Ridge')
    %b0=ridge(y,x,Lambdas,0);
    %[b0,FitInfo] = lasso(x,y,'Lambda', Lambdas);
    [H, H0, V, W, L, P_num] = myridge(y, x);
    
    
    % Cost functions WARNING!!! YOU IGHT WANT TO CALCULATE THE COST
    % FUNCTION ON THE CROSS VALIDATION DATASET!!!
    fprintf(1,'find Lambda with the minimum cost function\n')
    NbL=length(L);
    CostF = nan(NbL,1);
    for ll=1:NbL
        %YY = y- (b0(2:end,ll) .* x + b0(1,ll)); % Check the vector sizes
        YY = y- (H0(ll) + H(ll,:)*x')'; % Check the vector sizes
        YY2 = power(YY,2);
        CostF(ll)=mean(YY2);
    end
    DerivCostF = CostF(2:end)-CostF(1:end-1);
    Lambda = L(find(CostF==min(CostF)));
    if FIG==1
        figure()
        subplot(1,2,1)
        plot(L, CostF)
        xlabel('Lambdas ridge parameter')
        ylabel('Ridge Cost Function')
        vline(Lambda);
        subplot(1,2,2)
        plot(L(1:end-1), DerivCostF)
        xlabel('Lambdas ridge parameter')
        ylabel('derivative Ridge Cost Function')
        vline(Lambda);
        pause
        if Check==1
            fprintf(1, 'Calculate PC of spectro\n');
            [COEFF,SCORE,latent,tsquare]=princomp(x,'econ');
            nPC=60;
            ds=dataset();
            for ii=1:nPC
                ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
            end
            ds.y=y;
            mdl=LinearModel.fit(ds);
            PCSTRF=mdl.Coefficients.Estimate(2:end);
            STRF=COEFF(:,1:nPC)*PCSTRF;
            STRFM=reshape(STRF,Df, Dt);
            LongestStim = find(duration==max(duration));
            Fo_Indices=find(Spectro.fo{LongestStim}<=Flow);
            To_Indices=find((1000.*Spectro.to{LongestStim})<=Win);
            fprintf(1, 'Calculating STRF using the %d first PC of the spectro\n\n\n\n', nPC);
            figure(5)
            imagesc(Spectro.to{LongestStim}(To_Indices), Spectro.fo{LongestStim}(Fo_Indices), STRFM)
            axis xy
            title(sprintf('STRF obtained with linear model with %d PC of the spectro', nPC));
            pause
        end 
        for ll=1:NbL
            figure(4)
            imagesc(reshape(H(ll,:),Df,Dt))
            axis xy
            title(sprintf('STRF lambda=%f\n',L(ll)))
            pause
        end
          
    end
    RidgeLambdas(mm) = Lambda;
    RidgeB{mm} = b0;
end
    
%     jj=0;
%     ff=0;
%     R2A_temp=zeros(1,length(NBPC));
%     ModelPredict_temp = cell(1,length(NBPC));
%     PCSTRF_temp = cell(1,length(NBPC));
%     COEFF_temp = cell(1,length(NBPC));
%     LL_temp= zeros(1,length(NBPC));
%     Pvalue_temp=zeros(1,length(NBPC));
%     NEC_temp = zeros(1,length(NBPC));
%     SSres_temp = zeros(1,length(NBPC));
%     SSexp_temp = zeros(1,length(NBPC));
%     SStot_temp = zeros(1,length(NBPC));
%     
%     
%     
%     if FIG==1
%         figure(4)
%         subplot(2,1,2)
%         plot(NBPC, R2A_temp);
%         %hold on
%         ylabel('Adjusted R2 Acoustic Model')
%         xlabel('# PC')
%     end
%     %calculate the STRF
%     PCSTRF=PCSTRF_temp{jj-1};
%     COEFF = COEFF_temp{jj-1};
%     STRF=COEFF(:,1:NBPC(jj-1))*PCSTRF;
%     STRFM=reshape(STRF,NFreq, NTime);
%     LongestStim = find(duration==max(duration));
%     Fo_Indices=find(Spectro.fo{LongestStim}<=Flow);
%     To_Indices=find((1000.*Spectro.to{LongestStim})<=Win);
%     if FIG==1
%         ff = ff+1;
%         fprintf(1, 'Calculating STRF using the %d first PC of the spectro\n\n\n\n', NBPC(jj-1));
%         figure(ff)
%         imagesc(Spectro.to{LongestStim}(To_Indices), Spectro.fo{LongestStim}(Fo_Indices), STRFM)
%         axis xy
%         title(sprintf('STRF with %d PC of the spectro', NBPC(jj-1)));
%         pause
%     end
%     
%     %Store data for the acoustic model
%     NbOptPC(mm)=NBPC(jj-1);
%     R2A.Acoustic(mm) = R2A_temp(jj-1);
%     ModelPredict.Acoustic{mm} = ModelPredict_temp{jj-1};
%     NeuralResponse{mm} = y;
%     LL.Acoustic(mm) = LL_temp(jj-1);
%     NEC.Acoustic(mm) = NEC_temp(jj-1);
%     Pvalue.Acoustic(mm) = Pvalue_temp(jj-1);
%     SSres.Acoustic(mm) = SSres_temp(jj-1);
%     SSexp.Acoustic(mm) = SSexp_temp(jj-1);
%     SStot(mm) = SStot_temp(jj-1);
%     STRF_time{mm} = STRFM;
%     STRF_to{mm}=Spectro.to{LongestStim}(To_Indices);
%     STRF_fo{mm} = Spectro.fo{LongestStim}(Fo_Indices);
%     
%     % Now do the calculations for the semantic and AcSem models
%    
%     %Model with  VocType only
%     ds2=dataset();
%     ds2.Voctype=ordinal(VOC{mm});
%     ds2.y=y;
%     
%     mdl2=LinearModel.fit(ds2);  
%     R2A.Semantic(mm)=mdl2.Rsquared.Adjusted;
%     ModelPredict.Semantic{mm}=mdl2.predict;
%     LL.Semantic(mm) = mdl2.LogLikelihood;
%     NEC.Semantic(mm) = mdl2.NumEstimatedCoefficients;
%     tbl2=anova(mdl2,'summary');
%     Pvalue.Semantic(mm)=tbl2.pValue(2);
%     SSres.Semantic(mm) = mdl2.SSE;
%     SSexp.Semantic(mm) = mdl2.SSR;
%     
%     %Model with both PC of spectro and VocType
%     ds3=dataset();
%     for ii=1:NbOptPC(mm)
%             ds3.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
%     end
%     ds3.Voctype=ordinal(VOC{mm});
%     ds3.y=y;
% 
%     mdl3=LinearModel.fit(ds3);  
%     R2A.AcSem(mm)=mdl3.Rsquared.Adjusted;
%     ModelPredict.AcSem{mm}=mdl3.predict;
%     LL.AcSem(mm) = mdl3.LogLikelihood;
%     NEC.AcSem(mm) = mdl3.NumEstimatedCoefficients;
%     tbl3=anova(mdl3,'summary');
%     Pvalue.AcSem(mm)=tbl3.pValue(2);
%     SSres.AcSem(mm) = mdl3.SSE;
%     SSexp.AcSem(mm) =  mdl3.SSR;
%         
%     % Plot of the predicted spike rate given the voctype
%     CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
%     MeanValues=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
%     ModSem{mm}=MeanValues;
%     if FIG==1
%         ff = ff +1;
%         figure(ff);
%         plot(MeanValues, 1:length(MeanValues));
%         title('Predicted SR with semantic Model');
%         set(gca,'YTickLabel', unique(VOC{mm}));
%         ff=ff+1;
%         figure(ff)
%         gscatter(y, mdl2.predict, VOC{mm}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
%         ylabel('Predicted SR /ms with semantic model')
%         xlabel('Observed SR /ms')
%         pause
%     end
%     [HLRatio.AcAcSem(mm),PvalLRatio.AcAcSem(mm),stat,cValue] = lratiotest(LL.AcSem(mm),LL.Acoustic(mm),NEC.AcSem(mm) - NEC.Acoustic(mm));
%     [HLRatio.SemAcSem(mm),PvalLRatio.SemAcSem(mm),stat,cValue] = lratiotest(LL.AcSem(mm),LL.Semantic(mm),NEC.AcSem(mm) - NEC.Semantic(mm)); 
end 

