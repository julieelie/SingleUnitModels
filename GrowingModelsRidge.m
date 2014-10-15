function [R2, SSres, SSexp, SStot, ModelPredict, LL, NEC, PvalLRatio, HLRatio, NeuroRes, VOC, NbOptPC, Pvalue, Wins, NeuralResponse, STRF_time, STRF_to, STRF_fo, SemModel] = GrowingModelsRidge(Spectro, VocType, PSTH, Emitter, MinWin, MaxWin, Increment, ResDelay, NeuroRes)
FIG=1; % set to 1 for some debugging figures 2 for all debugging figures
Check=1;%set to 1 to compare ridge results with Linear model results
BootstrapSTRF=100; %set the number of time you want to estimate the STRF at each time window
if nargin<9
    NeuroRes = 'mean';
end
if nargin<8
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<7
    Increment = 5; %increase the size of the spectro window with a 5ms pace
end
if nargin<6
    MaxWin = 600; %maximum values the window of analysis can reach
end
if nargin<5
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
%Lambdas = 0:1e-5:1e-3;
%Lambdas = [1 10 20 50 100 10^3 10^4 10^5 10^6 10^7 10^8];


% Initialize a bunch of output variables
R2.Acoustic = nan(modNum,2);% R squared mean first column, sd second column
R2.Semantic = nan(modNum,2);
R2.AcSem = nan(modNum,2);
SSE.Acoustic = zeros(modNum,1);
SSE.Semantic = zeros(modNum,1);
SST.Acoustic = zeros(modNum,1);
SST.Semantic = zeros(modNum,1);
SSEdim.Acoustic = cell(modNum,1);
SSEdim.AcSem = cell(modNum,1);
SSTdim.Acoustic = cell(modNum,1);
SSTdim.AcSem = cell(modNum,1);
PropVal.Acoustic = nan(modNum,2);
PropVal.Sem = nan(modNum,2);
PropVal.AcSem = nan(modNum,2);
NDim.AcSem = nan(modNum,2);
Model.Acoustic.RidgeLambdas = nan(modNum,1);
Model.Acoustic.RidgeH = cell(modNum,1);
Model.Acoustic.RidgeH0 = cell(modNum,1);
Model.Acoustic.RidgeHPrime = cell(modNum,1);
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
Model.Sem.Coefficients = cell(modNum,1);
Model.Sem.CoefficientsNames = cell(modNum,1);
% NbOptPC = nan(1,modNum);
% PvalLRatio.AcAcSem = nan(modNum,1);
% PvalLRatio.SemAcSem = nan(modNum,1);
% HLRatio.AcAcSem = nan(modNum,1);
% HLRatio.SemAcSem = nan(modNum,1);
VOC = cell(modNum,1);
Model.AcSem.Coefficients = cell(modNum,1);
Model.AcSem.CoefficientsNames = cell(modNum,1);



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
    
    
    
    % Calculate STRF using ridge regression and cross-validation
    R2STRF=nan(BootstrapSTRF,1);
    ValProp=nan(BootstrapSTRF,1);
    LambdaIndex=nan(BootstrapSTRF,1);
    CostFSum = 0;
    MaxErrorSum=0;
    RidgeH=cell(BootstrapSTRF,1);
    RidgeH0=cell(BootstrapSTRF,1);
    RidgeHPrime=cell(BootstrapSTRF,1);
    SSEprog_sum = zeros(NbStim_local,1);
    SSTprog_sum = zeros(NbStim_local,1);
    P_num_min = NbStim_local;
    
    for bb=1:BootstrapSTRF
        % Construct a testing dataset and a validating dataset
        % remove one emitter per call category, if one category contain only
        % one emitter, don't use it in the model
        [ValSet, TrainSet] = create_cross_validation_sets(VOC{mm}, Emitter.Ename(Stim_local));
        ValProp(bb) = length(ValSet)./(length(ValSet)+length(TrainSet));

        % Run Ridge and find the best parameter (Lambdas) to
        % optimize the Costfunction of the ridge.
        fprintf(1,'Acoustic Model Ridge %d/%d\n', bb, BootstrapSTRF);
        if FIG>1
            pause
        end
        %b0=ridge(y,x,Lambdas,0);
        %[b0,FitInfo] = lasso(x,y,'Lambda', Lambdas);
        [RidgeLocalH, RidgeLocalH0, V, W, L, P_num,RidgeLocalHPrime,RidgeLocalH0Prime] = myridge(y(TrainSet), x(TrainSet,:));

        % Cost functions
        fprintf(1,'find Lambda with the minimum cost function\n')
        NbL=length(L);
        CostF = nan(NbL,1);
        MaxError = sum((y(ValSet)-mean(y(TrainSet))).^2);
        for ll=1:NbL
            %YY = y- (b0(2:end,ll) .* x + b0(1,ll)); % Check the vector sizes
            YY = y(ValSet)- (RidgeLocalH0(ll) + RidgeLocalH(ll,:)*x(ValSet,:)')'; % Check the vector sizes
            YY2 = YY.^2;
            CostF(ll)=sum(YY2);
        end
        MaxErrorSum = MaxErrorSum + MaxError;
        CostFSum = CostFSum + CostF;
        ll_min=find(CostF==min(CostF));
        ll_min=ll_min(end);
        LambdaIndex(bb)=ll_min;
        RidgeH{bb} = RidgeLocalH(ll_min,:);
        RidgeH0{bb} = RidgeLocalH0(ll_min);
        RidgeHPrime{bb} = RidgeLocalHPrime(ll_min);
        
        % Calculate R2 and errors with all dimensions of the best ridge model
        YY = y(ValSet)- (RidgeLocalH0(ll_min) + RidgeLocalH(ll_min,:)*x(ValSet,:)')';
        YY2 = YY.^2;
        SSE.Acoustic(mm) = SSE.Acoustic(mm) + sum(YY2);
        SST.Acoustic(mm) = SST.Acoustic(mm) + sum((y(ValSet)-mean(y(TrainSet))).^2);
        R2STRF(bb)= 1 - sum(YY2)./sum((y(ValSet)-mean(y(TrainSet))).^2);
        
        % Calculate the errror of the best ridge model as we increase the
        % number of dimensions
        WL =1./ (diag(W).^2 + L(ll_min));
        WL_h = diag(WL).^(0.5);
        x_prime =  (WL_h * V'* x(:,:)')';
        SSEprog=nan(P_num,1);
        for dimd=1:P_num
            %YYPrime = y(ValSet) - (mean(y(TrainSet))+RidgeLocalHPrime(ll_min,1:dimd)*x_prime(:,1:dimd)')';
            YYPrime = y(ValSet) - (mean(y(TrainSet))+RidgeLocalHPrime(ll_min,1:dimd)*x_prime(ValSet,1:dimd)' - RidgeLocalHPrime(ll_min,1:dimd)*mean(x_prime(:,1:dimd),1)')';
            SSEprog(dimd) = sum(YYPrime.^2);
        end
        SSEprog_sum = SSEprog_sum + [SSEprog; zeros(NbStim_local-P_num,1)];
        SSTprog=repmat(sum((y(ValSet)-mean(y(TrainSet))).^2),NbStim_local,1);
        SSTprog_sum = SSTprog_sum + SSTprog;
        P_num_min = min(P_num_min, P_num);
        
        % Plot Lambda values and STRF
        if FIG>1
            figure(3)
            plot(log(L), CostF)
            xlabel('log of Lambdas ridge parameter')
            ylabel('Ridge Cost Function')
            title(sprintf('training size %d validating size %d',length(TrainSet), length(ValSet)))
            vline(log(L(ll_min)));
            pause(1)
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
                figure(4)
                imagesc(Spectro.to{LongestStim}(To_Indices), Spectro.fo{LongestStim}(Fo_Indices), STRFM)
                axis xy
                title(sprintf('STRF obtained with linear model with %d PC of the spectro', nPC));
                pause(1)
            end 
            for ll=1:NbL
                if ll==find(CostF==min(CostF))
                    figure(5)
                    imagesc(reshape(RidgeLocalH(ll,:),Df,Dt))
                    axis xy
                    title(sprintf('THIS IS THE ONE!!!\nSTRF log of lambda=%f\n R2A=%f',log(L(ll)),R2STRF(bb)))
                else
                    figure(6)
                    imagesc(reshape(RidgeLocalH(ll,:),Df,Dt))
                    axis xy
                    title(sprintf('STRF log of lambda=%f\n',log(L(ll))))
                end
                pause(1)
            end

        end
    end
    ll_summin = find(CostFSum==min(CostFSum));
    ll_summin=ll_summin(end);
    PropVal.Acoustic(mm,1) = mean(ValProp);
    PropVal.Acoustic(mm,2) = std(ValProp);
    R2.Acoustic(mm,1) = mean(R2STRF);
    R2.Acoustic(mm,2) = std(R2STRF);
    SSEdim.Acoustic{mm} = SSEprog_sum(1:P_num_min);
    SSTdim.Acoustic{mm} = SSTprog_sum(1:P_num_min);
    
    if FIG>0
        fprintf(1,'R2A= %f +/- %f\n',R2.Acoustic(mm,1), R2.Acoustic(mm,2));
        figure(7)
        subplot(1,2,1)
        hist(LambdaIndex, max(LambdaIndex))
        subplot(1,2,2)
        plot(log(L),CostFSum)
        hline(MaxErrorSum);
        xlabel('log of Lambdas ridge parameter')
        ylabel('Sum of Ridge Cost Functions')
        vline(log(L(ll_summin)));
    end
    
    % Use the average Lambda to calculate the optimal STRF on all dataset
    [Model.Acoustic.RidgeH{mm}, Model.Acoustic.RidgeH0{mm}, V, W, L, P_num, Model.Acoustic.RidgeHPrime{mm}] = myridge(y, x,'Lambda Indices',ll_summin);
    Model.Acoustic.RidgeLambdas(mm) = L;
    if FIG>0
        figure(8)
        subplot(1,3,1)
        imagesc(reshape(Model.Acoustic.RidgeH{mm},Df,Dt))
        axis xy
        title(sprintf('OPTIMAL STRF calculated on whole dataset\nlog of lambda=%f\n R2A=%f+/-%f calculated with a cross vaidation\nincluding %f+/-%f %% of the data',log(L),R2.Acoustic(mm,1), R2.Acoustic(mm,2), PropVal.Acoustic(mm,1), PropVal.Acoustic(mm,2)))
    end
    
    
    % Calculate the semantic model using a linear model and
    % cross-validation
    R2Sem=nan(BootstrapSTRF,1);
    ValPropSem=nan(BootstrapSTRF,1);
    
    for bb=1:BootstrapSTRF
        % Construct a testing dataset and a validating dataset
        % remove one emitter per call category, if one category contain only
        % one emitter, don't use it in the model
        [ValSet, TrainSet] = create_cross_validation_sets(VOC{mm}, Emitter.Ename(Stim_local));
        ValPropSem(bb) = length(ValSet)./(length(ValSet)+length(TrainSet));
        
        % Linear Model with training set
        fprintf(1,'Semantic Model %d/%d\n', bb, BootstrapSTRF);
        ds2=dataset();
        ds2.Voctype=VOC{mm}(TrainSet);
        ds2.y=y(TrainSet);
        mdl2=fitlm(ds2,'CategoricalVars',1);
        % this is equivalent to:
        %mdl2=LinearModel.fit(ds2, 'CategoricalVars',1);
        %mdl2=fitlm(ds2,'CategoricalVars',1);
        %mdl2=LinearModel.fit(ds2);
        
        % R2 calculation with validating set
%         CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
%         Predicted_y=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
%         Predicted_y_voctype=unique(VOC{mm}(TrainSet));
%         Predicted_yval=nan(length(ValSet),1);
%         YY2_sum=0;
%         for vv=1:length(ValSet)
%             Val = ValSet(vv);
%             Val_type = VOC{mm}(Val);
%             YY2_sum = YY2_sum + (y(Val) - Predicted_y(find(strcmp(Predicted_y_voctype,Val_type)))).^2;
%             Predicted_yval(vv)=Predicted_y(find(strcmp(Predicted_y_voctype,Val_type)));
%         end
        % Instead of for loop above, predict the response using predict or feval function
        Predicted_yval = feval(mdl2,VOC{mm}(ValSet));
        SSE.Semantic(mm) = SSE.Semantic(mm) + sum((y(ValSet)-Predicted_yval).^2);
        SST.Semantic(mm) = SST.Semantic(mm) + sum((y(ValSet)-mean(y(TrainSet))).^2);
        R2Sem(bb)= 1 - sum((y(ValSet)-Predicted_yval).^2)./sum((y(ValSet)-mean(y(TrainSet))).^2);
    end
    
    R2.Semantic(mm,1)=mean(R2Sem);
    R2.Semantic(mm,2)=std(R2Sem);
    PropVal.Sem(mm,1) = mean(ValPropSem);
    PropVal.Sem(mm,2) = std(ValPropSem);
    
    % Optimal Semantic model using all data set
    ds2=dataset();
    ds2.Voctype=ordinal(VOC{mm});
    ds2.y=y;
    mdl2=LinearModel.fit(ds2);
    CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
    Predicted_y=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
    Model.Sem.Coefficients{mm}=Predicted_y;
    Model.Sem.CoefficientsNames{mm}=unique(VOC{mm});
    if FIG>0
        figure(8)
        ss=subplot(1,3,2);
        plot(Model.Sem.Coefficients{mm}, 1:length(Model.Sem.Coefficients{mm}));
        title('Predicted SR with Semantic Model');
        set(ss,'YTickLabel', Model.Sem.CoefficientsNames{mm});
    end
    
    
    % Calculate the model combining both Ridge STRF and Semantic labels
    % First rescale and rotate the spectrograms of the stimuli in the ridge
    % space
    WL =1./ (diag(W).^2 + Model.Acoustic.RidgeLambdas(mm));
    WL_h = diag(WL).^(0.5);
    x_ridge =  (WL_h * V'* x')';
    
    % Calculate the Full model using a linear model in the ridge space and
    % cross-validation
    R2AcSem=nan(BootstrapSTRF,1);
    ValPropAcSem=nan(BootstrapSTRF,1);
    NumSignifEdim=nan(BootstrapSTRF,1);
    SSEprog_sumAcSem = zeros(P_num,1);
    SSTprog_sumAcSem = zeros(P_num,1);
    
    for bb=1:BootstrapSTRF
        % Construct a testing dataset and a validating dataset
        % remove one emitter per call category, if one category contain only
        % one emitter, don't use it in the model
        [ValSet, TrainSet] = create_cross_validation_sets(VOC{mm}, Emitter.Ename(Stim_local));
        ValPropAcSem(bb) = length(ValSet)./(length(ValSet)+length(TrainSet));
        
        % Linear Model with training set
        fprintf(1,'Full Model %d/%d\n', bb, BootstrapSTRF);
%         ds3=dataset();
%         for ii=1:size(x_ridge, 2)
%             ds3.(sprintf('EIGENVAL%d',ii)) = x_ridge(TrainSet,ii);
%         end
%         ds3.Voctype=ordinal(VOC{mm}(TrainSet));
%         ds3.y=y(TrainSet);
%         mdl3=LinearModel.fit(ds3);
%         tbl3=anova(mdl3,'summary');
%         Pvalue.Semantic(mm)=tbl3.pValue(2);
        
        % Progressive Linear model with training set
        %[Coefficients, CoefficientsNames, NumSignifEdim(bb)] = myProgressiveLinearModel(y(TrainSet), x_ridge(TrainSet,:),VOC{mm}(TrainSet),0);
        ValidatingSet.x=x_ridge(ValSet,:);
        ValidatingSet.y=y(ValSet,:);
        ValidatingSet.Cat = VOC{mm}(ValSet);
        [Coefficients, CoefficientsNames, NumSignifEdim(bb),SSEAcSem, SSTAcSem] = myProgressiveLinearModel(y(TrainSet), x_ridge(TrainSet,:),VOC{mm}(TrainSet),1,size(x_ridge,2),'ValidatingSet',ValidatingSet);
        SSEprog_sumAcSem = SSEprog_sumAcSem + [SSEAcSem; zeros(P_num - length(SSEAcSem),1)];
        SSTprog_sumAcSem = SSTprog_sumAcSem + [SSTAcSem; zeros(P_num - length(SSTAcSem),1)];
        
        
        % Calculate predicted reponse of the validating dataset with this
        % model and R2
        x_ValSet = x_ridge(ValSet,1:NumSignifEdim(bb));
        Ind_strf = zeros(size(CoefficientsNames,1));
        for nn = 1:NumSignifEdim(bb)
            Ind_strf = Ind_strf + strcmp(CoefficientsNames, sprintf('DIM%d',nn));
        end
        y_predicted_STRF = x_ValSet * Coefficients(1,find(Ind_strf));
        y_predicted_Voc = nan(length(ValSet),1);
        y_predicted_Int = nan(length(ValSet),1);
        for ts=1:length(ValSet)
            Ind_vocType = strcmp(CoefficientsNames,sprintf('Cattype_%s',VOC{mm}{ValSet(ts)}));
            Ind_int = zeros(size(CoefficientsNames));
            for nn = 1:NumSignifEdim(bb)
                Ind_int = Ind_int + strcmp(CoefficientsNames, sprintf('DIM%d:Cattype_%s',nn,VOC{mm}{ValSet(ts)}));
            end
            y_predicted_Voc(ts) = Coefficients(1) + sum(Coefficients.*Ind_vocType);
            if sum(Ind_int)>0
                y_predicted_Int(ts) = x_ValSet(ts,:)*Coefficients(find(Ind_int));
            else
                y_predicted_Int(ts) = 0;
            end
        end
        y_predicted = y_predicted_STRF + y_predicted_Voc + y_predicted_Int;
        R2AcSem(bb)= 1 - sum((y(ValSet)-y_predicted).^2)./sum((y(ValSet)-mean(y(ValSet))).^2);
    end
    R2.AcSem(mm,1)=mean(R2AcSem);
    R2.AcSem(mm,2)=std(R2AcSem);
    PropVal.AcSem(mm,1) = mean(ValPropAcSem);
    PropVal.AcSem(mm,2) = std(ValPropAcSem);
    NDim.AcSem(mm,1)=mean(NumSignifEdim);
    NDim.AcSem(mm,2)=std(NumSignifEdim);
    SSEdim.AcSem{mm} = SSEprog_sumAcSem;
    SSTdim.AcSem{mm} = SSTprog_sumAcSem;
    
    figure(10)
    subplot(1,2,1)
    plot(SSEdim.Acoustic{mm})
    hold on
    plot(SSEdim.AcSem{mm},'r')
    ylabel('SSE')
    xlabel('Number of Acoustic Ridge Dimensions in the models');
    legend('SSE Acoustic Model','SSE Acoustic + Semantic Model');
    hold off
    subplot(1,2,2)
    plot(1-SSEdim.Acoustic{mm}./SSTdim.Acoustic{mm})
    hold on
    plot(1-SSEdim.AcSem{mm}./SSTdim.AcSem{mm},'r');
    ylabel('R2')
    xlabel('Number of Acoustic Ridge Dimensions in the models');
    legend('R2 Acoustic Model','R2 Acoustic + Semantic Model');
    hold off
    
    
    % Optimal Full model using all data set
    [Model.AcSem.Coefficients{mm}, Model.AcSem.CoefficientsNames{mm}] = myProgressiveLinearModel(y, x_ridge,VOC{mm},0,ceil(NDim.AcSem(mm,1)));
    
    
    % Prediction of the Acoustic Model using only the acoustic parameters in the eigenspace
    % used in the full model, the NDim.AcSem dimensions
    R2STRFbis = nan(BootstrapSTRF,1);
    for bb=1:BootstrapSTRF
        % Construct a testing dataset and a validating dataset
        % remove one emitter per call category, if one category contain only
        % one emitter, don't use it in the model
        [ValSet, TrainSet] = create_cross_validation_sets(VOC{mm}, Emitter.Ename(Stim_local));
        ValProp(bb) = length(ValSet)./(length(ValSet)+length(TrainSet));
        
        % Calculate R2
        YY = y(ValSet)- (RidgeH0{bb} + RidgeH{bb}(1:NDim.AcSem(mm,1))*x(ValSet,(1:NDim.AcSem(mm,1)))')';
        YY2 = YY.^2;
        R2STRFbis(bb)= 1 - sum(YY2)./sum((y(ValSet)-mean(y(ValSet))).^2);
    end
    
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

