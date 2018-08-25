function [Deviance, ParamModel, VOC, Wins] = GrowingModelsRidgeglmLambdaModBeta(Spectro, VocType, PSTH, Trials, Emitter, ParamWin, ParamModel, Cellname)
%resultsDirectory='/auto/tdrive/julie/NeuralData/SemanticGLMModel';
%resultsDirectory=pwd;
resultsDirectory='/Users/elie/Documents/CODE/SingleUnitModels/Resultsglmlasso/glmLasso_3Models_fitPerTrial/OptimalDevianceOffsetModel';

NbTrialStim=10;
SWITCH.FIG=2; % set to 1 for some debugging figures 2 for all debugging figures 3 for extra figures
SWITCH.Check=1;%set to 1 to compare ridge results with Linear model results
SWITCH.ZC=1;%Set to 0 if you don't want to z-score the values of the spectrograms before performing lasso glm
SWITCH.FitPerTrial = 1; %Set to 0 if you want to fit the model using the average response per stim





if nargin<7
    ParamModel.NeuroRes = 'mean';
    ParamModel.DISTR='poisson';%'normal'
    ParamModel.LINK='log';
    ParamModel.LINK_AR='identity';
    ParamModel.LAMBDA=53.3;%set to [] if you want to explore lamba values %53.3
    ParamModel.NUMLAMBDA=25;
    ParamModel.LAMBDARATIO=1e-4;
    ParamModel.BootstrapSTRF=15;
    ParamModel.ModelChoice=[1 1 0 1 1];% This logic vector indicates which models should
...be run (Acoustic, Semantic, AcSem without Offset,AcSem Semantic Offset, 
    ...AcSem Accoustic Offset) Note that the last two models need the 
    ...calculation of the first two models
end

if nargin<6
    ParamWin.MinWin = 40; %minimum size of the window of analysis from the begining and also size of analysis of spike rate
    ParamWin.MaxWin = 160; %maximum values the window of analysis can reach
    ParamWin.Increment = 40; %increase the size of the spectro window with a 5ms pace
    ParamWin.ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
Flow = 8000;%spectrograms are low passed at 8Khz for the calculations

%define the increasing size of the window of the spectrogram
Wins = ParamWin.MinWin:ParamWin.Increment:ParamWin.MaxWin;
%Wins=160;

% # of models to run on the data
modNum = length(Wins);

% Number of stims in the data set
NbStim = length(VocType);

% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
Alphas=[0.001 0.01 0.1 1]; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

%% Initialize a bunch of output variables
if ParamModel.ModelChoice(1)
    Deviance.Acoustic.bestvalues = cell(modNum);
    Deviance.Acoustic.DF = cell(modNum);
    Deviance.Acoustic.values = cell(modNum,length(Alphas));
    Deviance.Acoustic.lambda = cell(modNum,length(Alphas));
    Deviance.Acoustic.FitIndex=cell(modNum);
    Deviance.Acoustic.FitIndexAR=cell(modNum);
    Deviance.Acoustic.FitIndexARVal=cell(modNum);
    LL.Acoustic.values=cell(modNum,length(Alphas));
    LL.Acoustic.bestvalues= cell(modNum);
    Model.Acoustic.Lambdas = nan(modNum,length(Alphas));
    Model.Acoustic.B = cell(modNum,length(Alphas));
    Model.Acoustic.B0 = cell(modNum,length(Alphas));
    Model.Acoustic.ypredictVal = cell(modNum,length(Alphas));
    Model.Acoustic.Deviance=nan(modNum,length(Alphas));
    Model.Acoustic.DF=nan(modNum,length(Alphas));
end
if ParamModel.ModelChoice(2)
    Deviance.Sem.bestvalues = cell(modNum);
    Deviance.Sem.DF = cell(modNum);
    Deviance.Sem.FitIndex=cell(modNum);
    Deviance.Sem.FitIndexAR=cell(modNum);
    Deviance.Sem.FitIndexARVal=cell(modNum);
    LL.Sem.bestvalues=cell(modNum);
    Model.Semantic.B = cell(modNum,length(Alphas));
    Model.Semantic.B0 = cell(modNum,length(Alphas));
    Model.Semantic.ypredictVal = cell(modNum,1);
    Model.Semantic.Deviance=nan(modNum,length(Alphas));
    Model.Semantic.DF=nan(modNum,length(Alphas));
end
if ParamModel.ModelChoice(3)
    Deviance.AcSem.bestvalues = cell(modNum);
    Deviance.AcSem.DF = cell(modNum);
    Deviance.AcSem.values = cell(modNum,length(Alphas));
    Deviance.AcSem.lambda = cell(modNum,length(Alphas));
    Deviance.AcSem.FitIndex=cell(modNum);
    Deviance.AcSem.FitIndexAR=cell(modNum);
    Deviance.AcSem.FitIndexARVal=cell(modNum);
    LL.AcSem.values=cell(modNum,length(Alphas));
    LL.AcSem.bestvalues= cell(modNum);
    Model.AcSem.Lambdas = nan(modNum,length(Alphas));
    Model.AcSem.Bspectro = cell(modNum,length(Alphas));
    Model.AcSem.Bsem = cell(modNum,length(Alphas));
    Model.AcSem.B0 = cell(modNum,length(Alphas));
    Model.AcSem.ypredictVal = cell(modNum,length(Alphas));
    Model.AcSem.Deviance=nan(modNum,length(Alphas));
    Model.AcSem.DF=nan(modNum,length(Alphas));
end

if ParamModel.ModelChoice(4) && ParamModel.ModelChoice(1)
    Deviance.AcSemAc.bestvalues = cell(modNum);
    Deviance.AcSemAc.DF = cell(modNum);
    Deviance.AcSemAc.values = cell(modNum,length(Alphas));
    Deviance.AcSemAc.lambda = cell(modNum,length(Alphas));
    Deviance.AcSemAc.FitIndex=cell(modNum);
    Deviance.AcSemAc.FitIndexAR=cell(modNum);
    Deviance.AcSemAc.FitIndexARVal=cell(modNum);
    LL.AcSemAc.values=cell(modNum,length(Alphas));
    LL.AcSemAc.bestvalues= cell(modNum);
elseif ParamModel.ModelChoice(4) && ~ParamModel.ModelChoice(1)
    fprintf('WARNING: No way to calculate the AcSem model with Acoustic Offset if you do not calculate Acoustic model')
end

if ParamModel.ModelChoice(5) && ParamModel.ModelChoice(2)
    Deviance.AcSemSem.bestvalues = cell(modNum);
    Deviance.AcSemSem.DF = cell(modNum);
    Deviance.AcSemSem.values = cell(modNum,length(Alphas));
    Deviance.AcSemSem.lambda = cell(modNum,length(Alphas));
    Deviance.AcSemSem.FitIndex=cell(modNum);
    Deviance.AcSemSem.FitIndexAR=cell(modNum);
    Deviance.AcSemSem.FitIndexARVal=cell(modNum);
    LL.AcSemSem.values=cell(modNum,length(Alphas));
    LL.AcSemSem.bestvalues= cell(modNum);
elseif ParamModel.ModelChoice(5) && ~ParamModel.ModelChoice(2)
    fprintf('WARNING: No way to calculate the AcSem model with Semantic Offset if you do not calculate semantic model')
end

LL.Ceiling.bestvalues=cell(modNum);
LL.Floor.bestvalues=cell(modNum);
LL.AutoRegressive.bestvalues=cell(modNum);


PropVal.mean = nan(modNum,length(Alphas));
PropVal.std = nan(modNum,length(Alphas));
PropVal.values = cell(modNum);

Model.TickSTRFspectro.to = cell(modNum,1);
Model.TickSTRFspectro.fo = cell(modNum,1);
Model.MeanSpectroStim = cell(modNum,1);
Model.yobsCat = cell(modNum,length(Alphas));
Model.yobsVal = cell(modNum,length(Alphas));

                Model.Acoustic.SSTRresiduals = cell(modNum,length(Alphas));
                Model.Semantic.SSTRresiduals = cell(modNum,length(Alphas));
                Model.AcSem.SSTRresiduals = cell(modNum,length(Alphas));

Data.y_wholeset = cell(modNum,1);
Data.x_wholeset = cell(modNum,1);
Data.x_std_wholeset = cell(modNum,1);
Data.x_mean_wholeset = cell(modNum,1);
Data.X_voc_wholeset = cell(modNum,1);


VOC = cell(modNum,1);


% Initialize datasets to be used as passed information by the
% auto-regressive model
PreviousWin.y_dev_old=[];
PreviousWin.Stim_local_old=[];

%% Now loop through window sizes and calculate models
for mm = 1:modNum
    fprintf(1,'%d/%d models\n', mm, modNum);
    Win = Wins(mm);
    %% define new dataset depending on the size of the window of the model
    % loop through the stims and only keep the Win first ms of them when
    % they are longer than Win ms or disgard
    
    duration = nan(NbStim,1);
    for ss = 1:NbStim
        duration(ss)=Spectro.to{ss}(end)*1000; %converting s in ms here
    end
    Stim_local = find(duration >= (Win+ParamWin.ResDelay));% here we add ResDelay because we need to get sounds with corresponding psth that go ResDelay beyond the spectrogram of size Win
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
    if SWITCH.FitPerTrial
        y_dev=cell(NbStim_local,1);
        ymean=nan(NbStim_local,1);
        yvar=ymean;
    end
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
        
        % Values of max spike rate(y), mean spike rate (y) and exact number of spike per trial (y_dev) within the window
        % Check the distribution of responses (Gaussian or Poisson) for each stim
        if strcmp(ParamModel.NeuroRes, 'max')
            y(ss) = max(PSTH{dd}((Win-Param.MinWin+ParamWin.ResDelay):(Win+ParamWin.ResDelay)));
        elseif strcmp(ParamModel.NeuroRes, 'mean')
            y(ss) = mean(PSTH{dd}((Win-Param.MinWin+ParamWin.ResDelay):(Win+ParamWin.ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
            if SWITCH.FitPerTrial
                y_dev{ss}=nan(length(Trials{dd}),1);
                for tt=1:length(Trials{dd})
                    y_dev{ss}(tt)=sum((Trials{dd}{tt}>(Win-Param.MinWin+ParamWin.ResDelay)).*(Trials{dd}{tt}<(Win+ParamWin.ResDelay)));
                end
                ymean(ss)=mean(y_dev{ss});
                yvar(ss)=var(y_dev{ss});
            end
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', ParamModel.NeuroRes);
    
        end
    end
    
%     % Investigate how poisson or gaussian neural responses are for this
%     % neuron
%     MAX=max(max(yvar),max(ymean));
%     figure()
%     plot(yvar,ymean,'.', 'MarkerSize',10)
%     ylabel('Variance spike counts per stim')
%     xlabel('mean spike count per stim')
%     hold on
%     line([0 MAX], [0 MAX]);
%     R2 = 1- sum(power(yvar-ymean,2))/sum(power(yvar-repmat(mean(yvar),length(yvar),1),2));
%     fprintf(1,'Proportion of stims which response have both a normal and Poisson distribution: %f\n',sum(DistrPoisson(:,1).*DistrNormal(:,1))./NbStim_local);
%     fprintf(1,'Proportion of stims which responses are Poisson only:%f\n',sum((DistrPoisson(:,1)-DistrNormal(:,1))==1)./NbStim_local);
%     fprintf(1,'Proportion of stims which responses are Normal only:%f\n',sum((DistrNormal(:,1)-DistrPoisson(:,1))==1)./NbStim_local);
%     fprintf(1,'Proportion of stims which response have both a normal and Poisson distribution: %f\n',sum((DistrPoisson(:,1)+DistrNormal(:,1))==0)./NbStim_local);
%     
    % Take the log of the spectro and ground the output to supress -Inf
    % values
    x = 20*log10(abs(x));
    MAXI = max(max(x));
    x(x<(MAXI-80))=MAXI-80;
    
    % Calculate the z-score of x that will be used to construct final strf
    % with complete dataset
    x_mean = mean(x);
    x_std = std(x);
    x_zscore = (x-repmat(x_mean,size(x,1),1))./repmat(x_std,size(x,1),1);
    if sum(x_std==0)
        x_zscore(:,x_std==0)=0;
    end
    Data.x_std_wholeset{mm}=x_std;
    Data.x_mean_wholeset{mm}=x_mean;
    
    % Store the average spectro and ticks of the spectro and STRF
    Model.MeanSpectroStim{mm} = reshape(x_mean,Df,Dt);
    LongestStim = find(duration==max(duration));
    Fo_Indices=find(Spectro.fo{LongestStim}<=Flow);
    To_Indices=find((1000.*Spectro.to{LongestStim})<=Win);
    Model.TickSTRFspectro.to{mm} = Spectro.to{LongestStim}(To_Indices);
    Model.TickSTRFspectro.fo{mm} = Spectro.fo{LongestStim}(Fo_Indices);
    
    % Code the categorical data about call
    UVOC = unique(VOC{mm});
    X_voc = zeros(length(VOC{mm}), length(UVOC)-1);% Note that all zeros correspond to Ag call
    for vv=1:length(VOC{mm})
        for uv=2:length(UVOC)
            if strcmp(VOC{mm}(vv), UVOC{uv})
                X_voc(vv,uv-1)=1;
                break
            end
        end
    end
    
     % Get the whole dataset ready for trial fit
    if SWITCH.FitPerTrial
        y_wholeset = nan(length(y_dev)* NbTrialStim,1);
        y_wholeset_bestGuess = nan(length(y_dev)* NbTrialStim,1);
        x_wholeset = nan(length(y_dev)* NbTrialStim,size(x_zscore,2));
        X_voc_wholeset = nan(length(y_dev)* NbTrialStim,size(X_voc,2));
        yy=0;
        for TS=1:length(y_dev)
            YT=y_dev{TS};
            for tt=1:length(YT)
                yy=yy+1;
                y_wholeset(yy) = YT(tt);
                Y_temp=YT;
                Y_temp(tt)=[];
                y_wholeset_bestGuess(yy) = mean(Y_temp);
                x_wholeset(yy,:)=x_zscore(TS,:);
                X_voc_wholeset(yy,:) = X_voc(TS,:);
            end
        end
        x_wholeset=x_wholeset(1:yy,:);
        y_wholeset=y_wholeset(1:yy);
        y_wholeset_bestGuess=y_wholeset_bestGuess(1:yy);
        X_voc_wholeset=X_voc_wholeset(1:yy,:);
    else
        x_wholeset=x_zscore;
        y_wholeset=y;
        X_voc_wholeset = X_voc;
    end
    Data.x_wholeset{mm}=x_wholeset;
    Data.X_voc_wholeset{mm}=X_voc_wholeset;
    Data.y_wholeset{mm}=y_wholeset;
    Data.y_wholeset_bestGuess{mm}=y_wholeset_bestGuess;
    
    %% Calculate STRF Acoustic and AcSem using lasso glm and cross-validation
    
    if ParamModel.ModelChoice(1)
        Deviance_BestModel_Acoustic=nan(ParamModel.BootstrapSTRF,length(Alphas));
        LL_BestModel_Acoustic=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_Acoustic_local=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_Acoustic=nan(length(Alphas),1);
        FitIndexAR_Acoustic=nan(length(Alphas),1);
        FitIndexARVal_Acoustic=nan(length(Alphas),1);
        DF_BestModel_Acoustic=nan(ParamModel.BootstrapSTRF,length(Alphas));
        Lambda_BestModel_Acoustic=nan(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB_Ac=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB0_Ac=cell(ParamModel.BootstrapSTRF,length(Alphas));
    end
    
    if ParamModel.ModelChoice(2)
        Deviance_BestModel_Sem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        LL_BestModel_Sem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_Sem=nan(length(Alphas),1);
        FitIndex_Sem_local=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndexAR_Sem=nan(length(Alphas),1);
        FitIndexARVal_Sem=nan(length(Alphas),1);
        DF_BestModel_Sem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB_Sem=cell(ParamModel.BootstrapSTRF,length(Alphas));
    end
    
    if ParamModel.ModelChoice(3)
        Deviance_BestModel_AcSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        LL_BestModel_AcSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSem_local=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSem=nan(length(Alphas),1);
        FitIndexAR_AcSem=nan(length(Alphas),1);
        FitIndexARVal_AcSem=nan(length(Alphas),1);
        DF_BestModel_AcSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        Lambda_BestModel_AcSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBspectro_AcSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBsem_AcSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB0_AcSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
    end
    
    if  ParamModel.ModelChoice(4)
        DF_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,length(Alphas));
        Lambda_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB0_AcSemAc=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBsem_AcSemAc=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBspectro_AcSemAc=cell(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSemAc_local=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSemAc=nan(length(Alphas),1);
        FitIndexAR_AcSemAc=nan(length(Alphas),1);
        FitIndexARVal_AcSemAc=nan(length(Alphas),1);
        Deviance_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,length(Alphas));
        LL_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,length(Alphas));
    end
    if  ParamModel.ModelChoice(5)
        DF_BestModel_AcSemSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        Lambda_BestModel_AcSemSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        ModelB0_AcSemSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBsem_AcSemSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
        ModelBspectro_AcSemSem=cell(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSemSem_local=nan(ParamModel.BootstrapSTRF,length(Alphas));
        FitIndex_AcSemSem=nan(length(Alphas),1);
        FitIndexAR_AcSemSem=nan(length(Alphas),1);
        FitIndexARVal_AcSemSem=nan(length(Alphas),1);
        Deviance_BestModel_AcSemSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
        LL_BestModel_AcSemSem=nan(ParamModel.BootstrapSTRF,length(Alphas));
    end
        
    
    
    LL_BestGuess=nan(ParamModel.BootstrapSTRF,length(Alphas));
    LL_WorseGuess=nan(ParamModel.BootstrapSTRF,length(Alphas));
    LL_AR=nan(ParamModel.BootstrapSTRF,length(Alphas));
    LL_AR2=nan(ParamModel.BootstrapSTRF,length(Alphas));
    
    ValProp=nan(ParamModel.BootstrapSTRF,length(Alphas));
    ValSize=nan(ParamModel.BootstrapSTRF,length(Alphas));
   
    % Loop through alphas and cross validation bootstratps
    for aa=1:length(Alphas)
        Alpha=Alphas(aa);
        fprintf('Alpha=%f\n',Alpha);
        
        
        []=runLassoGlmModels(ParamModel.ModelChoice,ParamModel.BootstrapSTRF,VOC{mm},Emitter,StimLocal
        
        PropVal.mean(mm,aa) = mean(ValProp(:,aa));% Store the average proportion of stims in the validating set over bootstrap for that alpha
        PropVal.std(mm,aa) = std(ValProp(:,aa));% Store the std on proportion of stims in the validating set over bootstrap for that alpha
        
        if ParamModel.ModelChoice(1)
            Deviance.Acoustic.values{mm,aa}=Deviance_All_local_Ac;% Store all the deviance of all bootstraps for that alpha value
            Deviance.Acoustic.lambda{mm,aa}=Lambda_All_local_Ac;% Store all the lambdas of all bootstraps for that alpha value
            Model.Acoustic.ypredictVal{mm,aa} = Ypredict_Acoustic;
            LL.Acoustic.values{mm,aa}=LL_All_local_Ac;
            FitIndex_Acoustic(aa)=(sum(LL_BestModel_Acoustic(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_Acoustic(aa)=(sum(LL_BestModel_Acoustic(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexARVal_Acoustic(aa)=(sum(LL_BestModel_Acoustic(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR2(:,aa))-sum(LL_WorseGuess(:,aa)));
        end
        if ParamModel.ModelChoice(2)
            Model.Semantic.ypredictVal{mm,aa} = Ypredict_Semantic;
            FitIndex_Sem(aa)=(sum(LL_BestModel_Sem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_Sem(aa)=(sum(LL_BestModel_Sem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexARVal_Sem(aa)=(sum(LL_BestModel_Sem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR2(:,aa))-sum(LL_WorseGuess(:,aa)));
        end
        if ParamModel.ModelChoice(3)
            Deviance.AcSem.values{mm,aa}=Deviance_All_local_AcSem;% Store all the deviance of all bootstraps for that alpha value
            Deviance.AcSem.lambda{mm,aa}=Lambda_All_local_AcSem;% Store all the lambdas of all bootstraps for that alpha value
            Model.AcSem.ypredictVal{mm,aa} = Ypredict_AcSem;
            LL.AcSem.values{mm,aa}=LL_All_local_AcSem;
            FitIndex_AcSem(aa)=(sum(LL_BestModel_AcSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_AcSem(aa)=(sum(LL_BestModel_AcSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexARVal_AcSem(aa)=(sum(LL_BestModel_AcSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR2(:,aa))-sum(LL_WorseGuess(:,aa)));
        end
        if ParamModel.ModelChoice(4)
            Deviance.AcSemAc.values{mm,aa}=Deviance_All_local_AcSemAc;% Store all the deviance of all bootstraps for that alpha value
            Deviance.AcSemAc.lambda{mm,aa}=Lambda_All_local_AcSemAc;% Store all the lambdas of all bootstraps for that alpha value
            Model.AcSemAc.ypredictVal{mm,aa} = Ypredict_AcSemAc;
            LL.AcSemAc.values{mm,aa}=LL_All_local_AcSemAc;
            FitIndex_AcSemAc(aa)=(sum(LL_BestModel_AcSemAc(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_AcSemAc(aa)=(sum(LL_BestModel_AcSemAc(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexARVal_AcSemAc(aa)=(sum(LL_BestModel_AcSemAc(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR2(:,aa))-sum(LL_WorseGuess(:,aa)));
        end
        if ParamModel.ModelChoice(5)
            Deviance.AcSemSem.values{mm,aa}=Deviance_All_local_AcSemSem;% Store all the deviance of all bootstraps for that alpha value
            Deviance.AcSemSem.lambda{mm,aa}=Lambda_All_local_AcSemSem;% Store all the lambdas of all bootstraps for that alpha value
            Model.AcSemSem.ypredictVal{mm,aa} = Ypredict_AcSemSem;
            LL.AcSemSem.values{mm,aa}=LL_All_local_AcSemSem;
            FitIndex_AcSemSem(aa)=(sum(LL_BestModel_AcSemSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_AcSemSem(aa)=(sum(LL_BestModel_AcSemSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexARVal_AcSemSem(aa)=(sum(LL_BestModel_AcSemSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR2(:,aa))-sum(LL_WorseGuess(:,aa)));
        end
        
        Model.yobsVal{mm,aa} = Yobs;
        Model.yobsCat{mm,aa} = YVoc;
        
        
        
        
        
        if SWITCH.FIG>0
            % plot LogLikelihood of all models over bootstraps
            models=fieldnames(LL);
            %Colors='brmgcky';
            close(figure(14))
            figure(14)
            hold on
            if ParamModel.ModelChoice(1)
                plot(1:ParamModel.BootstrapSTRF,LL_BestModel_Acoustic(:,aa),'b.-')
                hline(mean(LL_BestModel_Acoustic(:,aa)),'b:')
            end
            if ParamModel.ModelChoice(2)
                plot(1:ParamModel.BootstrapSTRF,LL_BestModel_Sem(:,aa),'g.-')
                hline(mean(LL_BestModel_Sem(:,aa)),'g:')
            end
            if ParamModel.ModelChoice(3)
                plot(1:ParamModel.BootstrapSTRF,LL_BestModel_AcSem(:,aa),'m.-')
                hline(mean(LL_BestModel_AcSem(:,aa)),'m:')
            end
            if ParamModel.ModelChoice(4)
                plot(1:ParamModel.BootstrapSTRF,LL_BestModel_AcSemAc(:,aa),'.-','Color',[0.5 0 1])
                hline_color(mean(LL_BestModel_AcSemAc(:,aa)),':','','Color',[0.5 0 1])
            end
            if ParamModel.ModelChoice(5)
                plot(1:ParamModel.BootstrapSTRF,LL_BestModel_AcSemSem(:,aa),'r.-')
                hline(mean(LL_BestModel_AcSemSem(:,aa)),'r:')
            end
            plot(1:ParamModel.BootstrapSTRF,LL_BestGuess(:,aa),'y.-',1:ParamModel.BootstrapSTRF,LL_WorseGuess(:,aa),'k.-',1:ParamModel.BootstrapSTRF,LL_AR2(:,aa),'c.--');
            legend(models);
            hline(mean(LL_WorseGuess(:,aa)),'k:')
            hline(mean(LL_BestGuess(:,aa)),'y:')
            hline(mean(LL_AR2(:,aa)),'c:o')
            ylabel('LogLikelihood');
            xlabel('Bootstrap');
            hold off
        end

        %% Use the average Lambda to calculate the optimal models  
        if ParamModel.ModelChoice(1)
            % Use the average Lambda to calculate the optimal STRF Ac on all dataset
            [B_Ac_opt, FitInfo_Ac]=lassoglm(x_wholeset,y_wholeset,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK, 'Lambda',median(Lambda_BestModel_Acoustic(:,aa)));
            Model.Acoustic.B{mm,aa}=reshape(B_Ac_opt'.*x_std,Df,Dt);
            Model.Acoustic.Lambdas(mm,aa) = FitInfo_Ac.Lambda;
            Model.Acoustic.B0{mm,aa} = FitInfo_Ac.Intercept;

            if SWITCH.FIG>0
                figure(15)
                imagesc(Model.Acoustic.B{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AC lassoglm poisson whole dataset\nlog of lambda=%f',log(FitInfo_Ac.Lambda)))
            end
        end
        
        if ParamModel.ModelChoice(2)
            % Use all the dataset to calculate the optimal Semantic model
            MDL_Sem_opt=fitglm(X_voc_wholeset,y_wholeset,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK);
            Model.Semantic.B{mm,aa}=MDL_Sem_opt.Coefficients.Estimate(2:end);
            Model.Semantic.B0{mm,aa}=MDL_Sem_opt.Coefficients.Estimate(1);
            if SWITCH.FIG>0
                figure(16)
                ModelSemantic = repmat(Model.Semantic.B0{mm,aa},1,length(Model.Semantic.B{mm,aa})+1) + [0 Model.Semantic.B{mm,aa}'];
                plot(ModelSemantic,1:length(UVOC));
                title('Optimal Semantic Model');
                set(gca,'YTickLabel', UVOC);
            end
        end
        if ParamModel.ModelChoice(3)
            % Use the average Lambda to calculate the optimal STRF AcSem on all dataset
            [B_AcSem, FitInfo_AcSem]=lassoglm([X_voc_wholeset x_wholeset],y_wholeset,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK, 'Lambda',median(Lambda_BestModel_AcSem(:,aa)));
            Model.AcSem.Bspectro{mm,aa}=reshape(B_AcSem((size(X_voc,2)+1):end)'.*x_std,Df,Dt);
            Model.AcSem.Bsem{mm,aa}=B_AcSem(1:(size(X_voc,2)));
            Model.AcSem.Lambdas(mm,aa) = FitInfo_AcSem.Lambda;
            Model.AcSem.B0{mm,aa} = FitInfo_AcSem.Intercept;
            if SWITCH.FIG>0
                figure(17)
                subplot(1,2,1)
                imagesc(Model.AcSem.Bspectro{mm,aa})
                axis xy
                %title(sprintf('OPTIMAL STRF AcSem lassoglm poisson whole dataset\nlog of lambda=%f\n D=%f+/-%f (ratio of sum of LL for all bootstrap)\ncross-validation including %f+/-%f %% of the data',log10(FitInfo_AcSem.Lambda),Deviance.AcSem.mean(mm,aa), Deviance.AcSem.std(mm,aa), PropVal.mean(mm,aa), PropVal.std(mm,aa)))
                ss=subplot(1,2,2);
                ModelAcSem = repmat(Model.AcSem.B0{mm,aa},1,length(Model.AcSem.Bsem{mm,aa})+1) + [0 Model.AcSem.Bsem{mm,aa}'];
                plot(ModelAcSem,1:length(UVOC));
                title('Optimal AcSem Model');
                set(ss,'YTickLabel', UVOC);
            end
        end
        
        if ParamModel.ModelChoice(4)
            % Use the average Lambda to calculate the optimal STRF AcSemAc on all dataset
            Ac_Offset=x_wholeset*B_Ac_opt + repmat(Model.Acoustic.B0{mm,aa},size(x_wholeset,1),1);
            % Calculate The AcSemAc model
            [B_AcSemAc_local2, FitInfo_AcSemAc_local2]=lassoglm([X_voc_wholeset x_wholeset],y_wholeset,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'Standardize',0,'Lambda',median(Lambda_BestModel_AcSemAc(:,aa)), 'Offset', Ac_Offset);
            % then we need to add the coefficients of the acoustic
            % Model that were used to compute the offset
            B_AcSemAc_local2(1:length(MDL_Sem_opt.Coefficients.Estimate(2:end)))=B_AcSemAc_local2(1:length(MDL_Sem_opt.Coefficients.Estimate(2:end)));
            B_AcSemAc=B_AcSemAc_local2 + repmat([zeros(length(MDL_Sem_opt.Coefficients.Estimate(2:end)),1); B_Ac_opt],1,size(B_AcSemAc_local2,2));
            FitInfo_AcSemAc=FitInfo_AcSemAc_local2;
            FitInfo_AcSemAc.Intercept=FitInfo_AcSemAc.Intercept + Model.Acoustic.B0{mm,aa};
            Model.AcSemAc.Bspectro{mm,aa}=reshape(B_AcSemAc((size(X_voc,2)+1):end)'.*x_std,Df,Dt);
            Model.AcSemAc.Bsem{mm,aa}=B_AcSemAc(1:(size(X_voc,2)));
            Model.AcSemAc.Lambdas(mm,aa) = FitInfo_AcSemAc.Lambda;
            Model.AcSemAc.B0{mm,aa} = FitInfo_AcSemAc.Intercept;
            if SWITCH.FIG>0
                figure(18)
                subplot(1,2,1)
                imagesc(Model.AcSemAc.Bspectro{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AcSem AcOffset lassoglm poisson whole dataset\nlog of lambda=%f',log10(FitInfo_AcSemAc.Lambda)))
                ss=subplot(1,2,2);
                ModelAcSemAc = repmat(Model.AcSemAc.B0{mm,aa},1,length(Model.AcSemAc.Bsem{mm,aa})+1) + [0 Model.AcSemAc.Bsem{mm,aa}'];
                plot(ModelAcSemAc,1:length(UVOC));
                title('Optimal AcSem Model Ac Offset');
                set(ss,'YTickLabel', UVOC);
            end
        end
        
        if ParamModel.ModelChoice(5)
            % Calculate the Semantic Offset
            Sem_Offset=sum(repmat(MDL_Sem_opt.Coefficients.Estimate(2:end)',size(X_voc_wholeset,1),1).*X_voc_wholeset,2) + repmat(MDL_Sem_opt.Coefficients.Estimate(1),size(X_voc_wholeset,1),1);
            
            [B_AcSemSem_local, FitInfo_AcSemSem_local]=lassoglm([X_voc_wholeset x_wholeset],y_wholeset,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'Standardize',0,'Lambda',median(Lambda_BestModel_AcSemAc(:,aa)), 'Offset', Sem_Offset);
                
            % then we need to add the coefficients of the semantic
            % Model that were used to compute the offset
            B_AcSemSem_local(1:length(MDL_Sem_opt.Coefficients.Estimate(2:end)))=B_AcSemSem_local(1:length(MDL_Sem_opt.Coefficients.Estimate(2:end)));
            B_AcSemSem=B_AcSemSem_local + repmat([MDL_Sem_opt.Coefficients.Estimate(2:end)' zeros(1,size(B_AcSemSem_local,1)-length(MDL_Sem_opt.Coefficients.Estimate(2:end)))]',1,size(B_AcSemSem_local,2));
            FitInfo_AcSemSem=FitInfo_AcSemSem_local;
            FitInfo_AcSemSem.Intercept=FitInfo_AcSemSem.Intercept + MDL_Sem_opt.Coefficients.Estimate(1);

            Model.AcSemSem.Bspectro{mm,aa}=reshape(B_AcSemSem((size(X_voc,2)+1):end)'.*x_std,Df,Dt);
            Model.AcSemSem.Bsem{mm,aa}=B_AcSemSem(1:(size(X_voc,2)));
            Model.AcSemSem.Lambdas(mm,aa) = FitInfo_AcSemSem.Lambda;
            Model.AcSemSem.B0{mm,aa} = FitInfo_AcSemSem.Intercept;
            if SWITCH.FIG>0
                figure(19)
                subplot(1,2,1)
                imagesc(Model.AcSemSem.Bspectro{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AcSem SemOffset lassoglm poisson whole dataset\nlog of lambda=%f',log10(FitInfo_AcSemSem.Lambda)))
                ss=subplot(1,2,2);
                ModelAcSemSem = repmat(Model.AcSemSem.B0{mm,aa},1,length(Model.AcSemSem.Bsem{mm,aa})+1) + [0 Model.AcSemSem.Bsem{mm,aa}'];
                plot(ModelAcSemSem,1:length(UVOC));
                title('Optimal AcSem Model Ac Offset');
                set(ss,'YTickLabel', UVOC);
            end
        end
        
        
        %% I'm HERE Calculate the Deviances of optimal models
%         ypred_Acoustic = glmval([FitInfo_Ac.Intercept; B_Ac],x_wholeset,ParamModel.LINK);
%         ypred_AcSem = glmval([FitInfo_AcSem.Intercept; SemBoost.*B_AcSem],[X_voc_wholeset x_wholeset],ParamModel.LINK);
%         ypred_Semantic = glmval(MDL_Sem.Coefficients.Estimate,X_voc_wholeset,ParamModel.LINK);
%         if strcmp(ParamModel.DISTR,'poisson')
%             Model.Acoustic.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_Acoustic)) - (y_wholeset_bestGuess - ypred_Acoustic)));
%             Model.AcSem.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_AcSem)) - (y_wholeset_bestGuess - ypred_AcSem)));
%             Model.Semantic.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_Semantic)) - (y_wholeset_bestGuess - ypred_Semantic)));
%         elseif strcmp(ParamModel.DISTR,'normal')
%             Model.Acoustic.Deviance(mm,aa) =sum(power(y_wholeset - ypred_Acoustic,2));
%             Model.AcSem.Deviance(mm,aa) =sum(power(y_wholeset - ypred_AcSem,2));
%             Model.Semantic.Deviance(mm,aa) =sum(power(y_wholeset - ypred_Semantic,2));
%         else
%             fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
%         end
        
%         %% Calculate DF of optimal models
%         % Find the real degrees of freedom for each best model
%         % Note that in Matlab elastic net the parameters of Lasso (L1)
%         % and Ridge (L2) can be retrieved from the parameters of the
%         % elastic net A and L as: L2=L*(1-A)/2 and L1=L*A
%         % The real DF is calculated from the formula as indicated in
%         % (http://web.stanford.edu/~hastie/TALKS/enet_talk.pdf)
%         % AcSem
%         fprintf('Optimal AcSem model: Calculate exact number of parameters\n')
%         L2_AcSem = Model.AcSem.Lambdas(mm,aa)*(1-Alpha)/2;
%         %L1_AcSem = Lambda_BestModel_AcSem(bb,aa)*Alpha;
%         Active_Predictors_AcSem = find(SemBoost.*B_AcSem);
%         X_AcSem = [X_voc_wholeset x_wholeset];
%         Xact_AcSem = X_AcSem(:,Active_Predictors_AcSem);
%         Xact_AcSem_inv = (Xact_AcSem' * Xact_AcSem + L2_AcSem .* eye(length(Active_Predictors_AcSem)))^(-1);
%         Model.AcSem.DF(mm,aa) = trace(Xact_AcSem * Xact_AcSem_inv * Xact_AcSem');
%         
%         % Acoustic
%         fprintf('Acoustic model: Calculate exact number of parameters\n')
%         L2_Acoustic = Model.Acoustic.Lambdas(mm,aa)*(1-Alpha)/2;
%         %L1_Acoustic = Lambda_BestModel_Acoustic(bb,aa)*Alpha;
%         Active_Predictors_Ac = find(B_Ac);
%         Xact_Acoustic = x_wholeset(:,Active_Predictors_Ac);
%         Xact_Acoustic_inv = (Xact_Acoustic' * Xact_Acoustic + L2_Acoustic .* eye(length(Active_Predictors_Ac)))^(-1);
%         Model.Acoustic.DF(mm,aa) = trace(Xact_Acoustic * Xact_Acoustic_inv * Xact_Acoustic');
%         
%         % Semantic
%         fprintf('Semantic model: save number of parameters\n')
%         Model.Semantic.DF(mm,aa) = length(MDL_Sem.Coefficients.Estimate);
    end
    
    % Store all the deviance and DF of all bootstraps for all alpha values
    if ParamModel.ModelChoice(1)
        Deviance.Acoustic.bestvalues{mm} = Deviance_BestModel_Acoustic;
        Deviance.Acoustic.DF{mm} = DF_BestModel_Acoustic;
        Deviance.Acoustic.FitIndex{mm}=FitIndex_Acoustic;
        Deviance.Acoustic.FitIndexAR{mm}=FitIndexAR_Acoustic;
        Deviance.Acoustic.FitIndexARVal{mm}=FitIndexARVal_Acoustic;
        LL.Acoustic.bestvalues{mm} = LL_BestModel_Acoustic;
    end
    if ParamModel.ModelChoice(2)
        Deviance.Sem.DF{mm} = DF_BestModel_Sem;
        Deviance.Sem.bestvalues{mm}=Deviance_BestModel_Sem;
        Deviance.Sem.FitIndex{mm}=FitIndex_Sem;
        Deviance.Sem.FitIndexAR{mm}=FitIndexAR_Sem;
        Deviance.Sem.FitIndexARVal{mm}=FitIndexARVal_Sem;
        LL.Sem.bestvalues{mm} = LL_BestModel_Sem;
    end
    if ParamModel.ModelChoice(3)
        Deviance.AcSem.bestvalues{mm} = Deviance_BestModel_AcSem;
        Deviance.AcSem.DF{mm} = DF_BestModel_AcSem;
        Deviance.AcSem.FitIndex{mm}=FitIndex_AcSem;
        Deviance.AcSem.FitIndexAR{mm}=FitIndexAR_AcSem;
        Deviance.AcSem.FitIndexARVal{mm}=FitIndexARVal_AcSem;
        LL.AcSem.bestvalues{mm} = LL_BestModel_AcSem;
    end
    if ParamModel.ModelChoice(4)
        Deviance.AcSemAc.bestvalues{mm} = Deviance_BestModel_AcSemAc;
        Deviance.AcSemAc.DF{mm} = DF_BestModel_AcSemAc;
        Deviance.AcSemAc.FitIndex{mm}=FitIndex_AcSemAc;
        Deviance.AcSemAc.FitIndexAR{mm}=FitIndexAR_AcSemAc;
        Deviance.AcSemAc.FitIndexARVal{mm}=FitIndexARVal_AcSemAc;
        LL.AcSemAc.bestvalues{mm} = LL_BestModel_AcSemAc;
    end
    if ParamModel.ModelChoice(5)
        Deviance.AcSemSem.bestvalues{mm} = Deviance_BestModel_AcSemSem;
        Deviance.AcSemSem.DF{mm} = DF_BestModel_AcSemSem;
        Deviance.AcSemSem.FitIndex{mm}=FitIndex_AcSemSem;
        Deviance.AcSemSem.FitIndexAR{mm}=FitIndexAR_AcSemSem;
        Deviance.AcSemSem.FitIndexARVal{mm}=FitIndexARVal_AcSemSem;
        LL.AcSemSem.bestvalues{mm} = LL_BestModel_AcSemSem;
    end
    
    PropVal.values{mm} = ValSize;
    LL.Ceiling.bestvalues{mm} = LL_BestGuess;
    LL.Floor.bestvalues{mm} = LL_WorseGuess;
    LL.AutoRegressive.bestvalues{mm} = LL_AR;
    LL.AutoRegressiveVal.bestvalues{mm} = LL_AR2;
    
    
    
%     % Calculate the probability for Deviance and AIC  AcSem < Acoustic
%     % AcSem < Sem
%     DiffDeviance_AcSemAcoustic = Deviance_BestModel_AcSem - Deviance_BestModel_Acoustic;
%     DiffAIC_AcSemAcoustic = 2.*DF_BestModel_AcSem - 2.*DF_BestModel_Acoustic + DiffDeviance_AcSemAcoustic;
%     DiffDeviance_AcSemSemantic = Deviance_BestModel_AcSem - Deviance_BestModel_Sem;
%     DiffAIC_AcSemSemantic = 2.*DF_BestModel_AcSem - 2.*DF_BestModel_Sem + DiffDeviance_AcSemSemantic;
%     PDeviance_AcSemAcoustic{mm} = nan(length(Alphas),1);
%     PAIC_AcSemAcoustic{mm} = nan(length(Alphas),1);
%     PDeviance_AcSemSemantic{mm} = nan(length(Alphas),1);
%     PAIC_AcSemSemantic{mm}=nan(length(Alphas),1);
%     for aa=1:length(Alphas)
%         PDeviance_AcSemAcoustic{mm}(aa) = length(find(DiffDeviance_AcSemAcoustic(:,aa)<0))/ParamModel.BootstrapSTRF;
%         PAIC_AcSemAcoustic{mm}(aa) = length(find(DiffAIC_AcSemAcoustic(:,aa)<0))/ParamModel.BootstrapSTRF;
%         PDeviance_AcSemSemantic{mm}(aa) = length(find(DiffDeviance_AcSemSemantic(:,aa)<0))/ParamModel.BootstrapSTRF;
%         PAIC_AcSemSemantic{mm}(aa) = length(find(DiffAIC_AcSemSemantic(:,aa)<0))/ParamModel.BootstrapSTRF;
%     end
    
    %% keep track of the data and stim used for the next window
    PreviousWin.Stim_local_old = Stim_local;
    PreviousWin.y_dev_old = y_dev;
    
        
end


% Save data 
    save(fullfile(resultsDirectory,sprintf('%s_GLMPoisson.mat',Cellname)))
end

