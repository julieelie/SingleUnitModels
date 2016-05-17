function [Deviance, NeuroRes, VOC, Wins] = GrowingModelsRidgeglm(Spectro, VocType, PSTH, Trials, Emitter, MinWin, MaxWin, Increment, ResDelay, NeuroRes, DISTR, LINK, Cellname)
resultsDirectory='/auto/tdrive/julie/NeuralData/SemanticGLMModel';
%resultsDirectory='/Users/elie/Documents/CODE/SingleUnitModels/Resultsglmlasso/glmLasso_3Models_fitPerTrial/OptimalDevianceOffsetModel';
FIG=0; % set to 1 for some debugging figures 2 for all debugging figures
Check=1;%set to 1 to compare ridge results with Linear model results
BootstrapSTRF=25; %set the number of time you want to estimate the STRF at each time window
ZC=1;%Set to 0 if you don't want to z-score the values of the spectrograms before performing lasso glm
SemBoost = 1;
FitPerTrial = 1; %Set to 0 if you want to fit the model using the average response per stim
NbTrialStim=10;
LAMBDARATIO=1e-4;
OFFSET=2;%set 0 if you don't want the Acsem model to be calculated with an offset obtained from the semantic model
% set to 1 for AcSem calculated with Semantic Offset
%set to 2 for AcSem calculated with Semantic and Acoustic Offset

if nargin<12
    LINK='log'; %'identity'
end

if nargin<11
    DISTR='poisson';%'normal'
end

if nargin<10
    NeuroRes = 'mean';
end
if nargin<9
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<8
    Increment = 40; %increase the size of the spectro window with a 5ms pace
end
if nargin<7
    MaxWin = 160; %maximum values the window of analysis can reach
end
if nargin<6
    MinWin = 40; %minimum size of the window of analysis from the begining and also size of analysis of spike rate
end
Flow = 8000;%spectrograms are low passed at 8Khz for the calculations

%define the increasing size of the window of the spectrogram
%Wins = MinWin:Increment:MaxWin;
Wins=160;

% # of models to run on the data
modNum = length(Wins);

% Number of stims in the data set
NbStim = length(VocType);

% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
Alphas=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

%% Initialize a bunch of output variables
Deviance.Acoustic.bestvalues = cell(modNum);
Deviance.Acoustic.DF = cell(modNum);
Deviance.Acoustic.values = cell(modNum,length(Alphas));
Deviance.Acoustic.lambda = cell(modNum,length(Alphas));
Deviance.Acoustic.FitIndex=cell(modNum);
Deviance.Acoustic.FitIndexAR=cell(modNum);
Deviance.Sem.bestvalues = cell(modNum);
Deviance.Sem.DF = cell(modNum);
Deviance.Sem.FitIndex=cell(modNum);
Deviance.AcSem = Deviance.Acoustic;
LL.Acoustic.values=cell(modNum,length(Alphas));
LL.Acoustic.bestvalues= cell(modNum);
LL.AcSem=LL.Acoustic;
LL.AcSem2=LL.Acoustic;
LL.Ceiling.bestvalues=cell(modNum);
LL.Floor.bestvalues=cell(modNum);
LL.AutoRegressive.bestvalues=cell(modNum);

PropVal.mean = nan(modNum,length(Alphas));
PropVal.std = nan(modNum,length(Alphas));
PropVal.values = cell(modNum);

Model.Acoustic.Lambdas = nan(modNum,length(Alphas));
Model.Acoustic.B = cell(modNum,length(Alphas));
Model.Acoustic.B0 = cell(modNum,length(Alphas));
Model.Acoustic.ypredictVal = cell(modNum,length(Alphas));
Model.Acoustic.Deviance=nan(modNum,length(Alphas));
Model.Acoustic.DF=nan(modNum,length(Alphas));
Model.AcSem.Lambdas = nan(modNum,length(Alphas));
Model.AcSem.Bspectro = cell(modNum,length(Alphas));
Model.AcSem.Bsem = cell(modNum,length(Alphas));
Model.AcSem.B0 = cell(modNum,length(Alphas));
Model.AcSem.ypredictVal = cell(modNum,length(Alphas));
Model.AcSem.Deviance=nan(modNum,length(Alphas));
Model.AcSem.DF=nan(modNum,length(Alphas));
Model.Semantic.B = cell(modNum,length(Alphas));
Model.Semantic.B0 = cell(modNum,length(Alphas));
Model.Semantic.ypredictVal = cell(modNum,1);
Model.Semantic.Deviance=nan(modNum,length(Alphas));
Model.Semantic.DF=nan(modNum,length(Alphas));
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

PDeviance_AcSemAcoustic = cell(modNum,1);
PAIC_AcSemAcoustic = cell(modNum,1);
PDeviance_AcSemSemantic = cell(modNum,1);
PAIC_AcSemSemantic=cell(modNum,1);

VOC = cell(modNum,1);

if OFFSET>1
    Deviance.AcSem2 = Deviance.Acoustic;
    Model.AcSem2.ypredictVal = cell(modNum,length(Alphas));
end

% Initialize datasets to be used as passed information by the
% auto-regressive model
y_dev_old=[];
Stim_local_old=[];

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
    if FitPerTrial
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
        if strcmp(NeuroRes, 'max')
            y(ss) = max(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));
        elseif strcmp(NeuroRes, 'mean')
            y(ss) = mean(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
            if FitPerTrial
                y_dev{ss}=nan(length(Trials{dd}),1);
                for tt=1:length(Trials{dd})
                    y_dev{ss}(tt)=sum((Trials{dd}{tt}>(Win-MinWin+ResDelay)).*(Trials{dd}{tt}<(Win+ResDelay)));
                end
                ymean(ss)=mean(y_dev{ss});
                yvar(ss)=var(y_dev{ss});
            end
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', NeuroRes);
    
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
    if FitPerTrial
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
    Deviance_BestModel_Acoustic=nan(BootstrapSTRF,length(Alphas));
    Deviance_BestModel_AcSem=nan(BootstrapSTRF,length(Alphas));
    Deviance_BestModel_Sem=nan(BootstrapSTRF,length(Alphas));
    LL_BestModel_Acoustic=nan(BootstrapSTRF,length(Alphas));
    LL_BestModel_AcSem=nan(BootstrapSTRF,length(Alphas));
    LL_BestGuess=nan(BootstrapSTRF,length(Alphas));
    LL_WorseGuess=nan(BootstrapSTRF,length(Alphas));
    LL_AR=nan(BootstrapSTRF,length(Alphas));
    LL_BestModel_Sem=nan(BootstrapSTRF,length(Alphas));
    FitIndex_Acoustic_local=nan(BootstrapSTRF,length(Alphas));
    FitIndex_AcSem_local=nan(BootstrapSTRF,length(Alphas));
    FitIndex_Sem_local=nan(BootstrapSTRF,length(Alphas));
    FitIndex_Acoustic=nan(length(Alphas),1);
    FitIndex_AcSem=nan(length(Alphas),1);
    FitIndex_Sem=nan(length(Alphas),1);
    FitIndexAR_Acoustic=nan(length(Alphas),1);
    FitIndexAR_AcSem=nan(length(Alphas),1);
    FitIndexAR_Sem=nan(length(Alphas),1);
    DF_BestModel_Acoustic=nan(BootstrapSTRF,length(Alphas));
    DF_BestModel_AcSem=nan(BootstrapSTRF,length(Alphas));
    DF_BestModel_Sem=nan(BootstrapSTRF,length(Alphas));
    Lambda_BestModel_Acoustic=nan(BootstrapSTRF,length(Alphas));
    Lambda_BestModel_AcSem=nan(BootstrapSTRF,length(Alphas));
    Lambda_BestModel_AcSem2=nan(BootstrapSTRF,length(Alphas));
    ValProp=nan(BootstrapSTRF,length(Alphas));
    ValSize=nan(BootstrapSTRF,length(Alphas));
    ModelB_Ac=cell(BootstrapSTRF,length(Alphas));
    ModelB0_Ac=cell(BootstrapSTRF,length(Alphas));
    ModelBspectro_AcSem=cell(BootstrapSTRF,length(Alphas));
    ModelBsem_AcSem=cell(BootstrapSTRF,length(Alphas));
    ModelB0_AcSem=cell(BootstrapSTRF,length(Alphas));
    ModelB_Sem=cell(BootstrapSTRF,length(Alphas));
    
    if OFFSET>1
        DF_BestModel_AcSem2=nan(BootstrapSTRF,length(Alphas));
        ModelB0_AcSem2=cell(BootstrapSTRF,length(Alphas));
        ModelBsem_AcSem2=cell(BootstrapSTRF,length(Alphas));
        ModelBspectro_AcSem2=cell(BootstrapSTRF,length(Alphas));
        FitIndex_AcSem2_local=nan(BootstrapSTRF,length(Alphas));
        FitIndex_AcSem2=nan(length(Alphas),1);
        FitIndexAR_AcSem2=nan(length(Alphas),1);
        Deviance_BestModel_AcSem2=nan(BootstrapSTRF,length(Alphas));
        LL_BestModel_AcSem2=nan(BootstrapSTRF,length(Alphas));
   end
    
    
    % Loop through alphas and cross validation bootstratps
    for aa=1:length(Alphas)
        Alpha=Alphas(aa);
        fprintf('Alpha=%f\n',Alpha);
        Deviance_All_local_Ac=cell(BootstrapSTRF,1);
        Lambda_All_local_Ac = cell(BootstrapSTRF,1);
        LL_All_local_Ac=cell(BootstrapSTRF,1);
        Deviance_All_local_AcSem=cell(BootstrapSTRF,1);
        Lambda_All_local_AcSem = cell(BootstrapSTRF,1);
        LL_All_local_AcSem=cell(BootstrapSTRF,1);
        Ypredict_Acoustic = cell(BootstrapSTRF,1);
        Yobs = cell(BootstrapSTRF,1);
        YVoc = cell(BootstrapSTRF,1);
        Ypredict_Semantic = cell(BootstrapSTRF,1);
        Ypredict_AcSem = cell(BootstrapSTRF,1);
        if OFFSET>1
            Deviance_All_local_AcSem2=cell(BootstrapSTRF,1);
            Lambda_All_local_AcSem2 = cell(BootstrapSTRF,1);
            Ypredict_AcSem2 = cell(BootstrapSTRF,1);
            LL_All_local_AcSem2=cell(BootstrapSTRF,1);
        end
        
        for bb=1:BootstrapSTRF
            tic
            %% Construct a testing dataset and a validating dataset
            % remove one emitter per call category, if one category contain only
            % one emitter, don't use it in the model
            [ValSet, TrainSet] = create_cross_validation_sets(VOC{mm}, Emitter.Ename(Stim_local));
            ValProp(bb,aa) = length(ValSet)./(length(ValSet)+length(TrainSet));
            ValSize(bb,aa) = length(ValSet);
            
            if ZC==1
                    % z-score spectrogram data
                    x_Train_mean = mean(x(TrainSet,:));
                    x_Train_std = std(x(TrainSet,:));
                    x_Train = (x(TrainSet,:)-repmat(x_Train_mean,length(TrainSet),1))./repmat(x_Train_std,length(TrainSet),1);
                    x_Train(:,x_Train_std==0)=0;%set these parameters to zero as they do not cary any information and would give Inf or NaN
                    x_Val = (x(ValSet,:)-repmat(x_Train_mean,length(ValSet),1))./repmat(x_Train_std,length(ValSet),1);
                    x_Val(:,x_Train_std==0)=0;%set these parameters to zero as they do not cary any information and would give Inf or NaN
                else
                    x_Train = x(TrainSet,:);
                    x_Val = x(ValSet,:);
            end
            
            %% construct the vector of responses for the training dataset
            if FitPerTrial
                yy=0;
                x_Train_new = nan(NbTrialStim.*length(TrainSet),size(x_Train,2));
                X_voc_Train = nan(NbTrialStim.*length(TrainSet),size(X_voc,2));
                y_Train=nan(NbTrialStim.*length(TrainSet),1);
                if ~isempty(y_dev_old)
                    y_input_AR=nan(NbTrialStim.*length(TrainSet),3);
                end
                for TS=1:length(TrainSet) 
                    YT=y_dev{TrainSet(TS)};
                    for tt=1:length(YT)
                        yy=yy+1;
                        y_Train(yy) = YT(tt);
                        x_Train_new(yy,:)=x_Train(TS,:);
                        X_voc_Train(yy,:) = X_voc(TrainSet(TS),:);
                        if ~isempty(y_dev_old)
                            Local_trialset=YT;
                            Local_trialset(tt)=[];
                            y_input_AR(yy,1)=mean(Local_trialset);
                            Stim_num_real=Stim_local(TrainSet(TS));
                            Stim_num_old=find(Stim_local_old==Stim_num_real);
                            Local_trialset=y_dev_old{Stim_num_old};
                            y_input_AR(yy,3)=Local_trialset(tt);
                            Local_trialset(tt)=[];
                            y_input_AR(yy,2)=mean(Local_trialset);
                        end
                    end
                end
                x_Train=x_Train_new(1:yy,:);
                X_voc_Train=X_voc_Train(1:yy,:);
                y_Train=y_Train(1:yy);
                if ~isempty(y_dev_old)
                    y_input_AR=y_input_AR(1:yy,:);
                end
            else
                y_Train = y(TrainSet);
                X_voc_Train = X_voc(TrainSet,:)
            end
            
%             %% Set the value of the parameters for SemBoost
%             SemBoost=max(max(abs(x_Train)));
%             X_voc_Train_local=2*SemBoost.*X_voc_Train;

           
            
            %% First run Semantic glm Model and use the Beta.*x as an offset for the Acsem Model
            %Semantic Model
            fprintf(1,'Semantic Model poisson glm %d/%d\n', bb, BootstrapSTRF);
            
            MDL_Sem=fitglm(X_voc_Train,y_Train,'Distribution',DISTR,'Link',LINK);
            
            
            % Calculate the Semantic Offset
            if OFFSET
                Sem_Offset=sum(repmat(MDL_Sem.Coefficients.Estimate(2:end)',size(X_voc_Train,1),1).*X_voc_Train,2) + repmat(MDL_Sem.Coefficients.Estimate(1),size(X_voc_Train,1),1);
            end
            %% Run Acoustic Lassoglm model and find the best parameter (Lambdas) for that alpha to minimize the deviance.
            fprintf(1,'Acoustic Model lasso glm %d/%d\n', bb, BootstrapSTRF);
            if sum(sum(isnan(x_Train)))
                fprintf(1,'x_Train has %d NAN values and ZC=%d\n %d STD values of x(TrainSet,:) equal zero\n',sum(sum(isnan(x_Train))),ZC,sum(x_Train_std==0));
            end
            if sum(sum(isinf(x_Train)))
                fprintf(1,'x_Train has inf values and ZC=%d\n',ZC);
            end
                
            tic
            [B_Ac, FitInfo_Ac]=lassoglm(x_Train,y_Train,DISTR,'Alpha',Alpha,'Link',LINK,'NumLambda',25,'Standardize',0,'LambdaRatio',LAMBDARATIO);
            toc
            if FIG>2
                lassoPlot(B_Ac,FitInfo_Ac,'PlotType','Lambda','XScale','log')
                pause()
            end
            

            %% Run Acoustic +Semantic Lassoglm model with Semantic Offset and find the best parameter (Lambdas) for that alpha to minimize the deviance.
            
            fprintf(1,'Semantic + Acoustic Model lasso glm %d/%d\n', bb, BootstrapSTRF);
            if OFFSET
                
                [B_AcSem_local, FitInfo_AcSem_local]=lassoglm([SemBoost.*X_voc_Train x_Train],y_Train,DISTR,'Alpha',Alpha,'Link',LINK,'NumLambda',25,'Standardize',0,'LambdaRatio',LAMBDARATIO, 'Offset', Sem_Offset);
                
                % then we need to add the coefficients of the semantic
                % Model that were used to compute the offset
                B_AcSem_local(1:length(MDL_Sem.Coefficients.Estimate(2:end)))=B_AcSem_local(1:length(MDL_Sem.Coefficients.Estimate(2:end))).*SemBoost;
                B_AcSem=B_AcSem_local + repmat([MDL_Sem.Coefficients.Estimate(2:end)' zeros(1,size(B_AcSem_local,1)-length(MDL_Sem.Coefficients.Estimate(2:end)))]',1,size(B_AcSem_local,2));
                FitInfo_AcSem=FitInfo_AcSem_local;
                FitInfo_AcSem.Intercept=FitInfo_AcSem.Intercept + MDL_Sem.Coefficients.Estimate(1);
            else
                
                [B_AcSem_local, FitInfo_AcSem]=lassoglm([SemBoost.*X_voc_Train x_Train],y_Train,DISTR,'Alpha',Alpha,'Link',LINK,'NumLambda',25,'Standardize',0,'LambdaRatio',LAMBDARATIO);
                B_AcSem_local(1:length(MDL_Sem.Coefficients.Estimate(2:end)))=B_AcSem_local(1:length(MDL_Sem.Coefficients.Estimate(2:end))).*SemBoost;
                B_AcSem=B_AcSem_local;
            end
            if FIG>2
                lassoPlot(B_AcSem,FitInfo_AcSem,'PlotType','Lambda','XScale','log')
                pause()
            end
            
            %% Calculate predicted responses for Floor, ceiling model and autoregressive model
            if FitPerTrial
                % construct the vector of responses for the validating
                % dataset (actual values, best guess and worse guess)
                yy=0;
                y_Val=nan(length(ValSet).*NbTrialStim,1);
                y_Val_bestGuess = nan(length(ValSet).*NbTrialStim,1);
                y_Val_worseGuess = nan(length(ValSet).*NbTrialStim,1);
                if ~isempty(y_dev_old)
                    y_Val_input_AR = nan(length(ValSet).*NbTrialStim,3);
                end
                
                for vv=1:length(ValSet)
                        YT=y_dev{ValSet(vv)};
                        for tt=1:length(YT)
                            yy=yy+1;
                            y_Val(yy)=YT(tt);
                            YT_localmean = YT;
                            YT_localmean(tt)=[];
                            y_Val_bestGuess(yy)=mean(YT_localmean);
                            if ~isempty(y_dev_old)
                                y_Val_input_AR(yy,1)=y_Val_bestGuess(yy);
                                Stim_num_real=Stim_local(ValSet(vv));
                                Stim_num_old=find(Stim_local_old==Stim_num_real);
                                Local_trialset=y_dev_old{Stim_num_old};
                                y_Val_input_AR(yy,3)=Local_trialset(tt);
                                Local_trialset(tt)=[];
                                y_Val_input_AR(yy,2)=mean(Local_trialset);
                            end
                            Meansperstim_Cat_worseguess=cell(length(UVOC),1);
                            for ww=1:length(ValSet)
                                if ww~=vv
                                    YTw=y_dev{ValSet(ww)};
                                    MeanperStim_local=mean(YTw);
                                else
                                    MeanperStim_local=y_Val_bestGuess(yy);
                                end
                                for uu=1:length(UVOC)
                                    if strcmp(VOC{mm}(ValSet(ww)),UVOC{uu})
                                        Meansperstim_Cat_worseguess{uu}=[Meansperstim_Cat_worseguess{uu} MeanperStim_local];
                                    end
                                end
                            end
                            MeanperCat_worseguess=nan(length(UVOC),1);
                            EmptyCat=[];
                            for uu=1:length(UVOC)
                                if isempty(Meansperstim_Cat_worseguess{uu})
                                    EmptyCat=[EmptyCat uu];
                                else
                                    MeanperCat_worseguess(uu)=mean(Meansperstim_Cat_worseguess{uu});
                                end
                            end
                            if ~isempty(EmptyCat)
                                for ee=1:length(EmptyCat)
                                    MeanperCat_worseguess(EmptyCat(ee))=mean(MeanperCat_worseguess(setdiff(1:length(UVOC),EmptyCat)));% this line is to provide a mean for the worse guess for category not represented in the dataset and it is the average spike rate on other categories
                                end
                            end
                            y_Val_worseGuess(yy)=mean(MeanperCat_worseguess);
                        end
                end
                y_Val=y_Val(1:yy);% Actual values for the validating dataset
                y_Val_bestGuess=y_Val_bestGuess(1:yy);% Ceiling model: Best prediction (average of other trials for the same stim)
                y_Val_worseGuess=y_Val_worseGuess(1:yy);% Floor Model: worse prediction(average of the mean 
                
                if isempty(y_dev_old)
                    ypred_AR=y_Val_bestGuess; %For the first window no auto-regresssive model for now, just predicion from other trials
                else
                    % Auto-regressive model
                    fprintf(1,'Auto-regressive Model poisson glm %d/%d\n', bb, BootstrapSTRF);
                    y_Val_input_AR= y_Val_input_AR(1:yy,:);% Autoregressive model input y: validating input for the autoregressive model
                    MDL_AR=fitglm(y_input_AR,y_Train,'Distribution',DISTR,'Link',LINK);
                end
                    
            end
            

            %% Calculate Deviances, loglikelihood and predicted values for validating dataset
            % Deviance of Acoustic model for the validating dataset for each lambda
            fprintf(1,'Acoustic Model: Calculate deviance for each lambda on the validating dataset\n');
            NbL_Ac=length(FitInfo_Ac.Lambda);% number of lambdas tested
            Deviance_M_Ac= nan(NbL_Ac,1);
            LL_Ac= nan(NbL_Ac,1);
            if FitPerTrial
                % Construct the vector of the predicted response using the
                % Acoustic model
                ypred_Ac=nan(length(ValSet).*NbTrialStim,NbL_Ac);
                for ll=1:NbL_Ac
                    yy=0;
                    ypred_Ac_temp = glmval([FitInfo_Ac.Intercept(ll); B_Ac(:,ll)],x_Val,LINK);
                    for vv=1:length(ValSet)
                        YT=y_dev{ValSet(vv)};
                        for tt=1:length(YT)
                            yy=yy+1;
                            ypred_Ac(yy,ll)=ypred_Ac_temp(vv);
                        end
                    end
                    if strcmp(DISTR,'poisson')
                         Deviance_M_Ac(ll) = sum(2*(y_Val .* (log((y_Val_bestGuess+(y_Val_bestGuess==0)) ./ ypred_Ac(1:yy,ll))) - (y_Val_bestGuess - ypred_Ac(1:yy,ll))));
                         LL_Ac(ll) = sum(y_Val .* log(ypred_Ac(1:yy,ll)+(ypred_Ac(1:yy,ll)==0))-log(factorial(y_Val))-ypred_Ac(1:yy,ll));
                    elseif strcmp(DISTR,'normal')
                        Deviance_M_Ac(ll) = sum(power(y_Val - ypred_Ac(1:yy,ll),2));
                    else
                        fprintf(1,'your distribution is unrecognize');
                    end
                end
                LL_All_local_Ac{bb}=LL_Ac;
                LL_BestGuess(bb,aa) = sum(y_Val .* log(y_Val_bestGuess+(y_Val_bestGuess==0))-log(factorial(y_Val))-y_Val_bestGuess);
                LL_WorseGuess(bb,aa) = sum(y_Val .* log(y_Val_worseGuess+(y_Val_worseGuess==0))-log(factorial(y_Val))-y_Val_worseGuess);
                
            else
                y_Val=y(ValSet);
                for ll=1:NbL_Ac
                    ypred_Ac(:,ll) = glmval([FitInfo_Ac.Intercept(ll); B_Ac(:,ll)],x_Val,LINK);
                    if strcmp(DISTR,'poisson')
                         Deviance_M_Ac(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_Ac(:,ll))) - (y_Val - ypred_Ac(:,ll))));
                    elseif strcmp(DISTR,'normal')
                        Deviance_M_Ac(ll) = sum(power(y_Val - ypred_Ac(:,ll),2));
                    else
                        fprintf(1,'your distribution is unrecognize');
                    end
                end    
            end

            
            Deviance_All_local_Ac{bb} = Deviance_M_Ac;% Store All deviance values for that bootstrap
            Lambda_All_local_Ac{bb}=FitInfo_Ac.Lambda;% Store all Lambda tested for that booststrap

            % Deviance of Acoustic + Semantic model for the validating dataset for each lambda
            fprintf(1,'AcSem model: Calculate deviance for each lambda on the validating dataset\n');
            NbL_AcSem=length(FitInfo_AcSem.Lambda);% number of lambdas tested
            Deviance_M_AcSem= nan(NbL_AcSem,1);
            LL_AcSem=nan(NbL_AcSem,1);
            if FitPerTrial
                % Construct the vector of predicted responses using the
                % AcSem1 model
                ypred_AcSem=nan(length(ValSet).*NbTrialStim,NbL_AcSem);
                for ll=1:NbL_AcSem
                    ypred_AcSem_temp = glmval([FitInfo_AcSem.Intercept(ll); B_AcSem(:,ll)],[X_voc(ValSet,:) x_Val],LINK);
                    yy=0;
                    for vv=1:length(ValSet)
                        for tt=1:length(y_dev{ValSet(vv)})
                            yy=yy+1;
                            ypred_AcSem(yy,ll)=ypred_AcSem_temp(vv);
                        end
                    end
                    if strcmp(DISTR,'poisson')
                        Deviance_M_AcSem(ll) = sum(2*(y_Val .* (log((y_Val_bestGuess+(y_Val_bestGuess==0)) ./ ypred_AcSem(1:yy,ll))) - (y_Val_bestGuess - ypred_AcSem(1:yy,ll))));
                        LL_AcSem(ll) = sum(y_Val .* log(ypred_AcSem(1:yy,ll)+(ypred_AcSem(1:yy,ll)==0))-log(factorial(y_Val))-ypred_AcSem(1:yy,ll));
                    elseif strcmp(DISTR,'normal')
                        Deviance_M_AcSem(ll)=sum(power(y_Val - ypred_AcSem(1:yy,ll),2));
                    else
                        fprintf('unrecognize %s distribution for deviance calculation', DISTR);
                    end
                end
                ypred_AcSem=ypred_AcSem(1:yy,:);
                LL_All_local_AcSem{bb}=LL_AcSem;
            else
                ypred_AcSem(:,ll)=glmval([FitInfo_AcSem.Intercept(ll); B_AcSem(:,ll)],[X_voc(ValSet,:) x_Val],LINK);
                if strcmp(DISTR,'poisson')
                    Deviance_M_AcSem(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_AcSem(:,ll))) - (y_Val - ypred_AcSem(:,ll))));
                elseif strcmp(DISTR,'normal')
                    Deviance_M_AcSem(ll)=sum(power(y_Val - ypred_AcSem(:,ll),2));
                else
                    fprintf('unrecognize %s distribution for deviance calculation', DISTR);
                end
            end
                
            
            Deviance_All_local_AcSem{bb} = Deviance_M_AcSem;% Store All deviance values for that bootstrap
            Lambda_All_local_AcSem{bb}=FitInfo_AcSem.Lambda;% Store all Lambda tested for that booststrap

            % Deviance of Semantic model for the validating dataset for each lambda
            fprintf(1,'Semantic model: Calculate deviance on the validating dataset\n');
            ypred_Sem_temp = glmval(MDL_Sem.Coefficients.Estimate,X_voc(ValSet,:),LINK);
            if FitPerTrial
                % construc the vector of predicted responses using the semantic model 
                yy=0;
                ypred_Sem=nan(length(ValSet).*NbTrialStim,1);
                for vv=1:length(ValSet)
                    for tt=1:length(y_dev{ValSet(vv)})
                        yy=yy+1;
                        ypred_Sem(yy)=ypred_Sem_temp(vv);
                    end
                end
                ypred_Sem=ypred_Sem(1:yy);
            else
                ypred_Sem=ypred_Sem_temp;
            end
            if strcmp(DISTR,'poisson')
                Deviance_BestModel_Sem(bb,aa) = sum(2*(y_Val .* (log((y_Val_bestGuess+(y_Val_bestGuess==0)) ./ ypred_Sem)) - (y_Val_bestGuess - ypred_Sem)));
                LL_BestModel_Sem(bb,aa) = sum(y_Val .* log(ypred_Sem+(ypred_Sem==0))-log(factorial(y_Val))-ypred_Sem);
                FitIndex_Sem_local(bb,aa)=(LL_BestModel_Sem(bb,aa)-LL_WorseGuess(bb,aa))/(LL_BestGuess(bb,aa)-LL_WorseGuess(bb,aa));
            elseif strcmp(DISTR,'normal')
                Deviance_BestModel_Sem(bb,aa) =sum(power(y_Val - ypred_Sem,2));  
            else
                fprintf('unrecognize %s distribution for deviance calculation', DISTR);
            end   
            ModelB_Sem{bb,aa} = MDL_Sem.Coefficients.Estimate;% Parameters of the model at that bootstrap
            
            % LogLikelihood of the autoregressive model
            if ~isempty(y_dev_old)
                fprintf(1,'Auto-regressive model: Calculate deviance on the validating dataset\n');
                ypred_AR = glmval(MDL_AR.Coefficients.Estimate,y_Val_input_AR,LINK);
            end
            LL_AR(bb,aa) = sum(y_Val .* log(ypred_AR+(ypred_AR==0))-log(factorial(y_Val))-ypred_AR);
           %% Store the smallest Deviance Value, the corresponding best lambda and the corresponding best model
            % of best deviance, the parameters of the model for that best lambda for that booststrap
            fprintf(1,'Acoustic model: Find lambda with smallest value of Deviance\n');
            Dev_min_index_Ac=find(Deviance_M_Ac==min(Deviance_M_Ac));
            Dev_min_index_Ac=Dev_min_index_Ac(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
            Deviance_BestModel_Acoustic(bb,aa)= Deviance_M_Ac(Dev_min_index_Ac);
            LL_BestModel_Acoustic(bb,aa)=LL_Ac(Dev_min_index_Ac);
            FitIndex_Acoustic_local(bb,aa)=(LL_BestModel_Acoustic(bb,aa)-LL_WorseGuess(bb,aa))/(LL_BestGuess(bb,aa)-LL_WorseGuess(bb,aa));
            Lambda_BestModel_Acoustic(bb,aa)=FitInfo_Ac.Lambda(Dev_min_index_Ac);
            ModelB_Ac{bb,aa} = reshape(B_Ac(:,Dev_min_index_Ac)'.*x_Train_std,Df,Dt);% Parameters of the best model at that bootstrap
            ModelB0_Ac{bb,aa} = FitInfo_Ac.Intercept(Dev_min_index_Ac);% intercept of the best model
            fprintf(1,'AcSem model: Find lambda with smallest value of Deviance\n');
            Dev_min_index_AcSem=find(Deviance_M_AcSem==min(Deviance_M_AcSem));
            Dev_min_index_AcSem=Dev_min_index_AcSem(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
            Deviance_BestModel_AcSem(bb,aa)= Deviance_M_AcSem(Dev_min_index_AcSem);
            LL_BestModel_AcSem(bb,aa)=LL_AcSem(Dev_min_index_AcSem);
            FitIndex_AcSem_local(bb,aa)=(LL_BestModel_AcSem(bb,aa)-LL_WorseGuess(bb,aa))/(LL_BestGuess(bb,aa)-LL_WorseGuess(bb,aa));
            Lambda_BestModel_AcSem(bb,aa)=FitInfo_AcSem.Lambda(Dev_min_index_AcSem);
            ModelBspectro_AcSem{bb,aa} = reshape(B_AcSem((size(X_voc,2)+1):end,Dev_min_index_AcSem)'.*x_Train_std,Df,Dt);% Parameters of the best model at that bootstrap
            ModelBsem_AcSem{bb,aa} = B_AcSem(1:size(X_voc,2),Dev_min_index_AcSem);
            % identify the category for which there was no stim and
            % calculate our best guess for B for that category: the average
            % of the B of other categories
            noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
            if ~isempty(noModCat)
                if ~sum(B_AcSem(noModCat,Dev_min_index_AcSem))==0
                    fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                end
                ModelBsem_AcSem{bb,aa}(noModCat)=sum(B_AcSem(1:size(X_voc,2),Dev_min_index_AcSem))/(length(B_AcSem(1:size(X_voc,2),Dev_min_index_AcSem))-length(noModCat));
            end
            ModelB0_AcSem{bb,aa} = FitInfo_AcSem.Intercept(Dev_min_index_AcSem);% intercept of the best model

            %% Now that we now which Acoustic model is the best, calculate the offset of the best acoustic model to be used by the AcSem2 model
            
            if OFFSET>1
                % Calculate the Acoustic Offset
                Ac_Offset=x_Train*B_Ac(:,Dev_min_index_Ac) + repmat(ModelB0_Ac{bb,aa},size(x_Train,1),1);
                fprintf(1,'Semantic + Acoustic Model lasso glm with Acoustic Offset %d/%d\n', bb, BootstrapSTRF);
                % Calculate The AcSem2 model
                [B_AcSem_local2, FitInfo_AcSem_local2]=lassoglm([SemBoost.*X_voc_Train x_Train],y_Train,DISTR,'Alpha',Alpha,'Link',LINK,'NumLambda',25,'Standardize',0,'LambdaRatio',LAMBDARATIO, 'Offset', Ac_Offset);
                % then we need to add the coefficients of the acoustic
                % Model that were used to compute the offset
                B_AcSem_local2(1:length(MDL_Sem.Coefficients.Estimate(2:end)))=B_AcSem_local2(1:length(MDL_Sem.Coefficients.Estimate(2:end))).*SemBoost;
                B_AcSem2=B_AcSem_local2 + repmat([zeros(length(MDL_Sem.Coefficients.Estimate(2:end)),1); B_Ac(:,Dev_min_index_Ac)],1,size(B_AcSem_local,2));
                FitInfo_AcSem2=FitInfo_AcSem_local2;
                FitInfo_AcSem2.Intercept=FitInfo_AcSem2.Intercept + FitInfo_Ac.Intercept;
                
                % Deviance of Acoustic + Semantic model2 for the validating dataset for each lambda
                fprintf(1,'AcSem model2: Calculate deviance for each lambda on the validating dataset\n');
                NbL_AcSem2=length(FitInfo_AcSem2.Lambda);% number of lambdas tested
                Deviance_M_AcSem2= nan(NbL_AcSem2,1);
                LL_AcSem2 = nan(NbL_AcSem2,1);
                if FitPerTrial
                    ypred_AcSem2=nan(length(ValSet).*NbTrialStim,NbL_AcSem2);
                    for ll=1:NbL_AcSem2
                        ypred_AcSem2_temp = glmval([FitInfo_AcSem2.Intercept(ll); B_AcSem2(:,ll)],[X_voc(ValSet,:) x_Val],LINK);
                        yy=0;
                        for vv=1:length(ValSet)
                            for tt=1:length(y_dev{ValSet(vv)})
                                yy=yy+1;
                                ypred_AcSem2(yy,ll)=ypred_AcSem2_temp(vv);
                            end
                        end
                        if strcmp(DISTR,'poisson')
                            Deviance_M_AcSem2(ll) = sum(2*(y_Val .* (log((y_Val_bestGuess+(y_Val_bestGuess==0)) ./ ypred_AcSem2(1:yy,ll))) - (y_Val_bestGuess - ypred_AcSem2(1:yy,ll))));
                            LL_AcSem2(ll) = sum(y_Val .* log(ypred_AcSem2(1:yy,ll)+(ypred_AcSem2(1:yy,ll)==0))-log(factorial(y_Val))-ypred_AcSem2(1:yy,ll));
                        elseif strcmp(DISTR,'normal')
                            Deviance_M_AcSem2(ll)=sum(power(y_Val - ypred_AcSem2(1:yy,ll),2));
                        else
                            fprintf('unrecognize %s distribution for deviance calculation', DISTR);
                        end
                    end
                    ypred_AcSem2=ypred_AcSem2(1:yy,:);
                    LL_All_local_AcSem2{bb}=LL_AcSem2;
                else
                    ypred_AcSem2(:,ll)=glmval([FitInfo_AcSem2.Intercept(ll); B_AcSem2(:,ll)],[X_voc(ValSet,:) x_Val],LINK);
                    if strcmp(DISTR,'poisson')
                        Deviance_M_AcSem2(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_AcSem2(:,ll))) - (y_Val - ypred_AcSem2(:,ll))));
                    elseif strcmp(DISTR,'normal')
                        Deviance_M_AcSem2(ll)=sum(power(y_Val - ypred_AcSem2(:,ll),2));
                    else
                        fprintf('unrecognize %s distribution for deviance calculation', DISTR);
                    end
                end
                Deviance_All_local_AcSem2{bb} = Deviance_M_AcSem2;% Store All deviance values for that bootstrap
                Lambda_All_local_AcSem2{bb}=FitInfo_AcSem2.Lambda;% Store all Lambda tested for that booststrap
                % Calculate the smallest deviance
                fprintf(1,'AcSem2 model: Find lambda with smallest value of Deviance\n');
                Dev_min_index_AcSem2=find(Deviance_M_AcSem2==min(Deviance_M_AcSem2));
                Dev_min_index_AcSem2=Dev_min_index_AcSem2(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
                Deviance_BestModel_AcSem2(bb,aa)= Deviance_M_AcSem2(Dev_min_index_AcSem2);
                LL_BestModel_AcSem2(bb,aa)=LL_AcSem2(Dev_min_index_AcSem2);
                Lambda_BestModel_AcSem2(bb,aa)=FitInfo_AcSem2.Lambda(Dev_min_index_AcSem2);
                FitIndex_AcSem2_local(bb,aa)=(LL_BestModel_AcSem2(bb,aa)-LL_WorseGuess(bb,aa))/(LL_BestGuess(bb,aa)-LL_WorseGuess(bb,aa));
                ModelBspectro_AcSem2{bb,aa} = reshape(B_AcSem2((size(X_voc,2)+1):end,Dev_min_index_AcSem2)'.*x_Train_std,Df,Dt);% Parameters of the best model at that bootstrap
                ModelBsem_AcSem2{bb,aa} = B_AcSem2(1:size(X_voc,2),Dev_min_index_AcSem2);
                % identify the category for which there was no stim and
                % calculate our best guess for B for that category: the average
                % of the B of other categories
                noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSem2(noModCat,Dev_min_index_AcSem2))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    ModelBsem_AcSem2{bb,aa}(noModCat)=sum(B_AcSem2(1:size(X_voc,2),Dev_min_index_AcSem2))/(length(B_AcSem2(1:size(X_voc,2),Dev_min_index_AcSem2))-length(noModCat));
                end
                ModelB0_AcSem2{bb,aa} = FitInfo_AcSem2.Intercept(Dev_min_index_AcSem2);% intercept of the best model
                % Save predictions
                Ypredict_AcSem2{bb} = ypred_AcSem2(:,Dev_min_index_AcSem2);
                %Calculate degrees of freedom
                % AcSem
            fprintf('AcSem2 model: Calculate exact number of parameters\n')
            L2_AcSem2 = Lambda_BestModel_AcSem2(bb,aa)*(1-Alpha)/2;
            Active_Predictors_AcSem2 = find(B_AcSem2(:,Dev_min_index_AcSem2));
            X_AcSem = [X_voc_Train x_Train];
            Xact_AcSem2 = X_AcSem(:,Active_Predictors_AcSem2);
            Xact_AcSem_inv2 = (Xact_AcSem2' * Xact_AcSem2 + L2_AcSem2 .* eye(length(Active_Predictors_AcSem2)))^(-1);
            DF_BestModel_AcSem2(bb,aa) = trace(Xact_AcSem2 * Xact_AcSem_inv2 * Xact_AcSem2');
                
            end
            
            
            %% Store the predicted values of the best model on the validating dataset
            Yobs{bb} = y_Val;
            if FitPerTrial
                VocTypeTrial=cell(length(ValSet).*NbTrialStim,1);
                yy=0;
                for vv=1:length(ValSet)
                    for tt=1:length(y_dev{ValSet(vv)})
                        yy=yy+1;
                        VocTypeTrial{yy}=VOC{mm}{ValSet(vv)};
                    end
                end
                YVoc{bb}=VocTypeTrial(1:yy,:);
            else
                YVoc{bb} = VOC{mm}(ValSet);
            end
            Ypredict_Acoustic{bb} = ypred_Ac(:,Dev_min_index_Ac);
            Ypredict_Semantic{bb} = ypred_Sem;
            Ypredict_AcSem{bb} = ypred_AcSem(:,Dev_min_index_AcSem);
            
            %% Calculate degress of freedom for each model
            % Find the real degrees of freedom for each best model
            % Note that in Matlab elastic net the parameters of Lasso (L1)
            % and Ridge (L2) can be retrieved from the parameters of the
            % elastic net A and L as: L2=L*(1-A)/2 and L1=L*A
            % The real DF is calculated from the formula as indicated in
            % (http://web.stanford.edu/~hastie/TALKS/enet_talk.pdf)
            % AcSem
            fprintf('AcSem model: Calculate exact number of parameters\n')
            L2_AcSem = Lambda_BestModel_AcSem(bb,aa)*(1-Alpha)/2;
            %L1_AcSem = Lambda_BestModel_AcSem(bb,aa)*Alpha;
            Active_Predictors_AcSem = find(B_AcSem(:,Dev_min_index_AcSem));
            X_AcSem = [X_voc_Train x_Train];
            Xact_AcSem = X_AcSem(:,Active_Predictors_AcSem);
            Xact_AcSem_inv = (Xact_AcSem' * Xact_AcSem + L2_AcSem .* eye(length(Active_Predictors_AcSem)))^(-1);
            DF_BestModel_AcSem(bb,aa) = trace(Xact_AcSem * Xact_AcSem_inv * Xact_AcSem');
            
            % Acoustic
            fprintf('Acoustic model: Calculate exact number of parameters\n')
            L2_Acoustic = Lambda_BestModel_Acoustic(bb,aa)*(1-Alpha)/2;
            %L1_Acoustic = Lambda_BestModel_Acoustic(bb,aa)*Alpha;
            Active_Predictors_Ac = find(B_Ac(:,Dev_min_index_Ac));
            Xact_Acoustic = x_Train(:,Active_Predictors_Ac);
            Xact_Acoustic_inv = (Xact_Acoustic' * Xact_Acoustic + L2_Acoustic .* eye(length(Active_Predictors_Ac)))^(-1);
            DF_BestModel_Acoustic(bb,aa) = trace(Xact_Acoustic * Xact_Acoustic_inv * Xact_Acoustic');
            
            % Semantic
            fprintf('Semantic model: save number of parameters\n')
            DF_BestModel_Sem(bb,aa) = length(ModelB_Sem{bb,aa});
            
            %% Out put some results of deviance difference between models
            if FIG>0
                % Calculate the deviance difference for that bootstrap and
                % the AIC difference
                DiffDEV_AcSemAcoustic = Deviance_BestModel_AcSem(bb,aa)-Deviance_BestModel_Acoustic(bb,aa); %the more this value is negative the better the Acsem model is over the acoustic model
                DiffAIC_AcSemAcoustic = 2.*DF_BestModel_AcSem(bb,aa) - 2.*DF_BestModel_Acoustic(bb,aa) + DiffDEV_AcSemAcoustic;
                fprintf(1,'Deviance difference AcSem-Acoustic: %f (the more negative the better the AcSem model is\n',DiffDEV_AcSemAcoustic);
                fprintf(1,'AIC difference AcSem-Acoustic: %f (the more negative the better the AcSem model is\n',DiffAIC_AcSemAcoustic);
                DiffDEV_AcSemSemantic = Deviance_BestModel_AcSem(bb,aa)-Deviance_BestModel_Sem(bb,aa); %the more this value is negative the better the Acsem model is over the acoustic model
                DiffAIC_AcSemSemantic = 2.*DF_BestModel_AcSem(bb,aa) - 2.*DF_BestModel_Sem(bb,aa) + DiffDEV_AcSemSemantic;
                fprintf(1,'Deviance difference AcSem-Semantic: %f (the more negative the better the AcSem model is\n',DiffDEV_AcSemSemantic);
                fprintf(1,'AIC difference AcSem-Semantic: %f (the more negative the better the AcSem model is\n',DiffAIC_AcSemSemantic);
            end
            
            %% Plot Lambda values, STRFs = model coefficients, predictions
            if FIG>1
                % Plot lambda values as a function of deviance for each
                % model (AcSem and Acoustic)
                figure(2)
                subplot(1,2,1);
                plot(log10(FitInfo_Ac.Lambda), Deviance_M_Ac)
                xlabel('log of Lambdas')
                ylabel('Deviance of the Acoustic model')
                title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                vline(log10(FitInfo_Ac.Lambda(Dev_min_index_Ac)));
                subplot(1,2,2);
                plot(log10(FitInfo_AcSem.Lambda), Deviance_M_AcSem)
                xlabel('log of Lambdas')
                ylabel('Deviance of the AcSem model')
                title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                vline(log10(FitInfo_AcSem.Lambda(Dev_min_index_AcSem)));
                
                pause(1)
                if Check==1
                    fprintf(1, 'Calculate PC of spectro\n');
                    [COEFF,SCORE,latent,tsquare]=princomp(x_Train,'econ');
                    nPC=60;
                    ds=dataset();
                    for ii=1:nPC
                        ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
                    end
                    ds.y=y_Train;
                    mdl=LinearModel.fit(ds);
                    PCSTRF=mdl.Coefficients.Estimate(2:end);
                    STRF=COEFF(:,1:nPC)*PCSTRF;
                    STRFM=reshape(STRF,Df, Dt);
                    fprintf(1, 'Calculating STRF Ac using the %d first PC of the spectro\n\n\n\n', nPC);
                    figure(3)
                    imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, STRFM)
                    axis xy
                    title(sprintf('STRF Ac obtained with linear model with %d PC of the spectro', nPC));
                    pause(1)
                end
                
                % Plot coefficients Acoustic model
                for ll=1:NbL_Ac
                    if ll==Dev_min_index_Ac
                        figure(4)
                        if ZC==1
                            Lasso_STRF=reshape(B_Ac(:,ll)'.*x_Train_std,Df,Dt);
                        else
                            Lasso_STRF=reshape(B_Ac(:,ll),Df,Dt);
                        end
                        valc = max(abs(max(max(Lasso_STRF))), abs(min(min(Lasso_STRF))));
                        if valc==0
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF)
                        else
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF,[-valc valc])
                        end
                        axis xy
                        title(sprintf('THIS IS THE ONE!!!\nSTRF Ac log of lambda=%f\n Deviance=%f FitIndex=%f',log10(FitInfo_Ac.Lambda(ll)),Deviance_BestModel_Acoustic(bb,aa),FitIndex_Acoustic_local(bb,aa)))
                    else
                        figure(5)
                        if ZC==1
                            Lasso_STRF=reshape(B_Ac(:,ll)'.*x_Train_std,Df,Dt);
                        else
                            Lasso_STRF=reshape(B_Ac(:,ll),Df,Dt);
                        end
                        valc = max(abs(max(max(Lasso_STRF))), abs(min(min(Lasso_STRF))));
                        if valc==0
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF)
                        else
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF,[-valc valc])
                        end
                        axis xy
                        title(sprintf('STRF Ac log of lambda=%f\n',log10(FitInfo_Ac.Lambda(ll))))
                    end
                    pause(1)
                end
                
                % Plot Coefficients AcSem model
                for ll=1:NbL_AcSem
                    if ll==Dev_min_index_AcSem
                        figure(6)
                        if ZC==1
                            Lasso_STRF_AcSem=reshape(B_AcSem((size(X_voc,2)+1):end,ll)'.*x_Train_std,Df,Dt);
                        else
                            Lasso_STRF_AcSem=reshape(B_AcSem((size(X_voc,2)+1):end,ll),Df,Dt);
                        end
                        valc = max(abs(max(max(Lasso_STRF_AcSem))), abs(min(min(Lasso_STRF_AcSem))));
                        subplot(1,2,1)
                        if valc==0
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem)
                        else
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem,[-valc valc])
                        end
                        axis xy
                        title(sprintf('THIS IS THE ONE!!!\nSTRF AcSem log of lambda=%f\n Deviance=%f FitIndex=%f',log10(FitInfo_AcSem.Lambda(ll)),Deviance_BestModel_AcSem(bb,aa),FitIndex_AcSem_local(bb,aa)))
                        ss=subplot(1,2,2);
                        % Make sure we have one stim per category in the
                        % training set or set the coefficient value of the
                        % model for that category as the average of
                        % parameters for other categories instead of 0
                        B_AcSem_plot = [FitInfo_AcSem.Intercept(ll) ; B_AcSem(1:size(X_voc,2),ll)] + [0; repmat(FitInfo_AcSem.Intercept(ll),size(X_voc,2),1)];
                        noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                        if ~isempty(noModCat)
                            if ~sum(B_AcSem(noModCat,ll))==0
                                fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                            end
                            B_AcSem_plot(noModCat)=FitInfo_AcSem.Intercept(ll) + sum(B_AcSem(1:size(X_voc,2),ll))/(length(B_AcSem(1:size(X_voc,2),ll))-length(noModCat));
                        end
                        plot(B_AcSem_plot, 1:length(UVOC));
                        title('Coefficients Categories AcSem Model');
                        set(ss,'YTickLabel', UVOC);
                    else
                        figure(7)
                        if ZC==1
                            Lasso_STRF_AcSem=reshape(B_AcSem(9:end,ll)'.*x_Train_std,Df,Dt);
                        else
                            Lasso_STRF_AcSem=reshape(B_AcSem(9:end,ll),Df,Dt);
                        end
                        valc = max(abs(max(max(Lasso_STRF_AcSem))), abs(min(min(Lasso_STRF_AcSem))));
                        subplot(1,2,1)
                        if valc==0
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem)
                        else
                            imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem,[-valc valc])
                        end
                        axis xy
                        title(sprintf('STRF AcSem log of lambda=%f\n',log10(FitInfo_AcSem.Lambda(ll))));
                        ss=subplot(1,2,2);
                        % Make sure we have one stim per category in the
                        % training set or set the coefficient value of the
                        % model for that category as the average of
                        % parameters for other categories instead of 0
                        B_AcSem_plot = [FitInfo_AcSem.Intercept(ll) ; B_AcSem(1:size(X_voc,2),ll)] + [0; repmat(FitInfo_AcSem.Intercept(ll),size(X_voc,2),1)];
                        noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                        if ~isempty(noModCat)
                            if ~sum(B_AcSem(noModCat,ll))==0
                                fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                            end
                            B_AcSem_plot(noModCat)=FitInfo_AcSem.Intercept(ll) + sum(B_AcSem(1:size(X_voc,2),ll))/(length(B_AcSem(1:size(X_voc,2),ll)) -length(noModCat));
                        end
                        plot(B_AcSem_plot, 1:length(UVOC));
                        title('Coefficients Categories AcSem Model');
                        set(ss,'YTickLabel', UVOC);
                    end
                    pause(1)
                end
                if OFFSET>1
                    % Plot Coefficients AcSem2 model
                    for ll=1:NbL_AcSem2
                        if ll==Dev_min_index_AcSem2
                            figure(8)
                            if ZC==1
                                Lasso_STRF_AcSem2=reshape(B_AcSem2((size(X_voc,2)+1):end,ll)'.*x_Train_std,Df,Dt);
                            else
                                Lasso_STRF_AcSem2=reshape(B_AcSem2((size(X_voc,2)+1):end,ll),Df,Dt);
                            end
                            valc = max(abs(max(max(Lasso_STRF_AcSem2))), abs(min(min(Lasso_STRF_AcSem2))));
                            subplot(1,2,1)
                            if valc==0
                                imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem2)
                            else
                                imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem2,[-valc valc])
                            end
                            axis xy
                            title(sprintf('THIS IS THE ONE!!!\nSTRF AcSem2 log of lambda=%f\n Deviance=%f FitIndex=%f',log10(FitInfo_AcSem2.Lambda(ll)),Deviance_BestModel_AcSem2(bb,aa),FitIndex_AcSem2_local(bb,aa)))
                            ss=subplot(1,2,2);
                            % Make sure we have one stim per category in the
                            % training set or set the coefficient value of the
                            % model for that category as the average of
                            % parameters for other categories instead of 0
                            B_AcSem_plot2 = [FitInfo_AcSem2.Intercept(ll) ; B_AcSem2(1:size(X_voc,2),ll)] + [0; repmat(FitInfo_AcSem2.Intercept(ll),size(X_voc,2),1)];
                            noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                            if ~isempty(noModCat)
                                if ~sum(B_AcSem2(noModCat,ll))==0
                                    fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                                end
                                B_AcSem_plot2(noModCat)=FitInfo_AcSem2.Intercept(ll) + sum(B_AcSem2(1:size(X_voc,2),ll))/(length(B_AcSem2(1:size(X_voc,2),ll))-length(noModCat));
                            end
                            plot(B_AcSem_plot2, 1:length(UVOC));
                            title('Coefficients Categories AcSem2 Model');
                            set(ss,'YTickLabel', UVOC);
                        else
                            figure(9)
                            if ZC==1
                                Lasso_STRF_AcSem2=reshape(B_AcSem2(9:end,ll)'.*x_Train_std,Df,Dt);
                            else
                                Lasso_STRF_AcSem2=reshape(B_AcSem2(9:end,ll),Df,Dt);
                            end
                            valc = max(abs(max(max(Lasso_STRF_AcSem2))), abs(min(min(Lasso_STRF_AcSem2))));
                            subplot(1,2,1)
                            if valc==0
                                imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem2)
                            else
                                imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.fo{mm}, Lasso_STRF_AcSem2,[-valc valc])
                            end
                            axis xy
                            title(sprintf('STRF AcSem2 log of lambda=%f\n',log10(FitInfo_AcSem2.Lambda(ll))));
                            ss=subplot(1,2,2);
                            % Make sure we have one stim per category in the
                            % training set or set the coefficient value of the
                            % model for that category as the average of
                            % parameters for other categories instead of 0
                            B_AcSem_plot2 = [FitInfo_AcSem2.Intercept(ll) ; B_AcSem2(1:size(X_voc,2),ll)] + [0; repmat(FitInfo_AcSem2.Intercept(ll),size(X_voc,2),1)];
                            noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                            if ~isempty(noModCat)
                                if ~sum(B_AcSem2(noModCat,ll))==0
                                    fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                                end
                                B_AcSem_plot2(noModCat)=FitInfo_AcSem2.Intercept(ll) + sum(B_AcSem2(1:size(X_voc,2),ll))/(length(B_AcSem2(1:size(X_voc,2),ll)) -length(noModCat));
                            end
                            plot(B_AcSem_plot2, 1:length(UVOC));
                            title('Coefficients Categories AcSem2 Model');
                            set(ss,'YTickLabel', UVOC);
                        end
                        pause(1)
                    end
                end
                
                % Plot coefficient Semantic model
                figure(10)
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_Sem_plot = MDL_Sem.Coefficients.Estimate + [0; repmat(MDL_Sem.Coefficients.Estimate(1),length(MDL_Sem.Coefficients.Estimate)-1,1)];
                noModCat= find(sum(X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(MDL_Sem.Coefficients.Estimate(noModCat))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_Sem_plot(noModCat)=MDL_Sem.Coefficients.Estimate(1) + sum(MDL_Sem.Coefficients.Estimate(2:end))/(length(MDL_Sem.Coefficients.Estimate)-length(noModCat));
                end
                plot(B_Sem_plot, 1:length(UVOC));
                title(sprintf('Coefficients Categories Sem Model Deviance=%f\nFitIndex=%f',Deviance_BestModel_Sem(bb,aa),FitIndex_Sem_local(bb,aa)));
                set(gca,'YTickLabel', UVOC);
                pause(1)
                
                % Plot prediction for each model
                LocalModelPredict=cell(3,1);
                LocalModelPredict{3}=Ypredict_AcSem{bb};
                LocalModelPredict{1}=Ypredict_Acoustic{bb};
                LocalModelPredict{2}=Ypredict_Semantic{bb};
                if OFFSET>1
                    LocalModelPredict{4}=Ypredict_AcSem2{bb};
                end
                Legend= cell(length(LocalModelPredict),1);
                Legend{1}='Acoustic only';
                Legend{2}='Semantic only';
                Legend{3}='Acoustic + Semantic';
                if OFFSET>1
                    Legend{4}='Acoustic + Semantic';
                    MAXY=max([max(Ypredict_Acoustic{bb}) max(Ypredict_Semantic{bb}) max(Ypredict_AcSem{bb}) max(Ypredict_AcSem2{bb})]);
                else
                    MAXY=max([max(Ypredict_Acoustic{bb}) max(Ypredict_Semantic{bb}) max(Ypredict_AcSem{bb})]);
                end
                MAXX=max(y_Val);
                MAX=max(MAXX,MAXY);
                
                for jj=1:length(LocalModelPredict)
                    figure(11)
                    subplot(length(LocalModelPredict),1,jj);
                    h=gscatter(y_Val, LocalModelPredict{jj}, YVoc{bb}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
                    title(sprintf('%s', Legend{jj}));
                    xlabel('observed spike rate (spike per ms)')
                    ylabel('predicted spike rate (spike per ms)')
                    axis([0 MAX+MAX/4 0 MAX]);
                    hold on
                    plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
                    hold off
                end
                pause(1)
            end
            toc
        end
        PropVal.mean(mm,aa) = mean(ValProp(:,aa));% Store the average proportion of stims in the validating set over bootstrap for that alpha
        PropVal.std(mm,aa) = std(ValProp(:,aa));% Store the std on proportion of stims in the validating set over bootstrap for that alpha
        
        Deviance.Acoustic.values{mm,aa}=Deviance_All_local_Ac;% Store all the deviance of all bootstraps for that alpha value
        Deviance.Acoustic.lambda{mm,aa}=Lambda_All_local_Ac;% Store all the lambdas of all bootstraps for that alpha value
        LL.Acoustic.values{mm,aa}=LL_All_local_Ac;
        FitIndex_Acoustic(aa)=(sum(LL_BestModel_Acoustic(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
        FitIndexAR_Acoustic(aa)=(sum(LL_BestModel_Acoustic(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
        Deviance.AcSem.values{mm,aa}=Deviance_All_local_AcSem;% Store all the deviance of all bootstraps for that alpha value
        Deviance.AcSem.lambda{mm,aa}=Lambda_All_local_AcSem;% Store all the lambdas of all bootstraps for that alpha value
        LL.AcSem.values{mm,aa}=LL_All_local_AcSem;
        FitIndex_AcSem(aa)=(sum(LL_BestModel_AcSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
        FitIndexAR_AcSem(aa)=(sum(LL_BestModel_AcSem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
        FitIndex_Sem(aa)=(sum(LL_BestModel_Sem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
        FitIndexAR_Sem(aa)=(sum(LL_BestModel_Sem(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
        Model.yobsVal{mm,aa} = Yobs;
        Model.yobsCat{mm,aa} = YVoc;
        Model.Acoustic.ypredictVal{mm,aa} = Ypredict_Acoustic;
        Model.Semantic.ypredictVal{mm,aa} = Ypredict_Semantic;
        Model.AcSem.ypredictVal{mm,aa} = Ypredict_AcSem;
        if OFFSET>1
            Deviance.AcSem2.values{mm,aa}=Deviance_All_local_AcSem2;% Store all the deviance of all bootstraps for that alpha value
            Deviance.AcSem2.lambda{mm,aa}=Lambda_All_local_AcSem2;% Store all the lambdas of all bootstraps for that alpha value
            Model.AcSem2.ypredictVal{mm,aa} = Ypredict_AcSem2;
            LL.AcSem2.values{mm,aa}=LL_All_local_AcSem2;
            FitIndex_AcSem2(aa)=(sum(LL_BestModel_AcSem2(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_BestGuess(:,aa))-sum(LL_WorseGuess(:,aa)));
            FitIndexAR_AcSem2(aa)=(sum(LL_BestModel_AcSem2(:,aa))-sum(LL_WorseGuess(:,aa)))/(sum(LL_AR(:,aa))-sum(LL_WorseGuess(:,aa)));
        end

        %% Use the average Lambda to calculate the optimal models  
        % Use the average Lambda to calculate the optimal STRF Ac on all dataset
        [B_Ac, FitInfo_Ac]=lassoglm(x_wholeset,y_wholeset,DISTR,'Alpha',Alpha,'Link',LINK, 'Lambda',median(Lambda_BestModel_Acoustic(:,aa)));
        Model.Acoustic.B{mm,aa}=reshape(B_Ac'.*x_std,Df,Dt);
        Model.Acoustic.Lambdas(mm,aa) = FitInfo_Ac.Lambda;
        Model.Acoustic.B0{mm,aa} = FitInfo_Ac.Intercept;
        
        if FIG>0
            figure(10)
            imagesc(Model.Acoustic.B{mm,aa})
            axis xy
            %PROBLEM HEREtitle(sprintf('OPTIMAL STRF AC lassoglm poisson whole dataset\nlog of lambda=%f\n D=%f+/-%f (ratio of sum of LL for all bootstrap)\ncross-validation including %f+/-%f %% of the data',log(FitInfo_Ac.Lambda),Deviance.Acoustic.mean(mm,aa), Deviance.Acoustic.std(mm,aa), PropVal.mean(mm,aa), PropVal.std(mm,aa)))
        end

        % Use the average Lambda to calculate the optimal STRF AcSem on all dataset
        [B_AcSem, FitInfo_AcSem]=lassoglm([SemBoost.*X_voc_wholeset x_wholeset],y_wholeset,DISTR,'Alpha',Alpha,'Link',LINK, 'Lambda',median(Lambda_BestModel_AcSem(:,aa)));
        Model.AcSem.Bspectro{mm,aa}=reshape(B_AcSem((size(X_voc,2)+1):end)'.*x_std,Df,Dt);
        Model.AcSem.Bsem{mm,aa}=SemBoost.*B_AcSem(1:(size(X_voc,2)));
        Model.AcSem.Lambdas(mm,aa) = FitInfo_AcSem.Lambda;
        Model.AcSem.B0{mm,aa} = FitInfo_AcSem.Intercept;
        if FIG>0
            figure(11)
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
        % Use all the dataset to calculate the optimal Semantic model
        MDL_Sem=fitglm(X_voc_wholeset,y_wholeset,'Distribution',DISTR,'Link',LINK);
        Model.Semantic.B{mm,aa}=MDL_Sem.Coefficients.Estimate(2:end);
        Model.Semantic.B0{mm,aa}=MDL_Sem.Coefficients.Estimate(1);
        if FIG>0
            figure(12)
            ModelSemantic = repmat(Model.Semantic.B0{mm,aa},1,length(Model.Semantic.B{mm,aa})+1) + [0 Model.Semantic.B{mm,aa}'];
            plot(ModelSemantic,1:length(UVOC));
            title('Optimal Semantic Model');
            set(gca,'YTickLabel', UVOC);
        end
        
        %% Calculate the Deviances of optimal models
        ypred_Acoustic = glmval([FitInfo_Ac.Intercept; B_Ac],x_wholeset,LINK);
        ypred_AcSem = glmval([FitInfo_AcSem.Intercept; SemBoost.*B_AcSem],[X_voc_wholeset x_wholeset],LINK);
        ypred_Semantic = glmval(MDL_Sem.Coefficients.Estimate,X_voc_wholeset,LINK);
        if strcmp(DISTR,'poisson')
            Model.Acoustic.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_Acoustic)) - (y_wholeset_bestGuess - ypred_Acoustic)));
            Model.AcSem.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_AcSem)) - (y_wholeset_bestGuess - ypred_AcSem)));
            Model.Semantic.Deviance(mm,aa) = sum(2*(y_wholeset .* (log((y_wholeset_bestGuess+(y_wholeset_bestGuess==0)) ./ ypred_Semantic)) - (y_wholeset_bestGuess - ypred_Semantic)));
        elseif strcmp(DISTR,'normal')
            Model.Acoustic.Deviance(mm,aa) =sum(power(y_wholeset - ypred_Acoustic,2));
            Model.AcSem.Deviance(mm,aa) =sum(power(y_wholeset - ypred_AcSem,2));
            Model.Semantic.Deviance(mm,aa) =sum(power(y_wholeset - ypred_Semantic,2));
        else
            fprintf('unrecognize %s distribution for deviance calculation', DISTR);
        end
        
        %% Calculate DF of optimal models
        % Find the real degrees of freedom for each best model
        % Note that in Matlab elastic net the parameters of Lasso (L1)
        % and Ridge (L2) can be retrieved from the parameters of the
        % elastic net A and L as: L2=L*(1-A)/2 and L1=L*A
        % The real DF is calculated from the formula as indicated in
        % (http://web.stanford.edu/~hastie/TALKS/enet_talk.pdf)
        % AcSem
        fprintf('Optimal AcSem model: Calculate exact number of parameters\n')
        L2_AcSem = Model.AcSem.Lambdas(mm,aa)*(1-Alpha)/2;
        %L1_AcSem = Lambda_BestModel_AcSem(bb,aa)*Alpha;
        Active_Predictors_AcSem = find(SemBoost.*B_AcSem);
        X_AcSem = [X_voc_wholeset x_wholeset];
        Xact_AcSem = X_AcSem(:,Active_Predictors_AcSem);
        Xact_AcSem_inv = (Xact_AcSem' * Xact_AcSem + L2_AcSem .* eye(length(Active_Predictors_AcSem)))^(-1);
        Model.AcSem.DF(mm,aa) = trace(Xact_AcSem * Xact_AcSem_inv * Xact_AcSem');
        
        % Acoustic
        fprintf('Acoustic model: Calculate exact number of parameters\n')
        L2_Acoustic = Model.Acoustic.Lambdas(mm,aa)*(1-Alpha)/2;
        %L1_Acoustic = Lambda_BestModel_Acoustic(bb,aa)*Alpha;
        Active_Predictors_Ac = find(B_Ac);
        Xact_Acoustic = x_wholeset(:,Active_Predictors_Ac);
        Xact_Acoustic_inv = (Xact_Acoustic' * Xact_Acoustic + L2_Acoustic .* eye(length(Active_Predictors_Ac)))^(-1);
        Model.Acoustic.DF(mm,aa) = trace(Xact_Acoustic * Xact_Acoustic_inv * Xact_Acoustic');
        
        % Semantic
        fprintf('Semantic model: save number of parameters\n')
        Model.Semantic.DF(mm,aa) = length(MDL_Sem.Coefficients.Estimate);
    end
    
    % Store all the deviance and DF of all bootstraps for all alpha values
    Deviance.Acoustic.bestvalues{mm} = Deviance_BestModel_Acoustic;
    Deviance.Acoustic.DF{mm} = DF_BestModel_Acoustic;
    Deviance.Acoustic.FitIndex{mm}=FitIndex_Acoustic;
    Deviance.Acoustic.FitIndexAR{mm}=FitIndexAR_Acoustic;
    LL.Acoustic.bestvalues{mm} = LL_BestModel_Acoustic;
    Deviance.AcSem.bestvalues{mm} = Deviance_BestModel_AcSem;
    Deviance.AcSem.DF{mm} = DF_BestModel_AcSem;
    Deviance.AcSem.FitIndex{mm}=FitIndex_AcSem;
    Deviance.AcSem.FitIndexAR{mm}=FitIndexAR_AcSem;
    LL.AcSem.bestvalues{mm} = LL_BestModel_AcSem;
    Deviance.Sem.DF{mm} = DF_BestModel_Sem;
    Deviance.Sem.bestvalues{mm}=Deviance_BestModel_Sem;
    Deviance.Sem.FitIndex{mm}=FitIndex_Sem;
    Deviance.Sem.FitIndexAR{mm}=FitIndexAR_Sem;
    LL.Sem.bestvalues{mm} = LL_BestModel_Sem;
    PropVal.values{mm} = ValSize;
    if OFFSET>1
        Deviance.AcSem2.bestvalues{mm} = Deviance_BestModel_AcSem2;
        Deviance.AcSem2.DF{mm} = DF_BestModel_AcSem2;
        Deviance.AcSem2.FitIndex{mm}=FitIndex_AcSem2;
        Deviance.AcSem2.FitIndexAR{mm}=FitIndexAR_AcSem2;
        LL.AcSem2.bestvalues{mm} = LL_BestModel_AcSem2;
    end
    LL.Ceiling.bestvalues{mm} = LL_BestGuess;
    LL.Floor.bestvalues{mm} = LL_WorseGuess;
    LL.AutoRegressive.bestvalues{mm} = LL_AR;
    % Calculate the probability for Deviance and AIC  AcSem < Acoustic
    % AcSem < Sem
    DiffDeviance_AcSemAcoustic = Deviance_BestModel_AcSem - Deviance_BestModel_Acoustic;
    DiffAIC_AcSemAcoustic = 2.*DF_BestModel_AcSem - 2.*DF_BestModel_Acoustic + DiffDeviance_AcSemAcoustic;
    DiffDeviance_AcSemSemantic = Deviance_BestModel_AcSem - Deviance_BestModel_Sem;
    DiffAIC_AcSemSemantic = 2.*DF_BestModel_AcSem - 2.*DF_BestModel_Sem + DiffDeviance_AcSemSemantic;
    PDeviance_AcSemAcoustic{mm} = nan(length(Alphas),1);
    PAIC_AcSemAcoustic{mm} = nan(length(Alphas),1);
    PDeviance_AcSemSemantic{mm} = nan(length(Alphas),1);
    PAIC_AcSemSemantic{mm}=nan(length(Alphas),1);
    for aa=1:length(Alphas)
        PDeviance_AcSemAcoustic{mm}(aa) = length(find(DiffDeviance_AcSemAcoustic(:,aa)<0))/BootstrapSTRF;
        PAIC_AcSemAcoustic{mm}(aa) = length(find(DiffAIC_AcSemAcoustic(:,aa)<0))/BootstrapSTRF;
        PDeviance_AcSemSemantic{mm}(aa) = length(find(DiffDeviance_AcSemSemantic(:,aa)<0))/BootstrapSTRF;
        PAIC_AcSemSemantic{mm}(aa) = length(find(DiffAIC_AcSemSemantic(:,aa)<0))/BootstrapSTRF;
    end
    
    %% keep track of the data and stim used for the next window
    Stim_local_old = Stim_local;
    y_dev_old = y_dev;
    
        
end


% Save data 
    save(fullfile(resultsDirectory,sprintf('%s_GLMPoisson.mat',Cellname)))
end

