function [MResults]=runLassoGlmModels1L(Alpha,VocType,Emitter,Stim_local,SWITCH,Data,ParamModel, PreviousWin,TickSTRFspectro)
if nargin<9
    TickSTRFspectro=[];
end
NbTrialStim=10;
VOC=VocType(Stim_local);
EEnames=Emitter.Ename(Stim_local);
bb=ParamModel.BootstrapSTRF;
if bb~=1
    error('!!!!Error: runLassoGlmModels1L only perform one iteration!!!\n !!!and the number of iterations requested =%d!!!!\n',bb);
end
    

% ParamModel.ModelChoice is a logic vector which indicates which models should
...be run (Acoustic, Semantic, AcSem without Offset,AcSem Semantic Offset, 
    ...AcSem Accoustic Offset) Note that the last two models need the 
    ...calculation of the first two models

%if ParamModel.ModelChoice(1)
    Deviance_All_local_Ac=cell(ParamModel.BootstrapSTRF,1);
    Lambda_All_local_Ac = cell(ParamModel.BootstrapSTRF,1);
    LL_All_local_Ac=cell(ParamModel.BootstrapSTRF,1);
    Ypredict_Acoustic = cell(ParamModel.BootstrapSTRF,1);
    Deviance_BestModel_Acoustic= nan(ParamModel.BootstrapSTRF,1);
    LL_BestModel_Acoustic =nan(ParamModel.BootstrapSTRF,1);
    FitIndex_Acoustic_local=nan(ParamModel.BootstrapSTRF,1);
    Lambda_BestModel_Acoustic=nan(ParamModel.BootstrapSTRF,1);
    ModelB_Ac = cell(ParamModel.BootstrapSTRF,1);% Parameters of the best model at that bootstrap
    ModelB0_Ac=nan(ParamModel.BootstrapSTRF,1);
    if SWITCH.DF
        DF_BestModel_Acoustic = nan(ParamModel.BootstrapSTRF,1);
    end
    Speed_Acoustic = nan(ParamModel.BootstrapSTRF,1);
%end
%if ParamModel.ModelChoice(2)
    Deviance_BestModel_Sem = nan(ParamModel.BootstrapSTRF,1);
    LL_BestModel_Sem = nan(ParamModel.BootstrapSTRF,1);
    Ypredict_Semantic = cell(ParamModel.BootstrapSTRF,1);
    ModelB_Sem = cell(ParamModel.BootstrapSTRF,1);
    FitIndex_Sem_local=nan(ParamModel.BootstrapSTRF,1);
    if SWITCH.DF
        DF_BestModel_Sem = nan(ParamModel.BootstrapSTRF,1);
    end
%end
%if ParamModel.ModelChoice(3)
    Deviance_All_local_AcSem=cell(ParamModel.BootstrapSTRF,1);
    Lambda_All_local_AcSem = cell(ParamModel.BootstrapSTRF,1);
    LL_All_local_AcSem=cell(ParamModel.BootstrapSTRF,1);
    Ypredict_AcSem = cell(ParamModel.BootstrapSTRF,1);
    Deviance_BestModel_AcSem = nan(ParamModel.BootstrapSTRF,1);
    LL_BestModel_AcSem = nan(ParamModel.BootstrapSTRF,1);
    FitIndex_AcSem_local = nan(ParamModel.BootstrapSTRF,1);
    Lambda_BestModel_AcSem = nan(ParamModel.BootstrapSTRF,1);
    ModelBspectro_AcSem = cell(ParamModel.BootstrapSTRF,1);% Parameters of the best model at that bootstrap
    ModelBsem_AcSem = cell(ParamModel.BootstrapSTRF,1);
    ModelB0_AcSem = nan(ParamModel.BootstrapSTRF,1);
    if SWITCH.DF
        DF_BestModel_AcSem = nan(ParamModel.BootstrapSTRF,1);
    end
%end
%if ParamModel.ModelChoice(4)
    Deviance_All_local_AcSemAc=cell(ParamModel.BootstrapSTRF,1);
    Lambda_All_local_AcSemAc = cell(ParamModel.BootstrapSTRF,1);
    LL_All_local_AcSemAc=cell(ParamModel.BootstrapSTRF,1);
    Ypredict_AcSemAc = cell(ParamModel.BootstrapSTRF,1);
    Deviance_BestModel_AcSemAc= nan(ParamModel.BootstrapSTRF,1);
    LL_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,1);
    Lambda_BestModel_AcSemAc=nan(ParamModel.BootstrapSTRF,1);
    FitIndex_AcSemAc_local=nan(ParamModel.BootstrapSTRF,1);
    ModelBspectro_AcSemAc = cell(ParamModel.BootstrapSTRF,1);% Parameters of the best model at that bootstrap
    ModelBsem_AcSemAc=cell(ParamModel.BootstrapSTRF,1);
    ModelB0_AcSemAc=nan(ParamModel.BootstrapSTRF,1);
    if SWITCH.DF
        DF_BestModel_AcSemAc = nan(ParamModel.BootstrapSTRF,1);
    end
%end
%if ParamModel.ModelChoice(5)
    Deviance_All_local_AcSemSem=cell(ParamModel.BootstrapSTRF,1);
    Lambda_All_local_AcSemSem = cell(ParamModel.BootstrapSTRF,1);
    Ypredict_AcSemSem = cell(ParamModel.BootstrapSTRF,1);
    LL_All_local_AcSemSem=cell(ParamModel.BootstrapSTRF,1);
    Deviance_BestModel_AcSemSem= nan(ParamModel.BootstrapSTRF,1);
    LL_BestModel_AcSemSem=nan(ParamModel.BootstrapSTRF,1);
    FitIndex_AcSemSem_local=nan(ParamModel.BootstrapSTRF,1);
    Lambda_BestModel_AcSemSem = nan(ParamModel.BootstrapSTRF,1);
    ModelBspectro_AcSemSem = cell(ParamModel.BootstrapSTRF,1);% Parameters of the best model at that bootstrap
    ModelBsem_AcSemSem = cell(ParamModel.BootstrapSTRF,1);
    ModelB0_AcSemSem = nan(ParamModel.BootstrapSTRF,1);
    if SWITCH.DF
        DF_BestModel_AcSemSem = nan(ParamModel.BootstrapSTRF,1);
    end
%end

if SWITCH.AR
    LL_ARVal = nan(ParamModel.BootstrapSTRF,1);
    %LL_AR = nan(ParamModel.BootstrapSTRF,1);
    Ypredict_ARVal = cell(ParamModel.BootstrapSTRF,1);
end
LL_BestGuess = nan(ParamModel.BootstrapSTRF,1);
LL_Saturated = nan(ParamModel.BootstrapSTRF,1);
LL_WorseGuess = nan(ParamModel.BootstrapSTRF,1);
YobsVal = cell(ParamModel.BootstrapSTRF,1);
ValSetFirstRep = cell(ParamModel.BootstrapSTRF,1);
NumStimTrials = nan(ParamModel.BootstrapSTRF,1);
XVocVal = cell(ParamModel.BootstrapSTRF,1);
YCeilModelVal = cell(ParamModel.BootstrapSTRF,1);
Ypredict_Floor = cell(ParamModel.BootstrapSTRF,1);
XSpecVal = cell(ParamModel.BootstrapSTRF,1);
ValProp = nan(ParamModel.BootstrapSTRF,1);
ValSize = nan(ParamModel.BootstrapSTRF,1);
ValSets = cell(ParamModel.BootstrapSTRF,1);
TrainSets = cell(ParamModel.BootstrapSTRF,1);
UVOC = unique(VOC);
if ParamModel.ZC==1
    X_meanZC = cell(ParamModel.BootstrapSTRF,1);
    X_stdZC = cell(ParamModel.BootstrapSTRF,1);
end

tic
%% Construct a testing dataset and a validating dataset if demanded
% remove one emitter per call category, if one category contain only
% one emitter, don't use it in the model
if ParamModel.CV
    [ValSet, TrainSet] = create_cross_validation_sets(VOC, EEnames);
else
    TrainSet=1:length(Stim_local);
    ValSet=TrainSet;
    if ParamModel.BootstrapSTRF>1
        fprintf(1,'Note that it does NOT make sense to run %d loops without Cross-validation!!!\n',ParamModel.BootstrapSTRF);
    end
end
ValProp(bb) = length(ValSet)./(length(ValSet)+length(TrainSet));
ValSize(bb) = length(ValSet);
ValSets{bb} = ValSet;
TrainSets{bb} = TrainSet;

%% z-score the spectrograms with the mean and STD of the training dataset if demanded
if any(ParamModel.ModelChoice([1,3:5]))
    if ParamModel.ZC==1
        % z-score spectrogram data
        x_Train_mean = mean(Data.x(TrainSet,:));
        x_Train_std = std(Data.x(TrainSet,:));
        x_Train = (Data.x(TrainSet,:)-repmat(x_Train_mean,length(TrainSet),1))./repmat(x_Train_std,length(TrainSet),1);
        x_Train(:,x_Train_std==0)=0;%set these parameters to zero as they do not cary any information and would give Inf or NaN
        x_Val = (Data.x(ValSet,:)-repmat(x_Train_mean,length(ValSet),1))./repmat(x_Train_std,length(ValSet),1);
        x_Val(:,x_Train_std==0)=0;%set these parameters to zero as they do not cary any information and would give Inf or NaN
        if SWITCH.FIG>1
            for ss = 1:length(TrainSet)
                figure(315)
                STIM=reshape(x_Train(ss,:), length(TickSTRFspectro.fo), length(TickSTRFspectro.to));
                imagesc(TickSTRFspectro.to*1000,TickSTRFspectro.fo, STIM)
                axis xy
                xlabel('Time (ms)')
                ylabel('Frequencies')
                title('TrainStim after z-score')
                pause(0.5)
            end
        end
    else
        x_Train = Data.x(TrainSet,:);
        x_Val = Data.x(ValSet,:);
    end
end
%% construct the vector of responses, the vectors of spectrograms, vocalization types and previous spike count for the training dataset
if strcmp(ParamModel.NeuroRes, 'count')
    yy=0;
    if any(ParamModel.ModelChoice([1,3:5]))
        x_Train_new = nan(NbTrialStim.*length(TrainSet),size(x_Train,2));
    end
    X_voc_Train = nan(NbTrialStim.*length(TrainSet),size(Data.X_voc,2));
    y_Train=nan(NbTrialStim.*length(TrainSet),1);
    if ~isempty(PreviousWin.y_old)
        y_input_AR=nan(NbTrialStim.*length(TrainSet),3);
    end
    for TS=1:length(TrainSet)
        YT=Data.y{TrainSet(TS)};
        for tt=1:length(YT)
            yy=yy+1;
            y_Train(yy) = YT(tt);
            if any(ParamModel.ModelChoice([1,3:5]))
                x_Train_new(yy,:)=x_Train(TS,:);
            end
            X_voc_Train(yy,:) = Data.X_voc(TrainSet(TS),:);
            if ~isempty(PreviousWin.y_old)
                Local_trialset=YT;
                Local_trialset(tt)=[];
                y_input_AR(yy,1)=mean(Local_trialset);
                Stim_num_real=Stim_local(TrainSet(TS));
                Stim_num_old=find(PreviousWin.Stim_local_old==Stim_num_real);
                Local_trialset=PreviousWin.y_old{Stim_num_old};
                y_input_AR(yy,3)=Local_trialset(tt);
                Local_trialset(tt)=[];
                y_input_AR(yy,2)=mean(Local_trialset);
            end
        end
    end
    if any(ParamModel.ModelChoice([1,3:5]))
        x_Train=x_Train_new(1:yy,:);
    end
    X_voc_Train=X_voc_Train(1:yy,:);
    y_Train=y_Train(1:yy);
    if ~isempty(PreviousWin.y_old)
        y_input_AR=y_input_AR(1:yy,:);
    end
else
    y_Train = Data.y(TrainSet);
    X_voc_Train = Data.X_voc(TrainSet,:);
end




%% First run Semantic glm Model and use the coefficients as an offset for the Acsem Model
%Semantic Model
if ParamModel.ModelChoice(2)
    fprintf(1,'Semantic Model poisson glm %d/%d\n', bb, ParamModel.BootstrapSTRF);
    MDL_Sem=fitglm(X_voc_Train,y_Train,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK);
    %MDL_Sem.Coefficients.Estimate=glmfit(X_voc_Train,y_Train,ParamModel.DISTR,'link',ParamModel.LINK);
    
    % Calculate the Semantic Offset
    if ParamModel.ModelChoice(5)
        Sem_Offset=sum(repmat(MDL_Sem.Coefficients.Estimate(2:end)',size(X_voc_Train,1),1).*X_voc_Train,2) + repmat(MDL_Sem.Coefficients.Estimate(1),size(X_voc_Train,1),1);
    end
end
%% Run Acoustic Lassoglm model
if ParamModel.ModelChoice(1)
    fprintf(1,'Acoustic Model lasso glm %d/%d\n', bb, ParamModel.BootstrapSTRF);
    if sum(sum(isnan(x_Train)))
        if ParamModel.ZC==1
            fprintf(1,'x_Train has %d NAN values and ParamModel.ZC=%d\n %d STD values of Data.x(TrainSet,:) equal zero\n',sum(sum(isnan(x_Train))),ParamModel.ZC,sum(x_Train_std==0));
        else
            fprintf(1,'x_Train has %d NAN values and ParamModel.ZC=%d\n',sum(sum(isnan(x_Train))),ParamModel.ZC);
        end
    end
    if sum(sum(isinf(x_Train)))
        fprintf(1,'x_Train has inf values and ParamModel.ZC=%d\n',ParamModel.ZC);
    end
    
    if isempty(ParamModel.LAMBDA{1})
        [B_Ac, FitInfo_Ac]=lassoglm_GDLim(x_Train,y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
    else
        [B_Ac, FitInfo_Ac]=lassoglm_GDLim(x_Train,y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'Standardize',0,'Lambda',ParamModel.LAMBDA{1});
    end
    Speed_Acoustic(bb)=FitInfo_Ac.elapsedTime;
    fprintf('Back to runLassoGlmModels after Acoustic Model\n')
    
    if SWITCH.FIG>2
        lassoPlot(B_Ac,FitInfo_Ac,'PlotType','Lambda','XScale','log')
        pause()
    end
end


%% Run Acoustic +Semantic Lassoglm model
if ParamModel.ModelChoice(3)
    fprintf(1,'Semantic + Acoustic Model lasso glm %d/%d\n', bb, ParamModel.BootstrapSTRF);
    if  isempty(ParamModel.LAMBDA{3})
        [B_AcSem_local, FitInfo_AcSem]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO);
        B_AcSem=B_AcSem_local;
    else
        [B_AcSem_local, FitInfo_AcSem]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link','Standardize',0,'Lambda',ParamModel.LAMBDA{3});
        B_AcSem=B_AcSem_local;
    end
    fprintf('Back to runLassoGlmModels After Semantic+Acoustic Model\n')
    if SWITCH.FIG>2
        lassoPlot(B_AcSem,FitInfo_AcSem,'PlotType','Lambda','XScale','log')
        pause()
    end
end

%% Run Acoustic +Semantic Lassoglm model with Semantic Offset
if ParamModel.ModelChoice(5)
    fprintf(1,'AcSem Model with Semantic offset lasso glm %d/%d\n', bb, ParamModel.BootstrapSTRF);
    if isempty(ParamModel.LAMBDA{5})
        [B_AcSemSem_local, FitInfo_AcSemSem_local]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO, 'Offset', Sem_Offset);
        fprintf('Back to runLassoGlmModels After AcSem model with semantic offset\n')
        % then we need to add the coefficients of the semantic
        % Model that were used to compute the offset
        B_AcSemSem=B_AcSemSem_local + repmat([MDL_Sem.Coefficients.Estimate(2:end)' zeros(1,size(B_AcSemSem_local,1)-length(MDL_Sem.Coefficients.Estimate(2:end)))]',1,size(B_AcSemSem_local,2));
        FitInfo_AcSemSem=FitInfo_AcSemSem_local;
        FitInfo_AcSemSem.Intercept=FitInfo_AcSemSem.Intercept + MDL_Sem.Coefficients.Estimate(1);
        
    else
        [B_AcSemSem_local, FitInfo_AcSemSem_local]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'Standardize',0,'Lambda',ParamModel.LAMBDA{5}, 'Offset', Sem_Offset);
        fprintf('Back to runLassoGlmModels After AcSem model with semantic offset\n')
        % then we need to add the coefficients of the semantic
        % Model that were used to compute the offset
        B_AcSemSem=B_AcSemSem_local + repmat([MDL_Sem.Coefficients.Estimate(2:end)' zeros(1,size(B_AcSemSem_local,1)-length(MDL_Sem.Coefficients.Estimate(2:end)))]',1,size(B_AcSemSem_local,2));
        FitInfo_AcSemSem=FitInfo_AcSemSem_local;
        FitInfo_AcSemSem.Intercept=FitInfo_AcSemSem.Intercept + MDL_Sem.Coefficients.Estimate(1);
    end
    if SWITCH.FIG>2
        lassoPlot(B_AcSemSem,FitInfo_AcSem,'PlotType','Lambda','XScale','log')
        pause()
    end
end



%% Create vector of obeserved values of the validating set, and predicted responses for Floor and ceiling model
if strcmp(ParamModel.NeuroRes, 'count')
    % construct the vector of responses for the validating
    % dataset (actual values, best guess and worse guess)
    yy=0;
    y_Val=nan(length(ValSet).*NbTrialStim,1);
    y_Val_bestGuess = nan(length(ValSet).*NbTrialStim,1);
    y_Val_worseGuess = nan(length(ValSet).*NbTrialStim,1);
    if ~isempty(PreviousWin.y_old)
        y_Val_input_AR = nan(length(ValSet).*NbTrialStim,3);
        y_Val_input_AR_size = nan(length(ValSet).*NbTrialStim,2);
    end
    
    for vv=1:length(ValSet)
        YT=Data.y{ValSet(vv)};
        for tt=1:length(YT)
            yy=yy+1;
            y_Val(yy)=YT(tt);
            YT_localmean = YT;
            YT_localmean(tt)=[];
            y_Val_bestGuess(yy)=mean(YT_localmean);
            if y_Val_bestGuess(yy)==0
                y_Val_bestGuess(yy)=1/(2*length(YT_localmean));% We can't predict a mean of 0. 1/2*(nb of trials-1) is our best guess of the minimum y
            end
            if ~isempty(PreviousWin.y_old)
                y_Val_input_AR(yy,1)=mean(YT_localmean);
                y_Val_input_AR_size(yy,1)=length(YT_localmean);
                Stim_num_real=Stim_local(ValSet(vv));
                Stim_num_old=find(PreviousWin.Stim_local_old==Stim_num_real);
                Local_trialset=PreviousWin.y_old{Stim_num_old};
                y_Val_input_AR(yy,3)=Local_trialset(tt);
                Local_trialset(tt)=[];
                y_Val_input_AR(yy,2)=mean(Local_trialset);
                y_Val_input_AR_size(yy,2)=length(Local_trialset);
            end
            Meansperstim_Cat_worseguess=cell(length(UVOC),1);
            for ww=1:length(ValSet)
                if ww~=vv
                    YTw=Data.y{ValSet(ww)};
                    MeanperStim_local=mean(YTw);
                else
                    MeanperStim_local=mean(YT_localmean);
                end
                for uu=1:length(UVOC)
                    if strcmp(VOC{ValSet(ww)},UVOC{uu})
                        Meansperstim_Cat_worseguess{uu}=[Meansperstim_Cat_worseguess{uu} MeanperStim_local];
                        break
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
    y_Val_worseGuess=y_Val_worseGuess(1:yy);% Floor Model: worse prediction(average of the mean spike counts per call category
    
    %% Run autoregressive model and calculate predicted values
    
    if ~isempty(PreviousWin.y_old) && SWITCH.AR
        %For the first window no auto-regresssive model for now
        % Auto-regressive model
        fprintf(1,'Auto-regressive Model poisson glm %d/%d\n', bb, ParamModel.BootstrapSTRF);
        y_Val_input_AR= y_Val_input_AR(1:yy,:);% Autoregressive model input y: validating input for the autoregressive model
        y_Val_input_AR_size= y_Val_input_AR_size(1:yy,:);
        %MDL_AR=fitglm(y_input_AR,y_Train,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK_AR,'intercept',false);
        %ypred_AR = glmval(MDL_AR.Coefficients.Estimate,y_Val_input_AR,ParamModel.LINK_AR, 'constant','off');
        
        %[MDL_AR.Coefficients.Estimate,MDL_AR.Deviance,MDL_AR.stats]=glmfit(y_input_AR,y_Train,ParamModel.DISTR,'link',ParamModel.LINK_AR,'constant','off');
        MDL_AR2=fitglm(y_Val_input_AR,y_Val,'Distribution',ParamModel.DISTR,'Link',ParamModel.LINK_AR,'intercept',false);
        %[MDL_AR2.Coefficients.Estimate,MDL_AR2.Deviance,MDL_AR2.stats]=glmfit(y_Val_input_AR,y_Val,ParamModel.DISTR,'link',ParamModel.LINK_AR,'constant','off');
        % Calculate the minimum spike count that the model should
        % give given the input variability to decide of a minimum
        % value for the predicted spike count with the
        % auto-regressive model
        y_predict_ARVal_Min=glmval(MDL_AR2.Coefficients.Estimate,[(1/2).^y_Val_input_AR_size(:,1) (1/2).^y_Val_input_AR_size(:,2) zeros(size(y_Val_input_AR_size,1),1)],ParamModel.LINK_AR, 'constant','off');
        y_predict_ARVal_Optimum=max(MDL_AR2.Fitted.Response,y_predict_ARVal_Min);
        
    end
    
end


%% Calculate Deviances, loglikelihood and predicted values for validating dataset
% Deviance of Ceiling model (best Guess), floor model (Worseguess) and Acoustic model for the validating dataset for each lambda

% LogLikelihood of ceiling and floor models
LL_Saturated(bb) = LL_Calculus(y_Val,y_Val);
LL_BestGuess(bb) = LL_Calculus(y_Val,y_Val_bestGuess);
LL_WorseGuess(bb) = LL_Calculus(y_Val,y_Val_worseGuess);
if ParamModel.ModelChoice(1)
    fprintf(1,'Acoustic Model: Calculate deviance for each lambda on the validating dataset\n');
    NbL_Ac=length(FitInfo_Ac.Lambda);% number of lambdas tested
    Deviance_M_Ac= nan(NbL_Ac,1);
    LL_Ac= nan(NbL_Ac,1);
    if strcmp(ParamModel.NeuroRes, 'count')
        % Construct the vector of the predicted response using the
        % Acoustic model
        ypred_Ac=nan(length(ValSet).*NbTrialStim,NbL_Ac);
        for ll=1:NbL_Ac
            yy=0;
            ypred_Ac_temp = glmval([FitInfo_Ac.Intercept(ll); B_Ac(:,ll)],x_Val,ParamModel.LINK);
            for vv=1:length(ValSet)
                YT=Data.y{ValSet(vv)};
                for tt=1:length(YT)
                    yy=yy+1;
                    ypred_Ac(yy,ll)=ypred_Ac_temp(vv);
                end
            end
            if strcmp(ParamModel.DISTR,'poisson')
                LL_Ac(ll) = LL_Calculus(y_Val,ypred_Ac(1:yy,ll));
                Deviance_M_Ac(ll) = -2*(LL_Ac(ll)-LL_Saturated(bb));
            elseif strcmp(ParamModel.DISTR,'normal')
                Deviance_M_Ac(ll) = sum(power(y_Val - ypred_Ac(1:yy,ll),2));
            else
                fprintf(1,'your distribution is unrecognize');
            end
        end
        ypred_Ac = ypred_Ac(1:yy,:);
        LL_All_local_Ac{bb}=LL_Ac;
    else
        y_Val=Data.y(ValSet);
        for ll=1:NbL_Ac
            ypred_Ac(:,ll) = glmval([FitInfo_Ac.Intercept(ll); B_Ac(:,ll)],x_Val,ParamModel.LINK);
            if strcmp(ParamModel.DISTR,'poisson')
                Deviance_M_Ac(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_Ac(:,ll))) - (y_Val - ypred_Ac(:,ll))));
            elseif strcmp(ParamModel.DISTR,'normal')
                Deviance_M_Ac(ll) = sum(power(y_Val - ypred_Ac(:,ll),2));
            else
                fprintf(1,'your distribution is unrecognize');
            end
        end
    end
    Deviance_All_local_Ac{bb} = Deviance_M_Ac;% Store All deviance values for that bootstrap
    Lambda_All_local_Ac{bb}=FitInfo_Ac.Lambda;% Store all Lambda tested for that booststrap
end

% Deviance of Acoustic + Semantic model for the validating dataset for each lambda
if ParamModel.ModelChoice(3)
    fprintf(1,'AcSem model: Calculate deviance for each lambda on the validating dataset\n');
    NbL_AcSem=length(FitInfo_AcSem.Lambda);% number of lambdas tested
    Deviance_M_AcSem= nan(NbL_AcSem,1);
    LL_AcSem=nan(NbL_AcSem,1);
    if strcmp(ParamModel.NeuroRes, 'count')
        % Construct the vector of predicted responses using the
        % AcSem1 model
        ypred_AcSem=nan(length(ValSet).*NbTrialStim,NbL_AcSem);
        for ll=1:NbL_AcSem
            ypred_AcSem_temp = glmval([FitInfo_AcSem.Intercept(ll); B_AcSem(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
            yy=0;
            for vv=1:length(ValSet)
                for tt=1:length(Data.y{ValSet(vv)})
                    yy=yy+1;
                    ypred_AcSem(yy,ll)=ypred_AcSem_temp(vv);
                end
            end
            if strcmp(ParamModel.DISTR,'poisson')
                LL_AcSem(ll) = LL_Calculus(y_Val,ypred_AcSem(1:yy,ll));
                Deviance_M_AcSem(ll) = -2*(LL_AcSem(ll)-LL_Saturated(bb));
            elseif strcmp(ParamModel.DISTR,'normal')
                Deviance_M_AcSem(ll)=sum(power(y_Val - ypred_AcSem(1:yy,ll),2));
            else
                fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
            end
        end
        ypred_AcSem=ypred_AcSem(1:yy,:);
        LL_All_local_AcSem{bb}=LL_AcSem;
    else
        ypred_AcSem(:,ll)=glmval([FitInfo_AcSem.Intercept(ll); B_AcSem(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
        if strcmp(ParamModel.DISTR,'poisson')
            Deviance_M_AcSem(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_AcSem(:,ll))) - (y_Val - ypred_AcSem(:,ll))));
        elseif strcmp(ParamModel.DISTR,'normal')
            Deviance_M_AcSem(ll)=sum(power(y_Val - ypred_AcSem(:,ll),2));
        else
            fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
        end
    end
    Deviance_All_local_AcSem{bb} = Deviance_M_AcSem;% Store All deviance values for that bootstrap
    Lambda_All_local_AcSem{bb}=FitInfo_AcSem.Lambda;% Store all Lambda tested for that booststrap
end

% Deviance of Semantic model for the validating dataset
if ParamModel.ModelChoice(2)
    fprintf(1,'Semantic model: Calculate deviance on the validating dataset\n');
    ypred_Sem_temp = glmval(MDL_Sem.Coefficients.Estimate,Data.X_voc(ValSet,:),ParamModel.LINK);
    if strcmp(ParamModel.NeuroRes, 'count')
        % construct the vector of predicted responses using the semantic model
        yy=0;
        ypred_Sem=nan(length(ValSet).*NbTrialStim,1);
        for vv=1:length(ValSet)
            for tt=1:length(Data.y{ValSet(vv)})
                yy=yy+1;
                ypred_Sem(yy)=ypred_Sem_temp(vv);
            end
        end
        ypred_Sem=ypred_Sem(1:yy);
    else
        ypred_Sem=ypred_Sem_temp;
    end
    if strcmp(ParamModel.DISTR,'poisson')
        LL_BestModel_Sem(bb) = LL_Calculus(y_Val, ypred_Sem);
        Deviance_BestModel_Sem(bb) = -2*(LL_BestModel_Sem(bb)-LL_Saturated(bb));
        FitIndex_Sem_local(bb)=(LL_BestModel_Sem(bb)-LL_WorseGuess(bb))/(LL_BestGuess(bb)-LL_WorseGuess(bb));
    elseif strcmp(ParamModel.DISTR,'normal')
        Deviance_BestModel_Sem(bb) =sum(power(y_Val - ypred_Sem,2));
    else
        fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
    end
    ModelB_Sem{bb} = MDL_Sem.Coefficients.Estimate;% Parameters of the model at that bootstrap
    Ypredict_Semantic{bb} = ypred_Sem; % Store the predicted values for that bootstrap
end

% Deviance of Acoustic + Semantic model with Semantic Offset for the validating dataset for each lambda
if ParamModel.ModelChoice(5)
    fprintf(1,'AcSem model Sem Offset: Calculate deviance for each lambda on the validating dataset\n');
    NbL_AcSemSem=length(FitInfo_AcSemSem.Lambda);% number of lambdas tested
    Deviance_M_AcSemSem= nan(NbL_AcSemSem,1);
    LL_AcSemSem=nan(NbL_AcSemSem,1);
    if strcmp(ParamModel.NeuroRes, 'count')
        % Construct the vector of predicted responses using the
        % AcSemSem model
        ypred_AcSemSem=nan(length(ValSet).*NbTrialStim,NbL_AcSemSem);
        for ll=1:NbL_AcSemSem
            ypred_AcSemSem_temp = glmval([FitInfo_AcSemSem.Intercept(ll); B_AcSemSem(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
            yy=0;
            for vv=1:length(ValSet)
                for tt=1:length(Data.y{ValSet(vv)})
                    yy=yy+1;
                    ypred_AcSemSem(yy,ll)=ypred_AcSemSem_temp(vv);
                end
            end
            if strcmp(ParamModel.DISTR,'poisson')
                LL_AcSemSem(ll) = LL_Calculus(y_Val,ypred_AcSemSem(1:yy,ll));
                Deviance_M_AcSemSem(ll) = -2*(LL_AcSemSem(ll)-LL_Saturated(bb));
            elseif strcmp(ParamModel.DISTR,'normal')
                Deviance_M_AcSemSem(ll)=sum(power(y_Val - ypred_AcSemSem(1:yy,ll),2));
            else
                fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
            end
        end
        ypred_AcSemSem=ypred_AcSemSem(1:yy,:);
        LL_All_local_AcSemSem{bb}=LL_AcSemSem;
    else
        ypred_AcSemSem(:,ll)=glmval([FitInfo_AcSemSem.Intercept(ll); B_AcSemSem(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
        if strcmp(ParamModel.DISTR,'poisson')
            Deviance_M_AcSemSem(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_AcSemSem(:,ll))) - (y_Val - ypred_AcSemSem(:,ll))));
        elseif strcmp(ParamModel.DISTR,'normal')
            Deviance_M_AcSemSem(ll)=sum(power(y_Val - ypred_AcSemSem(:,ll),2));
        else
            fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
        end
    end
    Deviance_All_local_AcSemSem{bb} = Deviance_M_AcSemSem;% Store All deviance values for that bootstrap
    Lambda_All_local_AcSemSem{bb}=FitInfo_AcSemSem.Lambda;% Store all Lambda tested for that booststrap
end

% LogLikelihood of the autoregressive model
if SWITCH.AR
    fprintf(1,'Auto-regressive model: Calculate LogLikelihood on the validating dataset\n');
    %LL_AR(bb) = LL_Calculus(y_Val,ypred_AR);
    if ~isempty(PreviousWin.y_old)
        %Calculate loglikelihood
        LL_ARVal(bb)=LL_Calculus(y_Val,y_predict_ARVal_Optimum);
    else
        LL_ARVal(bb)=LL_BestGuess(bb); % Set the value to BestGuess model value when mm=1
    end
    if SWITCH.FIG>1 && ~isempty(PreviousWin.y_old)
        figure(301)
        plot(y_Val,y_Val_bestGuess,'r.','MarkerSize',10)
        hold on
        plot(y_Val,y_predict_ARVal_Optimum,'b.','MarkerSize',10)
        plot(y_Val,MDL_AR2.Fitted.Response,'g.','MarkerSize',10)
        legend('Ceiling Model','Auto-regressive model training','Auto-regressive Model fitted values')
        title(sprintf('LL Ceiling Model=%f\nLL Auto-Regressive Model Val  glmeval fit=%f\nLL Auto-Regressive Model Val fit=%f',LL_BestGuess(bb),LL_ARVal(bb),LL_ARVal(bb)))
        xlabel('Observed spike count')
        ylabel('Predicted spike count')
    end
end
%% Store the smallest Deviance Value, the corresponding best lambda and the corresponding best model
% of best deviance, the parameters of the model for that best lambda for that booststrap

if ParamModel.ModelChoice(1)
    fprintf(1,'Acoustic model: Find lambda with smallest value of Deviance\n');
    Dev_min_index_Ac=find(Deviance_M_Ac==min(Deviance_M_Ac));
    Dev_min_index_Ac=Dev_min_index_Ac(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
    Deviance_BestModel_Acoustic(bb)= Deviance_M_Ac(Dev_min_index_Ac);
    LL_BestModel_Acoustic(bb)=LL_Ac(Dev_min_index_Ac);
    FitIndex_Acoustic_local(bb)=(LL_BestModel_Acoustic(bb)-LL_WorseGuess(bb))/(LL_BestGuess(bb)-LL_WorseGuess(bb));
    if SWITCH.AR
        % FitIndexAR_Acoustic_local=(LL_BestModel_Acoustic(bb)-LL_WorseGuess(bb))/(LL_AR(bb)-LL_WorseGuess(bb));
        FitIndexARVal_Acoustic_local=(LL_BestModel_Acoustic(bb)-LL_WorseGuess(bb))/(LL_ARVal(bb)-LL_WorseGuess(bb));
    end
    Lambda_BestModel_Acoustic(bb)=FitInfo_Ac.Lambda(Dev_min_index_Ac);
    if ParamModel.ZC==1
        ModelB_Ac{bb} = reshape(B_Ac(:,Dev_min_index_Ac)'.*x_Train_std,Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    else
        ModelB_Ac{bb} = reshape(B_Ac(:,Dev_min_index_Ac)',Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    end
    ModelB0_Ac(bb) = FitInfo_Ac.Intercept(Dev_min_index_Ac);% intercept of the best model
    Ypredict_Acoustic{bb} = ypred_Ac(:,Dev_min_index_Ac);% Store the predicted values of the model for that bootstrap   
end

if ParamModel.ModelChoice(2) && SWITCH.AR
    fprintf('Semantic Model: Calculate local FitIndex based on autoregressive model');
    %FitIndexAR_Sem_local=(LL_BestModel_Sem(bb)-LL_WorseGuess(bb))/(LL_AR(bb)-LL_WorseGuess(bb));
    FitIndexARVal_Sem_local=(LL_BestModel_Sem(bb)-LL_WorseGuess(bb))/(LL_ARVal(bb)-LL_WorseGuess(bb));
end

if ParamModel.ModelChoice(3)
    fprintf(1,'AcSem model: Find lambda with smallest value of Deviance\n');
    Dev_min_index_AcSem=find(Deviance_M_AcSem==min(Deviance_M_AcSem));
    Dev_min_index_AcSem=Dev_min_index_AcSem(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
    Deviance_BestModel_AcSem(bb)= Deviance_M_AcSem(Dev_min_index_AcSem);
    LL_BestModel_AcSem(bb)=LL_AcSem(Dev_min_index_AcSem);
    FitIndex_AcSem_local(bb)=(LL_BestModel_AcSem(bb)-LL_WorseGuess(bb))/(LL_BestGuess(bb)-LL_WorseGuess(bb));
    if SWITCH.AR
        %FitIndexAR_AcSem_local=(LL_BestModel_AcSem(bb)-LL_WorseGuess(bb))/(LL_AR(bb)-LL_WorseGuess(bb));
        FitIndexARVal_AcSem_local=(LL_BestModel_AcSem(bb)-LL_WorseGuess(bb))/(LL_ARVal(bb)-LL_WorseGuess(bb));
    end
    Lambda_BestModel_AcSem(bb)=FitInfo_AcSem.Lambda(Dev_min_index_AcSem);
    if ParamModel.ZC==1
        ModelBspectro_AcSem{bb} = reshape(B_AcSem((size(Data.X_voc,2)+1):end,Dev_min_index_AcSem)'.*x_Train_std,Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    else
        ModelBspectro_AcSem{bb} = reshape(B_AcSem((size(Data.X_voc,2)+1):end,Dev_min_index_AcSem)',Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    end
    ModelBsem_AcSem{bb} = B_AcSem(1:size(Data.X_voc,2),Dev_min_index_AcSem);
    % identify the category for which there was no stim and
    % calculate our best guess for B for that category: the average
    % of the B of other categories
    noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
    if ~isempty(noModCat)
        if ~sum(B_AcSem(noModCat,Dev_min_index_AcSem))==0
            fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
        end
        ModelBsem_AcSem{bb}(noModCat)=sum(B_AcSem(1:size(Data.X_voc,2),Dev_min_index_AcSem))/(length(B_AcSem(1:size(Data.X_voc,2),Dev_min_index_AcSem))-length(noModCat));
    end
    ModelB0_AcSem(bb) = FitInfo_AcSem.Intercept(Dev_min_index_AcSem);% intercept of the best model
    Ypredict_AcSem{bb} = ypred_AcSem(:,Dev_min_index_AcSem);%Store predicted values for that bootstrap
end

if ParamModel.ModelChoice(5)
    fprintf(1,'AcSem model Semantic Offset: Find lambda with smallest value of Deviance\n');
    Dev_min_index_AcSemSem=find(Deviance_M_AcSemSem==min(Deviance_M_AcSemSem));
    Dev_min_index_AcSemSem=Dev_min_index_AcSemSem(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
    Deviance_BestModel_AcSemSem(bb)= Deviance_M_AcSemSem(Dev_min_index_AcSemSem);
    LL_BestModel_AcSemSem(bb)=LL_AcSemSem(Dev_min_index_AcSemSem);
    FitIndex_AcSemSem_local(bb)=(LL_BestModel_AcSemSem(bb)-LL_WorseGuess(bb))/(LL_BestGuess(bb)-LL_WorseGuess(bb));
    if SWITCH.AR
        %FitIndexAR_AcSemSem_local=(LL_BestModel_AcSemSem(bb)-LL_WorseGuess(bb))/(LL_AR(bb)-LL_WorseGuess(bb));
        FitIndexARVal_AcSemSem_local=(LL_BestModel_AcSemSem(bb)-LL_WorseGuess(bb))/(LL_ARVal(bb)-LL_WorseGuess(bb));
    end
    Lambda_BestModel_AcSemSem(bb)=FitInfo_AcSemSem.Lambda(Dev_min_index_AcSemSem);
    if ParamModel.ZC==1
        ModelBspectro_AcSemSem{bb} = reshape(B_AcSemSem((size(Data.X_voc,2)+1):end,Dev_min_index_AcSemSem)'.*x_Train_std,Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    else
        ModelBspectro_AcSemSem{bb} = reshape(B_AcSemSem((size(Data.X_voc,2)+1):end,Dev_min_index_AcSemSem)',Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    end
    ModelBsem_AcSemSem{bb} = B_AcSemSem(1:size(Data.X_voc,2),Dev_min_index_AcSemSem);
    % identify the category for which there was no stim and
    % calculate our best guess for B for that category: the average
    % of the B of other categories
    noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
    if ~isempty(noModCat)
        if ~sum(B_AcSemSem(noModCat,Dev_min_index_AcSemSem))==0
            fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
        end
        ModelBsem_AcSemSem{bb}(noModCat)=sum(B_AcSemSem(1:size(Data.X_voc,2),Dev_min_index_AcSemSem))/(length(B_AcSemSem(1:size(Data.X_voc,2),Dev_min_index_AcSemSem))-length(noModCat));
    end
    ModelB0_AcSemSem(bb) = FitInfo_AcSemSem.Intercept(Dev_min_index_AcSemSem);% intercept of the best model
    Ypredict_AcSemSem{bb} = ypred_AcSemSem(:,Dev_min_index_AcSemSem);%Store predicted values for that bootstrap
end

%% Now that we now which Acoustic model is the best, calculate the offset of the best acoustic model to be used by the AcSemAc model
if ParamModel.ModelChoice(4)
    % Calculate the Acoustic Offset
    Ac_Offset=x_Train*B_Ac(:,Dev_min_index_Ac) + repmat(ModelB0_Ac(bb),size(x_Train,1),1);
    fprintf(1,'AcSem Model lasso glm with Acoustic Offset %d/%d\n', bb, ParamModel.BootstrapSTRF);
    % Calculate The AcSemAc model
    if isempty(ParamModel.LAMBDA{4})
        [B_AcSem_local2, FitInfo_AcSem_local2]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'NumLambda',ParamModel.NUMLAMBDA,'Standardize',0,'LambdaRatio',ParamModel.LAMBDARATIO, 'Offset', Ac_Offset);
    else
        [B_AcSem_local2, FitInfo_AcSem_local2]=lassoglm_GDLim([X_voc_Train x_Train],y_Train,ParamModel.DISTR,'Alpha',Alpha,'Link',ParamModel.LINK,'Standardize',0,'Lambda',ParamModel.LAMBDA{4}, 'Offset', Ac_Offset);
    end
    fprintf('Back to runLassoGlmModels After AcSem model with acoustic offset\n')
    % then we need to add the coefficients of the acoustic
    % Model that were used to compute the offset
    B_AcSemAc=B_AcSem_local2 + repmat([zeros(length(MDL_Sem.Coefficients.Estimate(2:end)),1); B_Ac(:,Dev_min_index_Ac)],1,size(B_AcSem_local2,2));
    FitInfo_AcSemAc=FitInfo_AcSem_local2;
    FitInfo_AcSemAc.Intercept=FitInfo_AcSemAc.Intercept + FitInfo_Ac.Intercept(Dev_min_index_Ac);
    
    % Deviance of Acoustic + Semantic model with Acoustic Offset for the validating dataset for each lambda
    fprintf(1,'AcSem model Ac Offset: Calculate deviance for each lambda on the validating dataset\n');
    NbL_AcSemAc=length(FitInfo_AcSemAc.Lambda);% number of lambdas tested
    Deviance_M_AcSemAc= nan(NbL_AcSemAc,1);
    LL_AcSemAc = nan(NbL_AcSemAc,1);
    if strcmp(ParamModel.NeuroRes, 'count')
        ypred_AcSemAc=nan(length(ValSet).*NbTrialStim,NbL_AcSemAc);
        for ll=1:NbL_AcSemAc
            ypred_AcSem2_temp = glmval([FitInfo_AcSemAc.Intercept(ll); B_AcSemAc(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
            yy=0;
            for vv=1:length(ValSet)
                for tt=1:length(Data.y{ValSet(vv)})
                    yy=yy+1;
                    ypred_AcSemAc(yy,ll)=ypred_AcSem2_temp(vv);
                end
            end
            if strcmp(ParamModel.DISTR,'poisson')
                LL_AcSemAc(ll) = LL_Calculus(y_Val,ypred_AcSemAc(1:yy,ll));
                Deviance_M_AcSemAc(ll) = -2*(LL_AcSemAc(ll)-LL_Saturated(bb));
            elseif strcmp(ParamModel.DISTR,'normal')
                Deviance_M_AcSemAc(ll)=sum(power(y_Val - ypred_AcSemAc(1:yy,ll),2));
            else
                fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
            end
        end
        ypred_AcSemAc=ypred_AcSemAc(1:yy,:);
        LL_All_local_AcSemAc{bb}=LL_AcSemAc;
    else
        ypred_AcSemAc(:,ll)=glmval([FitInfo_AcSemAc.Intercept(ll); B_AcSemAc(:,ll)],[Data.X_voc(ValSet,:) x_Val],ParamModel.LINK);
        if strcmp(ParamModel.DISTR,'poisson')
            Deviance_M_AcSemAc(ll) = sum(2*(y_Val .* (log((y_Val+(y_Val==0)) ./ ypred_AcSemAc(:,ll))) - (y_Val - ypred_AcSemAc(:,ll))));
        elseif strcmp(ParamModel.DISTR,'normal')
            Deviance_M_AcSemAc(ll)=sum(power(y_Val - ypred_AcSemAc(:,ll),2));
        else
            fprintf('unrecognize %s distribution for deviance calculation', ParamModel.DISTR);
        end
    end
    Deviance_All_local_AcSemAc{bb} = Deviance_M_AcSemAc;% Store All deviance values for that bootstrap
    Lambda_All_local_AcSemAc{bb}=FitInfo_AcSemAc.Lambda;% Store all Lambda tested for that booststrap
    % Calculate the smallest deviance
    fprintf(1,'AcSemAc model: Find lambda with smallest value of Deviance\n');
    Dev_min_index_AcSemAc=find(Deviance_M_AcSemAc==min(Deviance_M_AcSemAc));
    Dev_min_index_AcSemAc=Dev_min_index_AcSemAc(end);%choose the last lambda that gives the smallest deviance, which should be the biggest
    Deviance_BestModel_AcSemAc(bb)= Deviance_M_AcSemAc(Dev_min_index_AcSemAc);
    LL_BestModel_AcSemAc(bb)=LL_AcSemAc(Dev_min_index_AcSemAc);
    Lambda_BestModel_AcSemAc(bb)=FitInfo_AcSemAc.Lambda(Dev_min_index_AcSemAc);
    FitIndex_AcSemAc_local(bb)=(LL_BestModel_AcSemAc(bb)-LL_WorseGuess(bb))/(LL_BestGuess(bb)-LL_WorseGuess(bb));
    if SWITCH.AR
        %FitIndexAR_AcSemAc_local=(LL_BestModel_AcSemAc(bb)-LL_WorseGuess(bb))/(LL_AR(bb)-LL_WorseGuess(bb));
        FitIndexARVal_AcSemAc_local=(LL_BestModel_AcSemAc(bb)-LL_WorseGuess(bb))/(LL_ARVal(bb)-LL_WorseGuess(bb));
    end
    if ParamModel.ZC==1
        ModelBspectro_AcSemAc{bb} = reshape(B_AcSemAc((size(Data.X_voc,2)+1):end,Dev_min_index_AcSemAc)'.*x_Train_std,Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    else
        ModelBspectro_AcSemAc{bb} = reshape(B_AcSemAc((size(Data.X_voc,2)+1):end,Dev_min_index_AcSemAc)',Data.Df,Data.Dt);% Parameters of the best model at that bootstrap
    end
    ModelBsem_AcSemAc{bb} = B_AcSemAc(1:size(Data.X_voc,2),Dev_min_index_AcSemAc);
    % identify the category for which there was no stim and
    % calculate our best guess for B for that category: the average
    % of the B of other categories
    noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
    if ~isempty(noModCat)
        if ~sum(B_AcSemAc(noModCat,Dev_min_index_AcSemAc))==0
            fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
        end
        ModelBsem_AcSemAc{bb}(noModCat)=sum(B_AcSemAc(1:size(Data.X_voc,2),Dev_min_index_AcSemAc))/(length(B_AcSemAc(1:size(Data.X_voc,2),Dev_min_index_AcSemAc))-length(noModCat));
    end
    ModelB0_AcSemAc(bb) = FitInfo_AcSemAc.Intercept(Dev_min_index_AcSemAc);% intercept of the best model
    % Save predictions
    Ypredict_AcSemAc{bb} = ypred_AcSemAc(:,Dev_min_index_AcSemAc);
end


%% Store the data category and the values of the validating dataset
YobsVal{bb} = y_Val;
NumStimTrials(bb) = length(y_Val);
YCeilModelVal{bb} = y_Val_bestGuess;
if ~isempty(PreviousWin.y_old) && SWITCH.AR
    Ypredict_ARVal{bb} = y_predict_ARVal_Optimum;
end
Ypredict_Floor{bb} = y_Val_worseGuess;
if ParamModel.ZC==1
    X_meanZC{bb} = x_Train_mean;
    X_stdZC{bb} = x_Train_std;
end
if strcmp(ParamModel.NeuroRes, 'count')
    VocTypeTrial=cell(length(ValSet).*NbTrialStim,1);
    if any(ParamModel.ModelChoice([1,3:5]))
        xValperTrial=nan(length(ValSet).*NbTrialStim,size(x_Val,2));
    end
    yy=0;
    for vv=1:length(ValSet)
        for tt=1:length(Data.y{ValSet(vv)})
            yy=yy+1;
            VocTypeTrial{yy}=VOC{ValSet(vv)};
            if any(ParamModel.ModelChoice([1,3:5]))
                xValperTrial(yy,:)=x_Val(vv,:);
            end
        end
    end
    XVocVal{bb}=VocTypeTrial(1:yy,:);
    if any(ParamModel.ModelChoice([1,3:5]))
        XSpecVal{bb} = xValperTrial(1:yy,:);
    end
else
    XVocVal{bb} = VOC{ValSet};
    if any(ParamModel.ModelChoice([1,3:5]))
        XSpecVal{bb} = x_Val;
    end
end

% Determine what are the indices of the prediction of each of the first
% repetition in the vector of prediction that take into accounts all stims and all repetitions
ValSetFirstRep{bb} = nan(length(ValSet),1);
TT=0;
for vv=1:length(ValSet)
    if vv>1
        TT=length(Data.y{ValSet(vv-1)})+TT;
    end
    ValSetFirstRep{bb}(vv)= TT +1;
end
%% Calculate degress of freedom for each model
% Find the real degrees of freedom for each best model
% Note that in Matlab elastic net the parameters of Lasso (L1)
% and Ridge (L2) can be retrieved from the parameters of the
% elastic net A and L as: L2=L*(1-A)/2 and L1=L*A
% The real DF is calculated from the formula as indicated in
% (http://web.stanford.edu/~hastie/TALKS/enet_talk.pdf)
if ParamModel.ModelChoice(1) && SWITCH.DF
    % Acoustic
    fprintf('Acoustic model: Calculate exact number of parameters\n')
    L2_Acoustic = Lambda_BestModel_Acoustic(bb)*(1-Alpha)/2;
    %L1_Acoustic = Lambda_BestModel_Acoustic(bb,aa)*Alpha;
    Active_Predictors_Ac = find(B_Ac(:,Dev_min_index_Ac));
    Xact_Acoustic = x_Train(:,Active_Predictors_Ac);
    Xact_Acoustic_inv = (Xact_Acoustic' * Xact_Acoustic + L2_Acoustic .* eye(length(Active_Predictors_Ac)))^(-1);
    DF_BestModel_Acoustic(bb) = trace(Xact_Acoustic * Xact_Acoustic_inv * Xact_Acoustic');
end
if ParamModel.ModelChoice(2) && SWITCH.DF
    % Semantic
    fprintf('Semantic model: save number of parameters\n')
    DF_BestModel_Sem(bb) = length(ModelB_Sem{bb});
end
if ParamModel.ModelChoice(3) && SWITCH.DF
    % AcSem
    fprintf('AcSem model: Calculate exact number of parameters\n')
    L2_AcSem = Lambda_BestModel_AcSem(bb)*(1-Alpha)/2;
    %L1_AcSem = Lambda_BestModel_AcSem(bb,aa)*Alpha;
    Active_Predictors_AcSem = find(B_AcSem(:,Dev_min_index_AcSem));
    X_AcSem = [X_voc_Train x_Train];
    Xact_AcSem = X_AcSem(:,Active_Predictors_AcSem);
    Xact_AcSem_inv = (Xact_AcSem' * Xact_AcSem + L2_AcSem .* eye(length(Active_Predictors_AcSem)))^(-1);
    DF_BestModel_AcSem(bb) = trace(Xact_AcSem * Xact_AcSem_inv * Xact_AcSem');
end
if ParamModel.ModelChoice(4) && SWITCH.DF
    %Calculate degrees of freedom
    % AcSemAc
    fprintf('AcSem model Ac Offset: Calculate exact number of parameters\n')
    L2_AcSemAc = Lambda_BestModel_AcSemAc(bb)*(1-Alpha)/2;
    Active_Predictors_AcSemAc = find(B_AcSemAc(:,Dev_min_index_AcSemAc));
    X_AcSemAc = [X_voc_Train x_Train];
    Xact_AcSemAc = X_AcSemAc(:,Active_Predictors_AcSemAc);
    Xact_AcSem_inv2 = (Xact_AcSemAc' * Xact_AcSemAc + L2_AcSemAc .* eye(length(Active_Predictors_AcSemAc)))^(-1);
    DF_BestModel_AcSemAc(bb) = trace(Xact_AcSemAc * Xact_AcSem_inv2 * Xact_AcSemAc');
end
if ParamModel.ModelChoice(5) && SWITCH.DF
    % AcSemSem
    fprintf('AcSem model Sem Offset: Calculate exact number of parameters\n')
    L2_AcSemSem = Lambda_BestModel_AcSemSem(bb)*(1-Alpha)/2;
    %L1_AcSemSem(bb,aa)*Alpha;
    Active_Predictors_AcSemSem = find(B_AcSemSem(:,Dev_min_index_AcSemSem));
    X_AcSemSem = [X_voc_Train x_Train];
    Xact_AcSemSem = X_AcSemSem(:,Active_Predictors_AcSemSem);
    Xact_AcSemSem_inv = (Xact_AcSemSem' * Xact_AcSemSem + L2_AcSemSem .* eye(length(Active_Predictors_AcSemSem)))^(-1);
    DF_BestModel_AcSemSem(bb) = trace(Xact_AcSemSem * Xact_AcSemSem_inv * Xact_AcSemSem');
end


%% Out put some results of deviance difference between models
%             if SWITCH.FIG>1
%                 % Calculate the deviance difference for that bootstrap and
%                 % the AIC difference
%                 DiffDEV_AcSemAcoustic = Deviance_BestModel_AcSem(bb,aa)-Deviance_BestModel_Acoustic(bb,aa); %the more this value is negative the better the Acsem model is over the acoustic model
%                 DiffAIC_AcSemAcoustic = 2.*DF_BestModel_AcSem(bb,aa) - 2.*DF_BestModel_Acoustic(bb,aa) + DiffDEV_AcSemAcoustic;
%                 fprintf(1,'Deviance difference AcSem-Acoustic: %f (the more negative the better the AcSem model is\n',DiffDEV_AcSemAcoustic);
%                 fprintf(1,'AIC difference AcSem-Acoustic: %f (the more negative the better the AcSem model is\n',DiffAIC_AcSemAcoustic);
%                 DiffDEV_AcSemSemantic = Deviance_BestModel_AcSem(bb,aa)-Deviance_BestModel_Sem(bb,aa); %the more this value is negative the better the Acsem model is over the acoustic model
%                 DiffAIC_AcSemSemantic = 2.*DF_BestModel_AcSem(bb,aa) - 2.*DF_BestModel_Sem(bb,aa) + DiffDEV_AcSemSemantic;
%                 fprintf(1,'Deviance difference AcSem-Semantic: %f (the more negative the better the AcSem model is\n',DiffDEV_AcSemSemantic);
%                 fprintf(1,'AIC difference AcSem-Semantic: %f (the more negative the better the AcSem model is\n',DiffAIC_AcSemSemantic);
%             end
toc

%% Plot Lambda values, STRFs = model coefficients, predictions
if SWITCH.FIG>1
    if NbL_Ac>1
        NbMod=sum(ParamModel.ModelChoice([1 3:5]));
        Mod=find(ParamModel.ModelChoice);
        Mod(find(Mod==2))=[];
        if NbMod
            % Plot lambda values as a function of deviance for each
            % model (AcSem and Acoustic)
            figure(302)
            for ii=1:NbMod
                subplot(1,NbMod,ii);
                if Mod(ii)==1
                    plot(log10(FitInfo_Ac.Lambda), Deviance_M_Ac)
                    xlabel('log of Lambdas')
                    ylabel('Deviance of the Acoustic model')
                    title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                    vline(log10(FitInfo_Ac.Lambda(Dev_min_index_Ac)));
                elseif Mod(ii)==3
                    plot(log10(FitInfo_AcSem.Lambda), Deviance_M_AcSem)
                    xlabel('log of Lambdas')
                    ylabel('Deviance of the AcSem model')
                    title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                    vline(log10(FitInfo_AcSem.Lambda(Dev_min_index_AcSem)));
                elseif Mod(ii)==4
                    plot(log10(FitInfo_AcSemAc.Lambda), Deviance_M_AcSemAc)
                    xlabel('log of Lambdas')
                    ylabel('Deviance of the AcSem model Ac Offset')
                    title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                    vline(log10(FitInfo_AcSemAc.Lambda(Dev_min_index_AcSemAc)));
                elseif Mod(ii)==5
                    plot(log10(FitInfo_AcSemSem.Lambda), Deviance_M_AcSemSem)
                    xlabel('log of Lambdas')
                    ylabel('Deviance of the AcSem model Sem Offset')
                    title(sprintf('training size %d validating size %d, alpha %f',length(TrainSet), length(ValSet), Alpha))
                    vline(log10(FitInfo_AcSemSem.Lambda(Dev_min_index_AcSemSem)));
                end
            end
        end
    end
    
    pause(1)
    if SWITCH.Check==1
        fprintf(1, 'Calculate PC of spectro\n');
        [COEFF,SCORE,~,~]=princomp(x_Train,'econ');
        nPC=60;
        ds=dataset();
        for ii=1:nPC
            ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
        end
        ds.y=y_Train;
        mdl=LinearModel.fit(ds);
        PCSTRF=mdl.Coefficients.Estimate(2:end);
        STRF=COEFF(:,1:nPC)*PCSTRF;
        STRFM=reshape(STRF,Data.Df, Data.Dt);
        fprintf(1, 'Calculating STRF Ac using the %d first PC of the spectro\n\n\n\n', nPC);
        figure(303)
        imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, STRFM)
        axis xy
        title(sprintf('STRF Ac obtained with linear model with %d PC of the spectro', nPC));
        pause(1)
    end
    
    if ParamModel.ModelChoice(1)
        % Plot coefficients Acoustic model
        for ll=1:NbL_Ac
            if ll==Dev_min_index_Ac
                figure(304)
                if ParamModel.ZC==1
                    Lasso_STRF=reshape(B_Ac(:,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF=reshape(B_Ac(:,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF))), abs(min(min(Lasso_STRF))));
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF,[-valc valc])
                end
                axis xy
                title(sprintf('THIS IS THE ONE!!!\nSTRF Ac log of lambda=%f\n Deviance=%f FitIndex=%f\nFitIndexAR=%f',log10(FitInfo_Ac.Lambda(ll)),Deviance_BestModel_Acoustic(bb),FitIndex_Acoustic_local(bb),FitIndexARVal_Acoustic_local))
            else
                figure(305)
                if ParamModel.ZC==1
                    Lasso_STRF=reshape(B_Ac(:,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF=reshape(B_Ac(:,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF))), abs(min(min(Lasso_STRF))));
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF,[-valc valc])
                end
                axis xy
                title(sprintf('STRF Ac log of lambda=%f\n',log10(FitInfo_Ac.Lambda(ll))))
            end
            pause(1)
        end
    end
    
    if ParamModel.ModelChoice(3)
        % Plot Coefficients AcSem model
        for ll=1:NbL_AcSem
            if ll==Dev_min_index_AcSem
                figure(306)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSem=reshape(B_AcSem((size(Data.X_voc,2)+1):end,ll)'.*x_Train_std,Df,Dt);
                else
                    Lasso_STRF_AcSem=reshape(B_AcSem((size(Data.X_voc,2)+1):end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSem))), abs(min(min(Lasso_STRF_AcSem))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSem)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSem,[-valc valc])
                end
                axis xy
                title(sprintf('THIS IS THE ONE!!!\nSTRF AcSem log of lambda=%f\n Deviance=%f FitIndex=%f\nFitIndexAR=%f',log10(FitInfo_AcSem.Lambda(ll)),Deviance_BestModel_AcSem(bb),FitIndex_AcSem_local(bb),FitIndexARVal_AcSem_local))
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSem_plot = [FitInfo_AcSem.Intercept(ll) ; B_AcSem(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSem.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSem(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSem_plot(noModCat)=FitInfo_AcSem.Intercept(ll) + sum(B_AcSem(1:size(Data.X_voc,2),ll))/(length(B_AcSem(1:size(Data.X_voc,2),ll))-length(noModCat));
                end
                plot(B_AcSem_plot, 1:length(UVOC));
                title('Coefficients Categories AcSem Model');
                set(ss,'YTickLabel', UVOC);
            else
                figure(307)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSem=reshape(B_AcSem(9:end,ll)'.*x_Train_std,Df,Dt);
                else
                    Lasso_STRF_AcSem=reshape(B_AcSem(9:end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSem))), abs(min(min(Lasso_STRF_AcSem))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSem)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSem,[-valc valc])
                end
                axis xy
                title(sprintf('STRF AcSem log of lambda=%f\n',log10(FitInfo_AcSem.Lambda(ll))));
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSem_plot = [FitInfo_AcSem.Intercept(ll) ; B_AcSem(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSem.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSem(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSem_plot(noModCat)=FitInfo_AcSem.Intercept(ll) + sum(B_AcSem(1:size(Data.X_voc,2),ll))/(length(B_AcSem(1:size(Data.X_voc,2),ll)) -length(noModCat));
                end
                plot(B_AcSem_plot, 1:length(UVOC));
                title('Coefficients Categories AcSem Model');
                set(ss,'YTickLabel', UVOC);
            end
            pause(1)
        end
    end
    
    if ParamModel.ModelChoice(4)
        % Plot Coefficients AcSemAc model
        for ll=1:NbL_AcSemAc
            if ll==Dev_min_index_AcSemAc
                figure(308)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSemAc=reshape(B_AcSemAc((size(Data.X_voc,2)+1):end,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF_AcSemAc=reshape(B_AcSemAc((size(Data.X_voc,2)+1):end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSemAc))), abs(min(min(Lasso_STRF_AcSemAc))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemAc)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemAc,[-valc valc])
                end
                axis xy
                title(sprintf('THIS IS THE ONE!!!\nSTRF AcSem2 log of lambda=%f\n Deviance=%f FitIndex=%f\nFitIndexAR=%f',log10(FitInfo_AcSemAc.Lambda(ll)),Deviance_BestModel_AcSemAc(bb),FitIndex_AcSemAc_local(bb),FitIndexARVal_AcSemAc_local))
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSemAc_plot2 = [FitInfo_AcSemAc.Intercept(ll) ; B_AcSemAc(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSemAc.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSemAc(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSemAc_plot2(noModCat)=FitInfo_AcSemAc.Intercept(ll) + sum(B_AcSemAc(1:size(Data.X_voc,2),ll))/(length(B_AcSemAc(1:size(Data.X_voc,2),ll))-length(noModCat));
                end
                plot(B_AcSemAc_plot2, 1:length(UVOC));
                title('Coefficients Categories AcSemAc Model');
                set(ss,'YTickLabel', UVOC);
            else
                figure(309)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSemAc=reshape(B_AcSemAc(9:end,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF_AcSemAc=reshape(B_AcSemAc(9:end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSemAc))), abs(min(min(Lasso_STRF_AcSemAc))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemAc)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemAc,[-valc valc])
                end
                axis xy
                title(sprintf('STRF AcSemAc log of lambda=%f\n',log10(FitInfo_AcSemAc.Lambda(ll))));
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSemAc_plot2 = [FitInfo_AcSemAc.Intercept(ll) ; B_AcSemAc(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSemAc.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSemAc(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSemAc_plot2(noModCat)=FitInfo_AcSemAc.Intercept(ll) + sum(B_AcSemAc(1:size(Data.X_voc,2),ll))/(length(B_AcSemAc(1:size(Data.X_voc,2),ll)) -length(noModCat));
                end
                plot(B_AcSemAc_plot2, 1:length(UVOC));
                title('Coefficients Categories AcSemAc Model');
                set(ss,'YTickLabel', UVOC);
            end
            pause(1)
        end
    end
    
    if ParamModel.ModelChoice(5)
        % Plot Coefficients AcSemSem model
        for ll=1:NbL_AcSemSem
            if ll==Dev_min_index_AcSemSem
                figure(310)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSemSem=reshape(B_AcSemSem((size(Data.X_voc,2)+1):end,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF_AcSemSem=reshape(B_AcSemSem((size(Data.X_voc,2)+1):end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSemSem))), abs(min(min(Lasso_STRF_AcSemSem))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemSem)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemSem,[-valc valc])
                end
                axis xy
                title(sprintf('THIS IS THE ONE!!!\nSTRF AcSemSem log of lambda=%f\n Deviance=%f FitIndex=%f\nFitIndexAR=%f',log10(FitInfo_AcSemSem.Lambda(ll)),Deviance_BestModel_AcSemSem(bb),FitIndex_AcSemSem_local(bb),FitIndexARVal_AcSemSem_local))
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSemSem_plot2 = [FitInfo_AcSemSem.Intercept(ll) ; B_AcSemSem(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSemSem.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSemSem(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSemSem_plot2(noModCat)=FitInfo_AcSemSem.Intercept(ll) + sum(B_AcSemSem(1:size(Data.X_voc,2),ll))/(length(B_AcSemSem(1:size(Data.X_voc,2),ll))-length(noModCat));
                end
                plot(B_AcSemSem_plot2, 1:length(UVOC));
                title('Coefficients Categories AcSemSem Model');
                set(ss,'YTickLabel', UVOC);
            else
                figure(311)
                if ParamModel.ZC==1
                    Lasso_STRF_AcSemSem=reshape(B_AcSemSem(9:end,ll)'.*x_Train_std,Data.Df,Data.Dt);
                else
                    Lasso_STRF_AcSemSem=reshape(B_AcSemSem(9:end,ll),Data.Df,Data.Dt);
                end
                valc = max(abs(max(max(Lasso_STRF_AcSemSem))), abs(min(min(Lasso_STRF_AcSemSem))));
                subplot(1,2,1)
                if valc==0
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemSem)
                else
                    imagesc(TickSTRFspectro.to, TickSTRFspectro.fo, Lasso_STRF_AcSemSem,[-valc valc])
                end
                axis xy
                title(sprintf('STRF AcSemSem log of lambda=%f\n',log10(FitInfo_AcSemSem.Lambda(ll))));
                ss=subplot(1,2,2);
                % Make sure we have one stim per category in the
                % training set or set the coefficient value of the
                % model for that category as the average of
                % parameters for other categories instead of 0
                B_AcSemSem_plot2 = [FitInfo_AcSemSem.Intercept(ll) ; B_AcSemSem(1:size(Data.X_voc,2),ll)] + [0; repmat(FitInfo_AcSemSem.Intercept(ll),size(Data.X_voc,2),1)];
                noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
                if ~isempty(noModCat)
                    if ~sum(B_AcSemSem(noModCat,ll))==0
                        fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
                    end
                    B_AcSemSem_plot2(noModCat)=FitInfo_AcSemSem.Intercept(ll) + sum(B_AcSemSem(1:size(Data.X_voc,2),ll))/(length(B_AcSemSem(1:size(Data.X_voc,2),ll)) -length(noModCat));
                end
                plot(B_AcSemSem_plot2, 1:length(UVOC));
                title('Coefficients Categories AcSemSem Model');
                set(ss,'YTickLabel', UVOC);
            end
            pause(1)
        end
    end
    
    if ParamModel.ModelChoice(2)
        % Plot coefficient Semantic model
        figure(312)
        % Make sure we have one stim per category in the
        % training set or set the coefficient value of the
        % model for that category as the average of
        % parameters for other categories instead of 0
        B_Sem_plot = MDL_Sem.Coefficients.Estimate + [0; repmat(MDL_Sem.Coefficients.Estimate(1),length(MDL_Sem.Coefficients.Estimate)-1,1)];
        noModCat= find(sum(Data.X_voc(TrainSet,:),1)==0,1);
        if ~isempty(noModCat)
            if ~sum(MDL_Sem.Coefficients.Estimate(noModCat))==0
                fprintf('There is a problem here the coefficients of the model is expected to be zero when there is no stim in the training set for that category\n')
            end
            B_Sem_plot(noModCat)=MDL_Sem.Coefficients.Estimate(1) + sum(MDL_Sem.Coefficients.Estimate(2:end))/(length(MDL_Sem.Coefficients.Estimate)-length(noModCat));
        end
        plot(B_Sem_plot, 1:length(UVOC));
        title(sprintf('Coefficients Categories Sem Model Deviance=%f\nFitIndex=%f\nFitIndexAR=%f',Deviance_BestModel_Sem(bb),FitIndex_Sem_local(bb),FitIndexARVal_Sem_local));
        set(gca,'YTickLabel', UVOC);
        pause(1)
    end
    
    
    % Plot prediction for each model
    NbMod=sum(ParamModel.ModelChoice);
    Mod=find(ParamModel.ModelChoice);
    LocalModelPredict=cell(NbMod,1);
    Legend= cell(NbMod,1);
    MAXY_local=nan(NbMod,1);
    for ii =1:NbMod
        if Mod(ii)==1
            LocalModelPredict{ii}=Ypredict_Acoustic{bb};
            Legend{ii}='Acoustic';
        elseif Mod(ii)==2
            LocalModelPredict{ii}=Ypredict_Semantic{bb};
            Legend{ii}='Semantic';
        elseif Mod(ii)==3
            LocalModelPredict{ii}=Ypredict_AcSem{bb};
            Legend{ii}='AcSem';
        elseif Mod(ii)==4
            LocalModelPredict{ii}=Ypredict_AcSemAc{bb};
            Legend{ii}='AcSem Ac Offset';
        elseif Mod(ii)==5
            LocalModelPredict{ii}=Ypredict_AcSemSem{bb};
            Legend{ii}='AcSem Sem Offset';
        end
        MAXY_local(ii)=max(LocalModelPredict{ii});
    end
    
    MAXY=max(MAXY_local);
    MAXX=max(y_Val);
    MAX=max(MAXX,MAXY);
    
    for jj=1:length(LocalModelPredict)
        figure(313)
        subplot(length(LocalModelPredict),1,jj);
        gscatter(y_Val, LocalModelPredict{jj}, XVocVal{bb}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
        title(sprintf('%s', Legend{jj}));
        if strcmp(ParamModel.NeuroRes,'mean')
            xlabel('observed spike rate (spike per ms)')
            ylabel('predicted spike rate (spike per ms)')
        elseif strcmp(ParamModel.NeuroRes, 'count')
            xlabel('observed spike count (bin=20 ms)')
            ylabel('predicted spike count (bin=20 ms)')
        end
        axis([0 MAX+MAX/4 0 MAX]);
        hold on
        plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
        hold off
    end
    pause(1)
end

    
if SWITCH.FIG>1 && (length(Lambda_All_local_Ac{1})>1 || length(Lambda_All_local_AcSem{1})>1 || length(Lambda_All_local_AcSemAc{1})>1 || length(Lambda_All_local_AcSemSem{1})>1)
    ModelNames={'Acoustic','Semantic','AcSem', 'AcSemAc','AcSemSem'};
    fig314=figure(314);
    close fig314
    figure(314)
    if ParamModel.ModelChoice(1)==1
        subplot(1,3,1)
        plot(1:ParamModel.BootstrapSTRF,cumsum(LL_BestModel_Acoustic)'./(1:ParamModel.BootstrapSTRF),'b-');
        subplot(1,3,2)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(cumsum(Lambda_BestModel_Acoustic)'./(1:ParamModel.BootstrapSTRF)),'b-');
        Cummedian_Lambda_Acoustic=nan(ParamModel.BootstrapSTRF,1);
        for bb=1:ParamModel.BootstrapSTRF
            Cummedian_Lambda_Acoustic(bb)=median(Lambda_BestModel_Acoustic(1:bb));
        end
        subplot(1,3,3)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(Cummedian_Lambda_Acoustic),'b-');
    end
    if ParamModel.ModelChoice(3)==1
        subplot(1,3,1)
        hold on
        plot(1:ParamModel.BootstrapSTRF,cumsum(LL_BestModel_AcSem)'./(1:ParamModel.BootstrapSTRF),'m-');
        subplot(1,3,2)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(cumsum(Lambda_BestModel_AcSem)'./(1:ParamModel.BootstrapSTRF)),'m-');
        Cummedian_Lambda_AcSem=nan(ParamModel.BootstrapSTRF,1);
        for bb=1:ParamModel.BootstrapSTRF
            Cummedian_Lambda_AcSem(bb)=median(Lambda_BestModel_AcSem(1:bb));
        end
        subplot(1,3,3)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(Cummedian_Lambda_AcSem),'m-');
    end
    if ParamModel.ModelChoice(4)==1
        subplot(1,3,1)
        hold on
        plot(1:ParamModel.BootstrapSTRF,cumsum(LL_BestModel_AcSemAc)'./(1:ParamModel.BootstrapSTRF),'.-','Color',[0.5 0 1]);
        subplot(1,3,2)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(cumsum(Lambda_BestModel_AcSemAc)'./(1:ParamModel.BootstrapSTRF)),'.-','Color',[0.5 0 1]);
        Cummedian_Lambda_AcSemAc=nan(ParamModel.BootstrapSTRF,1);
        for bb=1:ParamModel.BootstrapSTRF
            Cummedian_Lambda_AcSemAc(bb)=median(Lambda_BestModel_AcSemAc(1:bb));
        end
        subplot(1,3,3)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(Cummedian_Lambda_AcSemAc),'.-','Color',[0.5 0 1]);
    end
    if ParamModel.ModelChoice(5)==1
        subplot(1,3,1)
        hold on
        plot(1:ParamModel.BootstrapSTRF,cumsum(LL_BestModel_AcSemSem)'./(1:ParamModel.BootstrapSTRF),'r-');
        subplot(1,3,2)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(cumsum(Lambda_BestModel_AcSemSem)'./(1:ParamModel.BootstrapSTRF)),'r-');
        Cummedian_Lambda_AcSemSem=nan(ParamModel.BootstrapSTRF,1);
        for bb=1:ParamModel.BootstrapSTRF
            Cummedian_Lambda_AcSemSem(bb)=median(Lambda_BestModel_AcSemSem(1:bb));
        end
        subplot(1,3,3)
        hold on
        plot(1:ParamModel.BootstrapSTRF,log10(Cummedian_Lambda_AcSemSem),'r-');
    end      
    Modelschoice=find(ParamModel.ModelChoice);
    Modelschoice(Modelschoice==2)=[];
    subplot(1,3,1)
    ylabel('Average LL')
    xlabel('Booststrap')
    title('Stabilization of LL mean over bootstrap')
    legend(ModelNames{Modelschoice});
    subplot(1,3,2)
    ylabel('Mean Lambda Model (log10)')
    xlabel('Booststrap')
    title('Stabilization of lambda mean over bootstrap')
    legend(ModelNames{Modelschoice});
    subplot(1,3,3)
    ylabel('Median Lambda Model (log10)')
    xlabel('Booststrap')
    title('Stabilization of lambda median over bootstrap')
    legend(ModelNames{Modelschoice});
    hold off
    pause(2)
end

%% Save output data
if ParamModel.ModelChoice(1)
    MResults.Deviance_All_local_Ac=Deviance_All_local_Ac;
    MResults.Lambda_All_local_Ac = Lambda_All_local_Ac;
    MResults.LL_All_local_Ac=LL_All_local_Ac;
    MResults.Ypredict_Acoustic = Ypredict_Acoustic;
    MResults.Deviance_BestModel_Acoustic= Deviance_BestModel_Acoustic;
    MResults.LL_BestModel_Acoustic =LL_BestModel_Acoustic;
    MResults.FitIndex_Acoustic_local=FitIndex_Acoustic_local;
    MResults.Lambda_BestModel_Acoustic=Lambda_BestModel_Acoustic;
    MResults.ModelB_Ac = ModelB_Ac;% Parameters of the best model at that bootstrap
    MResults.ModelB0_Ac=ModelB0_Ac;
    MResults.Speed_Acoustic=Speed_Acoustic;
   if SWITCH.DF
        MResults.DF_BestModel_Acoustic = DF_BestModel_Acoustic;
    end
end
if ParamModel.ModelChoice(2)
    MResults.Deviance_BestModel_Sem = Deviance_BestModel_Sem;
    MResults.LL_BestModel_Sem = LL_BestModel_Sem;
    MResults.Ypredict_Semantic = Ypredict_Semantic;
    MResults.ModelB_Sem = ModelB_Sem;
    MResults.FitIndex_Sem_local=FitIndex_Sem_local;
    if SWITCH.DF
        MResults.DF_BestModel_Sem = DF_BestModel_Sem;
    end
end
if ParamModel.ModelChoice(3)
    MResults.Deviance_All_local_AcSem=Deviance_All_local_AcSem;
    MResults.Lambda_All_local_AcSem = Lambda_All_local_AcSem;
    MResults.LL_All_local_AcSem=LL_All_local_AcSem;
    MResults.Ypredict_AcSem = Ypredict_AcSem;
    MResults.Deviance_BestModel_AcSem = Deviance_BestModel_AcSem;
    MResults.LL_BestModel_AcSem = LL_BestModel_AcSem;
    MResults.FitIndex_AcSem_local = FitIndex_AcSem_local;
    MResults.Lambda_BestModel_AcSem = Lambda_BestModel_AcSem;
    MResults.ModelBspectro_AcSem = ModelBspectro_AcSem;% Parameters of the best model at that bootstrap
    MResults.ModelBsem_AcSem = ModelBsem_AcSem;
    MResults.ModelB0_AcSem = ModelB0_AcSem;
    if SWITCH.DF
        MResults.DF_BestModel_AcSem = DF_BestModel_AcSem;
    end
end
if ParamModel.ModelChoice(4)
    MResults.Deviance_All_local_AcSemAc=Deviance_All_local_AcSemAc;
    MResults.Lambda_All_local_AcSemAc = Lambda_All_local_AcSemAc;
    MResults.LL_All_local_AcSemAc=LL_All_local_AcSemAc;
    MResults.Ypredict_AcSemAc = Ypredict_AcSemAc;
    MResults.Deviance_BestModel_AcSemAc= Deviance_BestModel_AcSemAc;
    MResults.LL_BestModel_AcSemAc=LL_BestModel_AcSemAc;
    MResults.Lambda_BestModel_AcSemAc=Lambda_BestModel_AcSemAc;
    MResults.FitIndex_AcSemAc_local=FitIndex_AcSemAc_local;
    MResults.ModelBspectro_AcSemAc = ModelBspectro_AcSemAc;% Parameters of the best model at that bootstrap
    MResults.ModelBsem_AcSemAc=ModelBsem_AcSemAc;
    MResults.ModelB0_AcSemAc=ModelB0_AcSemAc;
    if SWITCH.DF
        MResults.DF_BestModel_AcSemAc = DF_BestModel_AcSemAc;
    end
end
if ParamModel.ModelChoice(5)
    MResults.Deviance_All_local_AcSemSem=Deviance_All_local_AcSemSem;
    MResults.Lambda_All_local_AcSemSem = Lambda_All_local_AcSemSem;
    MResults.Ypredict_AcSemSem = Ypredict_AcSemSem;
    MResults.LL_All_local_AcSemSem=LL_All_local_AcSemSem;
    MResults.Deviance_BestModel_AcSemSem= Deviance_BestModel_AcSemSem;
    MResults.LL_BestModel_AcSemSem=LL_BestModel_AcSemSem;
    MResults.FitIndex_AcSemSem_local=FitIndex_AcSemSem_local;
    MResults.Lambda_BestModel_AcSemSem = Lambda_BestModel_AcSemSem;
    MResults.ModelBspectro_AcSemSem = ModelBspectro_AcSemSem;% Parameters of the best model at that bootstrap
    MResults.ModelBsem_AcSemSem = ModelBsem_AcSemSem;
    MResults.ModelB0_AcSemSem = ModelB0_AcSemSem;
    if SWITCH.DF
        MResults.DF_BestModel_AcSemSem = DF_BestModel_AcSemSem;
    end
end

if SWITCH.AR
    MResults.LL_ARVal = LL_ARVal;
    %MResults.LL_AR = LL_AR;
    MResults.Ypredict_ARVal = Ypredict_ARVal;
end
MResults.LL_Saturated = LL_Saturated;
MResults.LL_BestGuess = LL_BestGuess;
MResults.LL_WorseGuess = LL_WorseGuess;
MResults.YobsVal = YobsVal;
MResults.ValSetFirstRep = ValSetFirstRep;
MResults.NumStimTrials = NumStimTrials;
MResults.YCeilModelVal = YCeilModelVal;
MResults.Ypredict_Floor = Ypredict_Floor;
MResults.XVocVal = XVocVal;
MResults.ValProp = ValProp;
MResults.ValSize = ValSize;
MResults.ValSets = ValSets;
MResults.TrainSets = TrainSets;
if any(ParamModel.ModelChoice([1,3:5]))
    MResults.XSpecVal = XSpecVal;
    if ParamModel.ZC==1
        MResults.X_meanZC = X_meanZC;
        MResults.X_stdZC = X_stdZC;
    end
end

