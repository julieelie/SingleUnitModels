function [LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod(Spectro, VocType, PSTH, Trials, Emitter, MinWin, MaxWin, Increment, ResDelay, NeuroRes, DISTR, LINK, calfilename,DoneCalc)
FIG=0; % set to 1 for some debugging figures 2 for all debugging figures 3 for extra figures
SWITCH.AllAlpha=0; % Set to 1 if ou want to calculate cumulative information for all alphas
SWITCH.FIG=0; % figure set up for runLassoGlmModels.m set to 1 for some debugging figures 2 for all debugging figures 3 for extra figures
SWITCH.Check=0;%set to 1 to compare ridge results with Linear model results
SWITCH.DF=0; % Set to 1 to have the exact degrees of freedom of the models being calculated
SWITCH.AR=0; %Set to 1 to run Auto-regressive model
ParamModel.ZC=0;%Set to 1 if you want to z-score bin per bin the values of the spectrograms using the training set just before performing lasso glm
ParamModel.MeanSubstractSpec=0; % Set to 1 if you want to substract the average spectro bin per bin on all the dataset
ParamModel.LAMBDARATIO=1e-4;
ParamModel.NUMLAMBDA=10;%25?
ParamModel.Cum_Info = 1; % set to 1 to calculate cumulative information
ParamModel.NumSamples_MC_Cum_Info = [10^6 10^4 10^3 10^2]; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers

if ParamModel.Cum_Info
    getenv('HOSTNAME')
    if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))%savio Cluster
        FolderTempInfStorage = '/tmp/LocalTempStorageInfo';
    elseif ismac()
        FolderTempInfStorage = '/tmp/LocalTempStorageInfo';
    else %we are on strfinator or a cluster machine
        FolderTempInfStorage = '/auto/tdrive/julie/LocalTempStorageInfo';
    end
end

NbBootstrap_Lambda = 20;%20
NbBootstrap_Deviance = 10;%10

% Determine a list of alpha (parameter that range the regularization
% betweeen ridge (L2, alpha=0) and lasso (L1, alpha =1))
Alphas=0.001; % STRFs are easier to interpret using ridge than using Lasso and deviances are similar.

ParamModel.ModelChoice=[0 1 0 0 0];% This logic vector indicates which models should
...be run (Acoustic, Semantic, AcSem without Offset,AcSem Semantic Offset,
    ...AcSem Accoustic Offset) Note that the last two models need the
    ...calculation of the first two models
    
ParamModel.LINK_AR='identity';

if nargin<14
    PrevData = false;
else
    PrevData = true;
end

if nargin<13
    saveonline = false;
else
    saveonline = true;
end

if nargin<12
    ParamModel.LINK='log'; %'identity'
else
    ParamModel.LINK=LINK;
end
if nargin<11
    ParamModel.DISTR='poisson';%'normal'
else
    ParamModel.DISTR=DISTR;
end

if nargin<10
    ParamModel.NeuroRes = 'count';
else
    ParamModel.NeuroRes =NeuroRes;
end
if nargin<9
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<8
    Increment = 20; %increase the size of the spectro window with a 5ms pace
end
if nargin<7
    MaxWin = 200; %maximum values the window of analysis can reach
end
if nargin<6
    MinWin = 20; %minimum size of the window of analysis from the begining and also size of analysis of spike rate/count (20 for STRF, 10 for Infocalulations)
end
Flow = 8000;%spectrograms are low passed at 8Khz for the calculations

%define the increasing size of the window of the spectrogram
if PrevData
    Wins = DoneCalc.Wins;
else
    Wins = MinWin:Increment:MaxWin;
end

% # of models to run on the data
modNum = length(Wins);

% Number of stims in the data set
NbStim = length(VocType);

%% Define the average spectro
if ParamModel.MeanSubstractSpec && any(ParamModel.ModelChoice([1,3:5]))
    % Find a stim larger than the max window of the dataset
    Duration = nan(NbStim,1);
    for st=1:NbStim
        Duration(st)=Spectro.to{st}(end)*1000;
    end
    LongStim = find(Duration >= MaxWin);
    LongStim = LongStim(1);
    NTime = sum(Spectro.to{LongStim}*1000 <= MaxWin);
    NFreq = sum(Spectro.fo{LongStim} <= Flow);
    SpecMat = nan(NFreq, NTime, NbStim);
    %SpecBin = zeros(NFreq, NTime);
    for st=1:NbStim
        NTime_local = sum(Spectro.to{st}*1000 <= MaxWin);
        NFreq_local = sum(Spectro.fo{st} <= Flow);
        Spec_local = reshape(Spectro.spec{st}, length(Spectro.fo{st}), length(Spectro.to{st}));
        SpecMat(1:NFreq_local, 1:NTime_local,st) = Spec_local(1:NFreq_local, 1:NTime_local);
        %SpecBin(1:NFreq_local, 1:NTime_local) = SpecBin(1:NFreq_local, 1:NTime_local) + ones(NFreq_local, NTime_local);
    end
    LogSpecMat = 20*log10(abs(SpecMat));
    % Ground the values to supress -Inf max 100 dB of difference is
    % autorized between highest and lowest values of power
    MAXI = max(max(max(LogSpecMat)));
    LogSpecMatcorrected = LogSpecMat;
    LogSpecMatcorrected(LogSpecMatcorrected< (MAXI-80)) = MAXI-80;
    AvSpec = nanmean(LogSpecMatcorrected,3);
    if FIG
        figure()
        subplot(1,4,1)
        imagesc(nanmean(abs(SpecMat),3))
        axis xy
        colorbar()
        subplot(1,4,2)
        imagesc(nanmean(LogSpecMat,3))
        axis xy
        colorbar()
        subplot(1,4,3)
        imagesc(nanmean(LogSpecMatcorrected,3))
        axis xy
        colorbar()
        subplot(1,4,4)
        imagesc(Spectro.to{LongStim}(1:NTime)*1000,Spectro.fo{LongStim}(1:NFreq), AvSpec)
        xlabel('Time (ms)')
        ylabel('Frequencies')
        title('Average spectro')
        axis xy
        colorbar()
    end
end

%% Initialize a bunch of output variables
if PrevData
    % Load Previous data
    Deviance = DoneCalc.Deviance;
    LL = DoneCalc.LL;
    LambdaChoice = DoneCalc.LambdaChoice;
    Model = DoneCalc.Model;
    PropVal = DoneCalc.PropVal;
    Data = DoneCalc.Data;
    clear DoneCalc
    
    % Determine where the calculus stopped
    UnrunWindows =[];
    for ww=1:modNum
        if isempty(Data.MeanSpectroStim{ww})
            UnrunWindows = [UnrunWindows ww];
        end
    end
    StartWin = min(UnrunWindows);
    
    % Check that there is no mess in the data (weird stuffs for some
    % files...)
    try
        PreviousWin.Stim_local_old=Data.x_stim_indices_wholeset{StartWin - 1};
    catch
        fprintf('Previous Data format is not coherent, prefer to discard and start from first window\n')
        PrevData = 0;
    end
end

if PrevData
    % Initialize datasets to be used as passed information by the
    % auto-regressive model\
    NbStim_local_old = length(PreviousWin.Stim_local_old);
    if strcmp(ParamModel.NeuroRes, 'count')
        PreviousWin.y_old=cell(NbStim_local_old,1);
    else
        PreviousWin.y_old = nan(NbStim_local_old,1);%this matrix will contain the average spike rate in spikes/ms at that precise position and for all the stims choosen for that run
    end
    Win_old = Wins(StartWin - 1);
    
    % Getting spike trains of previous window
    for ss = 1:NbStim_local_old
        dd=PreviousWin.Stim_local_old(ss);
        % Values of max spike rate(y), mean spike rate (y) and exact number of spike per trial (y) within the window
        % Check the distribution of responses (Gaussian or Poisson) for each stim
        if strcmp(ParamModel.NeuroRes, 'max')
            PreviousWin.y_old(ss) = max(PSTH{dd}((Win_old-MinWin+ResDelay):(Win_old+ResDelay)));
        elseif strcmp(ParamModel.NeuroRes, 'mean')
            PreviousWin.y_old(ss) = mean(PSTH{dd}((Win_old-MinWin+ResDelay):(Win_old+ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
        elseif strcmp(ParamModel.NeuroRes, 'count')
            PreviousWin.y_old{ss}=nan(length(Trials{dd}),1);
            for tt=1:length(Trials{dd})
                PreviousWin.y_old{ss}(tt)=sum((Trials{dd}{tt}>(Win_old-MinWin+ResDelay)).*(Trials{dd}{tt}<(Win_old+ResDelay)));
            end
 % REVISE FRON HERE!!       elseif strcmp(ParamModel.NeuroRes, 'count_gaussfiltered')
            PreviousWin.y_old{ss}=nan(length(Trials{dd}),1);
            for tt=1:length(Trials{dd})
                PreviousWin.y_old{ss}(tt)=sum((Trials{dd}{tt}>(Win_old-MinWin+ResDelay)).*(Trials{dd}{tt}<(Win_old+ResDelay)));
            end
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', ParamModel.NeuroRes);
            
        end
    end
else
    cum_info.ExactHist5 = nan(modNum,length(Alphas));
    cum_info.EstMonteCarlo = nan(modNum,length(Alphas));
    cum_info.EstMonteCarlo2 = nan(modNum,length(Alphas));
    cum_info.EstMonteCarlo3 = nan(modNum,length(Alphas));
    cum_info.EstMonteCarlo4 = nan(modNum,length(Alphas));
    cum_info.EstMarkov2 = nan(modNum,length(Alphas));
    cum_info.EstMarkov3 = nan(modNum,length(Alphas));
    cum_info.EstMarkov4 = nan(modNum,length(Alphas));
    cum_info.EstMarkov5 = nan(modNum,length(Alphas));
    if ParamModel.ModelChoice(1)
        Deviance.Acoustic.DF = cell(modNum,length(Alphas));
        Deviance.Acoustic.values = cell(modNum,length(Alphas));
        Deviance.Acoustic.lambda = cell(modNum,length(Alphas));
        Deviance.Acoustic.FitIndex=nan(modNum,length(Alphas));
        Deviance.Acoustic.FitIndexAR=nan(modNum,length(Alphas));
        Deviance.Acoustic.FitIndexARVal=nan(modNum,length(Alphas));
        Deviance.Acoustic.ypredictVal = cell(modNum,length(Alphas));
        Deviance.Acoustic.Speed = cell(modNum,length(Alphas));
        LL.Acoustic.values=cell(modNum,length(Alphas));
        LambdaChoice.Acoustic.Lambdas = cell(modNum,length(Alphas));
        LambdaChoice.Acoustic.Speed = cell(modNum,length(Alphas));
        Model.Acoustic.Lambdas = nan(modNum,length(Alphas));
        Model.Acoustic.B = cell(modNum,length(Alphas));
        Model.Acoustic.B0 = cell(modNum,length(Alphas));
        Model.Acoustic.Speed = cell(modNum,length(Alphas));
        Model.Acoustic.y_predict = cell(modNum,length(Alphas));
        Model.Acoustic.info = nan(modNum,length(Alphas));
        if ParamModel.Cum_Info
            Model.Acoustic.cum_info=cum_info;
        end
        Model.Acoustic.P_YgivenS_all1 = cell(modNum,length(Alphas));
        Model.Acoustic.P_YgivenS_all2 = cell(modNum,length(Alphas));
    end
    if ParamModel.ModelChoice(2)
        Deviance.Semantic.values = cell(modNum,length(Alphas));
        Deviance.Semantic.DF = cell(modNum,length(Alphas));
        Deviance.Semantic.FitIndex=nan(modNum,length(Alphas));
        Deviance.Semantic.FitIndexAR=nan(modNum,length(Alphas));
        Deviance.Semantic.FitIndexARVal=nan(modNum,length(Alphas));
        Deviance.Semantic.ypredictVal = cell(modNum,length(Alphas));
        LL.Semantic.values=cell(modNum,length(Alphas));
        Model.Semantic.B = cell(modNum,length(Alphas));
        Model.Semantic.B0 = cell(modNum,length(Alphas));
        Model.Semantic.y_predict = cell(modNum,length(Alphas));
        Model.Semantic.info = nan(modNum,length(Alphas));
%         Model.Semantic.info_enddataset = nan(modNum,length(Alphas));
        if ParamModel.Cum_Info
            Model.Semantic.cum_info.current_dataset=cum_info;
%             Model.Semantic.cum_info.end_dataset=cum_info;
        end
        Model.Semantic.P_YgivenS_all1 = cell(modNum,length(Alphas));
        Model.Semantic.P_YgivenS_all2 = cell(modNum,length(Alphas));
        Model.Semantic.P_YgivenS_all1_enddataset = cell(modNum,length(Alphas));
        Model.Semantic.P_YgivenS_all2_enddataset = cell(modNum,length(Alphas));
    end
    if ParamModel.ModelChoice(3)
        Deviance.AcSem.values = cell(modNum,length(Alphas));
        Deviance.AcSem.DF = cell(modNum,length(Alphas));
        Deviance.AcSem.lambda = cell(modNum,length(Alphas));
        Deviance.AcSem.FitIndex=nan(modNum,length(Alphas));
        Deviance.AcSem.FitIndexAR=nan(modNum,length(Alphas));
        Deviance.AcSem.FitIndexARVal=nan(modNum,length(Alphas));
        Deviance.AcSem.ypredictVal = cell(modNum,length(Alphas));
        LL.AcSem.values=cell(modNum,length(Alphas));
        LambdaChoice.AcSem.Lambdas = cell(modNum,length(Alphas));
        Model.AcSem.Lambdas = nan(modNum,length(Alphas));
        Model.AcSem.Bspectro = cell(modNum,length(Alphas));
        Model.AcSem.Bsem = cell(modNum,length(Alphas));
        Model.AcSem.B0 = cell(modNum,length(Alphas));
        Model.AcSem.y_predict = cell(modNum,length(Alphas));
        Model.AcSem.info = nan(modNum,length(Alphas));
    end
    
    if ParamModel.ModelChoice(4) && ParamModel.ModelChoice(1)
        Deviance.AcSemAc.values = cell(modNum,length(Alphas));
        Deviance.AcSemAc.DF = cell(modNum,length(Alphas));
        Deviance.AcSemAc.lambda = cell(modNum,length(Alphas));
        Deviance.AcSemAc.FitIndex=nan(modNum,length(Alphas));
        Deviance.AcSemAc.FitIndexAR=nan(modNum,length(Alphas));
        Deviance.AcSemAc.FitIndexARVal=nan(modNum,length(Alphas));
        Deviance.AcSemAc.ypredictVal = cell(modNum,length(Alphas));
        LL.AcSemAc.values=cell(modNum,length(Alphas));
        LambdaChoice.AcSemAc.Lambdas = cell(modNum,length(Alphas));
        Model.AcSemAc.Lambdas = nan(modNum,length(Alphas));
        Model.AcSemAc.Bspectro = cell(modNum,length(Alphas));
        Model.AcSemAc.Bsem = cell(modNum,length(Alphas));
        Model.AcSemAc.B0 = cell(modNum,length(Alphas));
        Model.AcSemAc.y_predict = cell(modNum,length(Alphas));
        Model.AcSemAc.info = nan(modNum,length(Alphas));
        if ParamModel.Cum_Info
            Model.AcSemAc.cum_info = cum_info;
        end
        Model.AcSemAc.P_YgivenS_all1 = cell(modNum,length(Alphas));
        Model.AcSemAc.P_YgivenS_all2 = cell(modNum,length(Alphas));
    elseif ParamModel.ModelChoice(4) && ~ParamModel.ModelChoice(1)
        fprintf('WARNING: No way to calculate the AcSem model with Acoustic Offset if you do not calculate Acoustic model')
    end
    
    if ParamModel.ModelChoice(5) && ParamModel.ModelChoice(2)
        Deviance.AcSemSem.values = cell(modNum,length(Alphas));
        Deviance.AcSemSem.DF = cell(modNum,length(Alphas));
        Deviance.AcSemSem.lambda = cell(modNum,length(Alphas));
        Deviance.AcSemSem.FitIndex=nan(modNum,length(Alphas));
        Deviance.AcSemSem.FitIndexAR=nan(modNum,length(Alphas));
        Deviance.AcSemSem.FitIndexARVal=nan(modNum,length(Alphas));
        Deviance.AcSemSem.ypredictVal = cell(modNum,length(Alphas));
        LL.AcSemSem.values=cell(modNum,length(Alphas));
        LambdaChoice.AcSemSem.Lambdas = cell(modNum,length(Alphas));
        Model.AcSemSem.Lambdas = nan(modNum,length(Alphas));
        Model.AcSemSem.Bspectro = cell(modNum,length(Alphas));
        Model.AcSemSem.Bsem = cell(modNum,length(Alphas));
        Model.AcSemSem.B0 = cell(modNum,length(Alphas));
        Model.AcSemSem.y_predict = cell(modNum,length(Alphas));
        Model.AcSemSem.info = nan(modNum,length(Alphas));
        if ParamModel.Cum_Info
            Model.AcSemSem.cum_info=cum_info;
        end
        Model.AcSemSem.P_YgivenS_all1 = cell(modNum,length(Alphas));
        Model.AcSemSem.P_YgivenS_all2 = cell(modNum,length(Alphas));
    elseif ParamModel.ModelChoice(5) && ~ParamModel.ModelChoice(2)
        fprintf('WARNING: No way to calculate the AcSem model with Semantic Offset if you do not calculate semantic model')
    end
    
    LambdaChoice.ValSets=cell(modNum,length(Alphas));
    LambdaChoice.TrainSets=cell(modNum,length(Alphas));
    Deviance.ValSets=cell(modNum,length(Alphas));
    Deviance.TrainSets=cell(modNum,length(Alphas));
    Deviance.yobsCat = cell(modNum,length(Alphas));
    Deviance.yobsVal = cell(modNum,length(Alphas));
    Deviance.YCeilModelVal= cell(modNum,length(Alphas));
    Deviance.Ypredict_ARVal = cell(modNum,length(Alphas));
    Deviance.Ypredict_Floor = cell(modNum,length(Alphas));
    LL.Ceiling.values=cell(modNum,length(Alphas));
    LL.Floor.values=cell(modNum,length(Alphas));
    %LL.AutoRegressiveTV.values=cell(modNum,length(Alphas));
    if SWITCH.AR
        LL.AutoRegressiveVal.values=cell(modNum,length(Alphas));
    end
    LL.Saturated.values = cell(modNum,length(Alphas));
    
    
    PropVal.mean = nan(modNum,length(Alphas));
    PropVal.std = nan(modNum,length(Alphas));
    PropVal.values = cell(modNum,length(Alphas));
    
    Model.TickSTRFspectro.to = cell(modNum,1);
    Model.TickSTRFspectro.fo = cell(modNum,1);
    
    
    Data.y_wholeset = cell(modNum,1);
    Data.y_wholeset_AvObsDataperStim = cell(modNum,1);
    Data.x_wholeset = cell(modNum,1);
    Data.x_std_wholeset = cell(modNum,1);
    Data.x_mean_wholeset = cell(modNum,1);
    Data.X_voc_wholeset = cell(modNum,1);
    Data.x_stim_indices_wholeset = cell(modNum,1);
    Data.x_stim_repetition = cell(modNum,1);
    Data.stim_entropy = nan(modNum,1);
    Data.semanticcategory_entropy = nan(modNum,1);
    Data.y_wholeset_bestGuess = cell(modNum,1);
    Data.y_wholeset_Floor = cell(modNum,1);
    Data.y_wholeset_AR = cell(modNum,1);
    Data.InputData = cell(modNum,1);
    Data.MeanSpectroStim = cell(modNum,1);
    Data.stim_entropy = nan(modNum,1);
    Model.Floor.info = nan(modNum,1);
    Model.Floor.cum_info = cum_info;
    Model.Floor.y_predict = cell(modNum,1);
    Model.Floor.P_YgivenS_all1 = cell(modNum,1);
    Model.Floor.P_YgivenS_all2 = cell(modNum,1);
    Model.Ceiling.info = nan(modNum,1);
%     Model.Ceiling.info_enddataset = nan(modNum,1);
    Model.Ceiling.cum_info.current_dataset=cum_info;
%     Model.Ceiling.cum_info.end_dataset=cum_info;
    Model.Ceiling.P_YgivenS_all1_enddataset = cell(modNum,1);
    Model.Ceiling.P_YgivenS_all2_enddataset = cell(modNum,1);
    Model.Ceiling.P_YgivenS_all1 = cell(modNum,1);
    Model.Ceiling.P_YgivenS_all2 = cell(modNum,1);
    if SWITCH.AR
        Model.AR.info = nan(modNum,1);
        Model.AR.cum_info = cum_info;
        Model.AR.P_YgivenS_all1 = cell(modNum,1);
        Model.AR.P_YgivenS_all2 = cell(modNum,1);
    end
    
    Data.VOC = cell(modNum,1);
    
    % Initialize datasets to be used as passed information by the
    % auto-regressive model
    PreviousWin.y_old=[];
    PreviousWin.Stim_local_old=[];
    StartWin=1;
end

%% Now loop through window sizes and calculate models
for mm = StartWin:modNum
    fprintf(1,'********************************%d/%d models********************\n', mm, modNum);
    WindowTime = tic;
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
    Stim_local_end = find(duration >= (Wins(end)+ResDelay));% here are the stims that are used until the last window
    Data.stim_entropy(mm) = log2(NbStim_local);
    if NbStim_local<20
        sprintf('Only %d stims long enough to run the model: no model is run with window size %dms\n', NbStim_local, Win);
        break
    end
    %NBPC = [1:9 10:5:(NbStim_local*0.8)]; % This is not a good solution since the R2A profile of most cells show stairs-like structure with increasing number of PC we need to apply a ridge regression after the PCA.
    if any(ParamModel.ModelChoice([1,3:5]))
        Dt = sum((1000.*Spectro.to{Stim_local(1)})<= Win);
        Df = sum(Spectro.fo{Stim_local(1)}<= Flow);
        InputData.Dt = Dt;
        InputData.Df = Df;
        x= nan(NbStim_local,Df*Dt);%this matrix will contain the vectors of spectrograms for all the stims for that window size

        % Store the ticks of the spectro and STRF
        LongestStim = find(duration==max(duration));
        Fo_Indices=find(Spectro.fo{LongestStim}<=Flow);
        To_Indices=find((1000.*Spectro.to{LongestStim})<=Win);
        Model.TickSTRFspectro.to{mm} = Spectro.to{LongestStim}(To_Indices);
        Model.TickSTRFspectro.fo{mm} = Spectro.fo{LongestStim}(Fo_Indices);
        TickSTRFspectro.to=Model.TickSTRFspectro.to{mm};
        TickSTRFspectro.fo=Model.TickSTRFspectro.fo{mm};
    end
    if FIG>=4
        % Check spectro of stims
        for ss = 1:length(Spectro.spec)
            figure(12)
            STIM=reshape(Spectro.spec{ss}, length(Spectro.fo{ss}), length(Spectro.to{ss}));
            valc = max(abs(max(max(STIM))), abs(min(min(STIM))));
            if valc==0
                imagesc(Spectro.to{ss}*1000,Spectro.fo{ss}, STIM)
                axis xy
            else
                imagesc(Spectro.to{ss}*1000, Spectro.fo{ss}, STIM,[-valc valc])
                axis xy
            end
            xlabel('Time (ms)')
            ylabel('Frequencies')
            title('Stim before cutting')
        end
    end
    
    % Initializing outputs
    Data.VOC{mm} = VocType(Stim_local);
    if strcmp(ParamModel.NeuroRes, 'count')
        y=cell(NbStim_local,1);
        ymean=nan(NbStim_local,1);
        StimRep = nan(NbStim_local,1);
        yvar=ymean;
    else
        y = nan(NbStim_local,1);%this matrix will contain the average spike rate in spikes/ms at that precise position and for all the stims choosen for that run
    end
    
    % Getting spectro and spike trains
    for ss = 1:NbStim_local
        dd=Stim_local(ss);
        if any(ParamModel.ModelChoice([1,3:5]))
            %new spectro
            MatSpec = reshape(Spectro.spec{dd}, length(Spectro.fo{dd}), length(Spectro.to{dd}));
            FreqBelowFlow = find(Spectro.fo{dd}<=Flow);
            EndFreq = FreqBelowFlow(end);
            NFreq_local=length(FreqBelowFlow);
            if NFreq_local~=Df
                sprintf('WARNING!! Trouble with the size of the spectros for stim %d\n', dd);
            end
            TimeBelowWin = find((1000.*Spectro.to{dd})<= Win);
            EndTime = TimeBelowWin(end);
            NTime_local = length(TimeBelowWin);
            if NTime_local~=Dt
                sprintf('WARNING!! Trouble with the size of the spectros for stim %d\n', dd);
            end
            Newspectro=MatSpec(1:EndFreq,1:EndTime);
            x(ss,:)=reshape(Newspectro, 1, NFreq_local*NTime_local);
        end
        % Values of max spike rate(y), mean spike rate (y) and exact number of spike per trial (y) within the window
        % Check the distribution of responses (Gaussian or Poisson) for each stim
        if strcmp(ParamModel.NeuroRes, 'max')
            y(ss) = max(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));
        elseif strcmp(ParamModel.NeuroRes, 'mean')
            y(ss) = mean(PSTH{dd}((Win-MinWin+ResDelay):(Win+ResDelay)));% here we get the average Spike Rate over bins of 1ms so the spike rate is in spike/ms
        elseif strcmp(ParamModel.NeuroRes, 'count')
            y{ss}=nan(length(Trials{dd}),1);
            StimRep(ss) = length(Trials{dd});
            for tt=1:length(Trials{dd})
                y{ss}(tt)=sum((Trials{dd}{tt}>(Win-MinWin+ResDelay)).*(Trials{dd}{tt}<(Win+ResDelay)));
            end
            ymean(ss)=mean(y{ss});
            yvar(ss)=var(y{ss});
        else
            fprintf('please correctly write what kind of neural response you want to predict\n %s does not make any sense!!\n', ParamModel.NeuroRes);
            
        end
    end
    InputData.y=y;
    
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
    
    if any(ParamModel.ModelChoice([1,3:5]))
        % Take the log of the spectros
        %and ground the output to supress -Inf values
        % Substract the average spectro if requested
        % Save x as input for models
        InputData.x = 20*log10(abs(x));
        MAXI = max(max(InputData.x));
        InputData.x(InputData.x<(MAXI-80))=MAXI-80;
        if ParamModel.MeanSubstractSpec
            AvSpec_local = repmat(reshape(AvSpec(1:Df, 1:Dt), 1, Df*Dt),NbStim_local,1);
            InputData.raw_x = InputData.x;
            InputData.x = InputData.x - AvSpec_local;
        end
    end
    
    % Code the categorical data about call
    UVOC = unique(Data.VOC{mm});
    X_voc = zeros(length(Data.VOC{mm}), length(UVOC)-1);% Note that all zeros correspond to Ag call
    for vv=1:length(Data.VOC{mm})
        for uv=2:length(UVOC)
            if strcmp(Data.VOC{mm}(vv), UVOC{uv})
                X_voc(vv,uv-1)=1;
                break
            end
        end
    end
    InputData.X_voc=X_voc;
    
    % Calculate the entropy of the categories
    NbCallperCat = sum(X_voc,1);
    NbCallperCatall = [NbCallperCat length(Data.VOC{mm})-sum(NbCallperCat)];
    PCallperCat = NbCallperCatall./length(Data.VOC{mm});
    Data.semanticcategory_entropy(mm) = -sum(PCallperCat.*log2(PCallperCat));
    
    
    if FIG>=3 && any(ParamModel.ModelChoice([1,3:5]))
        for ss = 1:NbStim_local
            figure(11)
            subplot(1,4,1);
            STIM=reshape(x(ss,:), length(Model.TickSTRFspectro.fo{mm}), length(Model.TickSTRFspectro.to{mm}));
            M=max(max(abs(STIM)));
            imagesc(Model.TickSTRFspectro.to{mm}*1000,Model.TickSTRFspectro.fo{mm}, STIM, [-M M])
            axis xy
            xlabel('Time (ms)')
            ylabel('Frequencies')
            title('Stim before log')
            subplot(1,4,2)
            STIM=reshape(InputData.x(ss,:), length(Model.TickSTRFspectro.fo{mm}), length(Model.TickSTRFspectro.to{mm}));
            M=max(max(abs(STIM)));
            imagesc(Model.TickSTRFspectro.to{mm}*1000,Model.TickSTRFspectro.fo{mm}, STIM,[-M M])
            axis xy
            xlabel('Time (ms)')
            ylabel('Frequencies')
            title('Stim after log')
            
            ss3=subplot(1,4,3);
            bar(InputData.X_voc(ss,:));
            set(ss3,'XtickLabel',UVOC(2:end))
            xlabel(sprintf('Vocalization Type: %s',Data.VOC{mm}{ss}))
            subplot(1,4,4)
            bar(InputData.y{ss});
            xlabel('Trials')
            ylabel('Spike Count')
            pause();
        end
    end
    
    %% Calculate all models using poisson lasso glm and cross-validation
    % Loop through alphas
    for aa=1:length(Alphas)
        Alpha=Alphas(aa);
        fprintf('Alpha=%f\n',Alpha);
        
        %% Determine the best Lambda value for each model based on NbBootstrap_Lambda bootstrap
        if any(ParamModel.ModelChoice([1,3:5]))
            ParamModel.CV=1;%Set the cross-validation argument to 1 here
            ParamModel.BootstrapSTRF=NbBootstrap_Lambda;
            ParamModel.LAMBDA=cell(length(ParamModel.ModelChoice),1);%set to [] if you want to explore lamba values
            [MResults]=runLassoGlmModels(Alpha,VocType,Emitter,Stim_local,SWITCH, InputData, ParamModel,PreviousWin,TickSTRFspectro);
            
            if ParamModel.ModelChoice(1)
                ParamModel.LAMBDA{1}=median(MResults.Lambda_BestModel_Acoustic);
                LambdaChoice.Acoustic.Lambdas{mm,aa} = MResults.Lambda_BestModel_Acoustic;
                LambdaChoice.Acoustic.Speed{mm,aa} = MResults.Speed_Acoustic;
            end
            if ParamModel.ModelChoice(3)
                ParamModel.LAMBDA{3}=median(MResults.Lambda_BestModel_AcSem);
                LambdaChoice.AcSem.Lambdas{mm,aa} = MResults.Lambda_BestModel_AcSem;
            end
            if ParamModel.ModelChoice(4)
                ParamModel.LAMBDA{4}=median(MResults.Lambda_BestModel_AcSemAc);
                LambdaChoice.AcSemAc.Lambdas{mm,aa} = MResults.Lambda_BestModel_AcSemAc;
            end
            if ParamModel.ModelChoice(5)
                ParamModel.LAMBDA{5}=median(MResults.Lambda_BestModel_AcSemSem);
                LambdaChoice.AcSemSem.Lambdas{mm,aa} = MResults.Lambda_BestModel_AcSemSem;
            end
            LambdaChoice.ValSets{mm,aa} = MResults.ValSets;
            LambdaChoice.TrainSets{mm,aa} = MResults.TrainSets;
        end
        
        %% Re-run the models with the best lambda for each model to estimate the variability of the LL over NbBootstrap_Deviance bootstraps
        ParamModel.CV=1;%Set the cross-validation argument to 1 here
        ParamModel.BootstrapSTRF=NbBootstrap_Deviance;
        if any(ParamModel.ModelChoice([1,3:5]))
            [MResults]=runLassoGlmModels(Alpha,VocType,Emitter,Stim_local,SWITCH, InputData, ParamModel,PreviousWin,TickSTRFspectro);
        else
            [MResults]=runLassoGlmModels(Alpha,VocType,Emitter,Stim_local,SWITCH, InputData, ParamModel,PreviousWin);
        end
        
        
        PropVal.mean(mm,aa) = mean(MResults.ValProp);% Store the average proportion of stims in the validating set over bootstrap for that alpha
        PropVal.std(mm,aa) = std(MResults.ValProp);% Store the std on proportion of stims in the validating set over bootstrap for that alpha
        PropVal.values{mm,aa} = MResults.ValSize;
        
        
        if ParamModel.ModelChoice(1)
            Deviance.Acoustic.values{mm,aa}=MResults.Deviance_BestModel_Acoustic;% Store the deviance of all bootstraps for that alpha value and the best lambda
            Deviance.Acoustic.lambda{mm,aa}=ParamModel.LAMBDA{1};% Store the lambda used at of all bootstraps for that alpha value
            Deviance.Acoustic.ypredictVal{mm,aa} = MResults.Ypredict_Acoustic;
            Deviance.Acoustic.Speed{mm,aa} = MResults.Speed_Acoustic;
            LL.Acoustic.values{mm,aa}=MResults.LL_BestModel_Acoustic;
            Deviance.Acoustic.FitIndex(mm,aa)=(sum(MResults.LL_BestModel_Acoustic)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_BestGuess)-sum(MResults.LL_WorseGuess));
            if SWITCH.AR
                %Deviance.Acoustic.FitIndexAR(mm,aa)=(sum(MResults.LL_BestModel_Acoustic)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_AR)-sum(MResults.LL_WorseGuess));
                Deviance.Acoustic.FitIndexARVal(mm,aa)=(sum(MResults.LL_BestModel_Acoustic)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_ARVal)-sum(MResults.LL_WorseGuess));
            end
            if SWITCH.DF
                Deviance.Acoustic.DF{mm,aa} = MResults.DF_BestModel_Acoustic;
            end
        end
        if ParamModel.ModelChoice(2)
            Deviance.Semantic.values{mm,aa}=MResults.Deviance_BestModel_Sem;% Store the deviance of all bootstraps for that alpha value and the best lambda
            Deviance.Semantic.ypredictVal{mm,aa} = MResults.Ypredict_Semantic;
            LL.Semantic.values{mm,aa}=MResults.LL_BestModel_Sem;
            Deviance.Semantic.FitIndex(mm,aa)=(sum(MResults.LL_BestModel_Sem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_BestGuess)-sum(MResults.LL_WorseGuess));
            if SWITCH.AR
                %Deviance.Semantic.FitIndexAR(mm,aa)=(sum(MResults.LL_BestModel_Sem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_AR)-sum(MResults.LL_WorseGuess));
                Deviance.Semantic.FitIndexARVal(mm,aa)=(sum(MResults.LL_BestModel_Sem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_ARVal)-sum(MResults.LL_WorseGuess));
            end
            if SWITCH.DF
                Deviance.Semantic.DF{mm,aa} = MResults.DF_BestModel_Sem;
            end
        end
        if ParamModel.ModelChoice(3)
            Deviance.AcSem.values{mm,aa}=MResults.Deviance_BestModel_AcSem;% Store the deviance of all bootstraps for that alpha value
            Deviance.AcSem.lambda{mm,aa}=ParamModel.LAMBDA{3};% Store the lambdas used at of all bootstraps for that alpha value
            Deviance.AcSem.ypredictVal{mm,aa} = MResults.Ypredict_AcSem;
            LL.AcSem.values{mm,aa}=MResults.LL_BestModel_AcSem;
            Deviance.AcSem.FitIndex(mm,aa)=(sum(MResults.LL_BestModel_AcSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_BestGuess)-sum(MResults.LL_WorseGuess));
            if SWITCH.AR
                %Deviance.AcSem.FitIndexAR(mm,aa)=(sum(MResults.LL_BestModel_AcSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_AR)-sum(MResults.LL_WorseGuess));
                Deviance.AcSem.FitIndexARVal(mm,aa)=(sum(MResults.LL_BestModel_AcSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_ARVal)-sum(MResults.LL_WorseGuess));
            end
            if SWITCH.DF
                Deviance.AcSem.DF{mm,aa} = MResults.DF_BestModel_AcSem;
            end
        end
        if ParamModel.ModelChoice(4)
            Deviance.AcSemAc.values{mm,aa}=MResults.Deviance_BestModel_AcSemAc;% Store the deviance of all bootstraps for that alpha value
            Deviance.AcSemAc.lambda{mm,aa}=ParamModel.LAMBDA{4};% Store the lambdas used at of all bootstraps for that alpha value
            Deviance.AcSemAc.ypredictVal{mm,aa} = MResults.Ypredict_AcSemAc;
            LL.AcSemAc.values{mm,aa}=MResults.LL_BestModel_AcSemAc;
            Deviance.AcSemAc.FitIndex(mm,aa)=(sum(MResults.LL_BestModel_AcSemAc)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_BestGuess)-sum(MResults.LL_WorseGuess));
            if SWITCH.AR
                %Deviance.AcSemAc.FitIndexAR(mm,aa)=(sum(MResults.LL_BestModel_AcSemAc)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_AR)-sum(MResults.LL_WorseGuess));
                Deviance.AcSemAc.FitIndexARVal(mm,aa)=(sum(MResults.LL_BestModel_AcSemAc)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_ARVal)-sum(MResults.LL_WorseGuess));
            end
            if SWITCH.DF
                Deviance.AcSemAc.DF{mm,aa} = MResults.DF_BestModel_AcSemAc;
            end
        end
        if ParamModel.ModelChoice(5)
            Deviance.AcSemSem.values{mm,aa}=MResults.Deviance_BestModel_AcSemSem;% Store the deviance of all bootstraps for that alpha value
            Deviance.AcSemSem.lambda{mm,aa}=ParamModel.LAMBDA{5};% Store the lambdas used at of all bootstraps for that alpha value
            Deviance.AcSemSem.ypredictVal{mm,aa} = MResults.Ypredict_AcSemSem;
            LL.AcSemSem.values{mm,aa}=MResults.LL_BestModel_AcSemSem;
            Deviance.AcSemSem.FitIndex(mm,aa)=(sum(MResults.LL_BestModel_AcSemSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_BestGuess)-sum(MResults.LL_WorseGuess));
            if SWITCH.AR
                %Deviance.AcSemSem.FitIndexAR(mm,aa)=(sum(MResults.LL_BestModel_AcSemSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_AR)-sum(MResults.LL_WorseGuess));
                Deviance.AcSemSem.FitIndexARVal(mm,aa)=(sum(MResults.LL_BestModel_AcSemSem)-sum(MResults.LL_WorseGuess))/(sum(MResults.LL_ARVal)-sum(MResults.LL_WorseGuess));
            end
            if SWITCH.DF
                Deviance.AcSemSem.DF{mm,aa} = MResults.DF_BestModel_AcSemSem;
            end
        end
        
        Deviance.yobsVal{mm,aa} = MResults.YobsVal;
        Deviance.yobsCat{mm,aa} = MResults.XVocVal;
        Deviance.YCeilModelVal{mm,aa} = MResults.YCeilModelVal;
        if SWITCH.AR
            Deviance.Ypredict_ARVal{mm,aa} = MResults.Ypredict_ARVal;
        end
        Deviance.Ypredict_Floor{mm,aa} = MResults.Ypredict_Floor;
        LL.Ceiling.values{mm,aa} = MResults.LL_BestGuess;
        LL.Floor.values{mm,aa} = MResults.LL_WorseGuess;
        if SWITCH.AR
            LL.AutoRegressiveVal.values{mm,aa}=MResults.LL_ARVal;
        end
        LL.Saturated.values{mm,aa} = MResults.LL_Saturated;
        LL.NumStimTrials{mm,aa} = MResults.NumStimTrials;
        %LL.AutoRegressiveTV.values{mm,aa}=MResults.LL_AR;
        Deviance.TrainSets{mm,aa} = MResults.TrainSets;
        Deviance.ValSets{mm,aa} = MResults.ValSets;
        
        
        %% plot LogLikelihood of all models over bootstraps
        if FIG>0
            models=fieldnames(LL);
            %Colors='brmgcky';
            close(figure(14))
            Ll=[];
            figure(14)
            hold on
            if ParamModel.ModelChoice(1)
                plot(1:ParamModel.BootstrapSTRF,LL.Acoustic.values{mm,aa},'b.-')
                hline(mean(LL.Acoustic.values{mm,aa}),'b:')
                Ll=[Ll LL.Acoustic.values{mm,aa}];
            end
            if ParamModel.ModelChoice(2)
                plot(1:ParamModel.BootstrapSTRF,LL.Semantic.values{mm,aa},'g.-')
                hline(mean(LL.Semantic.values{mm,aa}),'g:')
                Ll = [Ll LL.Semantic.values{mm,aa}];
            end
            if ParamModel.ModelChoice(3)
                plot(1:ParamModel.BootstrapSTRF,LL.AcSem.values{mm,aa},'m.-')
                hline(mean(LL.AcSem.values{mm,aa}),'m:')
                Ll = [Ll LL.AcSem.values{mm,aa}];
            end
            if ParamModel.ModelChoice(4)
                plot(1:ParamModel.BootstrapSTRF,LL.AcSemAc.values{mm,aa},'.-','Color',[0.5 0 1])
                hline_color(mean(LL.AcSemAc.values{mm,aa}),':','','Color',[0.5 0 1])
                Ll = [Ll LL.AcSemAc.values{mm,aa}];
            end
            if ParamModel.ModelChoice(5)
                plot(1:ParamModel.BootstrapSTRF,LL.AcSemSem.values{mm,aa},'r.-')
                hline(mean(LL.AcSemSem.values{mm,aa}),'r:')
                Ll = [Ll LL.AcSemSem.values{mm,aa}];
            end
            plot(1:ParamModel.BootstrapSTRF,LL.Ceiling.values{mm,aa},'y.-',1:ParamModel.BootstrapSTRF,LL.Floor.values{mm,aa},'k.-',1:ParamModel.BootstrapSTRF,LL.AutoRegressiveVal.values{mm,aa},'c.--');
            legend(models);
            hline(mean(LL.Floor.values{mm,aa}),'k:')
            hline(mean(LL.Ceiling.values{mm,aa}),'y:')
            if SWITCH.AR
                hline(mean(LL.AutoRegressiveVal.values{mm,aa}),'c:o')
            end
            ylabel('LogLikelihood');
            xlabel('Bootstrap');
            axis([1 ParamModel.BootstrapSTRF max(mean(Ll))-2*min(std(Ll)) 0])
            hold off
        end
        
        
        %% Use the best Lambda to calculate the optimal models on all dataset
        ParamModel.CV=0;%Set the cross-validation argument to 0 here
        ParamModel.BootstrapSTRF=1;
        if any(ParamModel.ModelChoice([1,3:5]))
            [MResultsOptimal]=runLassoGlmModels1L(Alpha,VocType,Emitter,Stim_local,SWITCH, InputData, ParamModel,PreviousWin,TickSTRFspectro);
        else
            [MResultsOptimal]=runLassoGlmModels1L(Alpha,VocType,Emitter,Stim_local,SWITCH, InputData, ParamModel,PreviousWin);
        end
        clear ff
        if ParamModel.ModelChoice(1)
            Model.Acoustic.B{mm,aa}=MResultsOptimal.ModelB_Ac{1};
            Model.Acoustic.Lambdas(mm,aa) = MResultsOptimal.Lambda_BestModel_Acoustic;
            Model.Acoustic.B0{mm,aa} = MResultsOptimal.ModelB0_Ac;
            Model.Acoustic.Speed{mm,aa}=MResultsOptimal.Speed_Acoustic;
            Model.Acoustic.y_predict{mm,aa}=MResultsOptimal.Ypredict_Acoustic{1};
            if FIG>0
                figure(15)
                imagesc(Model.Acoustic.B{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AC lassoglm poisson whole dataset\nlog of lambda=%f',log(Model.Acoustic.Lambdas(mm,aa))))
                SpectroDim=size(Model.Acoustic.B{mm,aa});
                if ParamModel.ZC==1
                    ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./MResultsOptimal.X_stdZC{1},size(MResultsOptimal.XSpecVal{1},1),1);
                else
                    ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm,aa}, 1, SpectroDim(1)*SpectroDim(2)),size(MResultsOptimal.XSpecVal{1},1),1);
                end
                Predict_Ac_wholeset = sum(ModelAcoustic_local.*MResultsOptimal.XSpecVal{1},2) + repmat(Model.Acoustic.B0{mm,aa},size(MResultsOptimal.XSpecVal{1},1),1);
                figure(20)
                if exist('ff','var')
                    ff=ff+1;
                else
                    ff=1;
                end
                if mod(sum(ParamModel.ModelChoice),2)==0
                    ss=subplot(2,sum(ParamModel.ModelChoice)/2,ff);
                else
                    ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,ff);
                end
                plot(Predict_Ac_wholeset,MResultsOptimal.YobsVal{1},'bo','MarkerSize',5)
                ylabel('Observed Spike count ')
                xlabel('Predicted Spike count wo output non-linearity')
                title(sprintf('Acoustic Model Alpha %f',Alphas(aa)))
                AXES=get(ss,'yLim');
                hold on
                plot(Predict_Ac_wholeset,exp(Predict_Ac_wholeset),'r.')
                set(ss,'yLim',AXES)
                legend('Data','exponential fit')
            end
        end
        
        if ParamModel.ModelChoice(2)
            Model.Semantic.B{mm,aa}=MResultsOptimal.ModelB_Sem{1}(2:end);
            Model.Semantic.B0{mm,aa}=MResultsOptimal.ModelB_Sem{1}(1);
            Model.Semantic.y_predict{mm,aa}=MResultsOptimal.Ypredict_Semantic{1};
            if FIG>0
                % Code the categorical data about call
                X_voc_local = zeros(length(MResultsOptimal.XVocVal{1}), length(UVOC)-1);% Note that all zeros correspond to Ag call
                for vv=1:length(MResultsOptimal.XVocVal{1})
                    for uv=2:length(UVOC)
                        if strcmp(MResultsOptimal.XVocVal{1}(vv), UVOC{uv})
                            X_voc_local(vv,uv-1)=1;
                            break
                        end
                    end
                end
                
                figure(16)
                ModelSemantic = repmat(Model.Semantic.B0{mm,aa},1,length(Model.Semantic.B{mm,aa})+1) + [0 Model.Semantic.B{mm,aa}'];
                plot(ModelSemantic,1:length(UVOC));
                title('Optimal Semantic Model');
                set(gca,'YTickLabel', UVOC);
                ModelSem_local=repmat(Model.Semantic.B{mm,aa}',size(MResultsOptimal.XSpecVal{1},1),1);
                Predict_Sem_wholeset = sum(ModelSem_local.*X_voc_local,2) + repmat(Model.Semantic.B0{mm,aa},size(MResultsOptimal.XSpecVal{1},1),1);
                figure(20)
                if exist('ff','var')
                    ff=ff+1;
                else
                    ff=1;
                end
                if mod(sum(ParamModel.ModelChoice),2)==0
                    ss=subplot(2,sum(ParamModel.ModelChoice)/2,ff);
                else
                    ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,ff);
                end
                plot(Predict_Sem_wholeset,MResultsOptimal.YobsVal{1},'bo','MarkerSize',10)
                ylabel('Observed Spike count wo output non-linearity')
                xlabel('Predicted Spike count')
                title(sprintf('Semantic Model Alpha %f',Alphas(aa)))
                AXES=get(ss,'yLim');
                hold on
                plot(Predict_Sem_wholeset,exp(Predict_Sem_wholeset),'r.')
                set(ss,'yLim',AXES)
                legend('Data','exponential fit')
            end
        end
        if ParamModel.ModelChoice(3)
            Model.AcSem.Bspectro{mm,aa}=MResultsOptimal.ModelBspectro_AcSem{1};
            Model.AcSem.Bsem{mm,aa}=MResultsOptimal.ModelBsem_AcSem{1};
            Model.AcSem.Lambdas(mm,aa) = MResultsOptimal.Lambda_BestModel_AcSem;
            Model.AcSem.B0{mm,aa} = MResultsOptimal.ModelB0_AcSem;
            Model.AcSem.y_predict{mm,aa}=MResultsOptimal.Ypredict_AcSem{1};
            if FIG>0
                figure(17)
                subplot(1,2,1)
                imagesc(Model.AcSem.Bspectro{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AcSem lassoglm poisson whole dataset\nlog10 of lambda=%f\n',log10(Model.AcSem.Lambdas(mm,aa))))
                ss=subplot(1,2,2);
                ModelAcSem = repmat(Model.AcSem.B0{mm,aa},1,length(Model.AcSem.Bsem{mm,aa})+1) + [0 Model.AcSem.Bsem{mm,aa}'];
                plot(ModelAcSem,1:length(UVOC));
                title('Optimal AcSem Model');
                set(ss,'YTickLabel', UVOC);
                if ParamModel.ZC==1
                    ModelAcSem_local=repmat([Model.AcSem.Bsem{mm,aa}' reshape(Model.AcSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./MResultsOptimal.X_stdZC{1}],size(MResultsOptimal.XSpecVal{1},1),1);
                else
                    ModelAcSem_local=repmat([Model.AcSem.Bsem{mm,aa}' reshape(Model.AcSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))],size(MResultsOptimal.XSpecVal{1},1),1);
                end
                Predict_AcSem_wholeset = sum(ModelAcSem_local.*[MResultsOptimal.XVocVal{1} MResultsOptimal.XSpecVal{1}],2) + repmat(Model.AcSem.B0{mm,aa},size(MResultsOptimal.XSpecVal{1},1),1);
                figure(20)
                if exist('ff','var')
                    ff=ff+1;
                else
                    ff=1;
                end
                if mod(sum(ParamModel.ModelChoice),2)==0
                    ss=subplot(2,sum(ParamModel.ModelChoice)/2,ff);
                else
                    ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,ff);
                end
                plot(Predict_AcSem_wholeset,MResultsOptimal.YobsVal{1},'bo','MarkerSize',10)
                ylabel('Observed Spike count wo output non-linearity')
                xlabel('Predicted Spike count')
                title(sprintf('AcSem Model Alpha %f',Alphas(aa)))
                AXES=get(ss,'yLim');
                hold on
                plot(Predict_AcSem_wholeset,exp(Predict_AcSem_wholeset),'r.')
                set(ss,'yLim',AXES)
                legend('Data','exponential fit')
                
            end
        end
        
        if ParamModel.ModelChoice(4)
            Model.AcSemAc.Bspectro{mm,aa}=MResultsOptimal.ModelBspectro_AcSemAc{1};
            Model.AcSemAc.Bsem{mm,aa}=MResultsOptimal.ModelBsem_AcSemAc{1};
            Model.AcSemAc.Lambdas(mm,aa) = MResultsOptimal.Lambda_BestModel_AcSemAc;
            Model.AcSemAc.B0{mm,aa} = MResultsOptimal.ModelB0_AcSemAc;
            Model.AcSemAc.y_predict{mm,aa}=MResultsOptimal.Ypredict_AcSemAc{1};
            if FIG>0
                figure(18)
                subplot(1,2,1)
                imagesc(Model.AcSemAc.Bspectro{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AcSem AcOffset lassoglm poisson whole dataset\nlog10 of lambda=%f',log10(Model.AcSemAc.Lambdas(mm,aa))))
                ss=subplot(1,2,2);
                ModelAcSemAc = repmat(Model.AcSemAc.B0{mm,aa},1,length(Model.AcSemAc.Bsem{mm,aa})+1) + [0 Model.AcSemAc.Bsem{mm,aa}'];
                plot(ModelAcSemAc,1:length(UVOC));
                title('Optimal AcSem Model Ac Offset');
                set(ss,'YTickLabel', UVOC);
                if ParamModel.ZC==1
                    ModelAcSemAc_local=repmat([Model.AcSemAc.Bsem{mm,aa}' reshape(Model.AcSemAc.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./MResultsOptimal.X_stdZC{1}],size(MResultsOptimal.XSpecVal{1},1),1);
                else
                    ModelAcSemAc_local=repmat([Model.AcSemAc.Bsem{mm,aa}' reshape(Model.AcSemAc.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))],size(MResultsOptimal.XSpecVal{1},1),1);
                end
                Predict_AcSemAc_wholeset = sum(ModelAcSemAc_local.*[X_voc_local MResultsOptimal.XSpecVal{1}],2) + repmat(Model.AcSemAc.B0{mm,aa},size(MResultsOptimal.XSpecVal{1},1),1);
                figure(20)
                if exist('ff','var')
                    ff=ff+1;
                else
                    ff=1;
                end
                if mod(sum(ParamModel.ModelChoice),2)==0
                    ss=subplot(2,sum(ParamModel.ModelChoice)/2,ff);
                else
                    ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,ff);
                end
                plot(Predict_AcSemAc_wholeset,MResultsOptimal.YobsVal{1},'bo','MarkerSize',10)
                ylabel('Observed Spike count wo output non-linearity')
                xlabel('Predicted Spike count')
                title(sprintf('AcSemAc Model Alpha %f',Alphas(aa)))
                AXES=get(ss,'yLim');
                hold on
                plot(Predict_AcSemAc_wholeset,exp(Predict_AcSemAc_wholeset),'r.')
                set(ss,'yLim',AXES)
                legend('Data','exponential fit')
            end
        end
        
        if ParamModel.ModelChoice(5)
            Model.AcSemSem.Bspectro{mm,aa}=MResultsOptimal.ModelBspectro_AcSemSem{1};
            Model.AcSemSem.Bsem{mm,aa}=MResultsOptimal.ModelBsem_AcSemSem{1};
            Model.AcSemSem.Lambdas(mm,aa) = MResultsOptimal.Lambda_BestModel_AcSemSem;
            Model.AcSemSem.B0{mm,aa} = MResultsOptimal.ModelB0_AcSemSem;
            Model.AcSemSem.y_predict{mm,aa}=MResultsOptimal.Ypredict_AcSemSem{1};
            if FIG>0
                figure(19)
                subplot(1,2,1)
                imagesc(Model.AcSemSem.Bspectro{mm,aa})
                axis xy
                title(sprintf('OPTIMAL STRF AcSem SemOffset lassoglm poisson whole dataset\nlog10 of lambda=%f',log10(Model.AcSemSem.Lambdas(mm,aa))))
                ss=subplot(1,2,2);
                ModelAcSemSem = repmat(Model.AcSemSem.B0{mm,aa},1,length(Model.AcSemSem.Bsem{mm,aa})+1) + [0 Model.AcSemSem.Bsem{mm,aa}'];
                plot(ModelAcSemSem,1:length(UVOC));
                title('Optimal AcSem Model Ac Offset');
                set(ss,'YTickLabel', UVOC);
                if ParamModel.ZC==1
                    ModelAcSemSem_local=repmat([Model.AcSemSem.Bsem{mm,aa}' reshape(Model.AcSemSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./MResultsOptimal.X_stdZC{1}],size(MResultsOptimal.XSpecVal{1},1),1);
                else
                    ModelAcSemSem_local=repmat([Model.AcSemSem.Bsem{mm,aa}' reshape(Model.AcSemSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))],size(MResultsOptimal.XSpecVal{1},1),1);
                end
                Predict_AcSemSem_wholeset = sum(ModelAcSemSem_local.*[X_voc_local MResultsOptimal.XSpecVal{1}],2) + repmat(Model.AcSemSem.B0{mm,aa},size(MResultsOptimal.XSpecVal{1},1),1);
                figure(20)
                if exist('ff','var')
                    ff=ff+1;
                else
                    ff=1;
                end
                if mod(sum(ParamModel.ModelChoice),2)==0
                    ss=subplot(2,sum(ParamModel.ModelChoice)/2,ff);
                else
                    ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,ff);
                end
                plot(Predict_AcSemSem_wholeset,MResultsOptimal.YobsVal{1},'bo','MarkerSize',10)
                ylabel('Observed Spike count wo output non-linearity')
                xlabel('Predicted Spike count')
                title(sprintf('AcSem Model Alpha %f',Alphas(aa)))
                AXES=get(ss,'yLim');
                hold on
                plot(Predict_AcSemSem_wholeset,exp(Predict_AcSemSem_wholeset),'r.')
                set(ss,'yLim',AXES)
                legend('Data','exponential fit')
                pause()
            end
        end
        
        %% Calculate the information for each model
        y_wholeset_bestGuess_perstim = MResultsOptimal.YCeilModelVal{1}(MResultsOptimal.ValSetFirstRep{1});
        MaxYpredictInfo = max(y_wholeset_bestGuess_perstim);
        
        if ParamModel.ModelChoice(1)
            % Acoustic model
            fprintf('**Info on Acoustic**\n')
            [Model.Acoustic.info(mm,aa),Model.Acoustic.P_YgivenS_all1{mm,aa},Model.Acoustic.P_YgivenS_all2{mm,aa}] = info_model_Calculus(Model.Acoustic.y_predict{mm,aa}(MResultsOptimal.ValSetFirstRep{1}),MaxYpredictInfo,MinWin);
            if mm>1 && SWITCH.AllAlpha
                [Model.Acoustic.cum_info(mm,aa)]=info_cumulative_model_Calculus(Model.Acoustic.P_YgivenS_all1(1:mm,aa), mm,Data.x_stim_indices_wholeset, Stim_local);
            elseif SWITCH.AllAlpha
                Model.Acoustic.cum_info(mm,aa) = Model.Acoustic.info(mm,aa);
            end
        end
        
        if ParamModel.ModelChoice(2)
            % Semantic Model
            fprintf('**Info on Semantic**\n')
            AllStims_data=Model.Semantic.y_predict{mm,aa}(MResultsOptimal.ValSetFirstRep{1});
            [Model.Semantic.info(mm,aa),Model.Semantic.P_YgivenS_all1{mm,aa},Model.Semantic.P_YgivenS_all2{mm,aa}]  = info_model_Calculus(AllStims_data,MaxYpredictInfo,MinWin);
%             EndStims_data=nan(length(Stim_local_end),1);
%             for ss = 1:length(Stim_local_end)
%                 EndStims_data(ss) = AllStims_data(find(Stim_local==Stim_local_end(ss)));
%             end
%             [Model.Semantic.info_enddataset(mm,aa),Model.Semantic.P_YgivenS_all1_enddataset{mm,aa},Model.Semantic.P_YgivenS_all2_enddataset{mm,aa}]  = info_model_Calculus(EndStims_data,MaxYpredictInfo,MinWin);
%             
            if mm>1 && SWITCH.AllAlpha
                [Model.Semantic.cum_info(mm,aa)]=info_cumulative_model_Calculus(Model.Semantic.P_YgivenS_all1(1:mm,aa), mm,Data.x_stim_indices_wholeset, Stim_local);
            elseif SWITCH.AllAlpha
                Model.Semantic.cum_info(mm,aa) = Model.Semantic.info(mm,aa);
            end
        end
        
        if ParamModel.ModelChoice(4)
            % AcSemAc
            fprintf('**Info on AcSemAc**\n')
            [Model.AcSemAc.info(mm,aa),Model.AcSemAc.P_YgivenS_all1{mm,aa},Model.AcSemAc.P_YgivenS_all2{mm,aa}] = info_model_Calculus(Model.AcSemAc.y_predict{mm,aa}(MResultsOptimal.ValSetFirstRep{1}),MaxYpredictInfo,MinWin);
            if mm>1 && SWITCH.AllAlpha
                [Model.AcSemAc.cum_info(mm,aa)]=info_cumulative_model_Calculus(Model.AcSemAc.P_YgivenS_all1(1:mm,aa), mm,Data.x_stim_indices_wholeset, Stim_local);
            elseif SWITCH.AllAlpha
                Model.AcSemAc.cum_info(mm,aa) = Model.AcSemAc.info(mm,aa);
            end
        end
        
        if ParamModel.ModelChoice(5)
            % AcSemSem
            fprintf('**Info on AcSemSem**\n')
            [Model.AcSemSem.info(mm,aa),Model.AcSemSem.P_YgivenS_all1{mm,aa},Model.AcSemSem.P_YgivenS_all2{mm,aa}] = info_model_Calculus(Model.AcSemSem.y_predict{mm,aa}(MResultsOptimal.ValSetFirstRep{1}),MaxYpredictInfo,MinWin);
            if mm>1 && SWITCH.AllAlpha
                [Model.AcSemSem.cum_info(mm,aa)]=info_cumulative_model_Calculus(Model.AcSemSem.P_YgivenS_all1(1:mm,aa), mm,Data.x_stim_indices_wholeset, Stim_local);
            elseif SWITCH.AllAlpha
                Model.AcSemSem.cum_info(mm,aa) = Model.AcSemSem.info(mm,aa);
            end
        end
        
        
    end
    
    %% Save the whole dataset values on which optimal models are run (same for all alphas so just take the last one)
    %Data.x_std_wholeset{mm}=MResultsOptimal.X_stdZC{1};
    %Data.x_mean_wholeset{mm}=MResultsOptimal.X_meanZC{1};
    Data.InputData{mm} = InputData;%Maybe you want to supress this line as everything can be retrieve from the other saved data
    Data.x_stim_indices_wholeset{mm} = Stim_local;
    Data.x_stim_indices_wholeset{end} = Stim_local_end;
    Data.x_stim_repetition{mm} = StimRep;
    if any(ParamModel.ModelChoice([1,3:5]))
        Data.FLow = Flow;
    end
    %Data.x_wholeset{mm}=MResultsOptimal.XSpecVal{1};
    Data.X_voc_wholeset{mm}=MResultsOptimal.XVocVal{1};
    Data.y_wholeset{mm}=MResultsOptimal.YobsVal{1};
    Data.y_wholeset_bestGuess{mm}=MResultsOptimal.YCeilModelVal{1};
    if SWITCH.AR
        Data.y_wholeset_AR{mm}=MResultsOptimal.Ypredict_ARVal{1};
    end
    Data.y_wholeset_Floor{mm}=MResultsOptimal.Ypredict_Floor{1};
    Data.y_ValSetFirstRep{mm} = MResultsOptimal.ValSetFirstRep{1};
    
    %% Calculate the information for models that do not depend on alpha (Floor, Ceiling, AR, Semantic)
    % Ceiling Model
    fprintf('**Info on Ceiling Model**\n')
    %[Model.Ceiling.info(mm),Model.Ceiling.P_YgivenS_all1{mm},Model.Ceiling.P_YgivenS_all2{mm}] = info_model_Calculus(Data.y_wholeset_bestGuess{mm}(Data.y_ValSetFirstRep{mm}),MaxYpredictInfo,MinWin);
    % The input data above makes sense more or less to compare with the AR
    % model (I'm not sure to understand while I did not take the average of
    % the predictions over trials instead of only taking the first one here
    % though....)
    % As of today 08/17/2016 I propose to work with the average observed
    % spike count over trials
    Data.y_wholeset_AvObsDataperStim{mm} = ymean;
    [Model.Ceiling.info(mm),Model.Ceiling.P_YgivenS_all1{mm},Model.Ceiling.P_YgivenS_all2{mm}] = info_model_Calculus(Data.y_wholeset_AvObsDataperStim{mm},MaxYpredictInfo,MinWin);
    % Calculating the info using only the dataset that is present until the
    % last window
%     EndStims_data=nan(length(Stim_local_end),1);
%     for ss = 1:length(Stim_local_end)
%         EndStims_data(ss) = ymean(find(Stim_local==Stim_local_end(ss)));
%     end
%     [Model.Ceiling.info_enddataset(mm),Model.Ceiling.P_YgivenS_all1_enddataset{mm},Model.Ceiling.P_YgivenS_all2_enddataset{mm}] = info_model_Calculus(EndStims_data,MaxYpredictInfo,MinWin);
    
    % Floor
    fprintf('**Info on Null-model**\n')
    %Model.Floor.y_predict{mm} = repmat(mean(Data.y_wholeset{mm}(Data.y_ValSetFirstRep{mm})),length(Data.x_stim_repetition{mm}),1);
    % The input data above does not make much sense, why taking only the
    % first trial of each stim to calculate the average response??
    % As of today 08/17/2016 I propose something else:
    Data.y_wholeset_AvObsData(mm) = mean(ymean);
    [Model.Floor.info(mm),Model.Floor.P_YgivenS_all1{mm},Model.Floor.P_YgivenS_all2{mm}] = info_model_Calculus(repmat(Data.y_wholeset_AvObsData(mm),NbStim_local,1), MaxYpredictInfo, MinWin);
    
    
    % AR Model
    if mm>1 && SWITCH.AR
        fprintf('**Info on AR**\n')
        [Model.AR.info(mm),Model.AR.P_YgivenS_all1{mm},Model.AR.P_YgivenS_all2{mm}] = info_model_Calculus(Data.y_wholeset_AR{mm}(Data.y_ValSetFirstRep{mm}),MaxYpredictInfo,MinWin);
    elseif SWITCH.AR
        Model.AR.P_YgivenS_all1{mm}=Model.Ceiling.P_YgivenS_all1{mm}; %initializing values of probabilities of the AR models with the ones from the Ceiling model to calculate a cumulative at Win=2
        Model.AR.info(mm) = Model.Ceiling.info(mm);
    end
    
    %% Calculate cumulative information
    if ParamModel.Cum_Info
        if mm==1
            fprintf('set CumInfo = Info for the first window\n');
            infofields = fieldnames(cum_info);
            if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
                % Acoustic model
                fprintf('**CumInfo on Acoustic**\n')
                for fn=1:numel(infofields)
                    Model.Acoustic.cum_info.(infofields{fn})(mm,1) = Model.Acoustic.info(mm,1);
                end
            end
            
            if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
                % Semantic Model
                fprintf('**CumInfo on Semantic**\n')
                for fn=1:numel(infofields)
                    Model.Semantic.cum_info.current_dataset.(infofields{fn})(mm,1) = Model.Semantic.info(mm,1);
%                     Model.Semantic.cum_info.end_dataset.(infofields{fn})(mm,1) = Model.Semantic.info(mm,1);
                end
            end
            
            if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
                % AcSemAc
                fprintf('**CumInfo on AcSemAc**\n')
                for fn=1:numel(infofields)
                    Model.AcSemAc.cum_info.(infofields{fn})(mm,1) = Model.AcSemAc.info(mm,1);
                end
            end
            
            if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
                % AcSemSem
                fprintf('**CumInfo on AcSemSem**\n')
                for fn=1:numel(infofields)
                    Model.AcSemSem.cum_info.(infofields{fn})(mm,1) = Model.AcSemSem.info(mm,1);
                end
            end
            
            %Floor
            for fn=1:numel(infofields)
                Model.Floor.cum_info.(infofields{fn})(mm,1) = Model.Floor.info(mm,1);
            end
            
            %Ceiling
            for fn=1:numel(infofields)
                Model.Ceiling.cum_info.current_dataset.(infofields{fn})(mm,1) = Model.Ceiling.info(mm,1);
%                 Model.Ceiling.cum_info.end_dataset.(infofields{fn})(mm,1) = Model.Ceiling.info(mm,1);
            end
            if SWITCH.AR
                for fn=1:numel(infofields)
                    Model.AR.cum_info.(infofields{fn})(mm,1) = Model.AR.info(mm,1);
                end
            end
        end
        if mm>1
            [Model] = info_cumulative_wrapper(ParamModel,SWITCH,Model,mm,Data.x_stim_indices_wholeset);
        end
        %         if mm==modNum
        %             %get rid of temp folder and its content that might be
        %             %created by info_cumulative_model_Calculus
        %             [status]=system(sprintf('rm -rf %s',FolderTempInfStorage));
        %         end
        
    end
    %% Store the average spectro if we z-scored the data with a mean and STD
    % calculated for each window only
    if any(ParamModel.ModelChoice([1,3:5]))
        if ParamModel.ZC
            Data.MeanSpectroStim{mm} = reshape(MResultsOptimal.X_meanZC{1},Df,Dt);
        elseif ParamModel.MeanSubstractSpec
            Data.MeanSpectroStim{mm} = AvSpec(1:Df, 1:Dt);
        end
    else
        Data.MeanSpectroStim{mm} = 'NoMeanSpectroCalculated Ac model was likely not run';    
    end
    
    
    %% keep track of the data and stim used for the next window
    PreviousWin.Stim_local_old = Stim_local;
    PreviousWin.y_old = y;
    
    %% Save what we have for now
    if saveonline
        save(calfilename,'LambdaChoice','Deviance','LL','Model','ParamModel','Data','PropVal','Wins','-append');
    end
    WindowElapsed=toc(WindowTime);
    Days = floor(WindowElapsed/(60*60*24));
    ETRem = WindowElapsed - Days*60*60*24;
    Hours = floor(ETRem/(60*60));
    ETRem = ETRem - Hours*60*60;
    Minutes = floor(ETRem/60);
    ETRem = ETRem-60*Minutes;
    fprintf(1,'********************************* Done with window %d/%d after the code run for %d days %dh %dmin and %dsec\n',mm, modNum,Days,Hours,Minutes,ETRem);
end

end

