function [R2A, SSres, SSexp, SStot, ModelPredict, LL, NEC, PvalLRatio, HLRatio, NeuroRes, VOC, NbOptPC, Pvalue, Wins, NeuralResponse, STRF_time, STRF_to, STRF_fo, ModSem] = GrowingModels(Spectro, VocType, PSTH, MinWin, MaxWin, Increment, ResDelay, NeuroRes)
FIG=1; % set to 1 for debugging figures
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

% define the range of number of PC to investigate
NbStim = length(VocType);


% Initialize a bunch of output variables
R2A.Acoustic = nan(modNum,1);% adjsuted R squared
R2A.Semantic = nan(modNum,1);
R2A.AcSem = nan(modNum,1);
LL = R2A;%Loglikelihood
Pvalue = R2A;%pvalue of the anova on the model
NEC = R2A;%Number of estimated coeeficients in the model
ModelPredict.Acoustic = cell(modNum,1);
ModelPredict.Semantic = cell(modNum,1);
ModelPredict.AcSem = cell(modNum,1);
NeuralResponse = cell(modNum,1);
STRF_time = cell(modNum,1);
STRF_to = cell(modNum,1);
STRF_fo = cell(modNum,1);
ModSem = cell(modNum,1);
NbOptPC = nan(1,modNum);
PvalLRatio.AcAcSem = nan(modNum,1);
PvalLRatio.SemAcSem = nan(modNum,1);
HLRatio.AcAcSem = nan(modNum,1);
HLRatio.SemAcSem = nan(modNum,1);
VOC = cell(modNum,1);

SSres.Acoustic = nan(modNum,1);
SSres.Semantic = nan(modNum,1);
SSres.AcSem = nan(modNum,1);
SSexp = SSres;
SStot = nan(modNum,1);



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
    Stim_local = find(duration >= (Win+ResDelay));% pourquoi +ResDelay ici????
    NbStim_local = length(Stim_local);
    if NbStim_local<20
        sprintf('Only %d stims long enough to run the model: no model is run with window size %dms\n', NbStim_local, Win);
        break
    end
    NBPC = [1:9 10:5:(NbStim_local*0.8)]; % This is not a good solution since the R2A profile of most cells show stairs-like structure with increasing number of PC we need to apply a ridge regression after the PCA.
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
    
    % Calculate principal components of the spectros
    fprintf(1, 'Calculate PC of spectro\n');
    [COEFF,SCORE,latent,tsquare]=princomp(x,'econ');
    if FIG==1
        figure(4)
        subplot(2,1,1)
        plot(cumsum(latent/sum(latent)))
        xlabel('# PC')
        ylabel('% variance explained by the PC')
    end
    
    % Loop through increasing numbers of PC and find the optimal number of
    % PC to predict neural response based on spectro only
    pValue=0;
    jj=0;
    ff=0;
    R2A_temp=zeros(1,length(NBPC));
    ModelPredict_temp = cell(1,length(NBPC));
    PCSTRF_temp = cell(1,length(NBPC));
    COEFF_temp = cell(1,length(NBPC));
    LL_temp= zeros(1,length(NBPC));
    Pvalue_temp=zeros(1,length(NBPC));
    NEC_temp = zeros(1,length(NBPC));
    SSres_temp = zeros(1,length(NBPC));
    SSexp_temp = zeros(1,length(NBPC));
    SStot_temp = zeros(1,length(NBPC));
    
    while jj<length(NBPC) && pValue<0.05
        jj = jj + 1;
        nPC=NBPC(jj);
        % Fit the neural data with the PC of the spectro
        %% WARNING WE NEED TO Regress with a RIDGE here CODE IT MANUALLY!!!! THEN RERUN on ALL CELLS
        fprintf('Constructing models of the neural response with %d PC of spectro\n', nPC);
        ds=dataset();
        for ii=1:nPC
            ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
        end
        ds.y=y;

        mdl=LinearModel.fit(ds);    %Model with PC of spectro only
        R2A_temp(jj)=mdl.Rsquared.Adjusted;
        ModelPredict_temp{jj}=mdl.predict;
        LL_temp(jj) = mdl.LogLikelihood;
        NEC_temp(jj) = mdl.NumEstimatedCoefficients;
        tbl=anova(mdl,'summary');
        Pvalue_temp(jj)=tbl.pValue(2);
        SSres_temp(jj) = mdl.SSE;
        SSexp_temp(jj) = mdl.SSR;
        SStot_temp(jj) = mdl.SST;
        PCSTRF_temp{jj} = mdl.Coefficients.Estimate(2:end);
        COEFF_temp{jj} = COEFF;
        
        % Calculate the p-value of the lratiotest between acoustic models with
        % nPC PC and previous nPC value
        if jj>1
            dof = NEC_temp(1,jj) - NEC_temp(1,(jj-1));
            [h_temp,pValue,stat,cValue] = lratiotest(LL_temp(1,(jj)),LL_temp(1,jj-1),dof);
            if FIG==1
                sprintf('the loglikelihood test between the acoustic models calculated with %dPC and %dPC gives a pValue of:%f\n', nPC, NBPC(jj-1), pValue)
                sprintf('the adjusted R2 are %f with %dPC and %f with %dPC\n', R2A_temp(jj-1), NBPC(jj-1), R2A_temp(jj), NBPC(jj))
            end
        end
    end
    
    if FIG==1
        figure(4)
        subplot(2,1,2)
        plot(NBPC, R2A_temp);
        %hold on
        ylabel('Adjusted R2 Acoustic Model')
        xlabel('# PC')
    end
    %calculate the STRF
    PCSTRF=PCSTRF_temp{jj-1};
    COEFF = COEFF_temp{jj-1};
    STRF=COEFF(:,1:NBPC(jj-1))*PCSTRF;
    STRFM=reshape(STRF,NFreq, NTime);
    LongestStim = find(duration==max(duration));
    Fo_Indices=find(Spectro.fo{LongestStim}<=Flow);
    To_Indices=find((1000.*Spectro.to{LongestStim})<=Win);
    if FIG==1
        ff = ff+1;
        fprintf(1, 'Calculating STRF using the %d first PC of the spectro\n\n\n\n', NBPC(jj-1));
        figure(ff)
        imagesc(Spectro.to{LongestStim}(To_Indices), Spectro.fo{LongestStim}(Fo_Indices), STRFM)
        axis xy
        title(sprintf('STRF with %d PC of the spectro', NBPC(jj-1)));
        pause
    end
    
    %Store data for the acoustic model
    NbOptPC(mm)=NBPC(jj-1);
    R2A.Acoustic(mm) = R2A_temp(jj-1);
    ModelPredict.Acoustic{mm} = ModelPredict_temp{jj-1};
    NeuralResponse{mm} = y;
    LL.Acoustic(mm) = LL_temp(jj-1);
    NEC.Acoustic(mm) = NEC_temp(jj-1);
    Pvalue.Acoustic(mm) = Pvalue_temp(jj-1);
    SSres.Acoustic(mm) = SSres_temp(jj-1);
    SSexp.Acoustic(mm) = SSexp_temp(jj-1);
    SStot(mm) = SStot_temp(jj-1);
    STRF_time{mm} = STRFM;
    STRF_to{mm}=Spectro.to{LongestStim}(To_Indices);
    STRF_fo{mm} = Spectro.fo{LongestStim}(Fo_Indices);
    
    % Now do the calculations for the semantic and AcSem models
   
    %Model with  VocType only
    ds2=dataset();
    ds2.Voctype=ordinal(VOC{mm});
    ds2.y=y;
    
    mdl2=LinearModel.fit(ds2);  
    R2A.Semantic(mm)=mdl2.Rsquared.Adjusted;
    ModelPredict.Semantic{mm}=mdl2.predict;
    LL.Semantic(mm) = mdl2.LogLikelihood;
    NEC.Semantic(mm) = mdl2.NumEstimatedCoefficients;
    tbl2=anova(mdl2,'summary');
    Pvalue.Semantic(mm)=tbl2.pValue(2);
    SSres.Semantic(mm) = mdl2.SSE;
    SSexp.Semantic(mm) = mdl2.SSR;
    
    %Model with both PC of spectro and VocType
    ds3=dataset();
    for ii=1:NbOptPC(mm)
            ds3.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
    end
    ds3.Voctype=ordinal(VOC{mm});
    ds3.y=y;

    mdl3=LinearModel.fit(ds3);  
    R2A.AcSem(mm)=mdl3.Rsquared.Adjusted;
    ModelPredict.AcSem{mm}=mdl3.predict;
    LL.AcSem(mm) = mdl3.LogLikelihood;
    NEC.AcSem(mm) = mdl3.NumEstimatedCoefficients;
    tbl3=anova(mdl3,'summary');
    Pvalue.AcSem(mm)=tbl3.pValue(2);
    SSres.AcSem(mm) = mdl3.SSE;
    SSexp.AcSem(mm) =  mdl3.SSR;
        
    % Plot of the predicted spike rate given the voctype
    CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
    MeanValues=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
    ModSem{mm}=MeanValues;
    if FIG==1
        ff = ff +1;
        figure(ff);
        plot(MeanValues, 1:length(MeanValues));
        title('Predicted SR with semantic Model');
        set(gca,'YTickLabel', unique(VOC{mm}));
        ff=ff+1;
        figure(ff)
        gscatter(y, mdl2.predict, VOC{mm}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
        ylabel('Predicted SR /ms with semantic model')
        xlabel('Observed SR /ms')
        pause
    end
    [HLRatio.AcAcSem(mm),PvalLRatio.AcAcSem(mm),stat,cValue] = lratiotest(LL.AcSem(mm),LL.Acoustic(mm),NEC.AcSem(mm) - NEC.Acoustic(mm));
    [HLRatio.SemAcSem(mm),PvalLRatio.SemAcSem(mm),stat,cValue] = lratiotest(LL.AcSem(mm),LL.Semantic(mm),NEC.AcSem(mm) - NEC.Semantic(mm)); 
end 

