function [R2A, ModelPredict, LL, NEC, PValLRatio, h, y, voc, Best_nbPC, Pvalue]=Spectro_Neuro_model(MatfilePath)
FIG = 0;
%nPC=100;
NBPC=1:10:110; % to see how adjusted R2 evolve with NB of PC
%NBPC=38; %Mean Nb of PC that gives the best results of adjusted R2 as observed by running this code with the previous line with the script optimalPC_ModelNeuralResponse.m on 10/02/2012
R2A=zeros(3,1);
ModelPredict = cell(3,1);
LL= zeros(3,1);
Pvalue=zeros(3,1);              %results from the anova on each model
NEC = zeros(3,1);
PValLRatio = zeros(3,1);
h = zeros(3,1);
R2A_temp=zeros(3,length(NBPC));
ModelPredict_temp = cell(3,length(NBPC));
LL_temp= zeros(3,length(NBPC));
Pvalue_temp=zeros(3,length(NBPC));
NEC_temp = zeros(3,length(NBPC));

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%problem here because Wholevoc matfiles have stim that have different
%sizes.
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/FirstVocMat/FirstVoc_Site2_L1100R1450_e21_s0_ss1.mat')
Res = load(MatfilePath);

%% Need to get rid of mlnoise sections
DataSel=zeros(1,length(Res.VocType));
nvoc=0;
voctype=Res.VocType;
for dd=1:length(voctype);
    if strcmp(voctype{dd}, 'Ag')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Be')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'DC')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Di')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'LT')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Ne')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Te')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Th')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'song')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    elseif strcmp(voctype{dd}, 'Wh')
        nvoc=nvoc+1;
        DataSel(nvoc)=dd;
    end
end
DataSel=DataSel(1:nvoc);

%% Re-cut sections such that the sections begin with a sound
%% Loop through sections to find the optimal window of time (starting with the begining of each sound
% and ending as far as we can given the interstimulus interval
Spectros=Res.Spectro;
BegSoundInd = zeros(nvoc,1);
DurationSoundInd = zeros(nvoc,1);
DurationSectionInd = zeros(nvoc,1);
for ss = 1:nvoc
    sss = DataSel(ss);
    Spectro = reshape(Spectros{sss}, length(Res.Spectrofo), length(Res.Spectroto));
    TimeSpec = find(sum(Spectro,1)>0);
    BegSoundInd(ss) = TimeSpec(1);
    DurationSectionInd(ss) = size(Spectro,2)- TimeSpec(1)+1;
    DurationSoundInd(ss) = TimeSpec(end)- TimeSpec(1)+1;
end
DurationSpecInd = min(DurationSectionInd); % sounds should be analysed within this optimal window
DurationSpecms = DurationSpecInd*Res.Spectroto(end)*1000/length(Res.Spectroto);% size of this optimal window in ms

%Check the duration of the sections using the output  of Matfile construct
%WV to decide what should be the size of the window for spike rate
% find the new window in frequency to filter the sounds below 8000Hz (some
% sounds were low pass filtered at 8000Hz and some at 12000Hz)
FreqBelow8 = find(Res.Spectrofo<=8000);
EndFreq = FreqBelow8(end);
NFreq=length(FreqBelow8);

%% Loop through sections to isolate new psth, new MeanRate and new spectro
PSTH = nan(nvoc, ceil(DurationSpecms));
MeanRate = nan(nvoc,1);
NewSpectros = nan(nvoc, DurationSpecInd*NFreq);
for ss=1:nvoc
    dd = DataSel(ss);
    % new spectro
    Spectro = reshape(Spectros{dd}, length(Res.Spectrofo), length(Res.Spectroto));
    newspectro = Spectro(1:EndFreq,BegSoundInd(ss):(BegSoundInd(ss) + DurationSpecInd -1));%select the right time and frequency windows 
    NewSpectros(ss,:) = reshape(newspectro, 1,DurationSpecInd*NFreq);
    
    %new psth
    psth = Res.PSTH{dd};
    BegPSTHInd = floor(Res.Spectroto(BegSoundInd(ss))*1000);
    EndPSTHInd = ceil(Res.Spectroto(BegSoundInd(ss) + DurationSpecInd -1)*1000);
    PSTH(ss,:) = psth(BegPSTHInd:EndPSTHInd);
    
    %MeanRate
    Trials = Res.Trials{dd};
    SR = zeros(length(Trials),1);
    for tt = 1:length(Trials)
        SR(tt) = 1000*sum(Trials{tt}>BegPSTHInd && Trials{tt}<EndPSTHInd)/DurationSpecms; %calculate the spike rate per s
    end
    MeanRate(ss) = mean(SR);
end

%% Gettin ready for models

y= MeanRate;
x = cell2mat(Res.Spectro(DataSel));
voc = Res.VocType(DataSel);

x = 20*log10(abs(x)); % here we take the log of the spectro but pb: x
MAXI = max(max(x));
x(find(x<(MAXI-80)))=MAXI-80;
%contains 0 and log(x) gives inf and princomp does not support inf... :-(

%% Calculate the principal components of the spectro
fprintf(1, 'Calculate PC of spectro\n');
[COEFF,SCORE,latent,tsquare]=princomp(x,'econ');
if FIG==1
    figure(1)
    plot(cumsum(latent/sum(latent)))
end
pValue=0;
jj=0;
ff=0;
%while jj<length(NBPC)
while pValue<0.05 && jj<length(NBPC)
    jj = jj + 1;
    nPC=NBPC(jj);

    %% Fit the neural data with the PC of the spectro and/or the vocalization type and retrieve RMSE 
    %ds=dataset(SCORE(:,1:100),y);
    fprintf('Constructing three models of the neural response (spike rate) with the following variables:\n %d PC of spectro\nVocalization type\n%d PC of spectro + vocalization type\n', nPC);
    ds=dataset();

    for ii=1:nPC
        ds.(sprintf('SCORE%d',ii)) = SCORE(:,ii);
    end

    ds2=dataset();
    ds2.Voctype=ordinal(voc);
    ds2.y=y;
    
    
    
    %keep on going with other models
    ds3=ds;
    ds.y=y;
    ds3.Voctype=ordinal(voc);
    ds3.y=y;

    mdl=LinearModel.fit(ds);    %Model with PC of spectro only
    R2A_temp(1,jj)=mdl.Rsquared.Adjusted;
    ModelPredict_temp{1,jj}=mdl.predict;
    LL_temp(1,jj) = mdl.LogLikelihood;
    NEC_temp(1,jj) = mdl.NumEstimatedCoefficients;
    tbl=anova(mdl,'summary');
    Pvalue_temp(1,jj)=tbl.pValue(2);
    mdl2=LinearModel.fit(ds2);  %Model with  VocType only
    R2A_temp(2,jj)=mdl2.Rsquared.Adjusted;
    ModelPredict_temp{2,jj}=mdl2.predict;
    LL_temp(2,jj) = mdl2.LogLikelihood;
    NEC_temp(2,jj) = mdl2.NumEstimatedCoefficients;
    tbl2=anova(mdl2,'summary');
    Pvalue_temp(2,jj)=tbl2.pValue(2);
    mdl3=LinearModel.fit(ds3);  %Model with both PC of spectro and VocType
    R2A_temp(3,jj)=mdl3.Rsquared.Adjusted;
    ModelPredict_temp{3,jj}=mdl3.predict;
    LL_temp(3,jj) = mdl3.LogLikelihood;
    NEC_temp(3,jj) = mdl3.NumEstimatedCoefficients;
    tbl3=anova(mdl3,'summary');
    Pvalue_temp(3,jj)=tbl3.pValue(2);
    
    % Calculate the p-value of the lratiotest between acoustic models with
    % nPC PC and previous nPC value
    if jj>1
        dof = NEC_temp(1,jj) - NEC_temp(1,(jj-1));
        [h_temp,pValue,stat,cValue] = lratiotest(LL_temp(1,(jj)),LL_temp(1,jj-1),dof);
        sprintf('the loglikelihood test between the acoustic models calculated with %dPC and %dPC gives a pValue of:%f\n', nPC, NBPC(jj-1), pValue)
    end
    
    %calculate the STRF
    if FIG==1
        ff = ff+1;
        fprintf(1, 'Calculating STRF using the %d first PC of the spectro\n\n\n\n', nPC);
        PCSTRF=mdl.Coefficients.Estimate(2:end);
        STRF=COEFF(:,1:nPC)*PCSTRF;
        STRFM=reshape(STRF,length(Res.Spectrofo), length(Res.Spectroto));
        figure(ff)
        imagesc(Res.Spectroto, Res.Spectrofo, STRFM)
        axis xy
        title(sprintf('STRF with %d PC of the spectro', nPC));
        pause
    end
    
    % Plot of the predicted spike rate given the voctype
    if FIG==2
        ff = ff +1;
        CoeffEstimates=mdl2.Coefficients.Estimate(1:end);
        MeanValues=CoeffEstimates + [0 ; repmat(CoeffEstimates(1), (length(CoeffEstimates)-1),1)];
        plot(MeanValues, 1:length(MeanValues));
        title('Predicted SR with semantic Model');
        pause
    end
    
    
end
Best_nbPC=NBPC(jj-1); %Best number of PC calculated on the acoustic model and use for this unit for the calculations of the 3 models
fprintf('best nb of PC= %d\n', Best_nbPC);
R2A(1:3,1) = R2A_temp(1:3,(jj-1));
ModelPredict{1,1} = ModelPredict_temp{1,(jj-1)};
ModelPredict{2,1} = ModelPredict_temp{2,(jj-1)};
ModelPredict{3,1} = ModelPredict_temp{3,(jj-1)};
LL(1:3,1) = LL_temp(1:3,(jj-1));
NEC(1:3,1) =NEC_temp(1:3,(jj-1));
Pvalue(1:3,1)=Pvalue_temp(1:3,(jj-1));
if NEC(1,1)>NEC(2,1)
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(1,1),LL(2,1),NEC(1,1) - NEC(2,1));
elseif NEC(1,1)<NEC(2,1)
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(2,1),LL(1,1),NEC(2,1) - NEC(1,1));
else
    fprintf('same degree of freedom for Lratiotest between Acoustic and semantic');
    [h(1),PValLRatio(1,1),stat,cValue] = lratiotest(LL(1,1),LL(2,1),NEC(1,1) - NEC(2,1));
end
[h(2),PValLRatio(2,1),stat,cValue] = lratiotest(LL(3,1),LL(1,1),NEC(3,1) - NEC(1,1));
[h(3),PValLRatio(3,1),stat,cValue] = lratiotest(LL(3,1),LL(2,1),NEC(3,1) - NEC(2,1));

%fprintf('the RMSE of the models are:\n%d first PC of spectro: %f\nVocType only: %f\n%d first PC of spectro + VocType: %f\n', nPC, RMSE(1), RMSE(2), nPC, RMSE(3));
if FIG==2
    figure(12)
    plot(NBPC, R2A_temp(1,:), 'rs-', NBPC, R2A_temp(2,:), 'co-', NBPC, R2A_temp(3,:), 'g*-');
    axis([0 100 0 1])
end
end

