function [] = SpectroSemantic_Neuro_model_glm(MatfilePath, Cellname,DISTR, LINK,MinWin, MaxWin, Increment, ResDelay)
TimerVal=tic;
% restaure if SWITCH set to 1: function
%[FanoFactor_mean,MinWin,MaxWin,Increment,ResDelay,OptimalCoherenceWinsize]
%= SpectroSemantic_Neuro_model_glm(MatfilePath, Cellname,DISTR, LINK,MinWin, MaxWin, Increment, ResDelay)
SWITCH.FanoFactor=0;
SWITCH.BestBin=0;
if nargin<8
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<7
    Increment = 10; %increase the size of the spectro window with a 10ms pace
end
if nargin<6
    MaxWin = 150; %maximum values the window of analysis can reach
end
if nargin<5
    MinWin = 20; %minimum size of the window of analysis from the begining Choose 20ms according to coherence cutoffFrequency calculation
end

if nargin<4
    LINK = 'log';
end

if nargin<3
    DISTR = 'poisson';
end

if nargin<2
    [path,Cellname,ext]=fileparts(MatfilePath);
end

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
Res = load(MatfilePath);

%% Get the data ready
% Select first sections
Firsts = find(Res.Voc_orders == 1);
% Need to get rid of mlnoise sections and whine sections when they
% exist. I construct a vector of indices of the right sections
DataSel=zeros(1,length(Firsts));
nvoc=0;
voctype=Res.VocType;
for ii=1:length(Firsts);
    dd = Firsts(ii);
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
%     elseif strcmp(voctype{dd}, 'Wh')
%         nvoc=nvoc+1;
%         DataSel(nvoc)=dd;
    end
end
DataSel=DataSel(1:nvoc);

%% Select the spectrograms of the selected stims
Spectro.spec = Res.Spectro(DataSel);
Spectro.to = Res.Spectroto(DataSel);
Spectro.fo = Res.Spectrofo(DataSel);

%% Extract Emitter ID
Ename = cell(length(DataSel),1);
Esex = cell(length(DataSel),1);
Eage = cell(length(DataSel),1);
Erelated = cell(length(DataSel),1);
for ii=1:length(DataSel)
    dd=DataSel(ii);
    [Path,File,Ext] = fileparts(Res.Original_wavfiles{dd});
    Ename{ii} = File(1:10);
    Esex{ii} = File(12);
    Eage{ii} = File(13);
    Erelated{ii} = File(14);
end
Emitter.Ename = Ename;
Emitter.Esex = Esex;
Emitter.Eage = Eage;
Emitter.Erelated = Erelated;

%% Estimate Poisson assumption for data
if SWITCH.FanoFactor
    FanoFactor_mean=nan(length(MinWin),1);
    for MW=1:length(MinWin)
        [NeuroRes, PG_Index,FanoFactor_Index, Wins] = PoissonGaussianNeuralResponses(Spectro, Res.VocType(DataSel), Res.PSTH(DataSel), Res.Trials(DataSel),Cellname,MinWin(MW), MaxWin, Increment, ResDelay);
        % According to PoissonGaussianNeuralResponses, neural responses are more
        % poisson than gaussian.
        FanoFactor_mean(MW) = mean(FanoFactor_Index);
    end
end

%% Compute coherence to determine the optimal window to calculate psth
if SWITCH.BestBin
    %Data processing
    [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Res.Trials(DataSel),Spectro,MaxWin,ResDelay);
    % compute coherence
    SampleRate=1000; %bin size =1ms so sample Rate = 1000Hz
    [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,SampleRate);
    OptimalCoherenceWinsize = CoherenceStruct.freqCutoff;
    % According to this code 20ms is the best size for a majority of cells see
    % fig BestPSTHBin.fig
    % At that window size, the values of the FanoFactor over cells is very
    % close to 1. see fig PoissonFanoFactor.fig
end

%return
%% Inject the data in the models

[LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod(Spectro, Res.VocType(DataSel), Res.PSTH(DataSel),Res.Trials(DataSel),Emitter, MinWin, MaxWin, Increment, ResDelay,'count',DISTR,LINK);

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile','ModMat',['Models_GLMPoisson' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile','ModMat',['Models_GLMPoisson' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMat',['Models_GLMPoisson' Res.Site '.mat']);
end
ElapsedTime = toc(TimerVal);
Days = floor(ElapsedTime/(60*60*24));
ETRem = ElapsedTime - Days*60*60*24;
Hours = floor(ETRem/(60*60));
ETRem = ETRem - Hours*60*60;
Minutes = floor(ETRem/60);
ETRem = ETRem-60*Minutes;
fprintf(1,'Corrected Code run for %d days %dh %dmin and %dsec\n',Days,Hours,Minutes,ETRem);
fprintf(1,'Good calculation of AcSemAc');
save(calfilename,'MatfilePath', 'LambdaChoice', 'Spectro', 'Deviance','LL','Model','ParamModel','Data','PropVal','Wins','ElapsedTime');
end

