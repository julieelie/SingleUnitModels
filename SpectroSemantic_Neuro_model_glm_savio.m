function [OptimalCoherenceWinsize] = SpectroSemantic_Neuro_model_glm_savio(MatfilePath, SWITCH, ParamModel,Cellname)
%% Get the environment to figure out on which machine/cluster we are
getenv('HOSTNAME')
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))%savio Cluster
    Savio=1;
    addpath(genpath('/global/home/users/jelie/CODE'));
elseif ismac()
    Savio=0;
    Me = 1;
    addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'));
    addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'));
    addpath(genpath('/Users/elie/Documents/CODE/STRFLab/trunk'));
else %we are on strfinator or a cluster machine
    Savio = 0;
    Me = 0;
    addpath(genpath('/auto/fhome/julie/Code/SingleUnitModels'));
    addpath(genpath('/auto/fhome/julie/Code/GeneralCode'));
    addpath(genpath('/auto/fhome/julie/Code/strflab'));
end

%% Start a timer for the function
TimerVal=tic;

%% Deal with input parameters
if nargin<2
    SWITCH = struct();
end
if ~isfield(SWITCH,'FanoFactor') || isempty(ParamModel.FanoFactor)
    SWITCH.FanoFactor=0;
end
if ~isfield(SWITCH,'BestBin') || isempty(ParamModel.BestBin)
    SWITCH.BestBin=1;
end
if ~isfield(SWITCH,'Models') || isempty(ParamModel.Models)
    SWITCH.Models=0;
end
if ~isfield(SWITCH,'InfoCal') || isempty(ParamModel.InfoCal)
    SWITCH.InfoCal=0;%Set to 1 if you want to calculate information on spike trains and change the name of the output file so they indicate "Info"
end

if nargin<3
    ParamModel = struct();
end
if ~isfield(ParamModel,'LINK') || isempty(ParamModel.LINK)
    ParamModel.LINK='log'; %'identity'
end
if ~isfield(ParamModel,'DISTR') || isempty(ParamModel.DISTR)
    ParamModel.DISTR='poisson';%'normal'
end
if ~isfield(ParamModel,'NeuroRes') || isempty(ParamModel.NeuroRes)
    ParamModel.NeuroRes = 'count_gaussfiltered';
end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 10; % end point of the first analysis window (spectrogram and neural response)
end
if ~isfield(ParamModel,'MaxWin') || isempty(ParamModel.MaxWin)
    ParamModel.MaxWin = 1000; %end point of the last anaysis window for...
    ... neural response and end point of the largest analysis window for...
        ... spectrogram
end
if ~isfield(ParamModel,'Increment') || isempty(ParamModel.Increment)
    ParamModel.Increment = 10; %increase the size of the spectro window with a Xms pace
end
if ~isfield(ParamModel,'NeuroBin') || isempty(ParamModel.NeuroBin)
    ParamModel.NeuroBin = 10; % size of the window (ms) within which the neural response is analyzed
                               % The end of the window of analysis is
                               % determined by the Increment and ResDelay (see below).
end
if ~isfield(ParamModel,'ResDelay') || isempty(ParamModel.ResDelay)
    ParamModel.ResDelay = 0; % Delay in ms between the end of the...
    ... spectrogram window and the end of the neural response window
end


if nargin<4
    [path,Cellname,ext]=fileparts(MatfilePath);
end

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
if Savio %savio Cluster
    Dir_local='/global/scratch/jelie/MatFiles/';
    Res=loadfromTdrive_savio(MatfilePath, Dir_local);
elseif Me
    Dir_local='/Users/elie/Documents/CODE/data/matfile/FirstVoc1sMat/';
    if ~exist('ext','var')
        [path,Cellname,ext]=fileparts(MatfilePath);
    end
    Res = load([Dir_local Cellname ext]);
else
    Res = load(MatfilePath);
end



%% Get ready saving files and directories
OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');

if Savio
    OutputDir_local='/global/scratch/jelie/MatFiles/ModMat';
    OutputDirEx_local='/global/home/users/jelie/MatFiles/ModMat';
    OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatSavio');
elseif Me
    if SWITCH.InfoCal || SWITCH.BestBin
        OutputDir_local='/users/elie/Documents/CODE/data/matfile/ModMatInfo';
    else
        OutputDir_local='/users/elie/Documents/CODE/data/matfile/ModMatAcOnly';
    end
    OutputDir_final=OutputDir_local;
else
    if SWITCH.InfoCal || SWITCH.BestBin
        OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');
    else
        OutputDir_final=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatAcOnly');
    end
    OutputDir_local=OutputDir_final;
end
if SWITCH.InfoCal || SWITCH.BestBin
    calfilename_local=fullfile(OutputDir_local,['InfoPoissonGF_' Res.Site '.mat']);
    calfilename_final=fullfile(OutputDir_final,['InfoPoissonGF_' Res.Site '.mat']);
else
    calfilename_local=fullfile(OutputDir_local,['Models_GLMPoisson_' Res.Site '.mat']);
    calfilename_final=fullfile(OutputDir_final,['Models_GLMPoisson_' Res.Site '.mat']);
end
outfilename_local=fullfile(OutputDir_local,['slurm_out*' Res.Site '*.txt']);
outfilename_final=fullfile(OutputDir_final,['slurm_out*' Res.Site '*.txt']);

PrevData=0;
if Savio
    try
        DoneCalc=loadfromTdrive_savio(calfilename_final, OutputDir_local,1);
        fprintf('Found some data for this unit\n')
            if isfield(DoneCalc, 'Deviance') && isfield(DoneCalc, 'LL') && isfield(DoneCalc, 'LambdaChoice') && isfield(DoneCalc, 'Model') && isfield(DoneCalc, 'PropVal') && isfield(DoneCalc, 'Data') && isfield(DoneCalc, 'Wins') && ~isempty(DoneCalc.Model.MeanSpectroStim{1})
                PrevData = 1;
            else
                frpintf('Data are not complete enough to be used\n')
                system(['rm ' calfilename_local]);
                clear 'DoneCalc'
                PrevData = 0;
            end

        
    catch err
        fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
        PrevData = 0;
    end
%         if strcmp(err.identifier, 'MATLAB:load:couldNotReadFile')
%             fprintf('No previous Data working from the first window\n');
%             PrevData = 0;
%         elseif strcmp(err.identifier, 'MATLAB:load:cantReadFile')
%             fprintf('Previous Data File corrupted working from the first window\n');
%             PrevData = 0;
%         else
%             fprintf('Error loading previous Data: %s\nworking from the first window\n',err.identifier);
%             PrevData = 0;
%         end
else
    try
        DoneCalc=load(calfilename_final);
        fprintf('Found some data for this unit\n')
            if isfield(DoneCalc, 'Deviance') && isfield(DoneCalc, 'LL') && isfield(DoneCalc, 'LambdaChoice') && isfield(DoneCalc, 'Model') && isfield(DoneCalc, 'PropVal') && isfield(DoneCalc, 'Data') && isfield(DoneCalc, 'Wins') && ~isempty(DoneCalc.Data.MeanSpectroStim{1})
                PrevData = 1;
            else
                fprintf('Data are not complete enough to be used\n')
                system(['rm ' calfilename_local]);
                clear 'DoneCalc'
                PrevData = 0;
            end

        
    catch err
        fprintf('No Previous Data available or complete, working from the first window\nThe error is %s\n',err.identifier, err.message);
        PrevData = 0;
    end
end
    
% save the path for now if no previous file
if ~ PrevData
    save(calfilename_local,'MatfilePath', '-append')
end

%% Get the data ready
if SWITCH.Models % For models we use vocalization sections, only the first element of each vocalization sequence
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
else
    DataSel = 1:length(Res.VocType);
end


%% Select the spectrograms of the selected stims
if SWITCH.Models
    Spectro.spec = Res.Spectro(DataSel);
    Spectro.to = Res.Spectroto(DataSel);
    Spectro.fo = Res.Spectrofo(DataSel);

    if ~PrevData
        % save the Stimuli Spectrograms for now if not done previously
        save(calfilename_local,'Spectro')
    end
end

%% Extract Emitter ID
if SWITCH.Models
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
end

%% Compute coherence to determine the optimal window size, the scale at which neural response information is maximized
if SWITCH.BestBin
    
    %Data processing
    [HalfTrain1, HalfTrain2, NumTrials]=organiz_data4coherence(Res.Trials_GaussFiltered(DataSel),Res.PSTH_GaussFiltered(DataSel),ParamModel);
    % compute coherence
    SampleRate=1000; %bin size =1ms so sample Rate = 1000Hz
    [CoherenceStruct]=compute_coherence_mean(HalfTrain1, HalfTrain2,SampleRate);
    OptimalCoherenceWinsize = CoherenceStruct.freqCutoff;
    % According to this code 20ms is the best size for a majority of cells see
    % fig BestPSTHBin.fig
    % At that window size, the values of the FanoFactor over cells is very
    % close to 1. see fig PoissonFanoFactor.fig
    save(calfilename_local,'MatfilePath', 'OptimalCoherenceWinsize','-append');
end

%% Estimate Poisson assumption for data at the choosen bining
if SWITCH.FanoFactor
    FanoFactor_mean=nan(length(Bins),1);
    for MW=1:length(Bins)
        [NeuroRes, PG_Index,FanoFactor_Index, Wins] = PoissonGaussianNeuralResponses(Spectro, Res.VocType(DataSel), Res.PSTH(DataSel), Res.Trials(DataSel),Cellname,Bins(MW), MaxWin, Increment, ResDelay);
        % According to PoissonGaussianNeuralResponses, neural responses are more
        % poisson than gaussian.
        FanoFactor_mean(MW) = mean(FanoFactor_Index);
    end
    save(calfilename_local,'MatfilePath', 'FanoFactor_mean','-append');
end

%return
%% Inject the data in the models
if SWITCH.Models
    if PrevData
        [LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod( Spectro, Res.VocType(DataSel), Res.PSTH_GaussFiltered(DataSel),Res.Trials_GaussFiltered(DataSel),Emitter, ParamModel,calfilename_local, DoneCalc);
    else
        [LambdaChoice, Deviance, LL, Model, ParamModel, Data, PropVal, Wins] = GrowingModelsRidgeglmLambdaMod( Spectro, Res.VocType(DataSel), Res.PSTH_GaussFiltered(DataSel),Res.Trials_GaussFiltered(DataSel),Emitter, ParamModel,calfilename_local);
    end
ElapsedTime = toc(TimerVal);
Days = floor(ElapsedTime/(60*60*24));
ETRem = ElapsedTime - Days*60*60*24;
Hours = floor(ETRem/(60*60));
ETRem = ETRem - Hours*60*60;
Minutes = floor(ETRem/60);
ETRem = ETRem-60*Minutes;
fprintf(1,'Corrected Code run for %d days %dh %dmin and %dsec\n',Days,Hours,Minutes,ETRem);
fprintf(1,'Good calculation of AcSemAc\n');
fprintf(1,'Threshold for coordinate descent corrected: relying on L2Norm of the parameters''vector\n');
save(calfilename_local,'MatfilePath', 'LambdaChoice', 'Deviance','LL','Model','ParamModel','Data','PropVal','Wins','ElapsedTime','-append');
end

if Savio
    [Status1]=transfertoTdrive_savio(calfilename_local,calfilename_final);
    [Status2]=transfertoTdrive_savio(outfilename_local,[OutputDir_final '/']);
    if ~(Status1 || Status2)
        system(['mv ' OutputDirEx_local '/JobToDoSavio/Ex*' Res.Site '*.txt ' OutputDirEx_local '/JobDoneSavio/'])
    end
end

fprintf(1,'Ready to quit');
quit

end

