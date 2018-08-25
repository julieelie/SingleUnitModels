
resultsDirectory='/auto/k8/julie';
%addpath(genpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/'));
addpath('/auto/k1/queued');
%addpath('/Users/elie/Documents/MATLAB/tlab/trunk/src/h5analysis/Julie_neuralcode/cubehelix')

%% Retrieve some INput variables
% cell containing the path of each unit
load('/auto/k8/julie/SemanticAnalysis.mat','List_matfilepath');
LM = length(List_matfilepath);
% delay over which you want to calculate Global R2 in ms
R2Delay = 300;

%% Create some output variables
% Vectors for Global R2A Models resluts 
Acoustic = zeros(LM,1);
Semantic = zeros(LM,1);
AcSem = zeros(LM,1);


%Matrices for R2A values
RAcoustic = zeros(LM,140);
RSemantic = RAcoustic;
RAcSem = RAcoustic;

% Matrices for ll tests
Pval_AS_A = RAcoustic;
Pval_AS_S = RAcoustic;

% Vectors for sum of squares of residuals (error) of each model and F
% statistic
SSAcoustic = zeros(LM,1);
SSSemantic = zeros(LM,1);
SSAcSem = zeros(LM,1);
SStot=zeros(LM,1);


%MAtrices for nb of stim and number of cat on which predicting models are
%calculated
RNb_Stim = RAcoustic;
RNb_Cat = RAcoustic;
RNb_PC = RAcoustic;

% Vector of the F value for global models comparison (Global R2A
% comparisions)
FAS_Ac = nan(LM,1);
FAS_Sem = nan(LM,1);
FAcoustic=nan(LM,1);
FSemantic=nan(LM,1);
FAcSem=nan(LM,1);

% Vector of the p-value for global models comparison (Global R2A
% comparisions)
PAS_Ac = nan(LM,1);
PAS_Sem = nan(LM,1);
PAc = nan(LM,1);
PSem = nan(LM,1);
PAcSem = nan(LM,1);

%% Start to loop through files and extract information
cd /auto/k6/julie/matfile
input_dir=pwd;
ii = 0;

for hh = 1:LM
    [Idir, Matfile, ext] = fileparts(List_matfilepath{hh});
    
    %retrieve Model files
    MDmatfiles=dir(fullfile(Idir, 'Models*.mat'));
    rm=length(MDmatfiles);
    SS_Ind=zeros(rm,1);
    for ff=1:rm
        if ~isempty(strfind(MDmatfiles(ff).name, 'ss'))
            SS_Ind(ff)=1;
        end
    end
    IndicesMD=find(SS_Ind);
    MD=length(IndicesMD);
     
            
        

 %% Upload the Model file for that unit and calculate global R2A
    sprintf('Upload Models results for that unit and caluclate Global R2A\n')
    for mm=1:MD
        MDMatfile = MDmatfiles(IndicesMD(mm)).name;
        if strcmp(Matfile(9:end), MDMatfile(8:end-4))
            MMAT = load(fullfile(Idir, MDMatfile));
            ii=ii+1;
            break
        end
    end

    %% Calculate global adjusted R2 again (wrongly calculated in
    % SpectroSemantic_Neuro_model as for 02.24.2014. Calculate them on the
    % first 300ms window.
    % Calculate p-value of global R2A comparisons based on F calculated here
    % (wrongly calculated in FindSemanticNeuronsFinal as of 02.24.2014

    

    if ii==1
        Wins=MMAT.Wins;
        RAcoustic = RAcoustic(:,1:length(Wins));
        RSemantic = RSemantic(:,1:length(Wins));
        RAcSem = RAcSem(:, 1:length(Wins));
        Pval_AS_A = Pval_AS_A(:, 1:length(Wins));
        Pval_AS_S = Pval_AS_S(:, 1:length(Wins));
        RNb_Stim = RNb_Stim(:,1:sum(Wins<=R2Delay));
        RNb_Cat = RNb_Cat(:,1:sum(Wins<=R2Delay));
        RNb_PC = RNb_PC(:,1:sum(Wins<=R2Delay));
    else
        if find(MMAT.Wins==R2Delay)~=IndWin
            fprintf('We have a problem of windows for the models of %s there are not the same as the previous unit\n', MDmatfile);
            break
        end

    end
    % Extract R2 of the three models for each time point
    RAcoustic(ii,:) = MMAT.R2A.Acoustic;
    RSemantic(ii,:) = MMAT.R2A.Semantic;
    RAcSem(ii,:) = MMAT.R2A.AcSem;
    
    % Extract the results of the log likelihood ratio test at each time
    % point
    Pval_AS_A(ii,:) = MMAT.PValLRatio.AcAcSem;
    Pval_AS_S(ii,:) = MMAT.PValLRatio.SemAcSem;
    
    %Index of window equal to 300 ms so we can calculate global
    %R2 within that window size
    IndWin = find(MMAT.Wins==R2Delay);

    % extract the Sum of squares Error (residuals) and total for
    % the 3 models within the zone of interest from 0 up to R2Delay ms
    SSAcoustic(ii) = nansum(MMAT.SSres.Acoustic(1:IndWin));
    SSSemantic(ii) = nansum(MMAT.SSres.Semantic(1:IndWin));
    SSAcSem(ii) = nansum(MMAT.SSres.AcSem(1:IndWin));
    SStot(ii) = nansum(MMAT.SStot(1:IndWin));

    % extract the number of categories and nb of stims on wich the
    % model are calculated for each time point
    for vv=1:IndWin
        RNb_Cat(ii,vv)=length(unique(MMAT.voc{vv}));
        RNb_Stim(ii,vv)=length(MMAT.voc{vv});
        RNb_PC(ii,vv)= MMAT.Best_nbPC(vv);
    end

    % Calculate the Global R2A for that unit for that period of time from 0
    % to R2Delay
    NT = sum(RNb_Stim(ii,1:IndWin)); % # observations = sum of the number fo stims over all the time points
     NDL = nansum(RNb_Cat(ii,1:IndWin)); % # parameters semantic models = sum of number of categories used in the semantic model over all time points
     NPC = nansum(RNb_PC(ii,1:IndWin)); % # parameters acoustic models = sum of all differents PC used over time
     NMOD = sum(~isnan(RNb_PC(ii, 1:IndWin))); % # of models (or time points) that run for that cell
     DFtot = NT - NMOD;%The null model is constructed with NMOD mean values, one mean for each model
     DFAc = NT - NPC - NMOD;%NMOD here represents the bo term (or bias term) of the regression for eah model, the term that take into account the average spike rate of the cell for these stims
     DFSem = NT - NDL - NMOD;
     DFAcSem = NT -  NPC - NDL - NMOD;
     
     Acoustic(ii) = 1- (SSAcoustic(ii)/DFAc)/(SStot(ii)/DFtot);
     Semantic(ii) = 1 - (SSSemantic(ii)/DFSem)/(SStot(ii)/DFtot);
     AcSem(ii) = 1 - (SSAcSem(ii)/DFAcSem)/(SStot(ii)/DFtot);
     
     % Calculate the F statistics for Global models for that unit for that period of time from 0
    % to R2Delay
     FAS_Ac(ii) = ((SSAcoustic(ii) - SSAcSem(ii))/(DFAc - DFAcSem))/(SSAcSem(ii)/DFAcSem);
     FAS_Sem(ii) = ((SSSemantic(ii) - SSAcSem(ii))/(DFSem - DFAcSem))/(SSAcSem(ii)/DFAcSem);
     FAcoustic(ii) = ((SStot(ii) - SSAcoustic(ii))/(DFtot - DFAc))/(SSAcoustic(ii)/DFAc);
     FSemantic(ii) = ((SStot(ii) - SSSemantic(ii))/(DFtot - DFSem))/(SSSemantic(ii)/DFSem);
     FAcSem(ii) = ((SStot(ii) - SSAcSem(ii))/(DFtot - DFAcSem))/(SSAcSem(ii)/DFAcSem);
     
     %Calculate p-values of global model comparisions based on F statictics
     PAS_Ac(ii) = 1-fcdf(FAS_Ac(ii), NDL, DFAcSem);
     PAS_Sem(ii) = 1-fcdf(FAS_Sem(ii),NPC, DFAcSem);
     PAc(ii) = 1-fcdf(FAcoustic(ii), NPC, DFAc);
     PSem(ii) = 1-fcdf(FSemantic(ii), NDL, DFSem);
     PAcSem(ii) = 1-fcdf(FAcSem(ii), NPC + NDL, DFAcSem);
end

%% Compile results into structures
sprintf('Compile results\n')
GlobalR2A.Acoustic = Acoustic(1:ii);
GlobalR2A.Semantic = Semantic(1:ii);
GlobalR2A.AcSem = AcSem(1:ii);

R2A.Acoustic = RAcoustic(1:ii,:);
R2A.Semantic = RSemantic(1:ii,:);
R2A.AcSem = RAcSem(1:ii,:);
R2A.Wins = Wins;
R2A.RNb_Cat = RNb_Cat(1:ii,:);
R2A.RNb_Stim = RNb_Stim(1:ii,:);
R2A.RNb_PC = RNb_PC(1:ii,:);
R2A.Pval_AS_A = Pval_AS_A(1:ii,:);
R2A.Pval_AS_S = Pval_AS_S(1:ii,:);
GlobalR2A.SSAcoustic = SSAcoustic(1:ii);
GlobalR2A.SSSemantic = SSSemantic(1:ii);
GlobalR2A.SSAcSem = SSAcSem(1:ii);
GlobalR2A.SStot = SStot(1:ii);
GlobalR2A.FAS_Ac = FAS_Ac(1:ii);
GlobalR2A.FAS_Sem = FAS_Sem(1:ii);
GlobalR2A.FAcoustic = FAcoustic(1:ii);
GlobalR2A.FSemantic = FSemantic(1:ii);
GlobalR2A.FAcSem = FAcSem(1:ii);
GlobalR2A.PAS_Ac = PAS_Ac(1:ii);
GlobalR2A.PAS_Sem = PAS_Sem(1:ii);
GlobalR2A.PAc = PAc(1:ii);
GlobalR2A.PSem = PSem(1:ii);
GlobalR2A.PAcSem = PAcSem(1:ii);


%% Save data
save(fullfile(resultsDirectory, 'SemanticAnalysis_Encod.mat'))
sprintf('Data saved under %s\n', fullfile(resultsDirectory, 'SemanticAnalysis_Encod.mat'))



