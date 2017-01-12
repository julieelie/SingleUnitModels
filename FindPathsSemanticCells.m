addpath(genpath('/Users/elie/Documents/CODE/SingleUnitModels'))
addpath(genpath('/Users/elie/Documents/CODE/GeneralCode'))
load('/Users/elie/Documents/CODE/data/matfile/SemanticAnalysis_ShuffPredict.mat', 'SemanticIndex')
load('/Users/elie/Documents/CODE/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance3.mat', 'List_matfilepath', 'List_anat', 'ZONES', 'ZONES_List', 'SUBJ', 'Spike_shape')
SU=find(Spike_shape>0);
NU = length(List_matfilepath);
SemCellPV = find(SemanticIndex.quantileDiagUniCatMaxInv<0.01);
SemSU = intersect(SU, SemCellPV);

% % ID Cell selective for DC in EJN 2015
DC_cell=418; % other /auto/k6/julie/matfile/GreBlu9508M/ConfMat_Site4_L1500R1900_e23_s0_ss2.mat
TheFileID=List_matfilepath(DC_cell);
[P,TheFile,ext]=fileparts(TheFileID{1});
Matfilepath= fullfile('/Users/elie/Documents/CODE/data/matfile/WholeVocMat',['WholeVoc' TheFile(8:end) ext]);
SpectroSemantic_Neuro_model_glm(Matfilepath,'DC128','poisson','log',110,110)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'DC128','poisson','identity',110,110)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'DC128','normal','identity',110,110)
% 
% ID Cell selective for Ag in EJN 2015
Ag_cell=127;%/auto/k6/julie/matfile/BlaBro09xxF/ConfMat_Site3_L2500R2300_e22_s1_ss1.mat
TheFileID=List_matfilepath(Ag_cell);
[P,TheFile,ext]=fileparts(TheFileID{1});
Matfilepath= fullfile('/Users/elie/Documents/CODE/data/matfile/WholeVocMat',['WholeVoc' TheFile(8:end) ext]);
SpectroSemantic_Neuro_model_glm(Matfilepath,'Ag127','poisson','log',110,110)
% %SpectroSemantic_Neuro_model_glm(Matfilepath,'Ag127','poisson','identity',110,110)
% %SpectroSemantic_Neuro_model_glm(Matfilepath,'Ag127','normal','identity',110,110)

%Synergistic cell Poster Models STRF at 125ms
Matfilepath= fullfile('/Users/elie/Documents/CODE/data/matfile/WholeVocMat','WholeVoc_Site2_L1000R900_e13_s0_ss1.mat');
SpectroSemantic_Neuro_model_glm(Matfilepath,'SynCell','poisson','log',125,125)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'SynCell','poisson','identity',125,125)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'SynCell','normal','identity',125,125)

% % Acoustic cell Poster Models STRF at 105ms
Matfilepath= fullfile('/Users/elie/Documents/CODE/data/matfile/WholeVocMat','WholeVoc_Site2_L1100R1450_e14_s0_ss1.mat');
SpectroSemantic_Neuro_model_glm(Matfilepath,'ACCell','poisson','log',105,105)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'ACCell','poisson','identity',105,105)
%SpectroSemantic_Neuro_model_glm(Matfilepath,'ACCell','normal','identity',105,105)

% %AS cell poster Models STRF at 65ms
% Matfilepath= fullfile('/Users/elie/Documents/CODE/data/matfile/WholeVocMat','WholeVoc_Site3_L1200R1200_e20_s0_ss3.mat');
% SpectroSemantic_Neuro_model_glm(Matfilepath,'ASCell','poisson','identity',65,65)
% SpectroSemantic_Neuro_model_glm(Matfilepath,'ASCell','normal','identity',65,65)


%% SLURM
addpath('/auto/k1/queued/')
addpath('/auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode')
resultsDirectory='/auto/tdrive/julie/NeuralData/SemanticGLMModel';
MatfileDC128='/auto/tdrive/julie/k6/julie/matfile/GreBlu9508M/WholeVoc_Site4_L1500R1900_e23_s0_ss2.mat';
MatfileAg127='/auto/tdrive/julie/k6/julie/matfile/BlaBro09xxF/WholeVoc_Site3_L2500R2300_e22_s1_ss1.mat';
MatfileSYN='/auto/tdrive/julie/k6/julie/matfile/YelBlu6903F/WholeVoc_Site2_L1000R900_e13_s0_ss1.mat';
MatfileAc='/auto/tdrive/julie/k6/julie/matfile/GreBlu9508M/WholeVoc_Site2_L1100R1450_e14_s0_ss1.mat';
ListMatSlurm={MatfileDC128,MatfileAg127,MatfileSYN,MatfileAc};
cmd = 'SpectroSemantic_Neuro_model_glm(''%s'');';
for ff=1:length(ListMatSlurm)
    Matfile_local=ListMatSlurm{ff};
    [path,file,ext]=fileparts(Matfile_local);
    MATName=[file ext];
    jobParams = struct;
    jobParams.partition = 'all';
    jobParams.cpus = 2;
    jobParams.memory = 8000;
    jobParams.out = fullfile(resultsDirectory,sprintf('slurm_out_%s.txt', MATName));
    jobParams.err = jobParams.out;
    icmd = sprintf(cmd, Matfile_local);
    fprintf('Calling slurm_sbatch with command %s\n',icmd);
    slurm_sbatch(icmd,jobParams);
end
