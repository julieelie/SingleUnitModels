resultsDirectory='/auto/tdrive/julie/NeuralData/SemanticGLMModel';
addpath(genpath('/auto/fhome/julie/matlab/tlab'))
addpath(genpath('/auto/fhome/julie/matlab/STRFLab'))
load('/auto/tdrive/julie/k8/SemanticAnalysis_ShuffPredict.mat', 'SemanticIndex')
load('/auto/tdrive/julie/k8/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat', 'List_matfilepath', 'Spike_shape')
SU=find(Spike_shape>0);
SemCellPV = find(SemanticIndex.quantileDiagUniCatMaxInv<0.01);
SemSU = intersect(SU, SemCellPV);
NSU=length(SemSU);
MinWin = 10:5:50;
FanoFactor_allCells = nan(NSU,length(MinWin));
OptimalCoherenceWinsize = nan(NSU,1);
List_SemanticCellspath=cell(NSU,1);
for ss=1:NSU
    fprintf('File %d/%d\n',ss,NSU);
    List_SemanticCellspath{ss}=List_matfilepath{SemSU(ss)};
    [P,TheFile,ext]=fileparts(List_SemanticCellspath{ss});
    Matfilepath= fullfile(P,['WholeVoc' TheFile(8:end) ext]);
    [FanoFactor_allCells(ss,:),MinWin,MaxWin,Increment,ResDelay,OptimalCoherenceWinsize(ss)] = SpectroSemantic_Neuro_model_glm(Matfilepath);
end
filename='FanoFactor_CoherenceOptPSTHBin_SemCell.mat';
Windows_param.MinWin=MinWin;
Windows_param.MaxWin=MaxWin;
Windows_param.Increment=Increment;
Windows_param.ResDelay=ResDelay;
save(fullfile(resultsDirectory,filename),'FanoFactor_allCells','Windows_param','OptimalCoherenceWinsize','List_SemanticCellspath')
exit
