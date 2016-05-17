%% Set up the paths, retrieve list of cells already done, list of all cells
addpath('/auto/k1/queued/')
addpath('/auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode')
addpath('/auto/fhome/julie/matlab/tlab/src/slurmbot/matlab')

cd /auto/tdrive/julie/k6/julie/matfile/ModMatSavio/
DoneFile=dir('Models_GLMPoisson*');

%system('ssh TheunissenLab')
% cd /auto/tdrive/julie/k6/julie/matfile/ModMat
% matlab
% DoneFile=dir('Models_GLMPoisson*');
% save('/auto/tdrive/julie/k6/julie/matfile/DoneFile.mat','DoneFile')
% quit
% logout
% system('scp TheunissenLab:/auto/tdrive/julie/k6/julie/matfile/DoneFile.mat /global/home/users/jelie/MatFiles/DoneFile.mat')
% load('/global/home/users/jelie/MatFiles/DoneFile.mat');

load('/auto/tdrive/julie/NeuralData/SemanticGLMModel/FanoFactor_CoherenceOptPSTHBin_SemCell.mat','List_SemanticCellspath');

%% Set up the variables for slurm
JobParams = struct;
JobParams.Name = '3PSTHModels';
JobParams.Partition = 'savio';
JobParams.Account = 'fc_birdpow';
JobParams.Qos = 'savio_normal';
JobParams.NTasks = 1;
JobParams.CPU = 20;
SlurmParams.cmd = 'SpectroSemantic_Neuro_model_glm_savio(''%s'');';
SlurmParams.resultsDirectory='/global/scratch/jelie/MatFiles/ModMat';

%% Set up variables to identify cells left to run and order
MeanRateCellToDo = nan(length(List_SemanticCellspath),1);
CTD=0; %Counter for cell to do
MatfileToDo = cell(length(List_SemanticCellspath),1);
MatNameToDo = cell(length(List_SemanticCellspath),1);

%% Identify cells left to run and their average spike rate
for ff=1:length(List_SemanticCellspath)
    fprintf(1,'checking file %d/%d\n',ff,length(List_SemanticCellspath));
    [P,TheFile,ext]=fileparts(List_SemanticCellspath{ff});
    ThisFiledone=0;
    for df=1:length(DoneFile)
        DataFilename = ['Models_GLMPoisson_' TheFile(9:end) ext];
        if strcmp(DoneFile(df).name,DataFilename)
            try
                Res=load(DataFilename);
                % Check if all windows were calculated
                UnrunWindows = find(isnan(Res.Deviance.Acoustic.FitIndex));
                LastWinRun=~isempty(Res.Model.MeanSpectroStim{end});
                if isempty(UnrunWindows) && LastWinRun
                    fprintf('All windows calculated\n')
                    %all windows were calculated
                    ThisFiledone=1;
                    % move this file to the done folder
                    [status, ~] = system(sprintf('mv %s %s', DataFilename,fullfile('DoneCells', DataFilename)), '-echo')
                    [status, ~] = system(sprintf('mv slurm_out_WholeVoc_%s*.txt DoneCells/', DataFilename(19:end-4)), '-echo')
                    break
                else
                    fprintf('There are still some windows to run\n')
                    ThisFiledone=0;
                end
            catch
                fprintf('cannot open properly that file\n')
                ThisFiledone=0;
            end
        end
    end
    if ThisFiledone
        fprintf('File Done don''t recalculate it\n');
        % Don't recalculate that file
    else
        fprintf('Adding this file to the list to do\n');
        CTD=CTD+1;
        MatfileToDo{CTD}= fullfile(P,['WholeVoc' TheFile(8:end) ext]);
        
        % Calculate the mean SpikeRate for that file
        Res=load(MatfileToDo{CTD});
        %Select the stims we are working with
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
        MeanRateCellToDo(CTD) = mean(cell2mat(Res.MeanRate(DataSel)));
        MatNameToDo{CTD}=['WholeVoc' TheFile(8:end) ext];
    end
end
MeanRateCellToDo=MeanRateCellToDo(1:CTD);
MatNameToDo = MatNameToDo(1:CTD);
MatfileToDo = MatfileToDo(1:CTD);

% Find out the order in which jobs should be run
[~,IndOrd_asc]=sort(MeanRateCellToDo);
IndOrd=flip(IndOrd_asc);

%% Create Jobs'files
cd /auto/tdrive/julie/k6/julie/matfile/ModMatSavio/JobToDoSavio
for jj=1:length(IndOrd)
    IndexFile = IndOrd(jj);
    JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%%j.txt', MatNameToDo{IndexFile}));
    JobParams.err = JobParams.out;
    icmd = sprintf(SlurmParams.cmd, MatfileToDo{IndexFile});
    fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
    slurm_sbatch_savio(icmd,JobParams);
end
fprintf(1,'DONE Creating all Jobs'' files!!!\n');

%% Find the jobs of the interesting cells
IndicesHipCells = nan(4,1);
IndicesHipCells(1) = IndOrd(1); %Highest spike rate cell Site2_L1100R1450_e14_s0_ss1 %Ac Cell
IndicesHipCells(2) = IndOrd(end); % lowest spike rate cell
for ff=1:length(MatNameToDo)
%     if strcmp(MatNameToDo{ff}(10:end-4), 'Site2_L1000R900_e13_s0_ss1') %SYN Cell cannot be found in the dataset of single semantic units...
%         IndicesHipCells(5)=ff;
    if strcmp(MatNameToDo{ff}(10:end-4), 'Site3_L2500R2300_e22_s1_ss1')%Ag127
        IndicesHipCells(3)=ff;
    elseif strcmp(MatNameToDo{ff}(10:end-4), 'Site4_L1500R1900_e23_s0_ss2') %DC128
        IndicesHipCells(4)=ff;
    end
end



for jj=1:length(IndicesHipCells)
    IndexFile=IndicesHipCells(jj)
    JobParams.out = fullfile(SlurmParams.resultsDirectory,sprintf('slurm_out_%s_%%j.txt', MatNameToDo{IndexFile}));
    JobParams.err = JobParams.out;
    icmd = sprintf(SlurmParams.cmd, MatfileToDo{IndexFile});
    fprintf(1,'creating file slurm_sbatch with command %s\n',icmd);
    slurm_sbatch_savio(icmd,JobParams);
end
fprintf(1,'DONE Creating all Jobs'' files!!!\n');

   