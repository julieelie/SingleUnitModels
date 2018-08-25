%% Set up the paths, retrieve list of cells already done, list of all cells
load('/auto/tdrive/julie/NeuralData/SemanticGLMModel/FanoFactor_CoherenceOptPSTHBin_SemCell.mat','List_SemanticCellspath');
addpath('/auto/k1/queued/')
addpath('/auto/fhome/julie/matlab/tlab/src/h5analysis/Julie_neuralcode')
cd /auto/tdrive/julie/k6/julie/matfile/ModMatInfo
DoneFile=dir('Models_InfoPoisson*');

%% Set up the variables for slurm
global JobTimer
NJobs=10;
PeriodTimerCheck=60*60;% check every hour 60*60
JobParams = struct;
JobParams.partition = 'all';
JobParams.cpus = 4;
JobParams.memory = 16000;
SlurmParams.cmd = 'SpectroSemantic_Neuro_model_glm_savio(''%s'');';
SlurmParams.resultsDirectory=fullfile('/auto','tdrive','julie','k6','julie','matfile','ModMatInfo');

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
    for dd=1:length(DoneFile)
        if strcmp(DoneFile(dd).name,['Models_InfoPoisson' TheFile(9:end) ext])
            ThisFiledone=1;
            break
        end
    end
    if ThisFiledone
        % Don't recalculate that file
    else
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

%% Set up the timer that is going to check how many of my jobs are queued
% and start job if that number is below NJobs
JobTimer = timer;
JobTimer.BusyMode = 'queue';
JobTimer.ExecutionMode = 'fixedRate';
JobTimer.Name = 'JobTimerName';
JobTimer.Period = PeriodTimerCheck;
JobTimer.TimerFcn = {@jobchecker,NJobs, MatNameToDo, MatfileToDo, JobParams,SlurmParams};
JobTimer.UserData = IndOrd;

%Start the Timer
start(JobTimer);


   