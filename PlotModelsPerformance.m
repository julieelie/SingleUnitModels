Logistic=load('/Users/elie/Documents/CODE/data/matfile/LogisticForTyler.mat');
SWITCH.Info=0;
SWITCH.AcOnly=1;
STRF_show=1;
NL_show=1; % Old switch for the out put non-linearity code, commented now
LL_show=0; % To see loglikelihood plots
DiffDeviance_calc=0;
SignifThresh=0.01;
Pause_Switch=1;
Figure_Switch=1;
NumModeles=19;
%load('SynCell_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
%load('DC128_GLMPoisson.mat')
%load('Ag127_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
%load('ACCell_CellResultsAlphaLambdaDev_AllModels_poisson_log.mat')
% load('WholeVoc_Site2_L1000R900_e13_s0_ss1_GLMPoisson.mat')% SYN
% load('WholeVoc_Site2_L1100R1450_e14_s0_ss1_GLMPoisson.mat')% Ac
% load('WholeVoc_Site3_L2500R2300_e22_s1_ss1_GLMPoisson.mat') %Ag127
% load('WholeVoc_Site4_L1500R1900_e23_s0_ss2_GLMPoisson.mat')%DC128 Path:/auto/k6/julie/matfile/GreBlu9508M/WholeVoc_Site4_L1500R1900_e23_s0_ss2.mat
cd /Users/elie/Documents/CODE/data/matfile/ModMatAcOnly
%cd /Users/elie/Documents/CODE/data/matfile/ModMatSavio/DoneCells/DoneCellsFeb
PoissonFiles=dir('Models*Site4_L1500R1900_e23_s0_ss2.mat');

MinDevlog10Lambdas.Acoustic = [];
MinDevlog10Lambdas.AcSem = [];
MinDevlog10Lambdas.AcSem2 = [];

NumCells=length(PoissonFiles);
AC.FitIndexsumLL = nan(NumCells,6);
AC.FitIndex = nan(NumCells,6);
AC.MeanDiffDev = nan(NumCells,6);
AC.SignifDiffDev = nan(NumCells,7);

Signif.AcvsFl = zeros(NumCells,NumModeles);
Signif.AcSemvsFl = zeros(NumCells,NumModeles);
Signif.CeilvsFl = zeros(NumCells,NumModeles);
Signif.SemvsFl = zeros(NumCells,NumModeles);
Signif.AcSemvsAc = zeros(NumCells,NumModeles);


NumMod_Cells = nan(length(PoissonFiles),1);
NumCells=length(PoissonFiles);
LLAverage.Acoustic = nan(NumModeles,NumCells,2);
LLAverage.AcSemSem = nan(NumModeles,NumCells,2);
LLAverage.AcSemAc = nan(NumModeles,NumCells,2);
LLAverage.Semantic = nan(NumModeles,NumCells,2);
LLAverage.Ceiling = nan(NumModeles,NumCells,2);
LLAverage.Floor = nan(NumModeles,NumCells,2);
LLAverage.AutoRegressive = nan(NumModeles,NumCells,2);
LLAverage.AcSemOpt = nan(NumModeles,NumCells,2);
LLAverage.CeilingOpt = nan(NumModeles,NumCells,2);
LLAverage.Saturated = nan(NumModeles,NumCells,2);

ModelInfo.Acoustic = nan(NumModeles,NumCells,1);
ModelInfo.AcSemSem = nan(NumModeles,NumCells,1);
ModelInfo.AcSemAc = nan(NumModeles,NumCells,1);
ModelInfo.Semantic = nan(NumModeles,NumCells,1);
ModelInfo.Ceiling = nan(NumModeles,NumCells,1);
ModelInfo.Floor = nan(NumModeles,NumCells,1);
ModelInfo.AutoRegressive = nan(NumModeles,NumCells,1);
ModelInfo.AcSemOpt = nan(NumModeles,NumCells,1);
ModelInfo.CeilingOpt = nan(NumModeles,NumCells,1);
ModelInfo.Saturated = nan(NumModeles,NumCells,1);

SpikeCountAverage = nan(NumModeles,NumCells);
%ff=0;
for ff=1:length(PoissonFiles)
%for pp=1:10
    fprintf('File %d/%d %s\n',ff,length(PoissonFiles),PoissonFiles(ff).name);
    try
        load(PoissonFiles(ff).name)
        D=1;
    catch
        fprintf('problem with the file\n')
        D=0;
    end
    
    if D
        if SWITCH.AcOnly
            ParamModel.ModelChoice=[1 0 0 0 0];
        end
        NumMod=length(Wins);
        UnrunWindows =[];
        if ~isfield(Data, 'MeanSpectroStim')
            Data.MeanSpectroStim = Model.MeanSpectroStim;
            ParamModel.MeanSubstractSpec=0;
        end
        for ww=1:NumMod
            if isempty(Data.MeanSpectroStim{ww})
                UnrunWindows = [UnrunWindows ww];
            end
        end
        StartWin = min(UnrunWindows);
        if isempty(UnrunWindows)
            fprintf('All windows calculated\n')
        else
            fprintf('Looking only at some of the windows, there are missing ones\n')
            Wins=setdiff(Wins, Wins(UnrunWindows));
            NumMod = length(Wins);
            %ff=ff+1;
            
            %if Pause_Switch
        end
        NumMod_Cells(ff)=NumMod;
        close all
        %% Retrieve input dataset at each window
        Data.x_wholeset = cell(NumMod,1);
        Data.x_mean_wholeset = cell(NumMod,1);
        Data.x_std_wholeset = cell(NumMod,1);
        Data.y_wholeset_bestGuess_perstim = cell(NumMod,1);
        Data.y_wholeset_AR_perstim = cell(NumMod,1);
        Data.y_wholeset_perstim = cell(NumMod,1);
        Data.x_voc_wholeset_perstim = cell(NumMod,1);
        Data.X_voc=cell(NumMod,1);
        for mm=1:NumMod
            Win=Wins(mm);
            Data_local = getSpectroStims(Spectro, Data.x_stim_indices_wholeset{mm},Data.MeanSpectroStim{mm}, Data.x_stim_repetition{mm}, Data.FLow, Win, ParamModel.MeanSubstractSpec);
            if ParamModel.MeanSubstractSpec
                Data.x_wholeset{mm} = Data_local.x_minusMean;
            elseif ParamModel.ZC
                Data.x_wholeset{mm} = Data_local.x_ZC;% Note that here this dataset is per stim!
            end
            Data.x_mean_wholeset{mm} = Data_local.x_mean;
            Data.x_std_wholeset{mm} = Data_local.x_std;
            Data.y_wholeset_bestGuess_perstim{mm} = nan(length(Data.x_stim_repetition{mm}),1);
            %Data.y_wholeset_AR_perstim{mm} = nan(length(Data.x_stim_repetition{mm}),1);
            Data.y_wholeset_perstim{mm} = nan(length(Data.x_stim_repetition{mm}),1);
            Data.x_wholeset_perstim{mm} = nan(length(Data.x_stim_repetition{mm}),size(Data.x_wholeset{mm},2));
            Data.x_voc_wholeset_perstim{mm} = cell(length(Data.x_stim_repetition{mm}),1);
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
            Data.X_voc{mm}=X_voc;
            
            
            NT=0;% this count the total number of trials
            for ss=1:length(Data.x_stim_repetition{mm})
                Data.x_wholeset_perstim{mm}(ss,:) = Data.x_wholeset{mm}(NT+1,:);
                Data.y_wholeset_perstim{mm}(ss) = Data.y_wholeset{mm}(NT+1);
                Data.x_voc_wholeset_perstim{mm}{ss} = Data.X_voc_wholeset{mm}{NT+1};
                Data.y_wholeset_bestGuess_perstim{mm}(ss) = mean(Data.y_wholeset_bestGuess{mm}((NT+1) : (Data.x_stim_repetition{mm}(ss)+NT)));
                %Data.y_wholeset_AR_perstim{mm}(ss) = mean(Data.y_wholeset_AR{mm}((NT+1) : (Data.x_stim_repetition{mm}(ss)+NT)));
                NT = Data.x_stim_repetition{mm}(ss) +NT;
            end
            
        end
        
        
        %% Check the cell format of model parameters and correct their position in the structure
        for mm=1:NumMod
            if ParamModel.ModelChoice(1) && iscell(Model.Acoustic.B{mm})
                Model.Acoustic.B{mm} = Model.Acoustic.B{mm}{1};
            end
            if ParamModel.ModelChoice(2) && iscell(Model.Semantic.B0{mm})
                Model.Semantic.B{mm} = Model.Semantic.B0{mm}{1}(2:end);
                Model.Semantic.B0{mm} = Model.Semantic.B0{mm}{1}(1);
            end
            if ParamModel.ModelChoice(5) && iscell(Model.AcSemAc.Bspectro{mm})
                Model.AcSemAc.Bspectro{mm} = Model.AcSemAc.Bspectro{mm}{1};
            end
            if ParamModel.ModelChoice(4) && iscell(Model.AcSemSem.Bspectro{mm})
                Model.AcSemSem.Bspectro{mm} = Model.AcSemSem.Bspectro{mm}{1};
            end
            if ParamModel.ModelChoice(5) && iscell(Model.AcSemAc.Bsem{mm})
                Model.AcSemAc.Bsem{mm}=Model.AcSemAc.Bsem{mm}{1};
            end
            if ParamModel.ModelChoice(4) && iscell(Model.AcSemSem.Bsem{mm})
                Model.AcSemSem.Bsem{mm} = Model.AcSemSem.Bsem{mm}{1};
            end
        end
        
        %% Calculate the information of each model
        if SWITCH.Info%any(~isfield(Model, {'Acoustic.Info', 'Semantic.Info', 'AcSemAc.Info','AcSemSem.Info', 'Ceiling.Info','Floor.Info','AR.Info'}))
            Model.Ceiling.info = nan(NumMod,1);
            %Model.AR.info = nan(NumMod,1);
            Model.Acoustic.y_predict = cell(NumMod,1);
            Model.Acoustic.info = nan(NumMod,1);
            Model.Semantic.y_predict = cell(NumMod,1);
            Model.Semantic.info = nan(NumMod,1);
            Model.AcSemAc.y_predict = cell(NumMod,1);
            Model.AcSemAc.info = nan(NumMod,1);
            Model.AcSemSem.y_predict = cell(NumMod,1);
            Model.AcSemSem.info = nan(NumMod,1);
            Model.Floor.ypredict = cell(NumMod,1);
            Model.Floor.info = nan(NumMod,1);
            for mm=1:NumMod
                Wins(mm)
                fprintf('**Ceiling**\n')
                MaxYpredictInfo = max(Data.y_wholeset_bestGuess_perstim{mm});
                % Ceiling Model
                Model.Ceiling.info(mm) = info_model_Calculus(Data.y_wholeset_bestGuess_perstim{mm},MaxYpredictInfo);
                ModelInfo.Ceiling(mm,ff) = Model.Ceiling.info(mm);
                
                % AR Model
                %fprintf('**AR**\n')
                %Model.AR.info(mm) = info_model_Calculus(Data.y_wholeset_AR_perstim{mm},MaxYpredictInfo);
                % ModelInfo.AutoRegressive(mm,ff) = Model.AR.info(mm);
                
                % Acoustic model
                fprintf('**Acoustic**\n')
                b = [Model.Acoustic.B0{mm}; reshape(Model.Acoustic.B{mm},size(Model.Acoustic.B{mm},1)*size(Model.Acoustic.B{mm},2),1)];
                Model.Acoustic.y_predict{mm} = glmval(b, Data.x_wholeset_perstim{mm},ParamModel.LINK);
                Model.Acoustic.info(mm) = info_model_Calculus(Model.Acoustic.y_predict{mm},MaxYpredictInfo);
                ModelInfo.Acoustic(mm,ff) = Model.Acoustic.info(mm);
                
                % Semantic Model
                fprintf('**Semantic**\n')
                b = [Model.Semantic.B0{mm} ; Model.Semantic.B{mm}];
                Model.Semantic.y_predict{mm} = glmval(b, Data.X_voc{mm},ParamModel.LINK);
                Model.Semantic.info(mm) = info_model_Calculus(Model.Semantic.y_predict{mm},MaxYpredictInfo);
                ModelInfo.Semantic(mm,ff) = Model.Semantic.info(mm);
                
                % AcSemAc
                fprintf('**AcSemAc**\n')
                b = [Model.AcSemAc.B0{mm} ; Model.AcSemAc.Bsem{mm} ; reshape(Model.AcSemAc.Bspectro{mm}, size(Model.AcSemAc.Bspectro{mm},1)*size(Model.AcSemAc.Bspectro{mm},2),1)];
                Model.AcSemAc.y_predict{mm} = glmval(b, [Data.X_voc{mm} Data.x_wholeset_perstim{mm}],ParamModel.LINK);
                Model.AcSemAc.info(mm) = info_model_Calculus(Model.AcSemAc.y_predict{mm},MaxYpredictInfo);
                ModelInfo.AcSemAc(mm,ff) = Model.AcSemAc.info(mm);
                
                % AcSemSem
                fprintf('**AcSemSem**\n')
                b = [Model.AcSemSem.B0{mm} ; Model.AcSemSem.Bsem{mm} ; reshape(Model.AcSemSem.Bspectro{mm}, size(Model.AcSemSem.Bspectro{mm},1)*size(Model.AcSemSem.Bspectro{mm},2),1)];
                Model.AcSemSem.y_predict{mm} = glmval(b, [Data.X_voc{mm} Data.x_wholeset_perstim{mm}],ParamModel.LINK);
                Model.AcSemSem.info(mm) = info_model_Calculus(Model.AcSemSem.y_predict{mm},MaxYpredictInfo);
                ModelInfo.AcSemSem(mm,ff) = Model.AcSemSem.info(mm);
                
                % Floor
                fprintf('**Nul-model**\n')
                Model.Floor.y_predict{mm} = repmat(mean(Data.y_wholeset_perstim{mm}),length(Data.x_stim_repetition{mm}),1);
                Model.Floor.info(mm) = info_model_Calculus(Model.Floor.y_predict{mm},MaxYpredictInfo);
                ModelInfo.Floor(mm,ff) = Model.Floor.info(mm);
            end
        end
        %% Ploting loglikelihood as a function of windows if we have more than one window
        if LL_show
            if NumMod>1
                aa=0;
                LL_allVal.Acoustic = nan(20*NumMod,1);
                LL_allVal.AcSemSem = nan(20*NumMod,1);
                LL_allVal.AcSemAc = nan(20*NumMod,1);
                LL_allVal.AcSemOpt = nan(20*NumMod,1);
                LL_allVal.Semantic = nan(20*NumMod,1);
                LL_allVal.Ceiling = nan(20*NumMod,1);
                LL_allVal.Floor = nan(20*NumMod,1);
                LL_allVal.AutoRegressive = nan(20*NumMod,1);
                LLSUM.Acoustic = nan(NumMod,1);
                LLSUM.AcSemSem = nan(NumMod,1);
                LLSUM.AcSemAc = nan(NumMod,1);
                LLSUM.Semantic = nan(NumMod,1);
                LLSUM.Ceiling = nan(NumMod,1);
                LLSUM.Floor = nan(NumMod,1);
                LLSUM.AutoRegressive = nan(NumMod,1);
                %LL.NumStimTrials = cell(NumMod,1);
                %LL.Saturated.values = cell(NumMod,1);
                for mm=1:NumMod
                    %LL.NumStimTrials{mm} = nan(length(Deviance.yobsVal{mm}),1);
                    %LL.Saturated.values{mm}=nan(length(Deviance.yobsVal{mm}),1);
                    %for bb=1:length(Deviance.yobsVal{mm})
                    %LL.NumStimTrials{mm}(bb) = length(Deviance.yobsVal{mm}{bb});
                    %LL.Saturated.values{mm}(bb) = LL_Calculus(Deviance.yobsVal{mm}{bb}, Deviance.yobsVal{mm}{bb});
                    %end
                    if ParamModel.ModelChoice(1)
                        LLAverage.Acoustic(mm,ff,1) = mean(LL.Acoustic.values{mm}./LL.NumStimTrials{mm});
                        LLAverage.Acoustic(mm,ff,2) = 2*std(LL.Acoustic.values{mm}./LL.NumStimTrials{mm})/(length(LL.Acoustic.values{mm})^0.5);
                        LLSUM.Acoustic(mm) = sum(LL.Acoustic.values{mm}./LL.NumStimTrials{mm});
                        LL_allVal.Acoustic((aa+1):(aa+NIt))=LL.Acoustic.values{mm}./LL.NumStimTrials{mm};
                    end
                    if ParamModel.ModelChoice(2)
                        LLAverage.Semantic(mm,ff,1) = mean(LL.Semantic.values{mm}./LL.NumStimTrials{mm});
                        LLAverage.Semantic(mm,ff,2) = 2*std(LL.Semantic.values{mm}./LL.NumStimTrials{mm})/(length(LL.Semantic.values{mm})^0.5);
                        LLSUM.Semantic(mm) = sum(LL.Semantic.values{mm}./LL.NumStimTrials{mm});
                        LL_allVal.Semantic((aa+1):(aa+NIt)) = LL.Semantic.values{mm}./LL.NumStimTrials{mm};
                    end
                    if ParamModel.ModelChoice(4)
                        LLAverage.AcSemSem(mm,ff,1) =mean(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm});
                        LLAverage.AcSemSem(mm,ff,2) =2*std(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm})/(length(LL.AcSemSem.values{mm})^0.5);
                        LLSUM.AcSemSem(mm) =sum(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm});
                        LL_allVal.AcSemSem((aa+1):(aa+NIt)) = LL.AcSemSem.values{mm}./LL.NumStimTrials{mm};
                    end
                    if ParamModel.ModelChoice(5)
                        LLAverage.AcSemAc(mm,ff,1) = mean(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm});
                        LLAverage.AcSemAc(mm,ff,2) = 2*std(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm})/(length(LL.AcSemAc.values{mm})^0.5);
                        LLSUM.AcSemAc(mm) = sum(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm});
                        LL_allVal.AcSemAc((aa+1):(aa+NIt)) = LL.AcSemAc.values{mm}./LL.NumStimTrials{mm};
                        if ParamModel.Modelchoice(4)
                            if sum(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm})>=sum(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm})
                                LL_allVal.AcSemOpt((aa+1):(aa+NIt)) = LL.AcSemSem.values{mm}./LL.NumStimTrials{mm};
                            else
                                LL_allVal.AcSemOpt((aa+1):(aa+NIt)) = LL.AcSemAc.values{mm}./LL.NumStimTrials{mm};
                            end
                        end
                    end
                    LLAverage.Ceiling(mm,ff,1) = mean(LL.Ceiling.values{mm}./LL.NumStimTrials{mm});
                    LLAverage.Floor(mm,ff,1) = mean(LL.Floor.values{mm}./LL.NumStimTrials{mm});
                    LLAverage.AutoRegressive(mm,ff,1) = mean(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm});
                    LLAverage.Saturated(mm, ff, 1) = mean(LL.Saturated.values{mm}./LL.NumStimTrials{mm});
                    LLAverage.Ceiling(mm,ff,2) = 2*std(LL.Ceiling.values{mm}./LL.NumStimTrials{mm})/(length(LL.Ceiling.values{mm})^0.5);
                    LLAverage.Floor(mm,ff,2) = 2*std(LL.Floor.values{mm}./LL.NumStimTrials{mm})/(length(LL.Floor.values{mm})^0.5);
                    LLAverage.AutoRegressive(mm,ff,2) = 2*std(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm})/(length(LL.AutoRegressiveVal.values{mm})^0.5);
                    LLAverage.Saturated(mm, ff, 2) = 2*std(LL.Saturated.values{mm}./LL.NumStimTrials{mm})/(length(LL.AutoRegressiveVal.values{mm})^0.5);
                    
                    SpikeCountAverage(mm, ff) = mean(Data.y_wholeset{mm});
                    
                    LLSUM.Ceiling(mm) = sum(LL.Ceiling.values{mm}./LL.NumStimTrials{mm});
                    LLSUM.Floor(mm) = sum(LL.Floor.values{mm}./LL.NumStimTrials{mm});
                    LLSUM.AutoRegressive(mm) = sum(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm});
                    NIt=length(PropVal.values{mm});
                    
                    LL_allVal.Ceiling((aa+1):(aa+NIt)) = LL.Ceiling.values{mm}./LL.NumStimTrials{mm};
                    LL_allVal.Floor((aa+1):(aa+NIt)) = LL.Floor.values{mm}./LL.NumStimTrials{mm};
                    LL_allVal.AutoRegressive((aa+1):(aa+NIt)) = LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm};
                    aa=aa+NIt;
                end
                if ParamModel.ModelChoice(1)
                    LL_allVal.Acoustic = LL_allVal.Acoustic(1:aa);
                end
                if ParamModel.ModelChoice(4)
                    LL_allVal.AcSemSem = LL_allVal.AcSemSem(1:aa);
                end
                if ParamModel.ModelChoice(5)
                    LL_allVal.AcSemAc = LL_allVal.AcSemAc(1:aa);
                    if ParamModel.ModelChoice(4)
                        LL_allVal.AcSemOpt = LL_allVal.AcSemOpt(1:aa);
                    end
                end
                if ParamModel.ModelChoice(2)
                    LL_allVal.Semantic = LL_allVal.Semantic(1:aa);
                end
                LL_allVal.Ceiling = LL_allVal.Ceiling(1:aa);
                LL_allVal.Floor = LL_allVal.Floor(1:aa);
                LL_allVal.AutoRegressive = LL_allVal.AutoRegressive(1:aa);
                
                models=fieldnames(LLAverage);
                Colors='brmgckybr';
                L=nan(numel(models)-2,1);
                if Figure_Switch
                    figure(5)
                    for ii=1:(numel(models)-2)
                        LLAverage_local=LLAverage.(models{ii})(:,ff,1);
                        LLAverage_local_ste=LLAverage.(models{ii})(:,ff,2);
                        ss1=subplot(4,1,1);
                        shadedErrorBar(1:NumMod,LLAverage_local,LLAverage_local_ste,sprintf(Colors(ii)),1)
                        hold on
                        ss2=subplot(4,1,2);
                        L(ii)=plot(1:NumMod,LLAverage_local,sprintf(Colors(ii)));
                        hold on
                    end
                    legend(L,models);
                    subplot(4,1,1)
                    ylabel('LogLikelihood');
                    xlabel('Window position');
                    %set(ss1,'XLim',[0 NumMod+1]);
                    Xtickposition=get(ss1,'XTick');
                    set(ss1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
                    hold off
                    subplot(4,1,2)
                    ylabel('LogLikelihood');
                    xlabel('Window position');
                    %set(ss1,'XLim',[0 NumMod+1]);;
                    Xtickposition=get(ss2,'XTick');
                    set(ss2,'XTickLabel', [0 Wins(Xtickposition(2:end))]);
                    hold off
                end
                
                
                for mm=1:NumMod
                    if ParamModel.ModelChoice(4) && ParamModel.ModelChoice(5)
                        if sum(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm})>=sum(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm})
                            LLAverage.AcSemOpt(mm,ff,1) = LLAverage.AcSemSem(mm,ff,1);
                            LLAverage.AcSemOpt(mm,ff,2)=LLAverage.AcSemSem(mm,ff,2);
                            ModelInfo.ACSemOpt(mm,ff) = ModelInfo.AcSemSem(mm,ff);
                        else
                            LLAverage.AcSemOpt(mm,ff,1) = LLAverage.AcSemAc(mm,ff,1);
                            LLAverage.AcSemOpt(mm,ff,2) = LLAverage.AcSemAc(mm,ff,2);
                            ModelInfo.ACSemOpt(mm,ff) = ModelInfo.AcSemAc(mm,ff);
                        end
                    end
                    if sum(LL.Ceiling.values{mm}./LL.NumStimTrials{mm})>=sum(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm})
                        LLAverage.CeilingOpt(mm,ff,1)=LLAverage.Ceiling(mm,ff,1);
                        LLAverage.CeilingOpt(mm,ff,2)=LLAverage.Ceiling(mm,ff,2);
                        ModelInfo.CeilingOpt(mm,ff)=ModelInfo.Ceiling(mm,ff);
                    else
                        LLAverage.CeilingOpt(mm,ff,1)=LLAverage.AutoRegressive(mm,ff,1);
                        LLAverage.CeilingOpt(mm,ff,2)=LLAverage.AutoRegressive(mm,ff,2);
                        ModelInfo.CeilingOpt(mm,ff)=ModelInfo.AutoRegressive(mm,ff);
                    end
                end
                if Figure_Switch
                    Colors='brmgcky';
                    figure(5)
                    ss3=subplot(4,1,3);
                    shadedErrorBar(1:NumMod,LLAverage.Acoustic(:,ff,1),LLAverage.Acoustic(:,ff,2),'b',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.AcSemOpt(:,ff,1),LLAverage.AcSemOpt(:,ff,2),'r',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.Semantic(:,ff,1),LLAverage.Semantic(:,ff,2),'g',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.CeilingOpt(:,ff,1),LLAverage.CeilingOpt(:,ff,2),'c',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.Floor(:,ff,1),LLAverage.Floor(:,ff,2),'k',1)
                    ss4=subplot(4,1,4);
                    L=plot(1:NumMod,LLAverage.Acoustic(:,ff,1),'b',1:NumMod,LLAverage.AcSemOpt(:,ff,1),'r',1:NumMod,LLAverage.Semantic(:,ff,1),'g',1:NumMod,LLAverage.CeilingOpt(:,ff,1),'c',1:NumMod,LLAverage.Floor(:,ff,1),'k');
                    legend(L,'Acoustic','AcSem','Semantic','Ceiling','Floor');
                    subplot(4,1,3)
                    ylabel('LogLikelihood');
                    xlabel('Window position');
                    %set(ss1,'XLim',[0 NumMod+1]);
                    Xtickposition=get(ss3,'XTick');
                    set(ss3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
                    hold off
                    subplot(4,1,4)
                    ylabel('LogLikelihood');
                    xlabel('Window position');
                    %set(ss1,'XLim',[0 NumMod+1]);;
                    Xtickposition=get(ss4,'XTick');
                    set(ss4,'XTickLabel', [0 Wins(Xtickposition(2:end))]);
                    hold off
                    if Pause_Switch
                        pause()
                    end
                    set(ss1,'YLim',[-500 10])
                    set(ss2,'YLim',[-500 10])
                    set(ss3,'YLim',[-500 10])
                    set(ss4,'YLim',[-500 10])
                    
                    % Figure for grant
                    figure(16)
                    shadedErrorBar(1:NumMod,LLAverage.Acoustic(:,ff,1),LLAverage.Acoustic(:,ff,2),'b',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.AcSemOpt(:,ff,1),LLAverage.AcSemOpt(:,ff,2),'r',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.CeilingOpt(:,ff,1),LLAverage.CeilingOpt(:,ff,2),'c',1)
                    hold on
                    shadedErrorBar(1:NumMod,LLAverage.Floor(:,ff,1),LLAverage.Floor(:,ff,2),'k',1)
                    ylabel('LogLikelihood');
                    xlabel('Window position');
                    Xtickposition=get(gca,'XTick');
                    set(gca,'XTickLabel', [0 Wins(Xtickposition(2:end))])
                    figure(17)%for the legend....
                    L=plot(1:NumMod,LLAverage.Acoustic(:,ff,1),'b',1:NumMod,LLAverage.AcSemOpt(:,ff,1),'r',1:NumMod,LLAverage.CeilingOpt(:,ff,1),'c',1:NumMod,LLAverage.Floor(:,ff,1),'k');
                    legend(L,'Acoustic','AcSem','Ceiling','Floor');
                    if Pause_Switch
                        pause()
                    end
                end
            end
        end
        %% Testing differences between models based on Deviance calculated on validating dataset
        if DiffDeviance_calc
            DiffDevianceAcousticAcSem=nan(NumMod,2);
            DiffDevianceSemanticAcSem=nan(NumMod,2);
            DiffDevianceAcousticSemantic=nan(NumMod,2);
            DiffDevianceFloorCeiling=nan(NumMod,2);
            DiffDevianceAcousticCeiling=nan(NumMod,2);
            DiffDevianceSemanticCeiling=nan(NumMod,2);
            DiffDevianceAcSemCeiling=nan(NumMod,2);
            DiffDevianceAcousticFloor=nan(NumMod,2);
            DiffDevianceSemanticFloor=nan(NumMod,2);
            DiffDevianceAcSemFloor=nan(NumMod,2);
            pAcousticAcSem = nan(NumMod,1);
            pSemanticAcSem=nan(NumMod,1);
            pAcousticSemantic=nan(NumMod,1);
            pFloorCeiling=nan(NumMod,1);
            pAcousticCeiling=nan(NumMod,1);
            pSemanticCeiling=nan(NumMod,1);
            pAcSemCeiling=nan(NumMod,1);
            pAcousticFloor=nan(NumMod,1);
            pSemanticFloor=nan(NumMod,1);
            pAcSemFloor=nan(NumMod,1);
            AcSemOptId = cell(NumMod,1);
            
            for mm=1:NumMod
                %         % Find out which of the AcSem model is the best (smallest deviance)
                %         Deviance_AcSem_Local=nan(length(Deviance.AcSemAc.values{mm}),1);
                %         Ind=find(Deviance.AcSemAc.values{mm}>Deviance.AcSemSem.values{mm});
                %         Ind2=find(Deviance.AcSemAc.values{mm}<Deviance.AcSemSem.values{mm});
                %         Deviance_AcSem_Local(Ind)=Deviance.AcSemAc.values{mm}(Ind);
                %         Deviance_AcSem_Local(Ind2)=Deviance.AcSemAc.values{mm}(Ind2);
                % Find out which of the AcSem is the best over bootstrap (highest
                % sum of LL)
                if sum(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm})>=sum(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm})
                    AcSemOptId{mm} = 'AcSemSem';
                    Deviance_AcSem_Local=Deviance.AcSemSem.values{mm};
                    LL_AcSem_Local = LL.AcSemSem.values{mm};
                else
                    AcSemOptId{mm} = 'AcSemAc';
                    Deviance_AcSem_Local=Deviance.AcSemAc.values{mm};
                    LL_AcSem_Local = LL.AcSemAc.values{mm};
                end
                
                % Find out which of the Ceiling Models is the best over bootstrap (highest
                % sum of LL) Then change this in the next lines to have a better
                % code
                if sum(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm})>=sum(LL.Ceiling.values{mm}./LL.NumStimTrials{mm})
                    LL_OptBestModel_Local = LL.AutoRegressiveVal.values{mm};
                else
                    LL_OptBestModel_Local = LL.Ceiling.values{mm};
                end
                %
                
                % Floor vs ceiling model
                DiffDevianceFloorCeiling_local=(-2*(LL.Floor.values{mm} - LL_OptBestModel_Local))./LL.NumStimTrials{mm};
                DiffDevianceFloorCeiling(mm,1)=mean(DiffDevianceFloorCeiling_local);
                DiffDevianceFloorCeiling(mm,2)=2*std(DiffDevianceFloorCeiling_local)/(length(DiffDevianceFloorCeiling_local)^0.5);
                %one sample t-test
                [h,pFloorCeiling(mm),ci,stats]=ttest(DiffDevianceFloorCeiling_local);
                
                % Acoustic vs ceiling model
                DiffDevianceAcousticCeiling_local=(-2*(LL.Acoustic.values{mm} - LL_OptBestModel_Local))./LL.NumStimTrials{mm};
                DiffDevianceAcousticCeiling(mm,1)=mean(DiffDevianceAcousticCeiling_local);
                DiffDevianceAcousticCeiling(mm,2)=2*std(DiffDevianceAcousticCeiling_local)/(length(DiffDevianceAcousticCeiling_local)^0.5);
                %one sample t-test
                [h,pAcousticCeiling(mm),ci,stats]=ttest(DiffDevianceAcousticCeiling_local);
                
                % Semantic vs ceiling model
                DiffDevianceSemanticCeiling_local=(-2*(LL.Semantic.values{mm} - LL_OptBestModel_Local))./LL.NumStimTrials{mm};
                DiffDevianceSemanticCeiling(mm,1)=mean(DiffDevianceSemanticCeiling_local);
                DiffDevianceSemanticCeiling(mm,2)=2*std(DiffDevianceSemanticCeiling_local)/(length(DiffDevianceSemanticCeiling_local)^0.5);
                %one sample t-test
                [h,pSemanticCeiling(mm),ci,stats]=ttest(DiffDevianceSemanticCeiling_local);
                
                % AcSem vs ceiling model
                DiffDevianceAcSemCeiling_local=(-2*(LL_AcSem_Local - LL_OptBestModel_Local))./LL.NumStimTrials{mm};
                DiffDevianceAcSemCeiling(mm,1)=mean(DiffDevianceAcSemCeiling_local);
                DiffDevianceAcSemCeiling(mm,2)=2*std(DiffDevianceAcSemCeiling_local)/(length(DiffDevianceAcSemCeiling_local)^0.5);
                %one sample t-test
                [h,pAcSemCeiling(mm),ci,stats]=ttest(DiffDevianceAcSemCeiling_local);
                
                % Acoustic vs Floor model
                DiffDevianceAcousticFloor_local=(2*(LL.Acoustic.values{mm} - LL.Floor.values{mm}))./LL.NumStimTrials{mm};
                DiffDevianceAcousticFloor(mm,1)=mean(DiffDevianceAcousticFloor_local);
                DiffDevianceAcousticFloor(mm,2)=2*std(DiffDevianceAcousticFloor_local)/(length(DiffDevianceAcousticFloor_local)^0.5);
                %one sample t-test
                [h,pAcousticFloor(mm),ci,stats]=ttest(DiffDevianceAcousticFloor_local);
                
                % Semantic vs Floor model
                DiffDevianceSemanticFloor_local=(2*(LL.Semantic.values{mm} - LL.Floor.values{mm}))./LL.NumStimTrials{mm};
                DiffDevianceSemanticFloor(mm,1)=mean(DiffDevianceSemanticFloor_local);
                DiffDevianceSemanticFloor(mm,2)=2*std(DiffDevianceSemanticFloor_local)/(length(DiffDevianceSemanticFloor_local)^0.5);
                %one sample t-test
                [h,pSemanticFloor(mm),ci,stats]=ttest(DiffDevianceSemanticFloor_local);
                
                % AcSem vs Floor model
                DiffDevianceAcSemFloor_local=(2*(LL_AcSem_Local - LL.Floor.values{mm}))./LL.NumStimTrials{mm};
                DiffDevianceAcSemFloor(mm,1)=mean(DiffDevianceAcSemFloor_local);
                DiffDevianceAcSemFloor(mm,2)=2*std(DiffDevianceAcSemFloor_local)/(length(DiffDevianceAcSemFloor_local)^0.5);
                %one sample t-test
                [h,pAcSemFloor(mm),ci,stats]=ttest(DiffDevianceAcSemFloor_local);
                
                % Acoustic vs Acsem Model
                DiffDevianceAcousticAcSem_local=(Deviance.Acoustic.values{mm} - Deviance_AcSem_Local)./LL.NumStimTrials{mm};
                DiffDevianceAcousticAcSem(mm,1)=mean(DiffDevianceAcousticAcSem_local);
                DiffDevianceAcousticAcSem(mm,2)=2*std(DiffDevianceAcousticAcSem_local)/(length(DiffDevianceAcousticAcSem_local)^0.5);
                %one sample t-test
                [h,pAcousticAcSem(mm),ci,stats]=ttest(DiffDevianceAcousticAcSem_local);
                
                % Semantic vs Acsem Model
                DiffDevianceSemanticAcSem_local=(Deviance.Semantic.values{mm} - Deviance_AcSem_Local)./LL.NumStimTrials{mm};
                DiffDevianceSemanticAcSem(mm,1)=mean(DiffDevianceSemanticAcSem_local);
                DiffDevianceSemanticAcSem(mm,2)=2*std(DiffDevianceSemanticAcSem_local)/(length(DiffDevianceSemanticAcSem_local)^0.5);
                %one sample t-test
                [h,pSemanticAcSem(mm),ci,stats]=ttest(DiffDevianceSemanticAcSem_local);
                
                % Semantic vs Acoustic Model
                DiffDevianceAcousticSemantic_local=(Deviance.Acoustic.values{mm} - Deviance.Semantic.values{mm})./LL.NumStimTrials{mm};
                DiffDevianceAcousticSemantic(mm,1)=mean(DiffDevianceAcousticSemantic_local);
                DiffDevianceAcousticSemantic(mm,2)=2*std(DiffDevianceAcousticSemantic_local)/(length(DiffDevianceAcousticSemantic_local)^0.5);
                %one sample t-test
                [h,pAcousticSemantic(mm),ci,stats]=ttest(DiffDevianceAcousticSemantic_local);
                
            end
            
            %% Plot the difference of deviance over time
            %Select only windows for which models are similar or better than Floor
            % model
            GoodAc=intersect(find(DiffDevianceAcousticFloor(:,1)>=0),find(pAcousticFloor<SignifThresh));
            Signif.AcvsFl(ff,GoodAc)=1;
            GoodAcSem=intersect(find(DiffDevianceAcSemFloor(:,1)>=0), find(pAcSemFloor<SignifThresh));
            Signif.AcSemvsFl(ff,GoodAcSem)=1;
            GoodSem=intersect(find(DiffDevianceSemanticFloor(:,1)>=0), find(pSemanticFloor<SignifThresh));
            Signif.SemvsFl(ff,GoodSem)=1;
            Signif.AcSemvsAc(ff,intersect(find(DiffDevianceAcousticAcSem(:,1)>=0), find(pAcousticAcSem<SignifThresh)))=1;
            Signif.CeilvsFl(ff,intersect(find(DiffDevianceFloorCeiling(:,1)>=0), find(pFloorCeiling<SignifThresh)))=1;
            
            %     GoodAc=find(DiffDevianceAcousticFloor(:,1)>=0);
            %     GoodAcSem=find(DiffDevianceAcSemFloor(:,1)>=0);
            %     GoodSem=find(DiffDevianceSemanticFloor(:,1)>=0);
            GoodAcousticAcSem=union(GoodAc,GoodAcSem);
            GoodAcousticSem=union(GoodAc,GoodSem);
            GoodSemanticAcSem=union(GoodSem,GoodAcSem);
            if Figure_Switch
                if (length(GoodAcousticAcSem)<=1) || (length(GoodAcousticSem)<=1) || (length(GoodSemanticAcSem)<=1)
                    fprintf('There is al least one model for which there is no more that a single time point significantly performing better than the floor model/n->no figure of the difference of deviance/n');
                else
                    figure(6)
                    ss1=subplot(3,1,1);
                    shadedErrorBar(GoodAcousticAcSem,DiffDevianceAcousticAcSem(GoodAcousticAcSem,1),DiffDevianceAcousticAcSem(GoodAcousticAcSem,2),'-g',1)
                    hold on
                    shadedErrorBar(GoodSemanticAcSem,DiffDevianceSemanticAcSem(GoodSemanticAcSem,1),DiffDevianceSemanticAcSem(GoodSemanticAcSem,2),'.-b',1)
                    hold on
                    shadedErrorBar(GoodAcousticSem,DiffDevianceAcousticSemantic(GoodAcousticSem,1),DiffDevianceAcousticSemantic(GoodAcousticSem,2),'.-k',1)
                    set(ss1,'XLim',[0 NumMod])
                    xlabel('Time (ms)')
                    ylabel('Difference of Deviance (Saturated Model)')
                    Xtickposition=get(ss1,'XTick');
                    set(ss1,'XTickLabel', [0 Wins(Xtickposition(2:end))]);
                    hline(0,'k:')
                    hold off
                    
                    ss2=subplot(3,1,2);
                    plot(GoodAcousticAcSem,DiffDevianceAcousticAcSem(GoodAcousticAcSem,1),'.-g',GoodSemanticAcSem,DiffDevianceSemanticAcSem(GoodSemanticAcSem,1),'.-b',GoodAcousticSem,DiffDevianceAcousticSemantic(GoodAcousticSem,1),'.-k')
                    hold on
                    hline(0,'k:')
                    set(ss2,'XLim',[0 NumMod])
                    Xtickposition=get(ss2,'XTick');
                    set(ss2,'XTickLabel', [0 Wins(Xtickposition(2:end))]);
                    legend('Acoustic-AcSem','Semantic-AcSem','Acoustic-Semantic','Location','SouthEast')
                    hold off
                    if Pause_Switch
                        pause()
                    end
                    set(ss1,'YLim',[-10 10])
                    set(ss2,'YLim',[-10 10])
                end
            end
            %% Store data for population in terms of average difference of deviance
            AC.MeanDiffDev(ff,1)=mean(DiffDevianceSemanticAcSem(GoodSemanticAcSem,1));
            AC.MeanDiffDev(ff,2)=mean(DiffDevianceAcousticAcSem(GoodAcousticAcSem,1));
            AC.MeanDiffDev(ff,3)=mean(DiffDevianceAcousticSemantic(GoodAcousticSem,1));
            AC.MeanDiffDev(ff,4) = mean(-2.*(LL_allVal.Floor-LL_allVal.Ceiling));
            AC.MeanDiffDev(ff,5) = mean(-2.*(LL_allVal.Floor-LL_allVal.Acoustic));
            AC.MeanDiffDev(ff,6) = mean(-2.*(LL_allVal.Floor-LL_allVal.AcSemOpt));
            AC.MeanDiffDev(ff,7) = mean(-2.*(LL_allVal.Floor-LL_allVal.Semantic));
            %one sample t-test
            fprintf('Significant differences:\n')
            % Floor vs ceilings
            [h,pFloorCeil_whole,ci,stats]=ttest(-2.*(LL_allVal.Floor-LL_allVal.Ceiling));
            if (~isnan(h)) && h
                fprintf('***Significance Ceiling model over Floor p=%f, DiffDev=%f+/-%f\n',pFloorCeil_whole,mean(-2.*(LL_allVal.Floor-LL_allVal.Ceiling)),std(-2.*(LL_allVal.Floor-LL_allVal.Ceiling)));
                if mean(-2.*(LL_allVal.Floor-LL_allVal.Ceiling))>0
                    AC.SignifDiffDev(ff,4)=1;
                else
                    AC.SignifDiffDev(ff,4)=0;
                end
            else
                fprintf('Ceiling model:NS\n')
                AC.SignifDiffDev(ff,4)=h;
            end
            
            % Floor vs Acoustic
            [h,pFloorAc_whole,ci,stats]=ttest(-2.*(LL_allVal.Floor-LL_allVal.Acoustic));
            if (~isnan(h)) && h
                fprintf('***Significance Acoustic model over Floor p=%f, DiffDev=%f+/-%f\n',pFloorAc_whole,mean(-2.*(LL_allVal.Floor-LL_allVal.Acoustic)),std(-2.*(LL_allVal.Floor-LL_allVal.Acoustic)));
                if mean(-2.*(LL_allVal.Floor-LL_allVal.Acoustic))>0
                    AC.SignifDiffDev(ff,5)=1;
                else
                    AC.SignifDiffDev(ff,5)=0;
                end
            else
                fprintf('Acoustic model:NS\n')
                AC.SignifDiffDev(ff,5)=h;
            end
            
            % Floor vs AcSem
            [h,pFloorAcSem_whole,ci,stats]=ttest(-2.*(LL_allVal.Floor-LL_allVal.AcSemOpt));
            if (~isnan(h)) && h
                fprintf('***Significance AcSem Opt model over Floor p=%f, DiffDev=%f+/-%f\n',pFloorAcSem_whole,mean(-2.*(LL_allVal.Floor-LL_allVal.AcSemOpt)),std(-2.*(LL_allVal.Floor-LL_allVal.AcSemOpt)));
                if mean(-2.*(LL_allVal.Floor-LL_allVal.AcSemOpt))>0
                    AC.SignifDiffDev(ff,6)=1;
                else
                    AC.SignifDiffDev(ff,6)=0;
                end
            else
                fprintf('AcSem Opt model:NS\n')
                AC.SignifDiffDev(ff,6)=h;
            end
            
            % Floor vs Sem
            [h,pFloorSem_whole,ci,stats]=ttest(-2.*(LL_allVal.Floor-LL_allVal.Semantic));
            if (~isnan(h)) && h
                fprintf('***Significance Semantic model over Floor p=%f, DiffDev=%f+/-%f\n',pFloorAcSem_whole,mean(-2.*(LL_allVal.Floor-LL_allVal.Semantic)),std(-2.*(LL_allVal.Floor-LL_allVal.Semantic)));
                if mean(-2.*(LL_allVal.Floor-LL_allVal.Semantic))>0
                    AC.SignifDiffDev(ff,7)=1;
                else
                    AC.SignifDiffDev(ff,7)=0;
                end
            else
                fprintf('Semantic model:NS\n')
                AC.SignifDiffDev(ff,7)=h;
            end
            
            % Acoustic vs AcSem
            % for all windows
            [h,pAcousticAcSem_whole,ci,stats]=ttest(DiffDevianceAcousticAcSem(:,1));
            if (~isnan(h)) && h
                fprintf('***Significance Semantic parameters over all wins p=%f, DiffDev=%f+/-%f\n',pAcousticAcSem_whole,mean(DiffDevianceAcousticAcSem(:,1)),std(DiffDevianceAcousticAcSem(:,1)));
                if mean(DiffDevianceAcousticAcSem(:,1))>0
                    AC.SignifDiffDev(ff,2)=1;
                else
                    AC.SignifDiffDev(ff,2)=0;
                end
            else
                fprintf('Semantic parameters over all windows:NS\n')
                AC.SignifDiffDev(ff,2)=h;
            end
            %restricting the dataset to windows where both
            % models do better than the floor model
            [h,pAcousticAcSem_whole,ci,stats]=ttest(DiffDevianceAcousticAcSem(GoodAcousticAcSem,1));
            if (~isnan(h)) && h
                fprintf('***Significance Semantic parameters for wins where AcSem and Ac models better than Floor p=%f, DiffDev=%f+/-%f\n',pAcousticAcSem_whole,mean(DiffDevianceAcousticAcSem(GoodAcousticAcSem,1)),std(DiffDevianceAcousticAcSem(GoodAcousticAcSem,1)));
                %             if mean(DiffDevianceAcousticAcSem(GoodAcousticAcSem,1))>0
                %                 AC.SignifDiffDev(ff,2)=1;
                %             else
                %                 AC.SignifDiffDev(ff,2)=0;
                %             end
            else
                fprintf('Semantic parameters for wins where AcSem and Ac models better than Floor:NS\n')
                %             AC.SignifDiffDev(ff,2)=h;
            end
            
            % Semantic vs AcSem restricting the dataset to windows where both
            % models do better than the floor model
            [h,pSemanticAcSem_whole,ci,stats]=ttest(DiffDevianceSemanticAcSem(GoodSemanticAcSem,1));
            if (~isnan(h)) && h
                fprintf('***Significance Acoustic parameters for wins where AcSem and Sem models better than floor p=%f DiffDev=%f+/-%f\n',pSemanticAcSem_whole,mean(DiffDevianceSemanticAcSem(GoodSemanticAcSem,1)),std(DiffDevianceSemanticAcSem(GoodSemanticAcSem,1)));
                if mean(DiffDevianceSemanticAcSem(GoodSemanticAcSem,1))>0
                    %                 AC.SignifDiffDev(ff,1)=1;
                else
                    %                 AC.SignifDiffDev(ff,1)=0;
                end
            else
                fprintf('Acoustic parameters: NS\n')
                %             AC.SignifDiffDev(ff,1)=h;
            end
            
            % Semantic vs AcSem over all windows
            [h,pSemanticAcSem_whole,ci,stats]=ttest(DiffDevianceSemanticAcSem(:,1));
            if (~isnan(h)) && h
                fprintf('***Significance Acoustic parameters over all wins p=%f DiffDev=%f+/-%f\n',pSemanticAcSem_whole,mean(DiffDevianceSemanticAcSem(:,1)),std(DiffDevianceSemanticAcSem(:,1)));
                if mean(DiffDevianceSemanticAcSem(:,1))>0
                    AC.SignifDiffDev(ff,1)=1;
                else
                    AC.SignifDiffDev(ff,1)=0;
                end
            else
                fprintf('Acoustic parameters: NS\n')
                AC.SignifDiffDev(ff,1)=h;
            end
            
            % Acoustic vs Semantic restricting the dataset to windows where both
            % models do better than the floor model
            [h,pAcousticSemantic_whole,ci,stats]=ttest(DiffDevianceAcousticSemantic(GoodAcousticSem,1));
            if (~isnan(h)) && h
                fprintf('Performance Acoustic and Semantic models different p=%f DiffDev=%f+/-%f\n',pAcousticSemantic_whole,mean(DiffDevianceAcousticSemantic(:,1)),std(DiffDevianceAcousticSemantic(:,1)));
                if mean(DiffDevianceAcousticSemantic(:,1))>0
                    AC.SignifDiffDev(ff,3)=1;
                else
                    AC.SignifDiffDev(ff,3)=0;
                end
            else
                fprintf('Semantic vs Acoustic: NS\n')
                AC.SignifDiffDev(ff,3)=h;
            end
            if Pause_Switch
                pause()
            end
            
            % subplot(2,3,4)
            % boxplot(-DiffAIC_AcSemAcoustic)
            % xlabel('Alpha')
            % ylabel('Difference of AIC')
            % set(gca,'XTickLabel',Alphas)
            % title('Acoustic - AcSem');
            % subplot(2,3,5)
            % boxplot(-DiffAIC_AcSemSemantic)
            % xlabel('Alpha')
            % ylabel('Difference of AIC')
            % set(gca,'XTickLabel',Alphas)
            % title('Semantic- AcSem');
            
            
            %% Plot FitIndex values as a function of window position for each model
            % Go Through windows and findout which is the best reference at each
            % time point for the FitIndex. It could be either Ceiling model or
            % autoregressive one, decide depending on which of the two has the
            % highest loglikelihood.
            FitIndex_Acoustic=nan(NumMod,2);
            FitIndex_Semantic=nan(NumMod,2);
            FitIndex_AcSem=nan(NumMod,2);
            for mm=1:NumMod
                % Find out which of the Autoregressive or Ceiling is the best over bootstrap (highest
                % sum of LL)
                if sum(LL.AutoRegressiveVal.values{mm}./LL.NumStimTrials{mm})>=sum(LL.Ceiling.values{mm}./LL.NumStimTrials{mm})
                    LL_OptBestModel_Local = LL.AutoRegressiveVal.values{mm};
                else
                    LL_OptBestModel_Local = LL.Ceiling.values{mm};
                end
                % Find out which of the AcSem is the best over bootstrap (highest
                % sum of LL)
                if sum(LL.AcSemSem.values{mm}./LL.NumStimTrials{mm})>=sum(LL.AcSemAc.values{mm}./LL.NumStimTrials{mm})
                    Deviance_AcSem_Local=Deviance.AcSemSem.values{mm};
                    LL_AcSem_Local = LL.AcSemSem.values{mm};
                else
                    Deviance_AcSem_Local=Deviance.AcSemAc.values{mm};
                    LL_AcSem_Local = LL.AcSemAc.values{mm};
                end
                
                FitIndex_Acoustic_local=(LL.Acoustic.values{mm}-LL.Floor.values{mm})./(LL_OptBestModel_Local-LL.Floor.values{mm});
                FitIndex_Acoustic(mm,1)=mean(FitIndex_Acoustic_local);
                FitIndex_Acoustic(mm,2)=2*std(FitIndex_Acoustic_local)/(length(FitIndex_Acoustic_local)^0.5);
                FitIndex_Semantic_local=(LL.Semantic.values{mm}-LL.Floor.values{mm})./(LL_OptBestModel_Local-LL.Floor.values{mm});
                FitIndex_Semantic(mm,1)=mean(FitIndex_Semantic_local);
                FitIndex_Semantic(mm,2)=2*std(FitIndex_Semantic_local)/(length(FitIndex_Semantic_local)^0.5);
                FitIndex_AcSem_local=(LL_AcSem_Local-LL.Floor.values{mm})./(LL_OptBestModel_Local-LL.Floor.values{mm});
                FitIndex_AcSem(mm,1)=mean(FitIndex_AcSem_local);
                FitIndex_AcSem(mm,2)=2*std(FitIndex_AcSem_local)/(length(FitIndex_AcSem_local)^0.5);
                
                
                if (FitIndex_Acoustic(mm,1)>1) && (pAcousticCeiling(mm)>SignifThresh) % Acoustic model not different from ceiling
                    FitIndex_Acoustic(mm,1)=1;
                elseif FitIndex_Acoustic(mm,1)<0
                    FitIndex_Acoustic(mm,1)=0;
                end
                if (FitIndex_Semantic(mm,1)>1) && (pSemanticCeiling(mm)>SignifThresh) % Semantic model not different from ceiling
                    FitIndex_Semantic(mm,1)=1;
                elseif FitIndex_Semantic(mm,1)<0
                    FitIndex_Semantic(mm,1)=0;
                end
                if (FitIndex_AcSem(mm,1)>1) && (pAcSemCeiling(mm)>SignifThresh) % AcSem model not different from ceiling
                    FitIndex_AcSem(mm,1)=1;
                elseif FitIndex_AcSem(mm,1)<0
                    FitIndex_AcSem(mm,1)=0;
                end
            end
            % Only plot indices that makes sense (Ceiling>Floor)!
            SignifCeil=intersect(find(pFloorCeiling<SignifThresh),find(DiffDevianceFloorCeiling(:,1)>0));
            SignifAc=intersect(SignifCeil,intersect(find(pAcousticAcSem<SignifThresh),find(DiffDevianceAcousticAcSem(:,1)>0)));
            SignifSem=intersect(SignifCeil,intersect(find(pSemanticAcSem<SignifThresh),find(DiffDevianceSemanticAcSem(:,1)>0)));
            if Figure_Switch
                if length(SignifCeil)>1
                    figure(6)
                    ss3=subplot(3,1,3);
                    shadedErrorBar(SignifCeil,FitIndex_Acoustic(SignifCeil,1),FitIndex_Acoustic(SignifCeil,2),'b',1)
                    hold on
                    shadedErrorBar(SignifCeil,FitIndex_Semantic(SignifCeil,1),FitIndex_Semantic(SignifCeil,2),'r',1)
                    hold on
                    shadedErrorBar(SignifCeil,FitIndex_AcSem(SignifCeil,1),FitIndex_AcSem(SignifCeil,2),'c',1)
                    %plot(SignifCeil,FitIndex_Acoustic(SignifCeil),'-g',SignifCeil,FitIndex_Semantic(SignifCeil),'-r',SignifCeil,FitIndex_AcSem(SignifCeil),'-c')
                    hold on
                    plot(SignifAc, FitIndex_AcSem(SignifAc),'dc','MarkerSize',20)
                    hold on
                    plot(SignifSem, FitIndex_AcSem(SignifSem),'*c', 'MarkerSize',10)
                    xlabel('Time (ms)')
                    ylabel('Fit Index')
                    legend('Acoustic','Semantic','AcSem','AcSem > Ac','AcSem > Sem', 'Location','SouthEast')
                    set(ss3,'XLim',[0 14])
                    Xtickposition=get(ss3,'XTick');
                    set(ss3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
                    set(ss3,'YLim',[0 1])
                    hold off
                    
                    fprintf('Average FitIndices over windows:\nAcoustic: %f\nAcSem: %f\nSemantic: %f\n',mean(FitIndex_Acoustic(SignifCeil,1)),mean(FitIndex_AcSem(SignifCeil,1)), mean(FitIndex_Semantic(SignifCeil,1)))
                    fprintf('Medial FitIndices over windows:\nAcoustic: %f\nAcSem: %f\nSemantic: %f\n',median(FitIndex_Acoustic(SignifCeil,1)),median(FitIndex_AcSem(SignifCeil,1)), median(FitIndex_Semantic(SignifCeil,1)))
                    if Pause_Switch
                        pause()
                    end
                    
                    % Figure for the grant
                    figure(15)
                    shadedErrorBar(SignifCeil,FitIndex_Acoustic(SignifCeil,1),FitIndex_Acoustic(SignifCeil,2),'b',1)
                    hold on
                    shadedErrorBar(SignifCeil,FitIndex_AcSem(SignifCeil,1),FitIndex_AcSem(SignifCeil,2),'r',1)
                    %plot(SignifCeil,FitIndex_Acoustic(SignifCeil),'-b',SignifCeil,FitIndex_AcSem(SignifCeil),'-r')
                    hold on
                    plot(SignifAc, FitIndex_AcSem(SignifAc),'*c','MarkerSize',10)
                    hold on
                    plot(GoodAc,repmat(0.95,size(GoodAc)), '.b', 'MarkerSize', 30)
                    hold on
                    plot(GoodAcSem,repmat(0.925,size(GoodAcSem)), '.r', 'MarkerSize', 30)
                    hold on
                    plot(SignifCeil,repmat(0.9,size(SignifCeil)), '.k', 'MarkerSize', 30)
                    hold on
                    xlabel('Time (ms)')
                    ylabel('Fit Index')
                    legend('Acoustic','AcSem','AcSem > Ac', 'Ac>Floor', 'AcSem>Floor','Ceil>Floor', 'Location','SouthEast')
                    set(gca,'XLim',[0 14])
                    Xtickposition=get(gca,'XTick');
                    set(gca,'XTickLabel', [0 Wins(Xtickposition(2:end))])
                    set(gca,'YLim',[0 1])
                    if Pause_Switch
                        pause()
                    end
                    
                elseif isempty(SignifCeil)
                    fprintf('FitIndices do not make any sense since the ceiling model is always doing worse than the floor\n')
                else
                    fprintf('FitIndex makes sense only for one window since other than for that window the ceiling model is always doing worse than the floor\n')
                    fprintf('Average FitIndices over windows:\nAcoustic: %f\nAcSem: %f\nSemantic: %f\n',mean(FitIndex_Acoustic(SignifCeil,1)),mean(FitIndex_AcSem(SignifCeil,1)), mean(FitIndex_Semantic(SignifCeil,1)))
                    fprintf('Medial FitIndices over windows:\nAcoustic: %f\nAcSem: %f\nSemantic: %f\n',median(FitIndex_Acoustic(SignifCeil,1)),median(FitIndex_AcSem(SignifCeil,1)), median(FitIndex_Semantic(SignifCeil,1)))
                    if Pause_Switch
                        pause()
                    end
                end
            end
            % Store Fitindices values for the population
            AC.FitIndex(ff,1)=mean(FitIndex_Acoustic(SignifCeil,1));
            AC.FitIndex(ff,2)=mean(FitIndex_Semantic(SignifCeil,1));
            AC.FitIndex(ff,3)=mean(FitIndex_AcSem(SignifCeil,1));
            
            AC.FitIndexsumLL(ff,1)=(sum(LLAverage.Acoustic(:,ff,1))-sum(LLAverage.Floor(:,ff,1)))/(sum(LLAverage.CeilingOpt(:,ff,1))/sum(LLAverage.Floor(:,ff,1)));
            AC.FitIndexsumLL(ff,2)=(sum(LLAverage.Semantic(:,ff,1))-sum(LLAverage.Floor(:,ff,1)))/(sum(LLAverage.CeilingOpt(:,ff,1))/sum(LLAverage.Floor(:,ff,1)));
            AC.FitIndexsumLL(ff,3)=(sum(LLAverage.AcSemOpt(:,ff,1))-sum(LLAverage.Floor(:,ff,1)))/(sum(LLAverage.CeilingOpt(:,ff,1))/sum(LLAverage.Floor(:,ff,1)));
            
            
        end
        
        %% Plot weights of the projection of the Acoustic STRF on the logistic regression dimensions
        STRFsize='all';
        LRProjection=nan(NumMod,length(Logistic.vocTypes));
        for mm=1:NumMod
            if sum(Model.TickSTRFspectro.fo{mm}-Logistic.fo(1:length(Model.TickSTRFspectro.fo{mm})))
                fprintf('!!!!!WARNING Frequency resolution of STRF and logistic regression are different\n')
            end
            if strcmp(STRFsize,'all')
                PaceNb=length(Logistic.to)-length(Model.TickSTRFspectro.to{mm});
                LocalProjection=nan(PaceNb,length(Logistic.vocTypes));
                for pp=1:PaceNb
                    for vt =1:length(Logistic.vocTypes)
                        LR_Filter=reshape(Logistic.PC_LR(:,vt),Logistic.nf,Logistic.nt);
                        LR_Filter_scaled = LR_Filter(1:length(Model.TickSTRFspectro.fo{mm}),pp:(pp+length(Model.TickSTRFspectro.to{mm})-1));
                        ScalingTerm = ((sum(sum(LR_Filter_scaled.^2)))^0.5)*((sum(sum(Model.Acoustic.B{mm}.^2)))^0.5);
                        LocalProjection(pp,vt) = sum(sum(LR_Filter_scaled.*Model.Acoustic.B{mm}))/ScalingTerm;
                    end
                end
                LRProjection(mm,:)=max(LocalProjection);
            else
                PaceNb=length(Logistic.to)-sum(Model.TickSTRFspectro.to{mm}<=(STRFsize/1000));
                for vt=1:length(Logistic.vocTypes)
                    if Model.TickSTRFspectro.to{mm}(end)<=(STRFsize/1000)
                        LocalProjection=nan(PaceNb,1);
                        for pp=1:PaceNb
                            LR_Filter=reshape(Logistic.PC_LR(:,vt),Logistic.nf,Logistic.nt);
                            LR_Filter_scaled = LR_Filter(1:length(Model.TickSTRFspectro.fo{mm}),pp:(pp+length(Model.TickSTRFspectro.to{mm})-1));
                            ScalingTerm = ((sum(sum(LR_Filter_scaled.^2)))^0.5)*((sum(sum(Model.Acoustic.B{mm}{1}.^2)))^0.5);
                            LocalProjection(pp) = sum(sum(LR_Filter_scaled.*Model.Acoustic.B{mm}{1}))/ScalingTerm;
                        end
                        LRProjection(mm,vt)=max(LocalProjection);
                    else
                        NbTickBloTh=sum(Model.TickSTRFspectro.to{mm}<=(STRFsize/1000));
                        PaceNb2=length(Model.TickSTRFspectro.to{mm})-NbTickBloTh;
                        LocalProjection=nan(PaceNb,PaceNb2);
                        for pp=1:PaceNb
                            LR_Filter=reshape(Logistic.PC_LR(:,vt),Logistic.nf,Logistic.nt);
                            LR_Filter_scaled = LR_Filter(1:length(Model.TickSTRFspectro.fo{mm}),pp:(pp+NbTickBloTh-1));
                            for pp2=1:PaceNb2
                                Model_scaled=Model.Acoustic.B{mm}(:,pp2:pp2+NbTickBloTh-1);
                                ScalingTerm = ((sum(sum(LR_Filter_scaled.^2)))^0.5)*((sum(sum(Model_scaled.^2)))^0.5);
                                LocalProjection(pp,pp2) = sum(sum(LR_Filter_scaled.*Model_scaled))/ScalingTerm;
                            end
                        end
                        LRProjection(mm,vt)=max(max(LocalProjection));
                    end
                end
            end
            
        end
        if Figure_Switch
            figure(14)
            plot(LRProjection(:,1:7),'-')
            hold on
            plot(LRProjection(:,8:end),':')
            axis([0 8 0 0.5])
            legend(Logistic.vocTypes)
            Xtickposition=get(gca,'XTick');
            set(gca,'XTickLabel', (Xtickposition+1)*10)
            xlabel('Time (ms)')
            ylabel('Correlation to the LR filter')
            hold off
            if Pause_Switch
                pause()
            end
        end
        
        %% Plot of models at each window
        
        if STRF_show
            %% Plot STRF of given windows in the same colorscale
            % Find the extreme values in Acoustic STRFs
            MaxSTRFAc=nan(NumMod,1);
            MinSTRFAc = nan(NumMod,1);
            Gain = nan(NumMod,1);
            for mm=1:NumMod
                
                MaxSTRFAc(mm) = max(max(Model.Acoustic.B{mm}));
                MinSTRFAc(mm) = min(min(Model.Acoustic.B{mm}));
                Gain(mm)=max(abs(MaxSTRFAc(mm)), abs(MinSTRFAc(mm)));
            end
            %CLim = [min(MinSTRFAc(2:2:6)) max(MaxSTRFAc(2:2:6))];
            
            % Plot the gain
            if Figure_Switch
                figure(18)
                plot(Gain, 'k-')
                %axis([0 8 0 0.45])
                Xtickposition=get(gca,'XTick');
                set(gca,'XTickLabel', (Xtickposition+1)*10)
                xlabel('Time (ms)')
                ylabel('Gain of the Acoustic filters')
            end
            for mm=1:NumMod
                if ParamModel.ModelChoice(4) && ParamModel.ModelChoice(5) && DiffDeviance_calc
                    fprintf('The best combine model at %dms is %s\n', Wins(mm),AcSemOptId{mm})
                end
                
                UVOC=unique(Data.VOC{mm});
                X_voc = zeros(length(Data.X_voc_wholeset{mm}), length(UVOC)-1);% Note that all zeros correspond to Ag call
                for vv=1:length(Data.X_voc_wholeset{mm})
                    for uv=2:length(UVOC)
                        if strcmp(Data.X_voc_wholeset{mm}(vv), UVOC{uv})
                            X_voc(vv,uv-1)=1;
                            break
                        end
                    end
                end
                fignum=0;
                if ParamModel.ModelChoice(1)
                    figure(8)
                    imagesc(Model.TickSTRFspectro.to{mm},Model.TickSTRFspectro.fo{mm},Model.Acoustic.B{mm})
                    axis xy
                    title(sprintf('Acoutic lassoglm poisson whole dataset\nlog of lambda=%f',log(Model.Acoustic.Lambdas(mm))))
                    SpectroDim=size(Model.Acoustic.B{mm});
                    if ParamModel.CV
                        ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm}, 1, SpectroDim(1)*SpectroDim(2))./Data.x_std_wholeset{mm},size(Data.x_wholeset{mm},1),1);
                    elseif ParamModel.MeanSubstractSpec
                        ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm}, 1, SpectroDim(1)*SpectroDim(2)),size(Data.x_wholeset{mm},1),1);
                    end
                    Predict_Ac_wholeset = sum(ModelAcoustic_local.*Data.x_wholeset{mm},2) + repmat(Model.Acoustic.B0{mm},size(Data.x_wholeset{mm},1),1);
                    if mm>1
                        close figure 13
                    end
                    figure(13)
                    if exist('fignum','var')
                        fignum=fignum+1;
                    else
                        fignum=1;
                    end
                    if SWITCH.AcOnly==1
                        ss=subplot(1,2,1);
                    else
                        if mod(sum(ParamModel.ModelChoice),2)==0
                            ss=subplot(2,sum(ParamModel.ModelChoice)/2,fignum);
                        else
                            ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,fignum);
                        end
                    end
                    plot(Predict_Ac_wholeset,Data.y_wholeset{mm},'bo','MarkerSize',5)
                    ylabel('Observed Spike count ')
                    xlabel('Predicted Spike count wo output non-linearity')
                    title(sprintf('Acoustic Model Win=%dms',Wins(mm)))
                    AXES=get(ss,'yLim');
                    hold on
                    plot(Predict_Ac_wholeset,exp(Predict_Ac_wholeset),'r.')
                    set(ss,'yLim',AXES)
                    legend('Data','exponential fit')
                    hold off
                    if SWITCH.AcOnly==1
                        ss2=subplot(1,2,2);
                        plot(exp(Predict_Ac_wholeset),Data.y_wholeset{mm},'bo','MarkerSize',10)
                        ylabel('Observed Spike count')
                        xlabel('Predicted Spike count')
                        title(sprintf('Acoustic Model Win=%dms',Wins(mm)))
                        AXES=get(ss,'yLim');
                    end
                end
                
                if ParamModel.ModelChoice(2)
                    figure(9)
                    ModelSemantic = repmat(Model.Semantic.B0{mm},1,length(Model.Semantic.B{mm})+1) + [0 Model.Semantic.B{mm}'];
                    plot(ModelSemantic,1:length(ModelSemantic));
                    title('Semantic Model');
                    set(gca,'YTickLabel', UVOC);
                    ModelSem_local=repmat(Model.Semantic.B{mm}',size(Data.x_wholeset{mm},1),1);
                    Predict_Sem_wholeset = sum(ModelSem_local.*X_voc,2) + repmat(Model.Semantic.B0{mm},size(Data.x_wholeset{mm},1),1);
                    figure(13)
                    if exist('fignum','var')
                        fignum=fignum+1;
                    else
                        fignum=1;
                    end
                    
                    if mod(sum(ParamModel.ModelChoice),2)==0
                        ss=subplot(2,sum(ParamModel.ModelChoice)/2,fignum);
                    else
                        ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,fignum);
                    end
                    
                    plot(Predict_Sem_wholeset,Data.y_wholeset{mm},'bo','MarkerSize',10)
                    ylabel('Observed Spike count wo output non-linearity')
                    xlabel('Predicted Spike count')
                    title(sprintf('Semantic Model Win=%dms',Wins(mm)))
                    AXES=get(ss,'yLim');
                    hold on
                    plot(Predict_Sem_wholeset,exp(Predict_Sem_wholeset),'r.')
                    set(ss,'yLim',AXES)
                    legend('Data','exponential fit')
                    
                end
                if ParamModel.ModelChoice(3)
                    figure(10)
                    subplot(1,2,1)
                    imagesc(Model.TickSTRFspectro.to{mm},Model.TickSTRFspectro.fo{mm},Model.AcSem.Bspectro{mm})
                    axis xy
                    %title(sprintf('OPTIMAL STRF AcSem lassoglm poisson whole dataset\nlog of lambda=%f\n D=%f+/-%f (ratio of sum of LL for all bootstrap)\ncross-validation including %f+/-%f %% of the data',log10(FitInfo_AcSem.Lambda),Deviance.AcSem.mean(mm,aa), Deviance.AcSem.std(mm,aa), PropVal.mean(mm,aa), PropVal.std(mm,aa)))
                    ss=subplot(1,2,2);
                    ModelAcSem = repmat(Model.AcSem.B0{mm},1,length(Model.AcSem.Bsem{mm})+1) + [0 Model.AcSem.Bsem{mm}'];
                    plot(ModelAcSem,1:length(UVOC));
                    title('AcSem Model');
                    set(ss,'YTickLabel', UVOC);
                    if ParamModel.CV
                        ModelAcSem_local=repmat([Model.AcSem.Bsem{mm}' reshape(Model.AcSem.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))./Data.x_std_wholeset{mm}],size(Data.x_wholeset{mm},1),1);
                    elseif ParamModel.MeanSubstractSpec
                        ModelAcSem_local=repmat([Model.AcSem.Bsem{mm}' reshape(Model.AcSem.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))],size(Data.x_wholeset{mm},1),1);
                    end
                    Predict_AcSem_wholeset = sum(ModelAcSem_local.*[X_voc Data.x_wholeset{mm}],2) + repmat(Model.AcSem.B0{mm},size(Data.x_wholeset{mm},1),1);
                    figure(13)
                    if exist('fignum','var')
                        fignum=fignum+1;
                    else
                        fignum=1;
                    end
                    if mod(sum(ParamModel.ModelChoice),2)==0
                        ss=subplot(2,sum(ParamModel.ModelChoice)/2,fignum);
                    else
                        ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,fignum);
                    end
                    plot(Predict_AcSem_wholeset,Data.y_wholeset{mm},'bo','MarkerSize',10)
                    ylabel('Observed Spike count wo output non-linearity')
                    xlabel('Predicted Spike count')
                    title(sprintf('AcSem Model Win=%dms',Wins(mm)))
                    AXES=get(ss,'yLim');
                    hold on
                    plot(Predict_AcSem_wholeset,exp(Predict_AcSem_wholeset),'r.')
                    set(ss,'yLim',AXES)
                    legend('Data','exponential fit')
                    
                end
                if ParamModel.ModelChoice(4) && strcmp(AcSemOptId{mm}, 'AcSemAc')
                    figure(11)
                    subplot(1,2,1)
                    imagesc(Model.TickSTRFspectro.to{mm},Model.TickSTRFspectro.fo{mm},Model.AcSemAc.Bspectro{mm})
                    axis xy
                    %title(sprintf('OPTIMAL STRF AcSem AcOffset lassoglm poisson whole dataset\nlog of lambda=%f',log10(FitInfo_AcSemAc.Lambda)))
                    ss=subplot(1,2,2);
                    ModelAcSemAc = repmat(Model.AcSemAc.B0{mm},1,length(Model.AcSemAc.Bsem{mm})+1) + [0 Model.AcSemAc.Bsem{mm}'];
                    plot(ModelAcSemAc,1:length(UVOC));
                    title('AcSem Model Ac Offset');
                    set(ss,'YTickLabel', UVOC);
                    if ParamModel.CV
                        ModelAcSemAc_local=repmat([Model.AcSemAc.Bsem{mm}' reshape(Model.AcSemAc.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))./Data.x_std_wholeset{mm}],size(Data.x_wholeset{mm},1),1);
                    elseif ParamModel.MeanSubstractSpec
                        ModelAcSemAc_local=repmat([Model.AcSemAc.Bsem{mm}' reshape(Model.AcSemAc.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))],size(Data.x_wholeset{mm},1),1);
                    end
                    Predict_AcSemAc_wholeset = sum(ModelAcSemAc_local.*[X_voc Data.x_wholeset{mm}],2) + repmat(Model.AcSemAc.B0{mm},size(Data.x_wholeset{mm},1),1);
                    figure(13)
                    if exist('fignum','var')
                        fignum=fignum+1;
                    else
                        fignum=1;
                    end
                    if mod(sum(ParamModel.ModelChoice),2)==0
                        ss=subplot(2,sum(ParamModel.ModelChoice)/2,fignum);
                    else
                        ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,fignum);
                    end
                    plot(Predict_AcSemAc_wholeset,Data.y_wholeset{mm},'bo','MarkerSize',10)
                    ylabel('Observed Spike count wo output non-linearity')
                    xlabel('Predicted Spike count')
                    title(sprintf('AcSemAc Model Win=%dms',Wins(mm)))
                    AXES=get(ss,'yLim');
                    hold on
                    plot(Predict_AcSemAc_wholeset,exp(Predict_AcSemAc_wholeset),'r.')
                    set(ss,'yLim',AXES)
                    legend('Data','exponential fit')
                end
                if ParamModel.ModelChoice(5) && strcmp(AcSemOptId{mm}, 'AcSemSem')
                    figure(12)
                    subplot(1,2,1)
                    imagesc(Model.TickSTRFspectro.to{mm},Model.TickSTRFspectro.fo{mm},Model.AcSemSem.Bspectro{mm})
                    axis xy
                    %title(sprintf('OPTIMAL STRF AcSem SemOffset lassoglm poisson whole dataset\nlog of lambda=%f',log10(FitInfo_AcSemSem.Lambda)))
                    ss=subplot(1,2,2);
                    ModelAcSemSem = repmat(Model.AcSemSem.B0{mm},1,length(Model.AcSemSem.Bsem{mm})+1) + [0 Model.AcSemSem.Bsem{mm}'];
                    plot(ModelAcSemSem,1:length(UVOC));
                    title('AcSem Model Sem Offset');
                    set(ss,'YTickLabel', UVOC);
                    if ParamModel.CV
                        ModelAcSemSem_local=repmat([Model.AcSemSem.Bsem{mm}' reshape(Model.AcSemSem.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))./Data.x_std_wholeset{mm}],size(Data.x_wholeset{mm},1),1);
                    elseif ParamModel.MeanSubstractSpec
                        ModelAcSemSem_local=repmat([Model.AcSemSem.Bsem{mm}' reshape(Model.AcSemSem.Bspectro{mm}, 1, SpectroDim(1)*SpectroDim(2))],size(Data.x_wholeset{mm},1),1);
                    end
                    Predict_AcSemSem_wholeset = sum(ModelAcSemSem_local.*[X_voc Data.x_wholeset{mm}],2) + repmat(Model.AcSemSem.B0{mm},size(Data.x_wholeset{mm},1),1);
                    figure(13)
                    if exist('fignum','var')
                        fignum=fignum+1;
                    else
                        fignum=1;
                    end
                    if mod(sum(ParamModel.ModelChoice),2)==0
                        ss=subplot(2,sum(ParamModel.ModelChoice)/2,fignum);
                    else
                        ss=subplot(2,(sum(ParamModel.ModelChoice)+1)/2,fignum);
                    end
                    plot(Predict_AcSemSem_wholeset,Data.y_wholeset{mm},'bo','MarkerSize',10)
                    ylabel('Observed Spike count wo output non-linearity')
                    xlabel('Predicted Spike count')
                    title(sprintf('AcSemSem Model Win=%dms',Wins(mm)))
                    AXES=get(ss,'yLim');
                    hold on
                    plot(Predict_AcSemSem_wholeset,exp(Predict_AcSemSem_wholeset),'r.')
                    set(ss,'yLim',AXES)
                    legend('Data','exponential fit')
                end
                if Pause_Switch
                    pause()
                end
            end
            
            
        end
    end
end


%     %% Plot best STRFs for each alpha each model, each boostrap
%     % Also calculate an R2
%     R2.SS_Acoustic=nan(BootstrapSTRF,NAlphas);
%     R2.SS_AcSem = R2.SS_Acoustic;
%     R2.SS_Semantic = R2.SS_Acoustic;
%     R2.SS_T = R2.SS_Acoustic;
%     Av_Lambda.Acoustic = nan(BootstrapSTRF,NAlphas);
%     Av_Lambda.AcSem = Av_Lambda.Acoustic;
%     Av_Deviance.Acoustic = nan(BootstrapSTRF,NAlphas);
%     Av_Deviance.AcSem=Av_Deviance.Acoustic;
%     MinLambdas.AcSem=nan(BootstrapSTRF,NAlphas);
%     MinLambdas.Acoustic=nan(BootstrapSTRF,NAlphas);
%     
%     
%     for AA=1:1
%         for BB=1:BT
%             if STRF_show==1
%                 % plot Acoustic STRF
%                 figure(10)
%                 valc_Ac = max(abs(max(max(ModelB_Ac{BB,AA}))), abs(min(min(ModelB_Ac{BB,AA}))));
%                 if valc_Ac==0
%                     imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelB_Ac{BB,AA})
%                 else
%                     imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelB_Ac{BB,AA}, [-valc_Ac valc_Ac])
%                 end
%                 axis xy
%                 title(sprintf('Acoustic Alpha=%f %d/%d\nDeviance=%f Biais=%f',Alphas(AA),BB,BT,Deviance.Acoustic.bestvalues{mm}(BB,AA),ModelB0_Ac{BB,AA}));
% 
%                 % Plot Model coefficients AcSem
%                 figure(11)
%                 subplot(1,2,1)
%                 valc_AcSem = max(abs(max(max(ModelBspectro_AcSem{BB,AA}))), abs(min(min(ModelBspectro_AcSem{BB,AA}))));
%                 if valc_AcSem==0
%                     imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelBspectro_AcSem{BB,AA})
%                 else
%                     imagesc(Model.TickSTRFspectro.to{mm}, Model.TickSTRFspectro.to{mm},ModelBspectro_AcSem{BB,AA}, [-valc_AcSem valc_AcSem])
%                 end
%                 axis xy
%                 title(sprintf('AcSem Alpha=%f %d/%d\nDeviance=%f Biais=%f',Alphas(AA),BB,BT,Deviance.AcSem.bestvalues{mm}(BB,AA),ModelB0_AcSem{BB,AA}));
%                 ss=subplot(1,2,2);
%                 plot([ModelB0_AcSem{BB,AA}; ModelBsem_AcSem{BB,AA}] + [0; repmat(ModelB0_AcSem{BB,AA},length(ModelBsem_AcSem{BB,AA}),1)],1:length(UVOC));
%                 title('Coefficients Categories');
%                 set(ss,'YTickLabel', UVOC);
% 
%                 % Plot model coefficients Semantic
%                 figure(12)
%                 plot(ModelB_Sem{BB,AA} + [0 ; repmat(ModelB_Sem{BB,AA}(1),length(ModelB_Sem{BB,AA})-1,1)],1:length(UVOC));
%                 title(sprintf('Sem Model\ndeviance=%f', Deviance.Sem.bestvalues{mm}(BB,AA)));
%                 set(gca,'YTickLabel', UVOC);
% 
%                 % Plot prediction for each model
%                 LocalModelPredict=cell(3,1);
%                 LocalModelPredict{2}=Model.AcSem.ypredictVal{mm,AA}{BB};
%                 LocalModelPredict{1}=Model.Acoustic.ypredictVal{mm,AA}{BB};
%                 LocalModelPredict{3}=Model.Semantic.ypredictVal{mm,AA}{BB};
%                 R2.SS_T(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-repmat(mean(Model.yobsVal{mm,AA}{BB}),length(Model.yobsVal{mm,AA}{BB}),1),2));
%                 R2.SS_Acoustic(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.Acoustic.ypredictVal{mm,AA}{BB},2));
%                 R2.SS_Semantic(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.Semantic.ypredictVal{mm,AA}{BB},2));
%                 R2.SS_AcSem(BB,AA) = sum(power(Model.yobsVal{mm,AA}{BB}-Model.AcSem.ypredictVal{mm,AA}{BB},2));
%                 R2.Acoustic(BB,AA)=1-R2.SS_Acoustic(BB,AA)/R2.SS_T(BB,AA);
%                 R2.Semantic(BB,AA)=1-R2.SS_Semantic(BB,AA)/R2.SS_T(BB,AA);
%                 R2.AcSem(BB,AA)=1-R2.SS_AcSem(BB,AA)/R2.SS_T(BB,AA);
%                 Legend= cell(3,1);
%                 Legend{1}=sprintf('Acoustic only R2=%f',1-R2.SS_Acoustic(BB,AA)/R2.SS_T(BB,AA));
%                 Legend{3}=sprintf('Semantic only R2=%f',1-R2.SS_Semantic(BB,AA)/R2.SS_T(BB,AA));
%                 Legend{2}=sprintf('Acoustic + Semantic R2=%f', 1-R2.SS_AcSem(BB,AA)/R2.SS_T(BB,AA));
%                 MAXX=max(Model.yobsVal{AA}{BB});
%                 MAXY=max([max(LocalModelPredict{1}) max(LocalModelPredict{2}) max(LocalModelPredict{3})]);
%                 MAX=max(MAXX,MAXY);
% 
%                 for jj=1:3
%                     figure(13)
%                     subplot(1,3,jj);
%                     h=gscatter(Model.yobsVal{mm,AA}{BB}, LocalModelPredict{jj}, Model.yobsCat{mm,AA}{BB}, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
%                     title(sprintf('%s', Legend{jj}));
%                     xlabel('observed spike rate (spike per ms)')
%                     ylabel('predicted spike rate (spike per ms)')
%                     axis([0 MAX+MAX/4 0 MAX]);
%                     hold on
%                     plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
%                     hold off
%                 end
%             end
% 
%             % calculate the average best lambda and Deviance over boostrap
%             Av_Deviance.Acoustic(BB,AA)=mean(Deviance.Acoustic.bestvalues{mm}(1:BB,AA));
%             MinLambdas.Acoustic(BB,AA)=Deviance.Acoustic.lambda{mm,AA}{BB}(find(Deviance.Acoustic.values{mm,AA}{BB}==Deviance.Acoustic.bestvalues{mm}(BB,AA)));
%             Av_Lambda.Acoustic(BB,AA)=median(MinLambdas.Acoustic(1:BB,AA));
%             Av_Deviance.AcSem(BB,AA)=mean(Deviance.AcSem.bestvalues{mm}(1:BB,AA));
%             MinLambdas.AcSem(BB,AA)=Deviance.AcSem.lambda{mm,AA}{BB}(find(Deviance.AcSem.values{mm,AA}{BB}==Deviance.AcSem.bestvalues{mm}(BB,AA)));
%             Av_Lambda.AcSem(BB,AA)=median(MinLambdas.AcSem(1:BB,AA));
%             
%             if STRF_show==1
%                 figure(9)
%                 plot(log10(Deviance.Acoustic.lambda{mm,AA}{BB}),Deviance.Acoustic.values{mm,AA}{BB},'k-',log10(Deviance.AcSem.lambda{mm,AA}{BB}),Deviance.AcSem.values{mm,AA}{BB},'r-');
%                 hold on
%                 plot(log10(MinLambdas.Acoustic(BB,AA)),Deviance.Acoustic.bestvalues{mm}(BB,AA),'k*')
%                 hold on
%                 plot(log10(MinLambdas.AcSem(BB,AA)),Deviance.AcSem.bestvalues{mm}(BB,AA),'r*')
%                 hold on
%                 hline(Deviance.Sem.bestvalues{mm}(BB,AA),'g-')
%                 legend('Acoustic', 'AcSem', 'Min Dev Acoustic','Min Dev AcSem')
%                 xlabel('log10 lambda')
%                 ylabel('Deviance')
%                 hold off
%                 pause(1)
%             end
%         end
%     end
%     %% Plot the cumulative average best Deviance and best lambda along bootstrap
%     Color='mgcbrk';
% 
%     figure(14)
%     for AA=1:length(Alphas)
%         subplot(2,2,1)
%         hold on
%         plot(1:BootstrapSTRF,Av_Deviance.Acoustic(:,AA),sprintf('%s-',Color(AA)));
%         subplot(2,2,2)
%         hold on
%         plot(1:BootstrapSTRF,Av_Deviance.AcSem(:,AA),sprintf('%s-',Color(AA)));
%         subplot(2,2,3)
%         hold on
%         plot(1:BootstrapSTRF,log10(Av_Lambda.Acoustic(:,AA)),sprintf('%s-',Color(AA)));
%         subplot(2,2,4)
%         hold on
%         plot(1:BootstrapSTRF,log10(Av_Lambda.AcSem(:,AA)),sprintf('%s-',Color(AA)));
%     end
%     subplot(2,2,1)
%     ylabel('Average Deviance Acoustic Model')
%     xlabel('Booststrap')
%     subplot(2,2,2)
%     ylabel('Average Deviance AcSem Model')
%     xlabel('Booststrap')
%     subplot(2,2,3)
%     ylabel('Median Lambda Acoustic Model (log10)')
%     xlabel('Booststrap')
%     subplot(2,2,4)
%     ylabel('Median Lambda AcSem Model (log10)')
%     xlabel('Booststrap')
%     hold off
%     figure(14)
%     pause()
%     
%     if NAlphas>1
%         %% Plot deviance of Ac vs Acsem + diag and color per alpha see where it's
%         % most often above/below diag chose alpha
%         Color='mgcbrk';
%         extremeBestDeviance=nan(NAlphas,2);
%         legendcolor=cell(NAlphas,2);
%         H=nan(NAlphas,2);
%         for AA=1:NAlphas
%             figure(15)
%             subplot(1,2,1)
%             hold on
%             H(AA,1)=plot(Deviance.AcSem.bestvalues{mm}(:,AA),Deviance.Acoustic.bestvalues{mm}(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
%             extremeBestDeviance(AA,1)=min(min(Deviance.Acoustic.bestvalues{mm}(:,AA)),min(Deviance.AcSem.bestvalues{mm}(:,AA)));
%             extremeBestDeviance(AA,2)=max(max(Deviance.Acoustic.bestvalues{mm}(:,AA)),max(Deviance.AcSem.bestvalues{mm}(:,AA)));
%             legendcolor{AA,1}=num2str(Alphas(AA));
%             subplot(1,2,2)
%             hold on
%             %H(AA,2)=plot(mean(Deviance.AcSem.bestvalues{1}(:,AA)),mean(Deviance.Acoustic.bestvalues{1}(:,AA)),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
%             X=[mean(Deviance.AcSem.bestvalues{mm}(:,AA))-std(Deviance.AcSem.bestvalues{mm}(:,AA)) mean(Deviance.AcSem.bestvalues{mm}(:,AA))+std(Deviance.AcSem.bestvalues{mm}(:,AA))];
%             Y=[mean(Deviance.Acoustic.bestvalues{mm}(:,AA))-std(Deviance.Acoustic.bestvalues{mm}(:,AA)) mean(Deviance.Acoustic.bestvalues{mm}(:,AA))+std(Deviance.Acoustic.bestvalues{mm}(:,AA))];
%             H(AA,2)=line([mean(Deviance.AcSem.bestvalues{mm}(:,AA)) mean(Deviance.AcSem.bestvalues{mm}(:,AA))],Y,'Color',Color(AA));
%             H(AA,2)=line(X,[mean(Deviance.Acoustic.bestvalues{mm}(:,AA)) mean(Deviance.Acoustic.bestvalues{mm}(:,AA))],'Color',Color(AA));
%             legendcolor{AA,2}=sprintf('mean %s',num2str(Alphas(AA)));
%         end
%         subplot(1,2,1)
%         legend(H(:,1),legendcolor{:,1})
%         MAX = max(extremeBestDeviance(:,2));
%         MIN = min(extremeBestDeviance(:,1));
%         axis([MIN MAX+0.1 MIN MAX+0.1])
%         hold on
%         plot(MIN:MAX/10:(MAX+0.1),MIN:MAX/10:(MAX+0.1), 'k');
%         hold off
%         ylabel('Acoustic Model Deviance')
%         xlabel('AcSem Model Deviance')
%         subplot(1,2,2)
%         legend(H(:,2),legendcolor{:,2})
%         MAX = max(extremeBestDeviance(:,2));
%         MIN = min(extremeBestDeviance(:,1));
%         axis([MIN MAX+0.1 MIN MAX+0.1])
%         hold on
%         plot(MIN:MAX/10:(MAX+0.1),MIN:MAX/10:(MAX+0.1), 'k');
%         hold off
%         ylabel('Acoustic Model Deviance')
%         xlabel('AcSem Model Deviance')
%     end
%     
%     %% Correlation between deviance (AIC) and an r2??
%     if NAlphas>1
%         Color='mgcbrk';
%         extremeBestDeviance=nan(NAlphas,2);
%         legendcolor=cell(NAlphas,2);
%         H=nan(NAlphas,3);
%         for AA=1:NAlphas
%             figure(16)
%             subplot(2,2,1)
%             hold on
%             H(AA,1)=plot(Deviance.AcSem.bestvalues{mm}(:,AA),1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
%             %extremeBestDeviance(AA,1)=min(min(Deviance.AcSem.bestvalues{1}(:,AA)),min(Deviance.AcSem.bestvalues{1}(:,AA)));
%             %extremeBestDeviance(AA,2)=max(max(Deviance.AcSem.bestvalues{1}(:,AA)),max(Deviance.AcSem.bestvalues{1}(:,AA)));
%             legendcolor{AA,1}=num2str(Alphas(AA));
%             subplot(2,2,2)
%             hold on
%             H(AA,2)=plot(Deviance.Acoustic.bestvalues{mm}(:,AA),1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
%             %extremeBestDeviance(AA,1)=min(min(Deviance.Acoustic.bestvalues{1}(:,AA)),min(Deviance.AcSem.bestvalues{1}(:,AA)));
%             %extremeBestDeviance(AA,2)=max(max(Deviance.Acoustic.bestvalues{1}(:,AA)),max(Deviance.AcSem.bestvalues{1}(:,AA)));
%             legendcolor{AA,1}=num2str(Alphas(AA));
%             subplot(2,2,4)
%             hold on
%             H(AA,3)=plot(1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),sprintf('%s.',Color(AA)), 'MarkerSize', 20);
%         end
%         subplot(2,2,1)
%         legend(H(:,1),legendcolor{:,1})
%         ylabel('AcSem Model R2')
%         xlabel('AcSem Model Deviance')
%         hold off
%         subplot(2,2,2)
%         legend(H(:,2),legendcolor{:,1})
%         ylabel('Acoustic Model R2')
%         xlabel('Acoustic Model Deviance')
%         hold off
%         subplot(2,2,4)
%         legend(H(:,3),legendcolor{:,1})
%         ylabel('AcSem Model R2')
%         xlabel('Acoustic Model R2')
%         axis([0 1 0 1])
%         line([0 1],[0 1])
%         hold off
% 
%         ss=subplot(2,2,3);
%         MEAN_R2_AcSem = mean(1-R2.SS_AcSem./R2.SS_T,1);
%         MEAN_R2_Acoustic = mean(1-R2.SS_Acoustic./R2.SS_T,1);
%         MEAN_R2_Semantic = mean(1-R2.SS_Semantic./R2.SS_T,1);
%         plot(1:NAlphas,MEAN_R2_AcSem,'r.','MarkerSize',20)
%         hold on
%         plot(MEAN_R2_Acoustic,'k.','MarkerSize',20)
%         hold on
%         plot(MEAN_R2_Semantic,'g.','MarkerSize',20)
%         for AA=1:NAlphas
%             hold on
%             line([AA AA], [MEAN_R2_AcSem(AA)-std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_AcSem(AA)+std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1)],'Color','r');
%             hold on
%             line([AA AA], [MEAN_R2_Acoustic(AA)-std(1-R2.SS_AcSem(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_Acoustic(AA)+std(1-R2.SS_Acoustic(:,AA)./R2.SS_T(:,AA),1)],'Color','k');
%             hold on
%             line([AA AA], [MEAN_R2_Semantic(AA)-std(1-R2.SS_Semantic(:,AA)./R2.SS_T(:,AA),1) MEAN_R2_Semantic(AA)+std(1-R2.SS_Semantic(:,AA)./R2.SS_T(:,AA),1)],'Color','g');
%         end
%         hold on
%         Pl=nan(3,1);
%         Pl(1)=plot(1-sum(R2.SS_AcSem,1)./sum(R2.SS_T,1),'r*','MarkerSize',10);
%         hold on
%         Pl(2)=plot(1-sum(R2.SS_Acoustic,1)./sum(R2.SS_T,1),'k*','MarkerSize',10);
%         hold on
%         Pl(3)=plot(1-sum(R2.SS_Semantic,1)./sum(R2.SS_T,1),'g*','MarkerSize',10);
%         set(ss,'XTick',[1:NAlphas])
%         set(ss,'XTickLabel',Alphas)
%         legend(Pl,'AcSem','Acoustic','Semantic','Location','SouthEast')
%         ylabel('R-square')
%         xlabel('Alphas')
%         hold off
%     end
%     
%     %% Effect of the exp output non-linearity?
%     if NL_show
%         for aa=1:length(Alphas)
%             SpectroDim=size(Model.Acoustic.B{mm,aa});
%             ModelAcoustic_local=repmat(reshape(Model.Acoustic.B{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./x_std,size(Data.x_wholeset{mm},1),1);
%             Predict_Ac_wholeset = sum(ModelAcoustic_local.*Data.x_wholeset{mm},2) + repmat(Model.Acoustic.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
%             figure(17)
%             ss=subplot(2,3,aa);
%             plot(Predict_Ac_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
%             ylabel('Observed Spike count ')
%             xlabel('Predicted Spike count wo output non-linearity')
%             title(sprintf('Acoustic Model Alpha %f',Alphas(aa)))
%             AXES=get(ss,'yLim');
%             hold on
%             plot(Predict_Ac_wholeset,exp(Predict_Ac_wholeset),'r.')
%             set(ss,'yLim',AXES)
%             legend('Data','exponential fit')
% 
%             ModelAcSem_local=repmat([Model.AcSem.Bsem{mm,aa}' reshape(Model.AcSem.Bspectro{mm,aa}, 1, SpectroDim(1)*SpectroDim(2))./x_std],size(Data.x_wholeset{mm},1),1);
%             Predict_AcSem_wholeset = sum(ModelAcSem_local.*[SemBoost.*Data.X_voc_wholeset{mm} Data.x_wholeset{mm}],2) + repmat(Model.AcSem.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
%             figure(18)
%             ss=subplot(2,3,aa);
%             plot(Predict_AcSem_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
%             ylabel('Observed Spike count wo output non-linearity')
%             xlabel('Predicted Spike count')
%             title(sprintf('AcSem Model Alpha %f',Alphas(aa)))
%             AXES=get(ss,'yLim');
%             hold on
%             plot(Predict_AcSem_wholeset,exp(Predict_AcSem_wholeset),'r.')
%             set(ss,'yLim',AXES)
%             legend('Data','exponential fit')
% 
%             ModelSem_local=repmat(Model.Semantic.B{mm,aa}',size(Data.x_wholeset{mm},1),1);
%             Predict_Sem_wholeset = sum(ModelSem_local.*Data.X_voc_wholeset{mm},2) + repmat(Model.Semantic.B0{mm,aa},size(Data.x_wholeset{mm},1),1);
%             figure(19)
%             ss=subplot(2,3,aa);
%             plot(Predict_Sem_wholeset,Data.y_wholeset{mm},'b*','MarkerSize',10)
%             ylabel('Observed Spike count wo output non-linearity')
%             xlabel('Predicted Spike count')
%             title(sprintf('Semantic Model Alpha %f',Alphas(aa)))
%             AXES=get(ss,'yLim');
%             hold on
%             plot(Predict_Sem_wholeset,exp(Predict_Sem_wholeset),'r.')
%             set(ss,'yLim',AXES)
%             legend('Data','exponential fit')
%         end
%         pause()
%     end
% end
% 
% % According to this code, 53.3 would be a good value for lambda
% 10^mean(MinDevlog10Lambdas.Acoustic)
% 
% 

%Interesting Cells: Models_GLMPoissonSite1_L750R1100_e1_s0_ss1.mat
% to a less extent: File Models_GLMPoissonSite2_L2000R1600_e21_s1_ss1.mat
save('DataGatheredModelsPerformance.mat')
save('DataGatheredModelsPerformance.mat', 'ModelInfo', '-append')

SpikeRateAveragepers = SpikeCountAverage(:,:).*50; %Number of spikes per s instead of per window of 20ms
TotalNumCell = sum(~isnan(LLAverage.Acoustic(1,:,1)));

%% Find the cells which don't have weird values of LL
[~,col] = find(LLAverage.Acoustic(:,:,1)<(3*min(min(LLAverage.Floor(:,:,1)))));
[~,col2] = find(~isnan(LLAverage.Acoustic(:,:,1)));
%OkCells = setdiff(unique(col2), unique(col));% Sorting cells by bad LL of acoustic model
OkCells = unique(col2); %getting rid of non calculated cells

%% Figures

figure(20)
for cc=1:length(PoissonFiles)
    if isnan(AC.SignifDiffDev(cc,1)) || isnan(AC.SignifDiffDev(cc,2))
        plot3(AC.FitIndex(cc,1),AC.FitIndex(cc,2), AC.FitIndex(cc,3),'ko','Markersize',10,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor','k')
        % Most likely cells for which models are performing poorly compare
        % to floor
    elseif AC.SignifDiffDev(cc,1) && AC.SignifDiffDev(cc,2)
        plot3(AC.FitIndex(cc,1),AC.FitIndex(cc,2), AC.FitIndex(cc,3),'ko','Markersize',10,'MarkerFaceColor',[1 0.75 0])
        % Synergistic cells for which adding semantic and acoustic
        % parameters to imple models increase significantly the LL
    elseif AC.SignifDiffDev(cc,1) && ~AC.SignifDiffDev(cc,2)
        plot3(AC.FitIndex(cc,1),AC.FitIndex(cc,2), AC.FitIndex(cc,3),'ko','Markersize',10,'MarkerFaceColor','b')
        % Acoustic cells, Semantic parameters don't do much
    elseif ~AC.SignifDiffDev(cc,1) && AC.SignifDiffDev(cc,2)
        plot3(AC.FitIndex(cc,1),AC.FitIndex(cc,2), AC.FitIndex(cc,3),'ko','Markersize',10,'MarkerFaceColor','g')
        % Semantic cells, Acoustic parameters don't do much
    elseif ~AC.SignifDiffDev(cc,1) && ~AC.SignifDiffDev(cc,2)
        plot3(AC.FitIndex(cc,1),AC.FitIndex(cc,2), AC.FitIndex(cc,3),'ko','Markersize',10,'MarkerFaceColor','k')
        % as cells, Acoustic and Semantic parameters seem equally efficient
    end
    hold on
    %pause()
end
hold off
xlabel('FitIndex Acoustic Model')
ylabel('FitIndex Semantic Model')
zlabel('FitIndex AcSem Model')


%% Plot the average profile of LL over cells as a function of time windows
AcSemCells = find(AC.SignifDiffDev(:,6)==1); %Cells for which AcSem doing better than the Floor over windows
AcSemCells = intersect(AcSemCells, OkCells);
ff21=figure(21);
% Average LL of non informative cells (Ceil = Floor Model, over windows)
subplot(2,2,1);
NonInfCells = find(AC.SignifDiffDev(:,4)==0);
NonInfCells = intersect(NonInfCells, OkCells);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,NonInfCells,1),2), nanmean(LLAverage.Acoustic(:,NonInfCells,1),2), nanmean(LLAverage.AcSemOpt(:,NonInfCells,1),2), nanmean(LLAverage.Semantic(:,NonInfCells,1),2), nanmean(LLAverage.Floor(:,NonInfCells,1),2),nanmean(LLAverage.Saturated(:,NonInfCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,NonInfCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'NonInfCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(1,-0.05,sprintf('# NonInformative Cells=%d %0.1f%%',length(NonInfCells),100*length(NonInfCells)/TotalNumCell)) 
text(5,-0.1,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,NonInfCells)),100*length(intersect(AcSemCells,NonInfCells))/TotalNumCell)) 
hold off

% Average LL of Acoustic cells (Ceil > Floor Model, Ac>Floor, Sem=Floor over windows)
subplot(2,2,2);
AcousticCells = intersect(find(AC.SignifDiffDev(:,5)==1), find(AC.SignifDiffDev(:,7)==0));
AcousticCells = intersect(AcousticCells,OkCells);
AcousticCells = setdiff(AcousticCells, NonInfCells);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,AcousticCells,1),2), nanmean(LLAverage.Acoustic(:,AcousticCells,1),2), nanmean(LLAverage.AcSemOpt(:,AcousticCells,1),2), nanmean(LLAverage.Semantic(:,AcousticCells,1),2), nanmean(LLAverage.Floor(:,AcousticCells,1),2),nanmean(LLAverage.Saturated(:,AcousticCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,AcousticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'AcousticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(3,-0.1,sprintf('# Acoustic Cells=%d %0.1f%%',length(AcousticCells),100*length(AcousticCells)/TotalNumCell))
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,AcousticCells)),100*length(intersect(AcSemCells,AcousticCells))/TotalNumCell)) 
hold off

% Average LL of Semantic cells (Ceil > Floor Model, Ac=Floor, Sem>Floor over windows)
subplot(2,2,3);
SemanticCells = intersect(find(AC.SignifDiffDev(:,5)==0), find(AC.SignifDiffDev(:,7)==1));
SemanticCells = intersect(SemanticCells, OkCells);
SemanticCells = setdiff(SemanticCells, NonInfCells);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,SemanticCells,1),2), nanmean(LLAverage.Acoustic(:,SemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,SemanticCells,1),2), nanmean(LLAverage.Semantic(:,SemanticCells,1),2), nanmean(LLAverage.Floor(:,SemanticCells,1),2),nanmean(LLAverage.Saturated(:,SemanticCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SemanticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(3,-0.1,sprintf('# Semantic Cells=%d %0.1f%%',length(SemanticCells),100*length(SemanticCells)/TotalNumCell)) 
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,SemanticCells)),100*length(intersect(AcSemCells,SemanticCells))/TotalNumCell)) 
hold off

% Average LL of AcousticSemantic cells (Ceil > Floor Model, Ac>Floor, Sem>Floor over windows)
subplot(2,2,4);
AcousticSemanticCells = intersect(find(AC.SignifDiffDev(:,5)==1), find(AC.SignifDiffDev(:,7)==1));
AcousticSemanticCells = intersect(AcousticSemanticCells, OkCells);
AcousticSemanticCells = setdiff(AcousticSemanticCells, NonInfCells);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Acoustic(:,AcousticSemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Semantic(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Floor(:,AcousticSemanticCells,1),2),nanmean(LLAverage.Saturated(:,AcousticSemanticCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,AcousticSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'AcousticSemanticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)

text(3,-0.1,sprintf('# AcousticSemantic Cells=%d %0.1f%%',length(AcousticSemanticCells),100*length(AcousticSemanticCells)/TotalNumCell)) 
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,AcousticSemanticCells)),100*length(intersect(AcSemCells,AcousticSemanticCells))/TotalNumCell)) 
hold off

%% Study the relationship between LL and spike count
SpikeCount = repmat(1:10,10,1);
LL_test = nan(10,1);
for tt=1:10
    LL_test(tt)=LL_Calculus(SpikeCount(:,tt), SpikeCount(:,tt));
end
figure(50)
subplot(1,2,1)
plot(SpikeCount(1,:),LL_test)
xlabel('Spike Count (fake data)')
ylabel('LogLikelihood saturated (fake data)')
subplot(1,2,2)
SPA_vector=reshape(SpikeRateAveragepers(:,:,1),size(SpikeRateAveragepers,1) * size(SpikeRateAveragepers,2),1);
LLSat_vector = reshape(LLAverage.Saturated(:,:,1),size(LLAverage.Saturated,1) * size(LLAverage.Saturated,2),1);
plot(SPA_vector, LLSat_vector, '.')
xlabel('Data average spike count per cell per Win')
ylabel('Data LogLikelihood Saturated per Stim per trial')

%% Let's substract the LL of the saturated model to see what is left and plot 21 again
figure(22)
% Average LL of non informative cells (Ceil = Floor Model, over windows)
subplot(2,2,1);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,NonInfCells,1),2), nanmean(LLAverage.Acoustic(:,NonInfCells,1),2), nanmean(LLAverage.AcSemOpt(:,NonInfCells,1),2), nanmean(LLAverage.Semantic(:,NonInfCells,1),2), nanmean(LLAverage.Floor(:,NonInfCells,1),2)];
SatMatrix = repmat(nanmean(LLAverage.Saturated(:,NonInfCells,1),2),1,5);
LL_MatrixPlot_y_left = LL_MatrixPlot_y_left - SatMatrix;
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,5);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,NonInfCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'NonInfCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(1,-0.05,sprintf('# NonInformative Cells=%d %0.1f%%',length(NonInfCells),100*length(NonInfCells)/TotalNumCell)) 
text(5,-0.1,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,NonInfCells)),100*length(intersect(AcSemCells,NonInfCells))/TotalNumCell)) 
hold off

% Average LL of Acoustic cells (Ceil > Floor Model, Ac>Floor, Sem=Floor over windows)
subplot(2,2,2);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,AcousticCells,1),2), nanmean(LLAverage.Acoustic(:,AcousticCells,1),2), nanmean(LLAverage.AcSemOpt(:,AcousticCells,1),2), nanmean(LLAverage.Semantic(:,AcousticCells,1),2), nanmean(LLAverage.Floor(:,AcousticCells,1),2)];
SatMatrix = repmat(nanmean(LLAverage.Saturated(:,AcousticCells,1),2),1,5);
LL_MatrixPlot_y_left = LL_MatrixPlot_y_left - SatMatrix;
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,5);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,AcousticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'AcousticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(3,-0.1,sprintf('# Acoustic Cells=%d %0.1f%%',length(AcousticCells),100*length(AcousticCells)/TotalNumCell))
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,AcousticCells)),100*length(intersect(AcSemCells,AcousticCells))/TotalNumCell)) 
hold off

% Average LL of Semantic cells (Ceil > Floor Model, Ac=Floor, Sem>Floor over windows)
subplot(2,2,3);
LL_MatrixPlot_y_left = [ nanmean(LLAverage.CeilingOpt(:,SemanticCells,1),2), nanmean(LLAverage.Acoustic(:,SemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,SemanticCells,1),2), nanmean(LLAverage.Semantic(:,SemanticCells,1),2), nanmean(LLAverage.Floor(:,SemanticCells,1),2)];
SatMatrix = repmat(nanmean(LLAverage.Saturated(:,SemanticCells,1),2),1,5);
LL_MatrixPlot_y_left = LL_MatrixPlot_y_left - SatMatrix;
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,5);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SemanticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
text(3,-0.1,sprintf('# Semantic Cells=%d %0.1f%%',length(SemanticCells),100*length(SemanticCells)/TotalNumCell)) 
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,SemanticCells)),100*length(intersect(AcSemCells,SemanticCells))/TotalNumCell)) 
hold off

% Average LL of AcousticSemantic cells (Ceil > Floor Model, Ac>Floor, Sem>Floor over windows)
subplot(2,2,4);
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Acoustic(:,AcousticSemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Semantic(:,AcousticSemanticCells,1),2), nanmean(LLAverage.Floor(:,AcousticSemanticCells,1),2)];
SatMatrix = repmat(nanmean(LLAverage.Saturated(:,AcousticSemanticCells,1),2),1,5);
LL_MatrixPlot_y_left = LL_MatrixPlot_y_left - SatMatrix;
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,5);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,AcousticSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'AcousticSemanticCells';
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)

text(3,-0.1,sprintf('# AcousticSemantic Cells=%d %0.1f%%',length(AcousticSemanticCells),100*length(AcousticSemanticCells)/TotalNumCell)) 
text(3,-0.2,sprintf('# AcSem Cells among those=%d %0.1f%% of total',length(intersect(AcSemCells,AcousticSemanticCells)),100*length(intersect(AcSemCells,AcousticSemanticCells))/TotalNumCell)) 
hold off

%% Verify the average values of LL difference in each cell group
% High vs low diff of Ceil vs Floor cells
figure(23)
sp1= subplot(2,2,1);
InfCells = find(AC.SignifDiffDev(1:ff,4)==1);
InfCells = intersect(InfCells,OkCells);
plot(1,AC.MeanDiffDev(NonInfCells,4),'.k', 'MarkerSize',10)
hold on
plot(1, mean(AC.MeanDiffDev(NonInfCells,4)),'+r','MarkerSize',10)
hold on
plot(2, AC.MeanDiffDev(InfCells,4), '.k','MarkerSize',10)
hold on
plot(2, mean(AC.MeanDiffDev(InfCells,4)),'+r','MarkerSize',10)
hold off
xlim([0 3])
hline(0)
set(sp1,'XTickLabel', {' ', 'NonSignifCells', 'SignifCells', ' '})
ylabel('LLCeil - LLFloor')

sp2= subplot(2,2,2);
AcCells = intersect(find(AC.SignifDiffDev(1:ff,5)==1), InfCells);
NonAcCells = intersect(find(AC.SignifDiffDev(1:ff,5)==0),InfCells);
plot(1,AC.MeanDiffDev(NonAcCells,5),'.k','MarkerSize',10)
hold on
plot(1, mean(AC.MeanDiffDev(NonAcCells,5)),'+r','MarkerSize',10)
hold on
plot(2, AC.MeanDiffDev(AcCells,5), '.k','MarkerSize',10)
hold on
plot(2, mean(AC.MeanDiffDev(AcCells,5)),'+r','MarkerSize',10)
hold off
xlim([0 3])
%ylim([-2 8])
set(sp2,'XTickLabel', {' ', 'NonSignifCells', 'SignifCells', ' '})
ylabel('LLAc - LLFloor')
hline(0)

sp3= subplot(2,2,3);
SemCells = intersect(find(AC.SignifDiffDev(1:ff,7)==1), InfCells);
NonSemCells = intersect(find(AC.SignifDiffDev(1:ff,7)==0),InfCells);
plot(1,AC.MeanDiffDev(NonSemCells,7),'.k','MarkerSize',10)
hold on
plot(1, mean(AC.MeanDiffDev(NonSemCells,7)),'+r','MarkerSize',10)
hold on
plot(2, AC.MeanDiffDev(SemCells,7), '.k','MarkerSize',10)
hold on
plot(2, mean(AC.MeanDiffDev(SemCells,7)),'+r','MarkerSize',10)
hold off
xlim([0 3])
set(sp3,'XTickLabel', {' ', 'NonSignifCells', 'SignifCells', ' '})
ylabel('LLSem - LLFloor')
hline(0)

sp4= subplot(2,2,4);
AcSemSigCells = intersect(find(AC.SignifDiffDev(1:ff,6)==1), InfCells);
NonAcSemSigCells = intersect(find(AC.SignifDiffDev(1:ff,6)==0),InfCells);
plot(1,AC.MeanDiffDev(NonAcSemSigCells,6),'.k','MarkerSize',10)
hold on
plot(1, mean(AC.MeanDiffDev(NonAcSemSigCells,6)),'+r','MarkerSize',10)
hold on
plot(2, AC.MeanDiffDev(AcSemSigCells,6), '.k','MarkerSize',10)
hold on
plot(2, mean(AC.MeanDiffDev(AcSemSigCells,6)),'+r','MarkerSize',10)
hold off
xlim([0 3])
set(sp4,'XTickLabel', {' ', 'NonSignifCells', 'SignifCells', ' '})
ylabel('LLAcSem - LLFloor')
hline(0)

%% Plot the difference of log likelihood over time for each cell type
figure(24)
sp1=subplot(2,2,1);
LLAverage.CeilMinusFloor = LLAverage.CeilingOpt(:,:,1)-LLAverage.Floor(:,:,1);
for uu=1:size(NonInfCells)
    hold on
    plot(LLAverage.CeilMinusFloor(:,NonInfCells(uu)),'.-','Color',[0 0 0 0.2])
end
hold on
for uu=1:size(InfCells)
    hold on
    plot(LLAverage.CeilMinusFloor(:,InfCells(uu)),'.-','Color',[1 0 0 0.2])
end
hold on
plot(nanmean(LLAverage.CeilMinusFloor(:,NonInfCells,1),2),'.-k')
hold on
plot(nanmean(LLAverage.CeilMinusFloor(:,InfCells,1),2),'.-r')
xlabel('windows (ms)')
ylabel('Ceiling - Floor LogLikelihood')
fprintf('Please enlarge figure 24 so the X axis labels can plot correctly')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
hold off

sp2=subplot(2,2,2);
LLAverage.AcMinusFloor = LLAverage.Acoustic(:,:,1)-LLAverage.Floor(:,:,1);
for uu=1:size(NonAcCells)
    hold on
    plot(LLAverage.AcMinusFloor(:,NonAcCells(uu)),'.-','Color',[0.5 0.5 0.5])
end
hold on
for uu=1:size(AcCells)
    hold on
    plot(LLAverage.AcMinusFloor(:,AcCells(uu)),'.-','Color',[0.5 0 0])
end
hold on
plot(nanmean(LLAverage.AcMinusFloor(:,NonAcCells,1),2),'.-k')
hold on
plot(nanmean(LLAverage.AcMinusFloor(:,AcCells,1),2),'.-r')
xlabel('windows (ms)')
ylabel('Ac - Floor LogLikelihood')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
ylim([-20 15])
hold off


sp3=subplot(2,2,3);
LLAverage.SemMinusFloor = LLAverage.Semantic(:,:,1)-LLAverage.Floor(:,:,1);
for uu=1:size(NonSemCells)
    hold on
    plot(LLAverage.SemMinusFloor(:,NonSemCells(uu)),'.-','Color',[0.5 0.5 0.5])
end
hold on
for uu=1:size(SemCells)
    hold on
    plot(LLAverage.SemMinusFloor(:,SemCells(uu)),'.-','Color',[0.5 0 0])
end
hold on
plot(nanmean(LLAverage.SemMinusFloor(:,SemCells,1),2),'.-r')
hold on
plot(nanmean(LLAverage.SemMinusFloor(:,NonSemCells,1),2),'.-k')
xlabel('windows (ms)')
ylabel('Semantic - Floor LogLikelihood')
Xtickposition=get(sp3,'XTick');
set(sp3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
hold off

sp4=subplot(2,2,4);
LLAverage.AcSemMinusFloor = LLAverage.AcSemOpt(:,:,1)-LLAverage.Floor(:,:,1);
for uu=1:size(NonAcSemSigCells)
    hold on
    plot(LLAverage.AcSemMinusFloor(:,NonAcSemSigCells(uu)),'.-','Color',[0.5 0.5 0.5])
end
hold on
for uu=1:size(AcSemSigCells)
    hold on
    plot(LLAverage.AcSemMinusFloor(:,AcSemSigCells(uu)),'.-','Color',[0.5 0 0])
end
hold on
plot(nanmean(LLAverage.AcSemMinusFloor(:,NonAcSemSigCells,1),2),'.-k')
hold on
plot(nanmean(LLAverage.AcSemMinusFloor(:,AcSemSigCells,1),2),'.-r')
xlabel('windows (ms)')
ylabel('AcSem - Floor LogLikelihood')
Xtickposition=get(sp4,'XTick');
set(sp4,'XTickLabel', [0 Wins(Xtickposition(2:end))])
hold off
ylim([-20 15])

%% Piece of figure to work on
sp3=subplot(3,1,3);
LLAverage.AcSemMinusAc = LLAverage.AcSemOpt(:,:,1)-LLAverage.Acoustic(:,:,1);
for uu=1:size(LLAverage.AcSemMinusAc,2)
    hold on
    plot(LLAverage.AcSemMinusAc(:,uu),'.-','Color',[0 0.3 0.5])
    pause()
end
hold on
plot(nanmean(LLAverage.CeilMinusFloor(:,:,1),2),'.-k')
xlabel('windows (ms)')
ylabel('AcSem - Ac LogLikelihood')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
ylim([0 15])
hold off   

%% Plot the number of cells where ceiling>Floor (auditory cells for which the type of vocalization played actually somewhat predict neural activity)
figure(25)
plot(sum(Signif.CeilvsFl,1), '.-c')
hold on
plot(sum(Signif.AcvsFl,1), '.-b')
hold on
plot(sum(Signif.AcSemvsFl,1), '.-r')
hold on
plot(sum(Signif.AcSemvsAc,1), '.-g')
xlabel('windows (ms)')
ylabel('Number of cells')
Xtickposition=get(gca,'XTick');
set(gca,'XTickLabel', [0 Wins(Xtickposition(2:end))])
legend('Ceiling vs Floor', 'Acoustic vs Floor', 'AcSem vs Floor', 'AcSem vs Ac', 'Location','NorthEast')

%% Plot the evolution of cell best model over time
figure(26)
% NonInfCells first
Mat4NonInfCells = plotCellPerfAsMat(Signif, NonInfCells,NumMod);
sp1=subplot(2,2,1);
image(Mat4NonInfCells)
ylabel('Ceil=Floor Cells')
xlabel('window (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])


% Acoustic Cells
Mat4AcousticCells = plotCellPerfAsMat(Signif, AcousticCells,NumMod);
sp2=subplot(2,2,2);
image(Mat4AcousticCells)
ylabel('Acoustic Cells')
xlabel('window (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])

% Semantic Cells
Mat4SemanticCells = plotCellPerfAsMat(Signif, SemanticCells,NumMod);
sp3=subplot(2,2,3);
image(Mat4SemanticCells)
ylabel('Semantic Cells')
xlabel('window (ms)')
Xtickposition=get(sp3,'XTick');
set(sp3,'XTickLabel', [0 Wins(Xtickposition(2:end))])

% Acoustic and Semantic Cells
Mat4AcousticSemanticCells = plotCellPerfAsMat(Signif, AcousticSemanticCells,NumMod);
sp4=subplot(2,2,4);
image(Mat4AcousticSemanticCells)
ylabel('Acoustic Semantic Cells')
xlabel('window (ms)')
Xtickposition=get(sp4,'XTick');
set(sp4,'XTickLabel', [0 Wins(Xtickposition(2:end))])

%% Plot difference of LL between Floor and Ceiling over time for all cells odered in descending value of average difference of LL
[~,IndCeil] = sort(AC.MeanDiffDev(OkCells,4));
[~,IndSem] = sort(AC.MeanDiffDev(OkCells,7));
IndSem = OkCells(IndSem);
IndCeil = OkCells(IndCeil);
figure(27)
sp1=subplot(1,4,1)
imagesc(LLAverage.CeilMinusFloor(:,flip(IndCeil))')
colorbar()
ylabel('Cells  increasing Ceiling-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])

sp2=subplot(1,4,2)
imagesc(LLAverage.SemMinusFloor(:,flip(IndCeil))')
colorbar()
ylabel('Cells increasing Ceiling-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])


sp3=subplot(1,4,3)
imagesc(LLAverage.CeilMinusFloor(:,flip(IndSem))')
colorbar()
ylabel('Cells increasing Semantic-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp3,'XTick');
set(sp3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])


sp4=subplot(1,4,4)
imagesc(LLAverage.SemMinusFloor(:,flip(IndSem))')
colorbar()
ylabel('Cells increasing Semantic-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp4,'XTick');
set(sp4,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])

%% Plot difference of LL between Floor and Ceiling over time for all cells odered in descending value of average difference of LL
LLAverage.AcSemMinusAc=LLAverage.AcSemOpt(:,:,1)-LLAverage.Acoustic(:,:,1);
[~,IndCeil] = sort(AC.MeanDiffDev(OkCells,4));
IndCeil = OkCells(IndCeil);
figure(28)
sp1=subplot(3,4,1);
imagesc(LLAverage.CeilMinusFloor(:,flip(IndCeil))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])

sp2=subplot(3,4,2);
imagesc(LLAverage.SemMinusFloor(:,flip(IndCeil))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])

sp3=subplot(3,4,3);
imagesc(LLAverage.AcMinusFloor(:,flip(IndCeil))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp3,'XTick');
set(sp3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Acoustic - LL Null Model');
[cmin,cmax]=caxis;
caxis([-cmax cmax])


sp4=subplot(3,4,4);
LLAverage.AcSemMinusAc=LLAverage.AcSemOpt(:,:,1)-LLAverage.Acoustic(:,:,1);
imagesc(LLAverage.AcSemMinusAc(:,flip(IndCeil))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp4,'XTick');
set(sp4,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL AcSem - LL Acoustic');
caxis([-1 1])

% Semantic LL order
[~,IndSem] = sort(AC.MeanDiffDev(OkCells,7));
IndSem = OkCells(IndSem);
sp5=subplot(3,4,5);
imagesc(LLAverage.CeilMinusFloor(:,flip(IndSem))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp5,'XTick');
set(sp5,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])

sp6=subplot(3,4,6);
imagesc(LLAverage.SemMinusFloor(:,flip(IndSem))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp6,'XTick');
set(sp6,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])

sp7=subplot(3,4,7);
imagesc(LLAverage.AcMinusFloor(:,flip(IndSem))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp7,'XTick');
set(sp7,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Acoustic - LL Null Model');
[cmin,cmax]=caxis;
caxis([-cmax cmax])


sp8=subplot(3,4,8);
imagesc(LLAverage.AcSemMinusAc(:,flip(IndSem))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp8,'XTick');
set(sp8,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL AcSem - LL Acoustic');
caxis([-1 1])

% Acoustic order
[~,IndAc] = sort(AC.MeanDiffDev(OkCells,5));
IndAc = OkCells(IndAc);
sp9=subplot(3,4,9);
imagesc(LLAverage.CeilMinusFloor(:,flip(IndAc))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp9,'XTick');
set(sp9,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])

sp10=subplot(3,4,10);
imagesc(LLAverage.SemMinusFloor(:,flip(IndAc))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp10,'XTick');
set(sp10,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])

sp11=subplot(3,4,11);
imagesc(LLAverage.AcMinusFloor(:,flip(IndAc))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp11,'XTick');
set(sp11,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Acoustic - LL Null Model');
[cmin,cmax]=caxis;
caxis([-cmax cmax])


sp12=subplot(3,4,12);
imagesc(LLAverage.AcSemMinusAc(:,flip(IndAc))')
colorbar()
ylabel('Cells')
xlabel('windows (ms)')
Xtickposition=get(sp12,'XTick');
set(sp12,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL AcSem - LL Acoustic');
caxis([-1 1])

%% Plot the best Semantic cells and the worst semantic cells Average LL
% Average LL of super Semantic cells (Sem >> Floor Model, over windows)
figure(29)
subplot(1,3,1);
%LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Acoustic(:,SuperSemCells,1),2), nanmean(LLAverage.AcSemOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Semantic(:,SuperSemCells,1),2), nanmean(LLAverage.Floor(:,SuperSemCells,1),2),nanmean(LLAverage.Saturated(:,SuperSemCells,1),2)];
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,OkCells,1),2), nanmean(LLAverage.Semantic(:,OkCells,1),2), nanmean(LLAverage.Floor(:,OkCells,1),2),nanmean(LLAverage.Saturated(:,OkCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,4);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,OkCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'Cells';
myplotyyLLSPAsemOnly(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off

SuperSemCells = IndSem((end-19) : end);
subplot(1,3,2);
%LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Acoustic(:,SuperSemCells,1),2), nanmean(LLAverage.AcSemOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Semantic(:,SuperSemCells,1),2), nanmean(LLAverage.Floor(:,SuperSemCells,1),2),nanmean(LLAverage.Saturated(:,SuperSemCells,1),2)];
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Semantic(:,SuperSemCells,1),2), nanmean(LLAverage.Floor(:,SuperSemCells,1),2),nanmean(LLAverage.Saturated(:,SuperSemCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,4);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SuperSemCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SuperSemCells';
myplotyyLLSPAsemOnly(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off


% Average LL of Worst semantic cells (Semantic = Floor Model over windows)
subplot(1,3,3);
WorstSemanticCells = IndSem(1:20);
%LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,WorstSemanticCells,1),2), nanmean(LLAverage.Acoustic(:,WorstSemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,WorstSemanticCells,1),2), nanmean(LLAverage.Semantic(:,WorstSemanticCells,1),2), nanmean(LLAverage.Floor(:,WorstSemanticCells,1),2),nanmean(LLAverage.Saturated(:,WorstSemanticCells,1),2)];
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,WorstSemanticCells,1),2), nanmean(LLAverage.Semantic(:,WorstSemanticCells,1),2), nanmean(LLAverage.Floor(:,WorstSemanticCells,1),2),nanmean(LLAverage.Saturated(:,WorstSemanticCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,4);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,WorstSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'WorstSemanticCells';
myplotyyLLSPAsemOnly(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off


%% Plot the best Semantic cells and the worst semantic cells Info
% Average Info of all cells (Sem >> Floor Model, over windows)
figure(30)
subplot(1,3,1);
%Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,OkCells),2), nanmean(ModelInfo.Acoustic(:,OkCells),2), nanmean(ModelInfo.AcSemOpt(:,OkCells),2), nanmean(ModelInfo.Semantic(:,OkCells),2), nanmean(ModelInfo.Floor(:,OkCells),2)];
Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,OkCells),2), nanmean(ModelInfo.Semantic(:,OkCells),2), nanmean(ModelInfo.Floor(:,OkCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,OkCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'AllCells';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off


% Average Info of super Semantic cells (Sem >> Floor Model, over windows)
figure(30)
SuperSemCells = IndSem((end-19) : end);
subplot(1,3,2);
%Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,SuperSemCells),2), nanmean(ModelInfo.Acoustic(:,SuperSemCells),2), nanmean(ModelInfo.AcSemOpt(:,SuperSemCells),2), nanmean(ModelInfo.Semantic(:,SuperSemCells),2), nanmean(ModelInfo.Floor(:,SuperSemCells),2)];
Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,SuperSemCells),2), nanmean(ModelInfo.Semantic(:,SuperSemCells),2), nanmean(ModelInfo.Floor(:,SuperSemCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SuperSemCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SuperSemCells';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off


% Average Info of Worst semantic cells (Semantic = Floor Model over windows)
subplot(1,3,3);
WorstSemanticCells = IndSem(1:20);
%Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,WorstSemanticCells),2), nanmean(ModelInfo.Acoustic(:,WorstSemanticCells),2), nanmean(ModelInfo.AcSemOpt(:,WorstSemanticCells),2), nanmean(ModelInfo.Semantic(:,WorstSemanticCells),2), nanmean(ModelInfo.Floor(:,WorstSemanticCells),2)];
Info_MatrixPlot_y_left = [nanmean(ModelInfo.Ceiling(:,WorstSemanticCells),2), nanmean(ModelInfo.Semantic(:,WorstSemanticCells),2), nanmean(ModelInfo.Floor(:,WorstSemanticCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,WorstSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'WorstSemanticCells';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
hold off

%% Plot the info per spike for each individual cell for ceiling, floor and semantic models

[~,colzero]=find(~SpikeCountAverage(:,OkCells));
ZeroRateCells = OkCells(unique(colzero));
OkRateCells = setdiff(OkCells,ZeroRateCells);
ModelInfoperspike.Ceiling =  ModelInfo.Ceiling ./ SpikeCountAverage;
ModelInfoperspike.Floor =  ModelInfo.Floor ./ SpikeCountAverage;
ModelInfoperspike.Semantic =  ModelInfo.Semantic ./ SpikeCountAverage;

figure(31)
% All cells
subplot(1,3,1)
Info_MatrixPlot_y_left = [nanmean(ModelInfoperspike.Ceiling(:,OkRateCells),2), nanmean(ModelInfoperspike.Semantic(:,OkRateCells),2), nanmean(ModelInfoperspike.Floor(:,OkRateCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,OkRateCells,1),2);
SPA_Plot_x_right = 1:NumMod;
LegendInfo.CellType = 'AllCells';
LegendInfo.YleftAxis = 'Information per spike (bits/spike)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 0.7])
hold off


% SuperSemCells
subplot(1,3,2)
Info_MatrixPlot_y_left = [nanmean(ModelInfoperspike.Ceiling(:,SuperSemCells),2), nanmean(ModelInfoperspike.Semantic(:,SuperSemCells),2), nanmean(ModelInfoperspike.Floor(:,SuperSemCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SuperSemCells,1),2);
SPA_Plot_x_right = 1:NumMod;
LegendInfo.CellType = 'SuperSemCells';
LegendInfo.YleftAxis = 'Information per spike (bits/spike)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 0.7])
hold off


% WorstSemCells
subplot(1,3,3)
Info_MatrixPlot_y_left = [nanmean(ModelInfoperspike.Ceiling(:,WorstSemanticCells),2), nanmean(ModelInfoperspike.Semantic(:,WorstSemanticCells),2), nanmean(ModelInfoperspike.Floor(:,WorstSemanticCells),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,WorstSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
LegendInfo.CellType = 'WorstSemanticCells';
LegendInfo.YleftAxis = 'Information per spike (bits/spike)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 0.7])
hold off



%% Plot the proportion of info per spike for each individual cell for semantic model

PropInfoperspike.Semantic =  (ModelInfoperspike.Semantic -  ModelInfoperspike.Floor)./ (ModelInfoperspike.Ceiling -  ModelInfoperspike.Floor);

figure(32)
% All cells
sp1=subplot(1,3,1);
plot(1:NumMod, PropInfoperspike.Semantic(:,OkRateCells),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, median(PropInfoperspike.Semantic(:,OkRateCells),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('Proportion of Semantic Information')
xlabel('Time (ms)')
fprintf('Enlarge figure')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('All Cells')



% SuperSemCells
sp2=subplot(1,3,2);
plot(1:NumMod, PropInfoperspike.Semantic(:,SuperSemCells),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, median(PropInfoperspike.Semantic(:,SuperSemCells),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('Proportion of Semantic Information')
xlabel('Time (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Super Semantic Cells')

% WorstSemCells
sp3=subplot(1,3,3);
plot(1:NumMod, PropInfoperspike.Semantic(:,WorstSemanticCells),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, median(PropInfoperspike.Semantic(:,WorstSemanticCells),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('Proportion of Semantic Information')
xlabel('Time (ms)')
Xtickposition=get(sp3,'XTick');
set(sp3,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Worst Semantic Cells')

%% Matrices of information 
[~,IndCeilInf] = sort(nanmean(ModelInfo.Ceiling(:,OkCells),1));
[~,IndSemInf] = sort(nanmean(ModelInfo.Semantic(:,OkCells),1));
IndSemInf = OkCells(IndSemInf);
IndCeilInf = OkCells(IndCeilInf);
ModelInfo.CeilMinusFloor = ModelInfo.Ceiling - ModelInfo.Floor;
ModelInfo.SemMinusFloor = ModelInfo.Semantic - ModelInfo.Floor;
IndCeilOkRate = IndCeilInf;
IndCeilOkRate(find(IndCeilInf == ZeroRateCells))=[];
IndSemOkRate = IndSemInf;
IndSemOkRate(find(IndSemInf == ZeroRateCells))=[];
OkCells2 = setdiff(OkCells, ZeroRateCells);

figure(33)
sp1=subplot(1,3,1);
imagesc(ModelInfo.Ceiling(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells  increasing Info Ceiling Model')
xlabel('windows (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Ceiling model (bits)');
caxis([0 0.7])

sp2=subplot(1,3,2);
imagesc(ModelInfo.Floor(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells increasing Info Ceiling Model')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Null Model (bits)');
caxis([0 0.7])

sp2=subplot(1,3,3);
imagesc(SpikeRateAveragepers(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells increasing Info Ceiling Model')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Spike Rate (/s)');
%caxis([0 0.5])

%% Fig 34 Average values of info (both total and per spike) over cells along time and spike rate
figure(34)
subplot(1,2,1)
Info_MatrixPlot_y_left = [nanmedian(ModelInfo.Ceiling(:,OkCells2),2), nanmedian(ModelInfo.Semantic(:,OkCells2),2), nanmedian(ModelInfo.Floor(:,OkCells2),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,OkCells2,1),2);
SPA_Plot_x_right = 1:NumMod;
LegendInfo.CellType = 'AllCells';
LegendInfo.YleftAxis = 'Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins, [-0.03 0.25])
hold off

subplot(1,2,2)
Info_MatrixPlot_y_left = [nanmedian(ModelInfoperspike.Ceiling(:,OkCells2),2), nanmedian(ModelInfoperspike.Semantic(:,OkCells2),2), nanmedian(ModelInfoperspike.Floor(:,OkCells2),2)];
Info_MatrixPlot_x_left = repmat((1:NumMod)',1,3);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,OkCells2,1),2);
SPA_Plot_x_right = 1:NumMod;
LegendInfo.CellType = 'AllCells';
LegendInfo.YleftAxis = 'Information (bits/spike)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins, [-0.1 0.6])
hold off

%% Fig35 Matrix of Info for Semantic model
figure(35)
sp1=subplot(1,3,1);
imagesc(ModelInfo.Semantic(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells increasing Ceiling Information')
xlabel('windows (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Semantic model (bits)');
caxis([0 0.5])


sp2=subplot(1,3,2);
imagesc(ModelInfo.Ceiling(:,flip(IndSemOkRate))')
colorbar()
ylabel('Cells increasing Semantic Information')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Ceiling model (bits)');
caxis([0 0.7])


sp3=subplot(1,3,3);
imagesc(ModelInfo.Semantic(:,flip(IndSemOkRate))')
colorbar()
ylabel('Cells increasing Semantic Information')
xlabel('windows (ms)')
Xtickposition=get(sp3,'XTick');
set(sp3, 'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Semantic model(bits)');
caxis([0 0.5])

%% Fig 36: Check the LL of the semantic model compare to ceiling how confident we are in our models
figure(36)
sp1=subplot(1,2,1)
imagesc(LLAverage.CeilMinusFloor(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells  increasing Semantic Information')
xlabel('windows (ms)')
fprintf('enlarge figure\n')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Ceiling - LL Null Model');
caxis([0 0.8])

sp2=subplot(1,2,2)
imagesc(LLAverage.SemMinusFloor(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells increasing Semantic Information')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('LL Semantic - LL Null Model');
caxis([0 0.4])

%% Plot the average LL for the first 40 cells
figure(37)
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,IndCeilOkRate(end-39:end),1),2), nanmean(LLAverage.Semantic(:,IndCeilOkRate(end-39:end),1),2), nanmean(LLAverage.Floor(:,IndCeilOkRate(end-39:end),1),2),nanmean(LLAverage.Saturated(:,IndCeilOkRate(end-39:end),1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,4);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,IndCeilOkRate(end-39:end),1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'Cells';
myplotyyLLSPAsemOnly(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins)
ylim([-1.6 -0])
set(gca, 'YTick', -2:0.5:0)
set(gca, 'YTickLabel', -2:0.5:0)

hold off

%% Fig 38: Look at the proportion of info about semantic in the 15 first and then 20-35 first cells
figure(38)
sp1=subplot(1,2,1);
plot(1:NumMod, PropInfoperspike.Semantic(:,IndCeilOkRate(end-19:end)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmedian(PropInfoperspike.Semantic(:,IndCeilOkRate(end-19:end)),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('Proportion of Semantic Information')
xlabel('Time (ms)')
fprintf('enlarge figure\n')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Highly Semantic informative Cells')

sp2=subplot(1,2,2);
plot(1:NumMod, PropInfoperspike.Semantic(:,IndCeilOkRate(end-39:end-20)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmedian(PropInfoperspike.Semantic(:,IndCeilOkRate(end-39:end-20)),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('Proportion of Semantic Information')
xlabel('Time (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Moderately Semantic informative Cells')

%% Also look at the proportion of the response that we explain with the semantic model (LL Ratio)
LLSemMinusFloor = LLAverage.SemMinusFloor;
LLSemMinusFloor(LLSemMinusFloor<0)=0;
LLCeilMinusFloor = LLAverage.CeilMinusFloor;
LLCeilMinusFloor(LLCeilMinusFloor<0)=0;
GooFI.Semantic = LLSemMinusFloor ./ LLCeilMinusFloor;

figure(39)
sp1=subplot(1,2,1);
plot(1:NumMod, GooFI.Semantic(:,IndCeilOkRate(end-19:end)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmedian(GooFI.Semantic(:,IndCeilOkRate(end-19:end)),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('GooFI Semantic')
xlabel('Time (ms)')
fprintf('enlarge figure\n')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Highly Semantic informative Cells')

sp2=subplot(1,2,2);
plot(1:NumMod, GooFI.Semantic(:,IndCeilOkRate(end-39:end-20)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmedian(GooFI.Semantic(:,IndCeilOkRate(end-39:end-20)),2),'Color', [0 0.8 0 1]);
ylim([0 1])
ylabel('GooFI Semantic')
xlabel('Time (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Moderately Semantic informative Cells')


%% Figure 40 exact values of information about semantic
figure(40)
sp1=subplot(1,2,1);
plot(1:NumMod, ModelInfo.Semantic(:,IndCeilOkRate(end-19:end)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmean(ModelInfo.Semantic(:,IndCeilOkRate(end-19:end)),2),'Color', [0 0.8 0 1]);
ylim([0 0.4])
ylabel('Semantic Information (bits)')
xlabel('Time (ms)')
fprintf('enlarge figure\n')
pause()
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Highly Semantic informative Cells')

sp2=subplot(1,2,2);
plot(1:NumMod, ModelInfo.Semantic(:,IndCeilOkRate(end-39:end-20)),'Color', [0 0.8 0 0.2]);
hold on
plot(1:NumMod, nanmean(ModelInfo.Semantic(:,IndCeilOkRate(end-39:end-20)),2),'Color', [0 0.8 0 1]);
ylim([0 0.4])
ylabel('Semantic Information (bits)')
xlabel('Time (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Moderately Semantic informative Cells')

%%

ModelInfoperspike.CeilMinusFloor = ModelInfoperspike.Ceiling - ModelInfoperspike.Floor;
ModelInfoperspike.SemMinusFloor = ModelInfoperspike.Semantic - ModelInfoperspike.Floor;
sp5=subplot(2,4,5);
imagesc(ModelInfoperspike.CeilMinusFloor(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells  increasing Ceiling-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp5,'XTick');
set(sp5,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Ceiling model above Null Model (bits per spike)');
caxis([0 1])

sp6=subplot(2,4,6);
imagesc(ModelInfoperspike.SemMinusFloor(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells increasing Ceiling-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp6,'XTick');
set(sp6,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Semantic model above Null Model (bits/spike)');
caxis([0 1])


sp7=subplot(2,4,7);
imagesc(ModelInfoperspike.CeilMinusFloor(:,flip(IndSemOkRate))')
colorbar()
ylabel('Cells increasing Semantic-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp7,'XTick');
set(sp7,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Ceiling model above Null Model (bits/spike)');
caxis([0 1])


sp8=subplot(2,4,8);
imagesc(ModelInfoperspike.SemMinusFloor(:,flip(IndSemOkRate))')
colorbar()
ylabel('Cells increasing Semantic-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp4,'XTick');
set(sp4,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('Information of Semantic model above Null Model (bits/spike)');
caxis([0 1])
%% Proportion of Semantic Info

figure(34)
sp1=subplot(1,2,1);
imagesc(100.*PropInfoperspike.Semantic(:,flip(IndCeilOkRate))')
colorbar()
ylabel('Cells  increasing Ceiling-Null Model')
xlabel('windows (ms)')
Xtickposition=get(sp1,'XTick');
set(sp1,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('% Semantic Information');
caxis([0 100])

sp2=subplot(1,2,2);
imagesc(100.*PropInfoperspike.Semantic(:,flip(IndSemOkRate))')
colorbar()
ylabel('Cells increasing Sem Model')
xlabel('windows (ms)')
Xtickposition=get(sp2,'XTick');
set(sp2,'XTickLabel', [0 Wins(Xtickposition(2:end))])
title('% Semantic Information');
caxis([0 100])

%% To work on
sp3=subplot(2,2,3);
LL_MatrixPlot_y_left = nan(14,6,length(SuperSemCells));
for uu=1:length(SuperSemCells)
    LL_MatrixPlot_y_left(:,:,uu) = [LLAverage.CeilingOpt(:,uu,1), LLAverage.Acoustic(:,uu,1), LLAverage.AcSemOpt(:,uu,1), LLAverage.Semantic(:,uu,1), LLAverage.Floor(:,uu,1),LLAverage.Saturated(:,uu,1)];
end

LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = SpikeRateAveragepers(:,:,1);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SuperSemCells';
Transp = 0.2;
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins, Transp)

hold on
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Acoustic(:,SuperSemCells,1),2), nanmean(LLAverage.AcSemOpt(:,SuperSemCells,1),2), nanmean(LLAverage.Semantic(:,SuperSemCells,1),2), nanmean(LLAverage.Floor(:,SuperSemCells,1),2),nanmean(LLAverage.Saturated(:,SuperSemCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,SuperSemCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'SuperSemCells';
Transp = 1;
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins, Transp)
hold off



subplot(2,2,4)
for uu=1:size(WorstSemanticCells)
    hold on
    LL_MatrixPlot_y_left = [LLAverage.CeilingOpt(:,uu,1), LLAverage.Acoustic(:,uu,1), LLAverage.AcSemOpt(:,uu,1), LLAverage.Semantic(:,uu,1), LLAverage.Floor(:,uu,1),LLAverage.Saturated(:,uu,1)];
    LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
    SPA_Plot_y_right = SpikeRateAveragepers(:,uu,1);
    SPA_Plot_x_right = 1:NumMod;
    CellType = 'WorstSemCells';
    Transp = 0.2;
    myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins, Transp)
end
hold on
LL_MatrixPlot_y_left = [nanmean(LLAverage.CeilingOpt(:,WorstSemanticCells,1),2), nanmean(LLAverage.Acoustic(:,WorstSemanticCells,1),2), nanmean(LLAverage.AcSemOpt(:,WorstSemanticCells,1),2), nanmean(LLAverage.Semantic(:,WorstSemanticCells,1),2), nanmean(LLAverage.Floor(:,WorstSemanticCells,1),2),nanmean(LLAverage.Saturated(:,WorstSemanticCells,1),2)];
LL_MatrixPlot_x_left = repmat((1:NumMod)',1,6);
SPA_Plot_y_right = nanmean(SpikeRateAveragepers(:,WorstSemanticCells,1),2);
SPA_Plot_x_right = 1:NumMod;
CellType = 'WorstSemCells';
Transp = 1;
myplotyyLLSPA(LL_MatrixPlot_x_left,LL_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,CellType,Wins, Transp)
hold off


%% Piece of code to work on
subplot(1,2,2)
for mm=1:NumMod
    % Plot cells for which Ceil=Floor
    plot(mm,find(~Signif.CeilvsFl(:,mm)), 'ko','Markersize',10,'MarkerFaceColor',[1 1 1], 'MarkerEdgeColor','k');
    hold on
    % Plot Cells for which Acoustic model is better than the Floor and
    % AcSem either is not better than the Floor or/and is not doing better
    % than Ac (yellow)
    plot(mm,intersect(find(Signif.AcvsFl(:,mm)),find(~Signif.AcSemvsAc(:,mm))), 'ko','Markersize',10,'MarkerFaceColor',[1 1 0], 'MarkerEdgeColor','k');
    hold on
    % Plot Cells for which Ac and AcSem are better than floor and
    % AcSem>Ac (orange)
    Localcells = intersect(intersect(find(Signif.AcvsFl(:,mm)),find(Signif.AcSemvsFl(:,mm))),find(Signif.AcSemvsAc(:,mm)));
    plot(mm,Localcells, 'ko','Markersize',10,'MarkerFaceColor',[1 0.5 0], 'MarkerEdgeColor','k');
    hold on
    % Plot Cells for which only AcSem is better than floor and
    % AcSem>Ac (red)
    Localcells = intersect(intersect(find(~Signif.AcvsFl(:,mm)),find(Signif.AcSemvsFl(:,mm))),find(Signif.AcSemvsAc(:,mm)));   
    plot(mm,Localcells, 'ko','Markersize',10,'MarkerFaceColor',[1 0 0], 'MarkerEdgeColor','k');
    hold on
    % Plot Cells for which none of the two (Ac and AcSem) are better than floor
    plot(mm,intersect(find(~Signif.AcvsFl(:,mm)),find(~Signif.AcSemvsFl(:,mm))), 'ko','Markersize',10,'MarkerFaceColor',[0.8 0.8 0.8], 'MarkerEdgeColor','k');
    hold on
end
hold off


%% Population overview of model results over windows
% Number of cells that has a Ceiling model better than the floor
NumCeilFlCells = nansum(AC.SignifDiffDev(:,4));
fprintf('There are %d/%d cells for which the Ceiling model predict better than the floor: units have info about the various vocalization types\n',NumCeilFlCells, length(AC.SignifDiffDev(:,4)));
% How many of those have an acoustic model that predict better than the
% floor?
NumAcFlCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,5));
fprintf('There are %d/%d VT cells for which the acoustic model is doing better than floor\n',NumAcFlCells,NumCeilFlCells);
% How many of those have an AcSem model that predict better than the
% floor?
NumAcSemFlCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,6));
fprintf('There are %d/%d VT cells for which the AcSem model is doing better than floor\n',NumAcSemFlCells,NumCeilFlCells);
% How many of those have both an Ac and an AcSem model that predict better than the
% floor?
NumAcAcSemFlCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,5) .* AC.SignifDiffDev(:,6));
fprintf('There are %d/%d VT cells for which the AcSem and Acoustic models are doing better than floor\n',NumAcAcSemFlCells,NumCeilFlCells);

% How many of those have a semantic model that predict better than the
% floor?
NumSemFlCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,7));
fprintf('There are %d/%d VT cells for which the Semantic model is doing better than floor\n',NumSemFlCells,NumCeilFlCells);

% How many of those have both a semantic and an AcSem model that predict better than the
% floor?
NumSemAcSemFlCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,7).* AC.SignifDiffDev(:,6));
fprintf('There are %d/%d VT cells for which the Semantic and Acsem models are doing better than floor\n',NumSemAcSemFlCells,NumCeilFlCells);

% How many of those have aan AcSem model that is doing better than the
% Acoustic model?
NumAcAcSemCells = nansum(AC.SignifDiffDev(:,4) .* AC.SignifDiffDev(:,2));
fprintf('There are %d/%d VT cells for which the Acsem model is doing better than Acoustic\n',NumAcAcSemCells,NumCeilFlCells);