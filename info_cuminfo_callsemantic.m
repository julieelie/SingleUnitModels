function [ParamModel, Data, InputData, Wins]=info_cuminfo_callsemantic(Trials,VocType, ParamModel, Calfilename)
FIG=0;
if nargin<3
    ParamModel = struct();
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

% Number of bootstraps
if ~isfield(ParamModel, 'NbBoot_Info') || isempty(ParamModel.NbBoot_Info)
    ParamModel.NbBoot_Info = 100;
end

% Set parameters for the number of samples that should be tested in the
% MonteCarlo estimation of the cumulative information
if ~isfield(ParamModel, 'NumSamples_MC_Cum_Info') || isempty(ParamModel.NumSamples_MC_Cum_Info)
ParamModel.NumSamples_MC_Cum_Info = [10^2 10^3 10^4 10^5 10^6]; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers
end

% Set the Parameters of the Markov approximation of the cumulative
% information
if ~isfield(ParamModel, 'MarkovParameters_Cum_Info') || isempty(ParamModel.MarkovParameters_Cum_Info)
    ParamModel.MarkovParameters_Cum_Info = [2 3 4;1 1 1];
end

% Set the fix time history to calculate the exact cumulative information
% (above 4 will certainly bug the computer asking for too big matrices)
if ~isfield(ParamModel, 'ExactHist') || isempty(ParamModel.ExactHist)
    ParamModel.ExactHist = 4;
end

if nargin<4
    saveonline = 0;
else
    saveonline = 1;
end

% define the list of end points of spectrogram and neural responses windows
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;

% # of models to run on the data
WinNum = length(Wins);

% Number of stims in the data set
NbStims = length(Trials);

% Number of stimulus categories
IdCats = unique(VocType);
NbCat = length(IdCats);



%% Initialize output variables
InputData.InfoStim = nan(NbStims,WinNum);
InputData.InfoStim_trials = cell(1,NbStims);
InputData.InfoCat = nan(NbCat, WinNum);
Data.stim_entropy = nan(1,WinNum);
Data.category_entropy = nan(1,WinNum);
Data.stim_info = nan(1,WinNum);
Data.category_info = nan(1,WinNum);
Data.stim_info_biais = nan(1,WinNum);
Data.category_info_biais = nan(1,WinNum);
Data.stim_P_YgivenS = cell(1,WinNum);
Data.category_P_YgivenS = cell(1,WinNum);
Data.stim_P_YgivenS_Bootstim = cell(ParamModel.NbBoot_Info/10,WinNum);
Data.category_P_YgivenS_Bootstim = cell(ParamModel.NbBoot_Info/10,WinNum);
Data.stim_P_YgivenS_Boottrial = cell(ParamModel.NbBoot_Info/10,WinNum);
Data.category_P_YgivenS_Boottrial = cell(ParamModel.NbBoot_Info/10,WinNum);

%% Now loop through bins and calculate spike patterns and instantaneous information
for ww = 1:WinNum
    Tstart=tic;
    %fprintf(1,'%d/%d models\n', mm, modNum);
    Win = Wins(ww);
    FirstTimePoint = Win - ParamModel.NeuroBin+ ParamModel.ResDelay +1;
    LastTimePoint = Win + ParamModel.ResDelay;
     
    % Calculating info about the stims
    Boot_switch = 0;
    [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch);
    Data.stim_info(ww) = Local_Output.value;
    InputData.InfoStim(:,ww) = Local_Output.Inputdata;
    Data.stim_P_YgivenS{ww} = Local_Output.P_YgivenS;
    Data.stim_entropy(ww) = Local_Output.data_entropy;
    for ss=1:NbStims
        if ww==1
            InputData.InfoStim_trials{ss} = Local_Output.Inputdata_samples{ss};
        else
            InputData.InfoStim_trials{ss} = [InputData.InfoStim_trials{ss} Local_Output.Inputdata_samples{ss}];
        end
    end
    
    % Calculating info about the categories
    [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch, VocType);
    Data.category_info(ww) = Local_Output.value;
    InputData.InfoCat(:,ww) = Local_Output.Inputdata;
    Data.category_P_YgivenS{ww} = Local_Output.P_YgivenS;
    Data.category_entropy(ww) = Local_Output.data_entropy;
    
    % Bootstrapping the calculation of information accross stim/categories and trials
    % to calculate the positive bias
    Boot_switch = 2;
    Info_biais_stim = nan(1,ParamModel.NbBoot_Info);
    Info_biais_category = nan(1,ParamModel.NbBoot_Info);
    for bb=1:ParamModel.NbBoot_Info
        [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch);
        Info_biais_stim(bb) = Local_Output.value;
        if bb/10 == round(bb/10)
            Data.stim_P_YgivenS_Boottrial{bb/10,ww} = Local_Output.P_YgivenS;
        end
        [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch,VocType);
        Info_biais_category(bb) = Local_Output.value;
        if bb/10 == round(bb/10)
            Data.category_P_YgivenS_Boottrial{bb/10,ww} = Local_Output.P_YgivenS;
        end
    end
    Data.stim_info_biais(ww)=mean(Info_biais_stim);
    Data.category_info_biais(ww)=mean(Info_biais_category);
    
    % Bootstrapping the calculation of information accross trials within
    % stim/categories to calculate the incertainty of information due to
    % the incertainty on spike rates
    Boot_switch = 1;
    Info_boot_stim = nan(1,ParamModel.NbBoot_Info);
    Info_boot_category = nan(1,ParamModel.NbBoot_Info);
    for bb=1:ParamModel.NbBoot_Info
        [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch);
        Info_boot_stim(bb) = Local_Output.value;
        if bb/10 == round(bb/10)
            Data.stim_P_YgivenS_Bootstim{bb/10,ww} = Local_Output.P_YgivenS;
        end
        [Local_Output] = info_model_Calculus_wrapper(Trials, FirstTimePoint, LastTimePoint,Boot_switch,VocType);
        Info_boot_category(bb) = Local_Output.value;
        if bb/10 == round(bb/10)
            Data.category_P_YgivenS_Bootstim{bb/10,ww} = Local_Output.P_YgivenS;
        end
    end
    Data.stim_info_boot_var(ww)=var(Info_boot_stim);
    Data.category_info_boot_var(ww)=var(Info_boot_category);
    Data.stim_info_boot_mean(ww)=mean(Info_boot_stim);
    Data.category_info_boot_mean(ww)=mean(Info_boot_category);
    fprintf('Done bin %d/%d after %f sec\n', ww, WinNum, toc(Tstart));
end

 %% Save what we have for now
 if saveonline
     if exist(Calfilename, 'file')==2
        save(Calfilename,'Data','VocType','ParamModel','Wins','InputData','-append');
     else
         save(Calfilename,'Data','VocType','ParamModel','Wins','InputData');
     end
 end

 %% Plot the results if requested
 if FIG
     figure()
     subplot(2,1,1)
     plot(Data.stim_info,'LineWidth',2, 'Color',[0 0 0])
     hold on
     plot(Data.stim_info_biais,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(Data.stim_info_boot_mean,'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Bootstrap accross stims (biais)', 'Bootstrap within stims', 'Location','NorthEast');
     hold on
     plot(Data.stim_info_boot_mean + 2*Data.stim_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.stim_info_boot_mean - 2*Data.stim_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold on
     plot(Data.stim_entropy, 'LineStyle','-.','Color','r')
     hold off
     ylim([-0.5 max(Data.stim_entropy)+1])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Stimulus Information in bits')
     
     subplot(2,1,2)
     plot(Data.stim_info,'LineWidth',2, 'Color',[0 0 0])
     hold on
     plot(Data.stim_info_biais,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(Data.stim_info_boot_mean,'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Bootstrap accross stims (biais)', 'Bootstrap within stims', 'Location','NorthEast');
     hold on
     plot(Data.stim_info_boot_mean + 2*Data.stim_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.stim_info_boot_mean - 2*Data.stim_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold off
     YL = get(gca,'YLim');
     ylim([-0.2 YL(2)])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Stimulus Information in bits')
     
     figure()
     subplot(2,1,1)
     plot(Data.category_info,'LineWidth',2, 'Color',[0 0 0])
     hold on
     plot(Data.category_info_biais,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(Data.category_info_boot_mean,'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Bootstrap accross stims (biais)', 'Bootstrap within stims', 'Location','NorthEast');
     hold on
     plot(Data.category_info_boot_mean + 2*Data.category_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.category_info_boot_mean - 2*Data.category_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold on
     plot(Data.category_entropy, 'LineStyle','-.','Color','r')
     hold off
     ylim([-0.5 max(Data.category_entropy)+1])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Category Information in bits')
     
     subplot(2,1,2)
     plot(Data.category_info,'LineWidth',2, 'Color',[0 0 0])
     hold on
     plot(Data.category_info_biais,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(Data.category_info_boot_mean,'LineWidth',2, 'Color',[0 1 0])
     legend('Information','Bootstrap accross stims (biais)', 'Bootstrap within stims', 'Location','NorthEast');
     hold on
     plot(Data.category_info_boot_mean + 2*Data.category_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.category_info_boot_mean - 2*Data.category_info_boot_var, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold off
     YL = get(gca,'YLim');
     ylim([-0.2 YL(2)])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Category Information in bits')
 end
%% Now calculating cumulative information
% Getting output variables ready
Data.cum_info_stim = struct();
Data.cum_info_cat = struct();
HY_Markov_stim = struct();
HY_Markov_cat = struct();
for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
    ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
    Data.cum_info_stim.(sprintf('%s',ModelType))= nan(2,WinNum);
    Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
    Data.cum_info_cat.(sprintf('%s',ModelType))= nan(2,WinNum);
    Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
end
ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
Data.cum_info_stim.(sprintf('%s',ModelType)) = nan(1,WinNum);
Data.cum_info_cat.(sprintf('%s',ModelType)) = nan(1,WinNum);
Data.cum_info_stim.(sprintf('%s',ModelType))(1) = Data.stim_info(1);
Data.cum_info_cat.(sprintf('%s',ModelType))(1) = Data.category_info(1);
for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
    ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
    Data.cum_info_stim.(sprintf('%s',ModelType))= nan(1,WinNum);
    Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
    Data.cum_info_cat.(sprintf('%s',ModelType))= nan(1,WinNum);
    Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
    HY_Markov_stim.(sprintf('%s',ModelType)) = nan(1,WinNum);
    HY_Markov_cat.(sprintf('%s',ModelType)) = nan(1,WinNum);
end

%% Running through time bins and calculate cumulative information
for ww =2:WinNum
    Tstart2=tic;
    for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
        ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
         % Monte Carlo estimation with full memory cumulative information stimuli
        [Data.cum_info_stim.(sprintf('%s',ModelType))(:,ww),~]=info_cumulative_model_Calculus(Data.stim_P_YgivenS(1:ww),'Model#',1,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
        % Monte Carlo estimation with full memory cumulative information
        % categories
        [Data.cum_info_cat.(sprintf('%s',ModelType))(:,ww),~]=info_cumulative_model_Calculus(Data.category_P_YgivenS(1:ww),'Model#',2,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
    end
    
    ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
    % Exact calculation cumulative information on stims with ParamModel.ExactHist*10 ms memory
    [Data.cum_info_stim.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.stim_P_YgivenS(1:ww),'Model#',1,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
    % Exact calculation cumulative information on categories with ParamModel.ExactHist*10 ms memory
    [Data.cum_info_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.category_P_YgivenS(1:ww),'Model#',2,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
    
    for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
        ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
        if ww>2
            % Markov chain estimation of cumulative information on stims
            [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.stim_P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss), 'HY_old', HY_Markov_stim.(sprintf('%s',ModelType))(ww-1));
            % Markov chain estimation of cumulative information on
            % categories
            [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.category_P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss),'HY_old', HY_Markov_cat.(sprintf('%s',ModelType))(ww-1));
        else
            % Markov chain estimation of cumulative information on stims
            [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.stim_P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
            % Markov chain estimation of cumulative information on
            % categories
            [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.category_P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
        end
    end
    fprintf('Done cumulative information bin %d/%d after %f sec\n', ww, WinNum, toc(Tstart2));
end

 %% Save what we have for now
 if saveonline
     save(Calfilename,'Data','ParamModel', 'NbBoot','-append');
 end

%% Plot the cumulative information
if FIG
    figure()
    subplot(2,1,1)
    plot(1:WinNum, Data.stim_info, 1:WinNum,Data.cum_info_stim.MarkovEst4, 1:WinNum,Data.cum_info_stim.MarkovEst3, 1:WinNum,Data.cum_info_stim.MarkovEst2,1:WinNum,Data.cum_info_stim.MonteCarlo6(1,:), 1:WinNum,Data.cum_info_stim.MonteCarlo5(1,:), 1:WinNum,Data.cum_info_stim.MonteCarlo4(1,:), 1:WinNum,Data.cum_info_stim.MonteCarlo3(1,:),1:WinNum,Data.cum_info_stim.MonteCarlo2(1,:), 1:WinNum, Data.cum_info_stim.ExactMem0_5, 1:WinNum, cumsum(Data.stim_info), 'LineWidth',2)
    Color_Lines = get(gca,'ColorOrder');
    legend('Information', 'Cumulative Info Markov4', 'Cumulative Info Markov3','Cumulative Info Markov2','Cumulative Information MC 10^6','Cumulative Information MC 10^5','Cumulative Information MC 10^4','Cumulative Information MC 10^3', 'Cumulative Information MC 10^2', 'Exact Cumulative Information 40ms history', 'Cumulative sum of info', 'Location', 'NorthWest')
    Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Stimulus Information in bits')
    subplot(2,1,2)
    plot(1:WinNum, Data.category_info, 1:WinNum,Data.cum_info_cat.MarkovEst4, 1:WinNum,Data.cum_info_cat.MarkovEst3, 1:WinNum,Data.cum_info_cat.MarkovEst2,1:WinNum,Data.cum_info_cat.MonteCarlo6(1,:), 1:WinNum,Data.cum_info_cat.MonteCarlo5(1,:), 1:WinNum,Data.cum_info_cat.MonteCarlo4(1,:), 1:WinNum,Data.cum_info_cat.MonteCarlo3(1,:),1:WinNum,Data.cum_info_cat.MonteCarlo2(1,:), 1:WinNum, Data.cum_info_cat.ExactMem0_5, 1:WinNum, cumsum(Data.category_info), 'LineWidth',2)
    Color_Lines = get(gca,'ColorOrder');
    legend('Information', 'Cumulative Info Markov4', 'Cumulative Info Markov3','Cumulative Info Markov2','Cumulative Information MC 10^6','Cumulative Information MC 10^5','Cumulative Information MC 10^4','Cumulative Information MC 10^3', 'Cumulative Information MC 10^2', 'Exact Cumulative Information 40ms history', 'Cumulative sum of info', 'Location', 'NorthWest')
    Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Category Information in bits')
end
%% Calculating bootstrap values accross stims and across trials within stims
ModelTypeMCstim = sprintf('MonteCarlo%d_Bootstim',log10(ParamModel.NumSamples_MC_Cum_Info(end)));
Data.cum_info_stim.(sprintf('%s',ModelTypeMCstim)) = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.(sprintf('%s',ModelTypeMCstim)) = nan(ParamModel.NbBoot_Info/10,WinNum);
ModelTypeMCtrial = sprintf('MonteCarlo%d_Boottrial',log10(ParamModel.NumSamples_MC_Cum_Info(end)));
Data.cum_info_stim.(sprintf('%s',ModelTypeMCtrial)) = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.(sprintf('%s',ModelTypeMCtrial)) = nan(ParamModel.NbBoot_Info/10,WinNum);

Data.cum_info_stim.ExactMem0_5_Bootstim = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.ExactMem0_5_Bootstim = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_stim.ExactMem0_5_Boottrial = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.ExactMem0_5_Boottrial = nan(ParamModel.NbBoot_Info/10,WinNum);

ModelTypeMarstim = sprintf('MarkovEst%d_Bootstim',ParamModel.MarkovParameters_Cum_Info(1,end));
Data.cum_info_stim.(sprintf('%s',ModelTypeMarstim)) = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.(sprintf('%s',ModelTypeMarstim)) = nan(ParamModel.NbBoot_Info/10,WinNum);
ModelTypeMartrial = sprintf('MarkovEst%d_Boottrial',ParamModel.MarkovParameters_Cum_Info(1,end));
Data.cum_info_stim.(sprintf('%s',ModelTypeMartrial)) = nan(ParamModel.NbBoot_Info/10,WinNum);
Data.cum_info_cat.(sprintf('%s',ModelTypeMartrial)) = nan(ParamModel.NbBoot_Info/10,WinNum);

for bb=1:ParamModel.NbBoot_Info/10
    fprintf('Bootstrap %d/%d\n', bb, Nb_Boot/10);
    HY_Markov5.stim.bbtrial = nan(1,WinNum);
    HY_Markov5.stim.bbstim = nan(1,WinNum);
    HY_Markov5.cat.bbtrial = nan(1,WinNum);
    HY_Markov5.cat.bbstim = nan(1,WinNum);
    for ww=2:Nb_Win
        fprintf('Bootstrap %d/%d Time point %d/%d\n',bb, ParamModel.NbBoot_Info, ww, WinNum);
        
        % First for the cumulative information about stimuli
        P_YgivenS_local_trial = Data.stim_P_YgivenS_Boottrial(bb,1:ww);
        P_YgivenS_local_stim = Data.stim_P_YgivenS_Bootstim(bb,1:ww);
        
        % Monte Carlo estimation with full memory
        [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
        Data.cum_info_stim.(sprintf('%s',ModelTypeMCstim))(bb,ww)=Icum_EstMonteCarlo_temp(1);
        [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
        Data.cum_info_stim.(sprintf('%s',ModelTypeMCtrial))(bb,ww)=Icum_EstMonteCarlo_temp(1);
        
        % Exact calculation with 50 ms memory
        [Data.cum_info_stim.ExactMem0_5_Bootstim(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',5);
        [Data.cum_info_stim.ExactMem0_5_Boottrial(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',5);
        
        if ww==2
            % Markov chain estimation 50 ms
            [Data.cum_info_stim.(sprintf('%s',ModelTypeMarstim))(bb,ww),HY_Markov5.stim.bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
            [Data.cum_info_stim.(sprintf('%s',ModelTypeMartrial))(bb,ww),HY_Markov5.stim.bbtrial(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
        else
            % Markov chain estimation 50 ms
            [Data.cum_info_stim.(sprintf('%s',ModelTypeMarstim))(bb,ww),HY_Markov5.stim.bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov5.stim.bbstim(ww-1));
            [Data.cum_info_stim.(sprintf('%s',ModelTypeMartrial))(bb,ww),HY_Markov5.stim.bbtrial(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov5.stim.bbtrial(ww-1));
        end
        
        % Then same thing for the cumulative information about categories
        P_YgivenS_local_trial = Data.cat_P_YgivenS_Boottrial(bb,1:ww);
        P_YgivenS_local_stim = Data.cat_P_YgivenS_Bootstim(bb,1:ww);
        
        % Monte Carlo estimation with full memory
        [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
        Data.cum_info_cat.(sprintf('%s',ModelTypeMCstim))(bb,ww)=Icum_EstMonteCarlo_temp(1);
        [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
        Data.cum_info_cat.(sprintf('%s',ModelTypeMCtrial))(bb,ww)=Icum_EstMonteCarlo_temp(1);
        
        % Exact calculation with 50 ms memory
        [Data.cum_info_cat.ExactMem0_5_Bootstim(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',5);
        [Data.cum_info_cat.ExactMem0_5_Boottrial(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',5);
        
        if ww==2
            % Markov chain estimation 50 ms
            [Data.cum_info_cat.(sprintf('%s',ModelTypeMarstim))(bb,ww),HY_Markov5.cat.bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
            [Data.cum_info_cat.(sprintf('%s',ModelTypeMartrial))(bb,ww),HY_Markov5.cat.bbtrial(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
        else
            % Markov chain estimation 50 ms
            [Data.cum_info_cat.(sprintf('%s',ModelTypeMarstim))(bb,ww),HY_Markov5.cat.bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_stim,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov5.cat.bbstim(ww-1));
            [Data.cum_info_cat.(sprintf('%s',ModelTypeMartrial))(bb,ww),HY_Markov5.cat.bbtrial(ww),~]=info_cumulative_model_Calculus(P_YgivenS_local_trial,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov5.cat.bbtrial(ww-1));
        end
        
    end
end

%% Save what we have for now
 if saveonline
     save(Calfilename,'Data','ParamModel', 'NbBoot','-append');
 end
end