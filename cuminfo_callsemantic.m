function [Calfilename_local]=cuminfo_callsemantic(Data, ParamModel,  Calfilename)
FIG=0;

%% Deals with input parameters
if nargin<2
    ParamModel = struct();
end
if  ~isfield(ParamModel,'MinWin') || isempty(ParamModel.MinWin)
    ParamModel.MinWin = 10; % end point of the first analysis window (spectrogram and neural response)
end

if ~isfield(ParamModel,'MaxWin_cumInfo') || isempty(ParamModel.MaxWin_cumInfo)
    ParamModel.MaxWin_cumInfo = 600; %end point of the last anaysis window for...
    ... the calculation of cumulative information
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
if ~isfield(ParamModel, 'NbBoot_CumInfo') || isempty(ParamModel.NbBoot_CumInfo)
    ParamModel.NbBoot_CumInfo = 10;
end

% now use the predetermined sets of JK indices
Avail_NJK = length(ParamModel.SetIndices_JK);

if ParamModel.NbBoot_CumInfo > Avail_NJK
    fprintf('WARNING: Only %d possible calculations of sets of JK points while you are asking for %d in the calculation of cumulative information\n', Avail_NJK, ParamModel.NbBoot_CumInfo);
    ParamModel.NbBoot_CumInfo = Avail_NJK;
else
    fprintf('WARNING: %d possible calculations of sets of JK points. you are asking for %d in the calculation of cumulative information\n', Avail_NJK, ParamModel.NbBoot_CumInfo);
end

% Set parameters for the max number of samples that should be tested in the
% optimal MonteCarlo estimation of the cumulative information
if ~isfield(ParamModel, 'MaxNumSamples_MCopt_Cum_Info')
    ParamModel.MaxNumSamples_MCopt_Cum_Info = 5*10^6; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers
elseif isempty(ParamModel.MaxNumSamples_MCopt_Cum_Info)
    fprintf(1,'No Monte Carlo approximation of the cumulative information with an optimized number of samples will be calculated\n');
end

%% Set up parameters of the function
% define the list of end points of time bins at which information is
% calculated
Data.Wins_cumInfo = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin_cumInfo;

% # of bins in the data
WinNum_cumInfo = length(Data.Wins_cumInfo);

%% Configure Parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    MyParPool = parpool(min(ParamModel.NbBoot_Info,str2num(getenv('SLURM_CPUS_ON_NODE'))),'IdleTimeout', Inf);
    system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
    [~,JobID] = system('echo $SLURM_JOB_ID');
    parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];    
end


%% Now calculating cumulative information if requested by the presence of input parameters
fprintf(1,'Starting calculation of cumulative information\n')

%% Running optimal Monte Carlo sampling cumulative information
if strcmp(ParamModel.CIType, 'CIS')
    fprintf(1,'Cumulative information for stimuli\n')
    % Initializing output variables
    Data.cum_info_stim = struct();
    
    % Cumulative information for stimuli
    [Data.cum_info_stim.MonteCarloOpt_raw, Data.cum_info_stim.MonteCarloOpt_bcorr, Data.cum_info_stim.MonteCarloOpt_err, Data.cum_info_stim.MonteCarloOpt_Samples] = cumulative_info_poisson_model_MCJK_wrapper(Data.P_YgivenS, Data.P_YgivenS_Bootstrap, ParamModel.Mean_Ntrials_perstim, WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
    
    % Filling in the first value of cumulative info with information
    % value at bin 1
    Data.cum_info_stim.MonteCarloOpt_raw(1) = Data.stim_info(1);
    Data.cum_info_stim.MonteCarloOpt_bcorr(1) = Data.stim_info_bcorr(1);
    Data.cum_info_stim.MonteCarloOpt_err(1) = Data.stim_info_err(1);
     %% Save what we have for now
     Calfilename_local = [Calfilename(1:end-4) '_CIS.mat'];
     if exist(Calfilename_local, 'file')==2
         save(Calfilename_local,'Data','ParamModel','-append');
     else
         save(Calfilename_local,'Data','ParamModel');
     end
     
     %% Data for figure
     if FIG
         LocalData = Data.cum_info_stim;
     end
end

if strcmp(ParamModel.CIType, 'CIC')
    fprintf(1,'Cumulative information for semantic\n')
    % Initializing output variables
    Data.cum_info_cat = struct();
    % Cumulative information for categories
    [Data.cum_info_cat.MonteCarloOpt_raw, Data.cum_info_cat.MonteCarloOpt_bcorr, Data.cum_info_cat.MonteCarloOpt_err, Data.cum_info_cat.MonteCarloOpt_Samples] = cumulative_info_poisson_model_MCJK_wrapper(Data.P_YgivenC, Data.P_YgivenC_Bootstrap, ParamModel.Mean_Ntrials_perstim,WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
    
    % Filling in the first value of cumulative info with information
    % value at bin 1
    Data.cum_info_cat.MonteCarloOpt_raw(1) = Data.category_info(1);
    Data.cum_info_cat.MonteCarloOpt_bcorr(1) = Data.category_info_bcorr(1);
    Data.cum_info_cat.MonteCarloOpt_err(1) = Data.category_info_err(1);
    
    %% Save what we have for now
    Calfilename_local = [Calfilename(1:end-4) '_CIC.mat'];
    if exist(Calfilename_local, 'file')==2
        save(Calfilename_local,'Data','ParamModel','-append');
    else
        save(Calfilename_local,'Data','ParamModel');
    end
    %% Data for figure
     if FIG
         LocalData = Data.cum_info_cat;
     end
end

if strcmp(ParamModel.CIType, 'CICRand')
    fprintf(1,'Cumulative information for random categories\n')
    % Initializing output variables
    Data.cum_info_catRand = struct();
    % Cumulative information for categories
    [Data.cum_info_catRand.MonteCarloOpt_raw, Data.cum_info_catRand.MonteCarloOpt_bcorr, Data.cum_info_catRand.MonteCarloOpt_err, Data.cum_info_catRand.MonteCarloOpt_Samples] = cumulative_info_poisson_model_MCJK_wrapper(Data.P_YgivenCRand, Data.P_YgivenCRand_Bootstrap, ParamModel.Mean_Ntrials_perstim,WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
    
    % Filling in the first value of cumulative info with information
    % value at bin 1
    Data.cum_info_catRand.MonteCarloOpt_raw(1) = Data.category_info(1);
    Data.cum_info_catRand.MonteCarloOpt_bcorr(1) = Data.category_info_bcorr(1);
    Data.cum_info_catRand.MonteCarloOpt_err(1) = Data.category_info_err(1);
    
    %% Save what we have for now
    Calfilename_local = [Calfilename(1:end-4) '_CICRand.mat'];
    if exist(Calfilename_local, 'file')==2
        save(Calfilename_local,'Data','ParamModel','-append');
    else
        save(Calfilename_local,'Data','ParamModel');
    end
    
    %% Data for figure
     if FIG
         LocalData = Data.cum_info_catRand;
     end
end

%% Cumulative information if the neuron was using a constant rate coding
if strcmp(ParamModel.CIType, 'CISR')
    fprintf(1,'Cumulative information for stimuli with fixed rate\n')
    % Initializing output variables
    Data.cum_info_stim_csteRate = struct();

    P_YgivenS_local = cell(WinNum_cumInfo,1);
    P_YgivenS_Bootstrap_local = cell(1,ParamModel.NbBoot_CumInfo);
    for ww=1:WinNum_cumInfo
        P_YgivenS_local{ww} = Data.P_YgivenS_csteRate;
    end
    for bb=1:ParamModel.NbBoot_CumInfo
        NJKsets = size(ParamModel.SetIndices_JK{bb},1);
        P_YgivenS_Bootstrap_local{bb} = cell(NJKsets, WinNum_cumInfo);
        for jk=1:NJKsets
            for ww=1:WinNum_cumInfo
                P_YgivenS_Bootstrap_local{bb}{jk,ww} = Data.P_YgivenS_Bootstrap_csteRate{bb}{jk};
            end
        end
    end
    
    % Cumulative information for stimuli
    [Data.cum_info_stim_csteRate.MonteCarloOpt_raw, Data.cum_info_stim_csteRate.MonteCarloOpt_bcorr, Data.cum_info_stim_csteRate.MonteCarloOpt_err, ~] = cumulative_info_poisson_model_MCJK_wrapper(P_YgivenS_local, P_YgivenS_Bootstrap_local, ParamModel.Mean_Ntrials_perstim, WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
    
    % Filling in the first value of cumulative info with information
    % value at bin 1
    Data.cum_info_stim_csteRate.MonteCarloOpt_bcorr(1) = Data.stim_info_bcorr_csteRate(1);
    Data.cum_info_stim_csteRate.MonteCarloOpt_err(1) = Data.stim_info_err_csteRate(1);
    %% Save what we have for now
    Calfilename_local = [Calfilename(1:end-4) '_CISR.mat'];
    if exist(Calfilename_local, 'file')==2
        save(Calfilename_local,'Data','ParamModel','-append');
    else
        save(Calfilename_local,'Data','ParamModel');
    end
    
    %% Data for figure
     if FIG
         LocalData = Data.cum_info_stim_csteRate;
     end
end

if strcmp(ParamModel.CIType, 'CICR')
    fprintf(1,'Cumulative information for semantic with fixed rate\n')
    % Initializing output variables
    Data.cum_info_cat_csteRate = struct();

    P_YgivenC_local = cell(WinNum_cumInfo,1);
    P_YgivenC_Bootstrap_local = cell(1,ParamModel.NbBoot_CumInfo);
    for ww=1:WinNum_cumInfo
        P_YgivenC_local{ww} = Data.P_YgivenC_csteRate;
    end
    for bb=1:ParamModel.NbBoot_CumInfo
        NJKsets = size(ParamModel.SetIndices_JK{bb},1);
        P_YgivenC_Bootstrap_local{bb} = cell(NJKsets, WinNum_cumInfo);
        for jk=1:NJKsets
            for ww=1:WinNum_cumInfo
                P_YgivenC_Bootstrap_local{bb}{jk,ww} = Data.P_YgivenC_Bootstrap_csteRate{bb}{jk};
            end
        end
    end
    % Cumulative information for categories
    [Data.cum_info_cat_csteRate.MonteCarloOpt_raw, Data.cum_info_cat_csteRate.MonteCarloOpt_bcorr, Data.cum_info_cat_csteRate.MonteCarloOpt_err, ~] = cumulative_info_poisson_model_MCJK_wrapper(P_YgivenC_local, P_YgivenC_Bootstrap_local, ParamModel.Mean_Ntrials_perstim,WinNum_cumInfo, ParamModel.NbBoot_CumInfo, ParamModel.MaxNumSamples_MCopt_Cum_Info);
    
    % Filling in the first value of cumulative info with information
    % value at bin 1
    Data.cum_info_cat_csteRate.MonteCarloOpt_bcorr(1) = Data.category_info_bcorr_csteRate(1);
    Data.cum_info_cat_csteRate.MonteCarloOpt_err(1) = Data.category_info_err_csteRate(1);
    
    %% Save what we have for now
    Calfilename_local = [Calfilename(1:end-4) '_CICR.mat'];
    if exist(Calfilename_local, 'file')==2
        save(Calfilename_local,'Data','ParamModel','-append');
    else
        save(Calfilename_local,'Data','ParamModel');
    end
    
     %% Data for figure
     if FIG
         LocalData = Data.cum_info_cat_csteRate;
     end
end


%% Plot the cumulative information
if FIG
    ColorCode = get(groot,'DefaultAxesColorOrder');
    ColorCode = [ColorCode ; 0.85 0.6940 0.556; 0.301 0.078 0.741];
    figure()
    WinNum = length(Data.stim_info_bcorr);
    %subplot(3,1,1)
    Biais_MC = LocalData.MonteCarloOpt_raw - LocalData.MonteCarloOpt_bcorr;
    plot(1:WinNum,Data.stim_info_bcorr,'LineWidth',2, 'Color',ColorCode(5,:))
    hold on
    plot(1:WinNum,Data.category_info_bcorr,'LineWidth',2, 'Color',ColorCode(3,:))
    hold on
    plot(1:WinNum_cumInfo, LocalData.MonteCarloOpt_raw, 'LineWidth',2, 'Color',ColorCode(6,:))
    hold on
    plot(1:WinNum_cumInfo,LocalData.MonteCarloOpt_bcorr, 'LineWidth',2, 'Color',ColorCode(1,:))
    hold on
    plot(1:WinNum_cumInfo, Biais_MC, 'LineWidth',2, 'Color', [ColorCode(6,:) 0.5])
    hold on
    line(1:WinNum_cumInfo, Data.stim_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color',[0.8 0.5 0.3])
    hold on
    line(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color','r')
    legend('Stimulus Information','Semantic Information','Cumulative Information','Biais corrected Cumulative Information', 'Biais on cumulative information', 'Stimulus Entropy','Semantic category Entropy', 'Location','NorthEast');
    hold on
    line([0 WinNum_cumInfo], [0 0], 'LineStyle','-.','Color','k')
    hold on
    Xtickposition=get(gca,'XTick');
    set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
    xlabel('Time (ms)')
    ylabel('Information (bits)')
    shadedErrorBar([],Data.stim_info_bcorr, Data.stim_info_err,{'Color',ColorCode(5,:), 'LineStyle','-', 'LineWidth',1},1)
    hold on
    shadedErrorBar([],Data.category_info_bcorr, Data.category_info_err,{'Color',ColorCode(3,:), 'LineStyle','-', 'LineWidth',1},1)
    hold on
    shadedErrorBar([],LocalData.MonteCarloOpt_raw, LocalData.MonteCarloOpt_err,{'Color',ColorCode(6,:), 'LineStyle','-', 'LineWidth',1},1)
    hold on
    shadedErrorBar([],LocalData.MonteCarloOpt_bcorr, LocalData.MonteCarloOpt_err,{'Color',ColorCode(1,:), 'LineStyle','-', 'LineWidth',1},1)
    hold on
    plot(1:WinNum_cumInfo, Data.stim_entropy(1:WinNum_cumInfo), 'LineWidth',1, 'LineStyle','-.','Color', [0.8 0.5 0.3])
    hold on
    line(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color','r')
    hold off
    if strcmp(ParamModel.CIType, 'CIS')
        title('Stimulus cumulative information')
    elseif strcmp(ParamModel.CIType, 'CIC')
        title('Semantic cumulative information')
    elseif strcmp(ParamModel.CIType, 'CICRand')
        title('Cumulative information for random categories')
    elseif strcmp(ParamModel.CIType, 'CICR')
        title('Semantic Cumulative information with fixed rates')
    elseif strcmp(ParamModel.CIType, 'CISR')
        title('Stimulus  Cumulative information with fixed rates')
    end
        
        
%     subplot(3,1,2)
%     Biais_MC = Data.cum_info_cat.MonteCarloOpt_raw - Data.cum_info_cat.MonteCarloOpt_bcorr;
%     plot(1:WinNum,Data.category_info_bcorr,'LineWidth',2, 'Color',ColorCode(5,:))
%     hold on
%     plot(1:WinNum_cumInfo, Data.cum_info_cat.MonteCarloOpt_raw, 'LineWidth',2, 'Color',ColorCode(6,:))
%     hold on
%     plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_bcorr, 'LineWidth',2, 'Color',ColorCode(1,:))
%     hold on
%     plot(1:WinNum_cumInfo, Biais_MC, 'LineWidth',2, 'Color', [ColorCode(6,:) 0.5])
%     hold on
%     line(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineStyle','-.','Color','r')
%     legend('Information','Cumulative Information','Biais corrected Cumulative Information', 'Biais on cumulative information', 'Information upper-bound given dataset size', 'Location','NorthEast');
%     hold on
%     line([0 WinNum_cumInfo], [0 0], 'LineStyle','-.','Color','k')
%     hold on
%     ylim([-0.5 log2(NbStims)+2])
%     Xtickposition=get(gca,'XTick');
%     set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
%     xlabel('Time (ms)')
%     ylabel('Information (bits)')
%     shadedErrorBar([],Data.category_info_bcorr, Data.category_info_err,{'Color',ColorCode(5,:), 'LineStyle','-', 'LineWidth',1},1)
%     hold on
%     shadedErrorBar([],Data.cum_info_cat.MonteCarloOpt_raw, Data.cum_info_cat.MonteCarloOpt_err,{'Color',ColorCode(6,:), 'LineStyle','-', 'LineWidth',1},1)
%     hold on
%     shadedErrorBar([],Data.cum_info_cat.MonteCarloOpt_bcorr, Data.cum_info_cat.MonteCarloOpt_err,{'Color',ColorCode(1,:), 'LineStyle','-', 'LineWidth',1},1)
%     hold on
%     plot(1:WinNum_cumInfo, Data.category_entropy(1:WinNum_cumInfo), 'LineWidth',1, 'LineStyle','-.','Color', [0.8 0.5 0.3])
%     hold off
%     title('Semantic cumulative information')
%     
%     subplot(3,1,3)
%     plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_raw*100 ./ Data.cum_info_stim.MonteCarloOpt_raw,'LineWidth',2, 'Color',ColorCode(1,:))
%     hold on
%     plot(1:WinNum_cumInfo,Data.cum_info_cat.MonteCarloOpt_bcorr*100 ./ Data.cum_info_stim.MonteCarloOpt_bcorr,'LineWidth',2, 'Color',ColorCode(1,:), 'LineStyle', '--')
%     hold on
%     plot(Data.category_entropy*100 ./ Data.stim_entropy, 'LineStyle','-.','Color','r')
%     legend('% Semantic Information', ' % Semantic  Information biais corrected',' % Semantic Information upper-bound given dataset size', 'Location','NorthEast');
%     hold on
%     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
%     hold on
%     ylim([-0.5 log2(NbStims)+1])
%     Xtickposition=get(gca,'XTick');
%     set(gca,'XTickLabel', Xtickposition*ParamModel.NeuroBin)
%     xlabel('Time ms')
%     ylabel('% Semantic Cumulative Information')
%     title('Percentage of cumulative information about semantic categories')
%     hold off
end


%% get rid of temporary files for parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(MyParPool);
    system(['rm -r ' parcluster.JobStorageLocation])
end

end