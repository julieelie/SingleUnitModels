function [ParamModel, Data, InputData, Wins]=info_cuminfo_callsemantic(PSTH,JackKnifeTrials,VocType, ParamModel,  Calfilename)
FIG=0;
if nargin<4
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
if ~isfield(ParamModel, 'NbBoot_Info') || isempty(ParamModel.NbBoot_Info)
    ParamModel.NbBoot_Info = 16;
end

if ~isfield(ParamModel, 'NbBoot_CumInfo') || isempty(ParamModel.NbBoot_CumInfo)
    ParamModel.NbBoot_CumInfo = 16;
end

% Set parameters for the number of samples that should be tested in the
% MonteCarlo estimation of the cumulative information
if ~isfield(ParamModel, 'NumSamples_MC_Cum_Info')
    ParamModel.NumSamples_MC_Cum_Info = 10^6; %Set the number of samples for the Monte Carlo approximation of the cumulative information 10^7 takes too much memory prefers lower numbers
elseif isempty(ParamModel.NumSamples_MC_Cum_Info)
    fprintf(1,'No Monte Carlo approximation of the cumulative information will be calculated\n');
end

% Set the Parameters of the Markov approximation of the cumulative
% information
if ~isfield(ParamModel, 'MarkovParameters_Cum_Info')
    ParamModel.MarkovParameters_Cum_Info = [2 3 4;1 1 1];
elseif isempty(ParamModel.MarkovParameters_Cum_Info)
    fprintf(1,'No Markov approximation of the cumulative information will be calculated\n');
end

% Set the fix time history to calculate the exact cumulative information
% (above 4 will certainly bug the computer asking for too big matrices)
if ~isfield(ParamModel, 'ExactHist')
    ParamModel.ExactHist = 4;
elseif isempty(ParamModel.ExactHist)
    fprintf(1,'There will be no exact calculation of cumulative information\n');
end

if nargin<5
    saveonline = 0;
else
    saveonline = 1;
end

% define the list of end points of spectrogram and neural responses windows
Wins = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin;
Wins_cumInfo = ParamModel.MinWin:ParamModel.Increment:ParamModel.MaxWin_cumInfo;

% # of models to run on the data
WinNum = length(Wins);
WinNum_cumInfo = length(Wins_cumInfo);

% Number of stims in the data set
NbStims = length(PSTH);

% Number of stimulus categories
IdCats = unique(VocType);
NbCat = length(IdCats);

%% Configure Parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    MyParPool = parpool(str2num(getenv('SLURM_CPUS_ON_NODE')),'IdleTimeout', Inf);
    system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
    [~,JobID] = system('echo $SLURM_JOB_ID');
    parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];    
end

%% Initialize output variables
Rate4InfoStim = nan(NbStims,WinNum);
Rate4InfoStim_Boot = cell(1,ParamModel.NbBoot_Info);
InputData.VocType = VocType;
Stim_entropy = nan(1,WinNum);
Category_entropy = nan(1,WinNum);
Stim_info = nan(1,WinNum);
Stim_info_JKBoot = nan(ParamModel.NbBoot_Info,WinNum);
Category_info = nan(1,WinNum);
Category_info_JKBoot = nan(ParamModel.NbBoot_Info,WinNum);

P_YgivenS = cell(1,WinNum);
P_YgivenC = cell(1,WinNum);
P_YgivenS_Bootstrap = cell(ParamModel.NbBoot_CumInfo,WinNum);
P_YgivenC_Bootstrap = cell(ParamModel.NbBoot_CumInfo,WinNum);

%% Now loop through bins and calculate spike rates and instantaneous information
%parfor
for ww = 1:WinNum
    Tstart=tic;
    fprintf(1,'Instantaneous info exact spike patterns %d/%d window\n', ww, WinNum);
    Win = Wins(ww);
    FirstTimePoint = (Win - ParamModel.NeuroBin+ ParamModel.ResDelay)*ParamModel.Response_samprate/1000 +1;
    LastTimePoint = (Win + ParamModel.ResDelay)*ParamModel.Response_samprate/1000;
     
    % Calculating info about the stims and the categories for actual values of spike rates
    [Local_Output] = info_model_Calculus_wrapper2(PSTH, FirstTimePoint, LastTimePoint,VocType,ParamModel.Response_samprate);
    Rate4InfoStim(:,ww) = Local_Output.InputdataStim;
    Stim_info(ww) = Local_Output.stim_value;
    Category_info(ww) = Local_Output.cat_value;
    P_YgivenS{ww} = Local_Output.P_YgivenS;
    P_YgivenC{ww} = Local_Output.P_YgivenC;
    Stim_entropy(ww) = Local_Output.stim_entropy;
    Category_entropy(ww) = Local_Output.cat_entropy; 
   fprintf('Instantaneous Info: Done bin %d/%d after %f sec\n', ww, WinNum, toc(Tstart));
end
    
% Bootstrapping the calculation of information with Jackknife estimations of spike rate
%parfor
for bb=1:ParamModel.NbBoot_Info
    fprintf(1,'%d/%d bootstrap instantaneous info with Jackknife estimates of spike rates\n', bb, ParamModel.NbBoot_Info);
    
    % Choosing a different set of JK trials for the stims for each bootstrap
    NbStim = length(JackKnifeTrials);
    JackKnifePSTH = cell(1,NbStim);
    Rate4InfoStim_Boot{bb} = nan(NbStims,WinNum);
    for st = 1:NbStim
        PSTH_Local = JackKnifeTrials{st};
        NJK = size(PSTH_Local,1);
        JackKnifePSTH{st} = PSTH_Local(randperm(NJK,1),:);
    end

    % Then run the calculation of information on all windows
    for ww_in = 1:WinNum
        fprintf(1,'JK instantaneous info Boostrap %d %d/%d window\n',bb, ww_in, WinNum);
        Win = Wins(ww_in);
        FirstTimePoint = (Win - ParamModel.NeuroBin+ ParamModel.ResDelay)*ParamModel.Response_samprate/1000 +1;
        LastTimePoint = (Win + ParamModel.ResDelay)*ParamModel.Response_samprate/1000;        
        
        [Local_Output] = info_model_Calculus_wrapper2(JackKnifePSTH, FirstTimePoint, LastTimePoint, VocType, ParamModel.Response_samprate);
        
        Rate4InfoStim_Boot{bb}(:,ww_in) = Local_Output.InputdataStim;
        Stim_info_JKBoot(bb,ww_in) = Local_Output.stim_value;
        Category_info_JKBoot(bb,ww_in) = Local_Output.cat_value;
        if bb <= ParamModel.NbBoot_CumInfo
            P_YgivenS_Bootstrap{bb,ww_in} = Local_Output.P_YgivenS;
            P_YgivenC_Bootstrap{bb,ww_in} = Local_Output.P_YgivenC;
        end
    end
end

% Estimate information for infinite number of trials
% ordonnée à l'origine: b = (y1*x2 - y2*x1)/(x2-x1)
Data.stim_info_infT = (Stim_info./ParamModel.Mean_Ntrials_perstim(2) - mean(Stim_info_JKBoot,1)./ParamModel.Mean_Ntrials_perstim(1))./(1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
Data.category_info_infT = (Category_info./ParamModel.Mean_Ntrials_perstim(2) - mean(Category_info_JKBoot,1)./ParamModel.Mean_Ntrials_perstim(1))./(1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
Data.stim_info_JKBoot_std=std(Stim_info_JKBoot);
Data.category_info_JKBoot_std=std(Category_info_JKBoot);

Data.stim_info_JKBoot_infT = (repmat(Stim_info ./ ParamModel.Mean_Ntrials_perstim(2), ParamModel.NbBoot_Info,1) - Stim_info_JKBoot ./ ParamModel.Mean_Ntrials_perstim(1)) ./ (1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
Data.category_info_JKBoot_infT = (repmat(Category_info ./ ParamModel.Mean_Ntrials_perstim(2), ParamModel.NbBoot_Info,1) - Category_info_JKBoot ./ ParamModel.Mean_Ntrials_perstim(1)) ./ (1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
Data.stim_info_infT_std = std(Data.stim_info_JKBoot_infT);
Data.category_info_infT_std = std(Data.category_info_JKBoot_infT);


% Stuff in results in structure
InputData.VocType = VocType;
InputData.Rate4InfoStim_Boot = Rate4InfoStim_Boot;
InputData.Rate4InfoStim = Rate4InfoStim;

Data.P_YgivenS_Bootstrap = P_YgivenS_Bootstrap;
Data.P_YgivenC_Bootstrap = P_YgivenC_Bootstrap;

Data.P_YgivenS = P_YgivenS;
Data.P_YgivenC = P_YgivenC;

Data.stim_entropy = Stim_entropy;
Data.category_entropy = Category_entropy;
Data.stim_info = Stim_info;
Data.stim_info_JKBoot = Stim_info_JKBoot;
Data.category_info = Category_info;
Data.category_info_JKBoot = Category_info_JKBoot;

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
     for ss=1:NbStims
        plot(InputData.Rate4InfoStim(ss,:)./10,'LineWidth',2, 'Color','g')
        hold on
        Local_bootrate = nan(ParamModel.NbBoot_CumInfo, size(InputData.Rate4InfoStim_Boot{1},2));
        for bb=1:ParamModel.NbBoot_CumInfo
            plot(InputData.Rate4InfoStim_Boot{bb}(ss,:)./10, 'Color','k')
            hold on
            if bb==1
                plot(mean(Local_bootrate,1), 'LineWidth',2,'Color','r')%this is just for the legend
                hold on
                legend('Actual spike rate','individual bootstrap', 'Average bootstrapped spike rate')
            end
            Local_bootrate(bb,:) = InputData.Rate4InfoStim_Boot{bb}(ss,:);
        end
        plot(mean(Local_bootrate,1)./10, 'LineWidth',2,'Color','r')
        hold off
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
        xlabel('Time ms')
        ylabel('Spike rate /ms')
        title(sprintf('Stim %d/%d',ss, NbStims))
        pause(1)
     end
     
    figure()
    plot(var(InputData.Rate4InfoStim),'LineWidth',2)
    hold on
    for bb=1:ParamModel.NbBoot_CumInfo
            plot(var(InputData.Rate4InfoStim_Boot{bb}))
            if bb==1
                legend('Actual', 'Bootstrapped')
            end
            hold on
    end
    hold off
    Xtickposition=get(gca,'XTick');
    set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
    xlabel('Time ms')
    ylabel('Stimulus spike rate variance')
    pause()
    
         
     
     
     figure()
     subplot(2,1,1)
     plot(Data.stim_info,'LineWidth',2, 'Color',[0 0 0])
     hold on
     plot(Data.stim_info_infT,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(mean(Data.stim_info_JKBoot),'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Info 4 Infinite # Trials', 'Info 4 JK Trials', 'Location','NorthEast');
     hold on
     plot(Data.stim_info_infT + 2*Data.stim_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.stim_info_infT - 2*Data.stim_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
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
     plot(Data.stim_info_infT,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(mean(Data.stim_info_JKBoot),'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Info 4 Infinite # Trials', 'Info 4 JK Trials', 'Location','NorthEast');
     hold on
     plot(Data.stim_info_infT + 2*Data.stim_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.stim_info_infT - 2*Data.stim_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
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
     plot(Data.category_info_infT,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(mean(Data.category_info_JKBoot),'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Info 4 Infinite # Trials', 'Info 4 JK Trials', 'Location','NorthEast');
     hold on
     plot(Data.category_info_infT + 2*Data.category_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.category_info_infT - 2*Data.category_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
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
     plot(Data.category_info_infT,'LineWidth',2, 'Color',[0 0 1])
     hold on
     plot(mean(Data.category_info_JKBoot),'LineWidth',2, 'Color',[0 1 0])
     legend('Information', 'Info 4 Infinite # Trials', 'Info 4 JK Trials', 'Location','NorthEast');
     hold on
     plot(Data.category_info_infT + 2*Data.category_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     plot(Data.category_info_infT - 2*Data.category_info_infT_std, 'LineWidth',2,'Color',[0 1 0 0.4], 'LineStyle','--')
     hold on
     line([0 WinNum], [0 0], 'LineStyle','-.','Color','k')
     hold off
     YL = get(gca,'YLim');
     ylim([-0.2 YL(2)])
     Xtickposition=get(gca,'XTick');
     set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
     xlabel('Time ms')
     ylabel('Category Information in bits')
     pause()
 end
%% Now calculating cumulative information
% Getting output variables ready
if isempty(ParamModel.ExactHist) && isempty(ParamModel.MarkovParameters_Cum_Info) && isempty(ParamModel.NumSamples_MC_Cum_Info)
    % skip the next calculation we're not running the calculation of
    % cumulative information
else
    fprintf(1,'Starting calculation of cumulative information\n')
    Data.cum_info_stim = struct();
    Data.cum_info_cat = struct();
    
    if ~isempty(ParamModel.MarkovParameters_Cum_Info)
        HY_Markov_stim = struct();
        HY_Markov_cat = struct();
        for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
            ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
            Data.cum_info_stim.(sprintf('%s',ModelType))= nan(1,WinNum_cumInfo);
            Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
            Data.cum_info_cat.(sprintf('%s',ModelType))= nan(1,WinNum_cumInfo);
            Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
            HY_Markov_stim.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
            HY_Markov_cat.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
        end
    end
    if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
        for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
            ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
            Data.cum_info_stim.(sprintf('%s',ModelType))= nan(2,WinNum_cumInfo);
            Data.cum_info_stim.(sprintf('%s',ModelType))(1,1)= Data.stim_info(1);
            Data.cum_info_cat.(sprintf('%s',ModelType))= nan(2,WinNum_cumInfo);
            Data.cum_info_cat.(sprintf('%s',ModelType))(1,1)= Data.category_info(1);
        end
    end
    if ~isempty(ParamModel.ExactHist)
        ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
        Data.cum_info_stim.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
        Data.cum_info_cat.(sprintf('%s',ModelType)) = nan(1,WinNum_cumInfo);
        Data.cum_info_stim.(sprintf('%s',ModelType))(1) = Data.stim_info(1);
        Data.cum_info_cat.(sprintf('%s',ModelType))(1) = Data.category_info(1);
    end
    
    
    %% Running through time bins and calculate cumulative information
    for ww =2:WinNum_cumInfo
        Tstart2=tic;
        if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
            for ss=1:length(ParamModel.NumSamples_MC_Cum_Info)
                ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
                % Monte Carlo estimation with full memory cumulative information stimuli
                [Data.cum_info_stim.(sprintf('%s',ModelType))(:,ww),~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
                % Monte Carlo estimation with full memory cumulative information
                % categories
                [Data.cum_info_cat.(sprintf('%s',ModelType))(:,ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',2,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(ss));
            end
        end
        
        if ~isempty(ParamModel.ExactHist)
            ModelType = sprintf('ExactMem0_%d',ParamModel.ExactHist);
            % Exact calculation cumulative information on stims with ParamModel.ExactHist*10 ms memory
            [Data.cum_info_stim.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
            % Exact calculation cumulative information on categories with ParamModel.ExactHist*10 ms memory
            [Data.cum_info_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',2,'CalMode','Exact_Mem', 'Exact_history',ParamModel.ExactHist);
        end
        
        if ~isempty(ParamModel.MarkovParameters_Cum_Info)
            for ss=1:size(ParamModel.MarkovParameters_Cum_Info,2)
                ModelType = sprintf('MarkovEst%d',ParamModel.MarkovParameters_Cum_Info(1,ss));
                if ww>2
                    % Markov chain estimation of cumulative information on stims
                    [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss), 'HY_old', HY_Markov_stim.(sprintf('%s',ModelType))(ww-1));
                    % Markov chain estimation of cumulative information on
                    % categories
                    [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss),'HY_old', HY_Markov_cat.(sprintf('%s',ModelType))(ww-1));
                else
                    % Markov chain estimation of cumulative information on stims
                    [Data.cum_info_stim.(sprintf('%s',ModelType))(ww), HY_Markov_stim.(sprintf('%s',ModelType))(ww), ~]=info_cumulative_model_Calculus(Data.P_YgivenS(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
                    % Markov chain estimation of cumulative information on
                    % categories
                    [Data.cum_info_cat.(sprintf('%s',ModelType))(ww), HY_Markov_cat.(sprintf('%s',ModelType))(ww),~]=info_cumulative_model_Calculus(Data.P_YgivenC(1:ww),'Model#',1,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,ss));
                end
            end
        end
        fprintf('Done cumulative information bin %d/%d after %f sec\n', ww, WinNum_cumInfo, toc(Tstart2));
    end
    
    %% Save what we have for now
    if saveonline
        save(Calfilename,'Data','ParamModel','-append');
    end
    
    %% Plot the cumulative information
    if FIG
        figure()
        subplot(2,1,1)
        plot(1:WinNum_cumInfo, Data.stim_info(1:WinNum_cumInfo), 1:WinNum_cumInfo,Data.cum_info_stim.MarkovEst4, 1:WinNum_cumInfo,Data.cum_info_stim.MarkovEst3, 1:WinNum_cumInfo,Data.cum_info_stim.MarkovEst2,1:WinNum_cumInfo,Data.cum_info_stim.MonteCarlo6(1,:), 1:WinNum_cumInfo,Data.cum_info_stim.MonteCarlo5(1,:), 1:WinNum_cumInfo,Data.cum_info_stim.MonteCarlo4(1,:), 1:WinNum_cumInfo,Data.cum_info_stim.MonteCarlo3(1,:),1:WinNum_cumInfo,Data.cum_info_stim.MonteCarlo2(1,:), 1:WinNum_cumInfo, Data.cum_info_stim.ExactMem0_4, 1:WinNum_cumInfo, cumsum(Data.stim_info(1:WinNum_cumInfo)), 'LineWidth',2)
        Color_Lines = get(gca,'ColorOrder');
        legend('Information', 'Cumulative Info Markov4', 'Cumulative Info Markov3','Cumulative Info Markov2','Cumulative Information MC 10^6','Cumulative Information MC 10^5','Cumulative Information MC 10^4','Cumulative Information MC 10^3', 'Cumulative Information MC 10^2', 'Exact Cumulative Information 40ms history', 'Cumulative sum of info', 'Location', 'NorthWest')
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
        xlabel('Time ms')
        ylabel('Stimulus Information in bits')
        subplot(2,1,2)
        plot(1:WinNum_cumInfo, Data.category_info(1:WinNum_cumInfo), 1:WinNum_cumInfo,Data.cum_info_cat.MarkovEst4, 1:WinNum_cumInfo,Data.cum_info_cat.MarkovEst3, 1:WinNum_cumInfo,Data.cum_info_cat.MarkovEst2,1:WinNum_cumInfo,Data.cum_info_cat.MonteCarlo6(1,:), 1:WinNum_cumInfo,Data.cum_info_cat.MonteCarlo5(1,:), 1:WinNum_cumInfo,Data.cum_info_cat.MonteCarlo4(1,:), 1:WinNum_cumInfo,Data.cum_info_cat.MonteCarlo3(1,:),1:WinNum_cumInfo,Data.cum_info_cat.MonteCarlo2(1,:), 1:WinNum_cumInfo, Data.cum_info_cat.ExactMem0_4, 1:WinNum_cumInfo, cumsum(Data.category_info(1:WinNum_cumInfo)), 'LineWidth',2)
        Color_Lines = get(gca,'ColorOrder');
        legend('Information', 'Cumulative Info Markov4', 'Cumulative Info Markov3','Cumulative Info Markov2','Cumulative Information MC 10^6','Cumulative Information MC 10^5','Cumulative Information MC 10^4','Cumulative Information MC 10^3', 'Cumulative Information MC 10^2', 'Exact Cumulative Information 40ms history', 'Cumulative sum of info', 'Location', 'NorthWest')
        Xtickposition=get(gca,'XTick');
        set(gca,'XTickLabel', Xtickposition*ParamModel.Increment)
        xlabel('Time ms')
        ylabel('Category Information in bits')
    end
    %% Calculating bootstrap values accross stims and across trials within stims
    fprintf(1,'Bootstraping calculation of cumulative information\n')
    if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
        Cum_info_stim_LastMonteCarlo_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        Cum_info_cat_LastMonteCarlo_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        Cum_info_stim_LastMonteCarlo_Bootstrap(:,1) = repmat(Data.stim_info(1),ParamModel.NbBoot_CumInfo,1);
        Cum_info_cat_LastMonteCarlo_Bootstrap(:,1) = repmat(Data.category_info(1),ParamModel.NbBoot_CumInfo,1);
        %Cum_info_stim_LastMonteCarlo_Bootsample = nan(ParamModel.NbBoot_cumInfo,WinNum_cumInfo);
        %Cum_info_cat_LastMonteCarlo_Bootsample = nan(ParamModel.NbBoot_cumInfo,WinNum_cumInfo);
    end
    
    if ~isempty(ParamModel.ExactHist)
        Cum_info_stim_ExactMem0_4_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        Cum_info_cat_ExactMem0_4_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
    end
    
    if ~isempty(ParamModel.MarkovParameters_Cum_Info)
        Cum_info_stim_LastMarkov_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
        Cum_info_cat_LastMarkov_Bootstrap = nan(ParamModel.NbBoot_CumInfo,WinNum_cumInfo);
    end
    
%     if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
%         delete(gcp)
%         parpool(str2num(getenv('SLURM_CPUS_ON_NODE')));
%         system('mkdir -p /global/scratch/$USER/$SLURM_JOB_ID')
%         [~,JobID] = system('echo $SLURM_JOB_ID');
%         parcluster.JobStorageLocation = ['/global/scratch/jelie/' JobID];
%     end
    
    %parfor
    for bb=1:ParamModel.NbBoot_CumInfo
        Tstart3=tic;
        fprintf('Bootstrap CumInfo %d/%d\n', bb, ParamModel.NbBoot_CumInfo);
        if ~isempty(ParamModel.MarkovParameters_Cum_Info)
            HY_Markov4_stim_bbstim = nan(1,WinNum_cumInfo);
            HY_Markov4_cat_bbstim = nan(1,WinNum_cumInfo);
        end
        for ww=2:WinNum_cumInfo
            fprintf('Bootstrap CumInfo %d/%d Time point %d/%d\n',bb, ParamModel.NbBoot_CumInfo, ww, WinNum_cumInfo);
            
            % First for the cumulative information about stimuli
            P_YgivenS_local = P_YgivenS_Bootstrap(bb,1:ww);
            
            if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
                % Monte Carlo estimation with full memory
                [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenS_local,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
                Cum_info_stim_LastMonteCarlo_Bootstrap(bb,ww)=Icum_EstMonteCarlo_temp(1);
                % if you want to restore the bootstrap of MC on same
                % estimations of spike rate restore the following lines
                %[Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(Data.P_YgivenS(1,1:ww),'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
                %Cum_info_stim_LastMonteCarlo_Bootsample(bb,ww) = Icum_EstMonteCarlo_temp(1);
            end
            
            if ~isempty(ParamModel.ExactHist)
                % Exact calculation with 50 ms memory
                [Cum_info_stim_ExactMem0_4_Bootstrap(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',4);
            end
            
            if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                if ww==2
                    % Markov chain estimation 50 ms
                    [Cum_info_stim_LastMarkov_Bootstrap(bb,ww),HY_Markov4_stim_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
                else
                    % Markov chain estimation 50 ms
                    [Cum_info_stim_LastMarkov_Bootstrap(bb,ww),HY_Markov4_stim_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov4_stim_bbstim(ww-1));
                end
            end
            
            % Then same thing for the cumulative information about categories
            P_YgivenC_local = P_YgivenC_Bootstrap(bb,1:ww);
            
            if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
                % Monte Carlo estimation with full memory
                [Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(P_YgivenC_local,'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
                Cum_info_cat_LastMonteCarlo_Bootstrap(bb,ww)=Icum_EstMonteCarlo_temp(1);
                
                % if you want to restore the bootstrap of MC on same
                % estimations of spike rate restore the following lines
                %[Icum_EstMonteCarlo_temp,~]=info_cumulative_model_Calculus(Data.P_YgivenC(1,1:ww),'Model#',bb,'CalMode','MonteCarlo', 'MCParameter',ParamModel.NumSamples_MC_Cum_Info(end));
                %Cum_info_cat_LastMonteCarlo_Bootsample(bb,ww) = Icum_EstMonteCarlo_temp(1);
            end
            
            if ~isempty(ParamModel.ExactHist)
                % Exact calculation with 40 ms memory
                [Cum_info_cat_ExactMem0_4_Bootstrap(bb,ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'Model#',bb,'CalMode','Exact_Mem', 'Exact_history',4);
            end
            
            if ~isempty(ParamModel.MarkovParameters_Cum_Info)
                if ww==2
                    % Markov chain estimation 50 ms
                    [Cum_info_cat_LastMarkov_Bootstrap(bb,ww),HY_Markov4_cat_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end));
                else
                    % Markov chain estimation 50 ms
                    [Cum_info_cat_LastMarkov_Bootstrap(bb,ww),HY_Markov4_cat_bbstim(ww),~]=info_cumulative_model_Calculus(P_YgivenC_local,'CalMode','MarkovChain', 'MarkovParameters',ParamModel.MarkovParameters_Cum_Info(:,end), 'HY_old', HY_Markov4_cat_bbstim(ww-1));
                end
            end
        end
        fprintf('Done bootstrap %d on cumulative information after %f sec\n', bb, toc(Tstart3));
    end
    
    if ~isempty(ParamModel.NumSamples_MC_Cum_Info)
        ModelTypeMCstim = sprintf('MonteCarlo%d_JKBootstrap',log10(ParamModel.NumSamples_MC_Cum_Info(end)));
        Data.cum_info_stim.(sprintf('%s',ModelTypeMCstim)) = Cum_info_stim_LastMonteCarlo_Bootstrap;
        Data.cum_info_cat.(sprintf('%s',ModelTypeMCstim)) = Cum_info_cat_LastMonteCarlo_Bootstrap;

        % Estimate cumulative information for infinite number of trials and
        % stuff in results in structure
            % ordonnée à l'origine: b = (y1*x2 - y2*x1)/(x2-x1)
        ModelTypeMC = sprintf('MonteCarlo%d_',log10(ParamModel.NumSamples_MC_Cum_Info(end)));
        ModelType=sprintf('MonteCarlo%d',log10(ParamModel.NumSamples_MC_Cum_Info(ss)));
        
        Data.cum_info_stim.(sprintf('%sinfT',ModelTypeMC)) = (Data.cum_info_stim.(sprintf('%s',ModelType))(1,:) ./ ParamModel.Mean_Ntrials_perstim(2) - mean(Cum_info_stim_LastMonteCarlo_Bootstrap,1)./ParamModel.Mean_Ntrials_perstim(1))./(1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
        Data.cum_info_stim.(sprintf('%sJK_std',ModelTypeMC)) = std(Cum_info_stim_LastMonteCarlo_Bootstrap);
        Data.cum_info_cat.(sprintf('%sinfT',ModelTypeMC)) = (Data.cum_info_cat.(sprintf('%s',ModelType))(1,:)./ParamModel.Mean_Ntrials_perstim(2) - mean(Cum_info_cat_LastMonteCarlo_Bootstrap,1)./ParamModel.Mean_Ntrials_perstim(1))./(1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
        Data.cum_info_cat.(sprintf('%sJK_std',ModelTypeMC)) = std(Cum_info_cat_LastMonteCarlo_Bootstrap);
        Data.cum_info_stim.(sprintf('%sJKBoot_infT',ModelTypeMC)) = (repmat(Data.cum_info_stim.(sprintf('%s',ModelType))(1,:) ./ ParamModel.Mean_Ntrials_perstim(2), ParamModel.NbBoot_CumInfo,1) - Cum_info_stim_LastMonteCarlo_Bootstrap ./ ParamModel.Mean_Ntrials_perstim(1)) ./ (1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
        Data.cum_info_cat.(sprintf('%sJKBoot_infT',ModelTypeMC)) = (repmat(Data.cum_info_cat.(sprintf('%s',ModelType))(1,:) ./ ParamModel.Mean_Ntrials_perstim(2), ParamModel.NbBoot_CumInfo,1) - Cum_info_cat_LastMonteCarlo_Bootstrap ./ ParamModel.Mean_Ntrials_perstim(1)) ./ (1/ParamModel.Mean_Ntrials_perstim(2) - 1/ParamModel.Mean_Ntrials_perstim(1));
        Data.cum_info_stim.(sprintf('%sJKBoot_infT_std',ModelTypeMC)) = std(Data.cum_info_stim.(sprintf('%sJKBoot_infT',ModelTypeMC)));
        Data.cum_info_cat.(sprintf('%sJKBoot_infT_std',ModelTypeMC)) = std(Data.cum_info_cat.(sprintf('%sJKBoot_infT',ModelTypeMC)));
        
        % if you want to restore the bootstrap of MC on same
        % estimations of spike rate restore the following lines
        %ModelTypeMCsample = sprintf('MonteCarlo%d_Bootsample',log10(ParamModel.NumSamples_MC_Cum_Info(end)));
        %Data.cum_info_stim.(sprintf('%s',ModelTypeMCsample)) = Cum_info_stim_LastMonteCarlo_Bootsample;
        %Data.cum_info_cat.(sprintf('%s',ModelTypeMCsample)) = Cum_info_cat_LastMonteCarlo_Bootsample;
    end
    
    if ~isempty(ParamModel.ExactHist)
        Data.cum_info_stim.ExactMem0_4_Bootstim = Cum_info_stim_ExactMem0_4_Bootstrap;
        Data.cum_info_cat.ExactMem0_4_Bootstim = Cum_info_cat_ExactMem0_4_Bootstrap;
    end
    
    if ~isempty(ParamModel.MarkovParameters_Cum_Info)
        ModelTypeMarstim = sprintf('MarkovEst%d_Bootstim',ParamModel.MarkovParameters_Cum_Info(1,end));
        Data.cum_info_stim.(sprintf('%s',ModelTypeMarstim)) = Cum_info_stim_LastMarkov_Bootstrap;
        Data.cum_info_cat.(sprintf('%s',ModelTypeMarstim)) = Cum_info_cat_LastMarkov_Bootstrap;

    end
end




%% get rid of temporary files for parallel computing
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))
    delete(MyParPool);
    system(['rm -r ' parcluster.JobStorageLocation])
end


%% Save what we have for now
 if saveonline
     save(Calfilename,'Data','ParamModel','-append');
 end
end