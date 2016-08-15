function [Icum, HY, HYgivenS]=info_cumulative_model_Calculus(P_YgivenS,varargin)
tic
%% This function calculates the cumulative information of the neural response
% P_YgivenS is a cell array of length equals to the number of time points t in ...
... the neural responses. Each cell contains a c*N matrix of probabilities,
    ... the conditional probability of a given neural spike count c at the
    ... particular time point t for the stimulus N. The cumulative
    ... information is calculated/estimated for the last time point from
    ... the first time point.
    
%   [...] =
%   info_cumulative_model_Calculus(...,'PARAM1',VAL1,'PARAM2',VAL2,...)
%   specifies various parameters of the calculation.
% Valid parameters are the following:
%  Parameter      Value
%  'CalMode'    The type of calculation used to calculate the
%         exact entropy of neural responses in the cumulative information.
%         There is currently 3 exact calculations
%         that each have their own inconveniences: 'Exact_Mem' uses the
%         computer memory and will crash if the dataset is too large;
%         'Exact_Paths' do all the calculations one by one and will take
%         infinite time if the dataset is too large; 'Exact_HardDrive' will
%         also perfom all the calculations, but saving temporary
%         calculations on the hard drive, eventualy crashing for large
%         dataset. There are also 2 approximations: 'MarkovChain' uses the
%         parameters enterd as MarkovParameters to estimate the entropy of
%         the neural response with a history of MarkovParameters(1) and a
%         pace for the Markov chain of MarkovParameters(2); 'MonteCarlo'
%         uses a MonteCarlo estimation with a weight correction, the number
%         of samples used for the estimation is given by MCParameter.
%  'FirstTimePoint' is the identity/name of the first time point from which the 
%                   cumulative information is estimated. Its value
%                   is defaulted to 1.
%  'StimIndicesAll' is a cell array containing for each time point the ID of
%                   the stims (integers) for which a spike count could be calculated
%                   at the corresponding time point.
%  'StimIndicesLast' is a vector of the stims' ID (integers) for the last time point
%  'Model#'         is the index/name of the neuron/model for which the
%                   cumulative information is calculated. Its value is
%                   defaulted to 1.
%  'TempStoragePath' is the path to the folder where the temporary files
%                    generated by the 'Exact_HardDrive' calculations should
%                    be stored
%   'MarkovParameters' is a 2 column vector giving the time history of the
%                      Markov Chain (number of time points in the past used
%                      to calculate the chain as first element and time
%                      point pace as a second element, both are integers).
%   'MCParameter'    is the number of samples in the Monte Carlo estimation
%                    of the entropy of the response
%   'MinProbThresh'  set to 1 to discard probabilities lower than 1/life of
%                    a zebra finch in the lab, set to 0 to take into
%                    acccount all probabilities. Default 0.
%   'DebugFig'      set to 1 to see some debugging figures
%   'MinProb'       threshold of probability under which values are set to
%                   zero. Default is 1/number of 20ms bins in a life of a
%                   zebra finch in the lab (8 years)
%   'HY_old'        The entropy of the response of the previous time point 
%                   can be given for the Markov chain approximation to 
%                   enhance speed of calculation.
%   'ScaleY'        Set to 1 to scale the distribution of conditional
%                   probabilities given the stimulus the same way for both
%                   entropy calculations
%   'Exact_history' is the number of time point in the past that should be
%   taken into account to calculate the cumulative information at time t
%   for the exact calculations Exact_Mem and Exact_HardDrive. It is set to the maximum by default 

%% Sorting input arguments
pnames = {'CalMode', 'FirstTimePoint','StimIndicesAll','StimIndicesLast','Model#','TempStoragePath','MarkovParameters', 'MCParameter','MinProbThresh', 'DebugFig','MinProb','HY_old','ScaleY','Exact_history'};

% Calculating default values of input arguments
Stim_local = 1:size(P_YgivenS{end},2);
StimIndices_AllWin = cell(length(P_YgivenS),1);
for tt=1:length(P_YgivenS)
    StimIndices_AllWin{tt} = Stim_local;
end
getenv('HOSTNAME')
if ~isempty(strfind(getenv('HOSTNAME'),'.savio')) || ~isempty(strfind(getenv('HOSTNAME'),'.brc'))%savio Cluster
    FolderTempInfStorage = '/tmp/LocalTempStorageInfo';
elseif ismac()
    FolderTempInfStorage = '/tmp/LocalTempStorageInfo';
else %we are on strfinator or a cluster machine
    [~, User]=system('whoami');
    FolderTempInfStorage = fullfile('/auto/tdrive', User(1:end-1), 'LocalTempStorageInfo');
end
MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab

% Get input arguments
dflts  = {'Exact_Mem' 1 StimIndices_AllWin Stim_local 1 FolderTempInfStorage [5 1] 1000 0 0 MinProb [] 1 length(P_YgivenS)};
[CalculMode, Firstwin,StimIndices_AllWin, Stim_local,modnb,FolderTempInfStorage, Markov_history, N_MC, MinProbThresh, DebugFig, MinProb,HY_old, ScaleY,ExactMem_history] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%% Set some parameters
Lastwin = length(P_YgivenS);
win=Lastwin-Firstwin+1;

%% Construct the data set of individual p
P_YgivenS_local = cell(win,1);
P_YgivenS_local_rescaled1 = nan(size(P_YgivenS{1},1),length(Stim_local),win);
YIndMax = nan(win,1);

% Gather for each window before the current window win the probabilities of
% the stimuli that are still part of the calculations at the window win
for ww=1:(win-1)
    P_YgivenS_local_w = nan(size(P_YgivenS{1},1),length(Stim_local));
    for ss=1:length(Stim_local)
        st = Stim_local(ss);
        st_local = find(StimIndices_AllWin{ww}==st);
        P_YgivenS_local_rescaled1(:,ss,ww)= P_YgivenS{ww}(:,st_local) ./sum(P_YgivenS{ww}(:,st_local));
        if ScaleY
            P_YgivenS_local_w(:,ss) = P_YgivenS_local_rescaled1(:,ss,ww);
        else
            P_YgivenS_local_w(:,ss) = P_YgivenS{ww}(:,st_local);
        end
        
    end
    
    % Don't keep the probabilities corresponding to least probable spike
    % count values for all stimuli
    if MinProbThresh
        Badyy = find(sum((P_YgivenS_local_w < MinProb),2)==length(Stim_local));
        Goodyy=setdiff(1:size(P_YgivenS{1},1), Badyy);
        YIndMax(ww) = length(Goodyy);
        P_YgivenS_local{ww} = P_YgivenS_local_w(Goodyy,:);
    else
        YIndMax(ww) = size(P_YgivenS_local_w,1);
        P_YgivenS_local{ww} = P_YgivenS_local_w;
    end
        
        
end

P_YgivenS_local_rescaled1(:,:,win)=P_YgivenS{win}(:,:)./repmat(sum(P_YgivenS{win}(:,:),1),size(P_YgivenS{win},1),1);
% Don't keep the probabilities corresponding to least probable spike
% count values for all stimuli at the time point win
if MinProbThresh
    Badyy = find(sum((P_YgivenS{win} < MinProb),2)==length(Stim_local));
    Goodyy=setdiff(1:size(P_YgivenS{1},1), Badyy);
    YIndMax(win) = length(Goodyy);
    if ScaleY
        P_YgivenS_local{win}=P_YgivenS{win}(Goodyy,:)./repmat(sum(P_YgivenS{win}(Goodyy,:),1),size(P_YgivenS{win},1),1);
    else
        P_YgivenS_local{win}=P_YgivenS{win}(Goodyy,:);
    end
else
    YIndMax(win) = size(P_YgivenS{win},1);
    if ScaleY
        P_YgivenS_local{win}=P_YgivenS_local_rescaled1(:,:,win);
    else
        P_YgivenS_local{win}=P_YgivenS{win};
    end
end


%% Calculating conditional entropy
if strcmp(CalculMode, 'Exact_Mem') || strcmp(CalculMode, 'Exact_HardDrive')
    EarliestPastPoint = min(size(P_YgivenS_local_rescaled1,3)-1,ExactMem_history);
    P_YgivenS_local_rescaled1 = P_YgivenS_local_rescaled1(:,:,(end-EarliestPastPoint) : end);
end
HYgivenS = sum(sum(sum(-P_YgivenS_local_rescaled1.*log2(P_YgivenS_local_rescaled1 + (P_YgivenS_local_rescaled1==0)))))/length(Stim_local);
fprintf(1, 'Conditional entropy (entropy of the neural response given the stimulus): %f\n', HYgivenS);

%% First strategy for coding the exact entropy
if strcmp(CalculMode, 'Exact_Mem')
    fprintf(1, 'Calulation of Exact entropy(entropy of the neural response)\nusing the exact calculation and computer memory (crash down if too big)\n');
    if ExactMem_history==length(P_YgivenS_local)
        fprintf(1,'The history of the calculation is set to its max, all info redundancies are investigated\n');
    else
        fprintf(1,'The history of the calculation is set to %d time points or %d ms\nOnly the information redundancies present in that window are taken into account\n',ExactMem_history, ExactMem_history*10);
    end
     % Construct the dataset of p products (all the different possible combinations of possible responses for each stim)
    P_YgivenS_localchunked = P_YgivenS_local((end-EarliestPastPoint) : end);
    P_YgivenS_perstim_wins = insideMultiMat(P_YgivenS_localchunked,MinProbThresh, MinProb);%The dimension of this should be size(P_YgivenS{1},1),length(Stim_local)
    
    % Calculating exact entropy
    P_Y_wins = sum(P_YgivenS_perstim_wins,2)/length(Stim_local);%The dimension of this should be size(P_YgivenS{1},1)^win,1
    
    % rescale p-values so they sum to 1
    P_Y_wins_rescaled = P_Y_wins./sum(P_Y_wins);
    if round(sum(P_Y_wins),10)~=1 && ScaleY && ~MinProbThresh
        fprintf(1, 'WARNING: the distribution of joint probabilities is significantly rescaled (sum=%d) when it should not be given that the probability conditional to the stimulus wee already scales\n',sum(P_Y_wins));
    elseif round(sum(P_Y_wins),10)~=1 && ~ScaleY && ~MinProbThresh
        fprintf(1, 'The distribution of joint probabilities is rescaled (sum=%d), this is expected given that the data are not scaled the same way for the response entropy and the conditional entropy\n',sum(P_Y_wins));
    elseif round(sum(P_Y_wins),10)~=1 && ~ScaleY && MinProbThresh
        fprintf(1, 'The distribution of joint probabilities is rescaled (sum=%d), this is expected given that a threshold was applied on joint probablities to be taken into account in the calculation of the response entropy and distributions of proba are scaled differently for the calculation of the 2 entropies\n',sum(P_Y_wins));
    elseif round(sum(P_Y_wins),10)~=1 && ScaleY && MinProbThresh
        fprintf(1, 'The distribution of joint probabilities is rescaled (sum=%d), this is expected given that a threshold was applied on joint probablities to be taken into account in the calculation of the response entropy\n',sum(P_Y_wins));
    end
    % entropy
    HY = sum(-P_Y_wins_rescaled.*log2(P_Y_wins_rescaled + (P_Y_wins_rescaled==0)));
    fprintf(1, 'Exact entropy(entropy of the neural response)\nusing the exact calculation and computer memory (crash down if too big) : %f\n', HY);
end

%% Second strategy for coding exact entropy
% Summing up the probabilities on the fly
% Define a minimum probability under which the proba should be considered 0
if strcmp(CalculMode, 'Exact_Paths')
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the exact calculation and paths strategy (calculation can be very long)\n');
    % initialize the first path
    Path = ones(win,1);

    % calculate the proba of the path and sum them online until you get to the
    % last path
    NotLastPath = 1;
    HY = 0;
    NbP=0;
    while NotLastPath
        NbP = NbP + 1;
        %fprintf('The path tested is %d %d %d\n', Path);
        P_path = P_YgivenS_local{1}(Path(1),:);
        for ww=2:win
            P_path = P_path .* P_YgivenS_local{ww}(Path(ww),:);
            if MinProbThresh && mean(P_path)<MinProb
                Aborted_path = 1;
                break % No need to pursue this path it's too improbable
%             elseif ~MinProbThresh && mean(P_path)==0
%                 Aborted_path = 1;
%                 break % No need to pursue this path it's improbable
            else
                Aborted_path = 0;
            end
        end
        if ~Aborted_path
            % cumulate entropy
            Mean_P_path = mean(P_path);
            HY = HY - Mean_P_path*log2(Mean_P_path + (Mean_P_path==0));
            % update path
            [Path, NotLastPath] = updatePath(Path,YIndMax);
        else
            % No need to cumulate entropy
            % update path at a given window
            [Path, NotLastPath] = updatePath(Path,YIndMax,ww);
        end
    end
    fprintf(1, 'Exact entropy (entropy of the neural response)\nusing the exact calculation and paths strategy (calculation can be very long) : %f\n', HY);
end


%% Third Strategy to calculate the exact entropy: save the matrices of probabilities on hard drive 
% Construct the dataset of p products (all the different possible
% combinations of possible responses for each stim) and save this to a file
if strcmp(CalculMode, 'Exact_HardDrive')
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the exact calculation saving temporary data on hard drive (calculation can take large space)\n');
    if ExactMem_history==length(P_YgivenS_local)
        fprintf(1,'The history of the calculation is set to its max, all info redundancies are investigated\n');
    else
        fprintf(1,'The history of the calculation is set to %d time points or %d ms\nOnly the information redundancies present in that window are taken into account\n',ExactMem_history, ExactMem_history*10);
    end
     % Construct the dataset of p products (all the different possible combinations of possible responses for each stim)
    P_YgivenS_localchunked = P_YgivenS_local((end-EarliestPastPoint) : end);
    
    system(sprintf('mkdir %s',FolderTempInfStorage))
    [MatDim, FileMat] = insideMultiMat_F(P_YgivenS_localchunked,Firstwin,StimIndices_AllWin, Stim_local,modnb,FolderTempInfStorage,MinProbThresh, MinProb);%The dimension of this should be size(P_YgivenS{1},1),length(Stim_local) 

    % Calculating exact entropy
    fid_final=fopen(fullfile(FileMat.path,FileMat.name));
    P_Y_wins = nan(1,MatDim(2));

    % read the file matrix in chuncks of 10^6 columns or smaller if smaller
    N = floor(MatDim(2)/10^6);
    NbCol_chunks = [repmat(10^6, 1,N) MatDim(2) - N*10^6];
    Offsets = [0 cumsum(NbCol_chunks(1:end-1))].*MatDim(1).*8; % each double number is coded by 8 bytes

    for cc=1:length(NbCol_chunks)
        % correctly position into the old file and read out chunk of
        % matrix
        fseek(fid_final,Offsets(cc),'bof');
        Mat_local = fread(fid_final,[MatDim(1) NbCol_chunks(cc)], 'double');
        if cc==length(NbCol_chunks)
            % check that we read all the file
            Current_p = ftell(fid_final);
            fseek(fid_final,-2*8,'eof');
            Endoffile_p = ftell(fid_final);
            if Current_p ~=Endoffile_p
                fprintf('Issue here!! the columns if the matrix %s where not all used for the calculation\n',FileMat.name);
            end
        end
        if length(Stim_local)~=MatDim(1)
            fprintf('!!! Problem of dimension here some stims (rows in the file) have been lost??!!\n');
        end
        %fprintf('Info_cumulative_model_Calculus: column chunck %d\nIn P_Y_wins from %d to %d with NbCol Mat_local=%d\n', cc,(Offsets(cc)/(MatDim(1).*8)+1),sum(NbCol_chunks(1:cc)),size(Mat_local,2));
        P_Y_wins((Offsets(cc)/(MatDim(1).*8)+1):sum(NbCol_chunks(1:cc)))=mean(Mat_local,1);
    end
    % rescale p-values so they sum to 1
    P_Y_wins_rescaled = P_Y_wins./sum(P_Y_wins);
    % entropy
    HY = sum(-P_Y_wins_rescaled.*log2(P_Y_wins_rescaled + (P_Y_wins_rescaled==0)));

    % Close file and delete temp directory with its content if we are at the
    % last window
    fclose(fid_final);
    if win>1 
        system(sprintf('rm %s/Temp_CumInf%d_%d_%d', FolderTempInfStorage,modnb,Firstwin,win-1));
    end
    fprintf(1, 'Exact entropy (entropy of the neural response)\nusing the exact calculation saving temporary data on hard drive (calculation can take large space): %f\n', HY);
end


%% Fourth Strategy: Calculate the entropy of the response with a Markov chain approach
if strcmp(CalculMode, 'MarkovChain')
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the Markov Chain estimation with a history of %d and a pace of %d\n',Markov_history(1), Markov_history(2));
    %design a T_chain dimensions matrix  
    HY = Markov_entropy(P_YgivenS_local,Markov_history,MinProbThresh, MinProb,HY_old);
    fprintf(1, 'Exact entropy (entropy of the neural response)\nusing the Markov Chain estimation: %f\n', HY);
end

%% Fifth Strategy: Monte Carlo Estimate of the entropy of the response
if strcmp(CalculMode, 'MonteCarlo')
    fprintf(1, 'Calculation of Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples\n',N_MC);
    % Calculate the cdf of the marginal probabilities (probabilities of responses at each time point)
    P_Yt=cell(win,1);
    for ww=1:win
        P_Yt{ww} = cumsum(mean(P_YgivenS_local{ww},2)./sum(mean(P_YgivenS_local{ww},2)));
    end
    % Set the number of MonteCarlo samples N_MC
    PY_MC = nan(N_MC,1); %exact joint probability of MC sequences of responses
    QY_MC = nan(N_MC,1); %estimation of the joint probability of MC sequences of responses using the marginals
    PYgivenS_MC = nan(N_MC,length(Stim_local)); %exact joint probability of MC sequences of responses for each stimulus
    
    % Loop through the number of samples
    for ss=1:N_MC
        % Probabilistically determine a sequence of responses using the marginal
        % probabilities
        % for each time point, take a random number from the uniform
        % distribution 0-1 and find the neural response that correspond to
        % that probability
        P_YgivenS_local_Resp = nan(win,length(Stim_local));
        for ww=1:win
            u=rand(1);
            Resp_MC = sum(P_Yt{ww}<u)+1;
            P_YgivenS_local_Resp(ww,:) = P_YgivenS_local{ww}(Resp_MC,:);
        end
        
       % Calculate the exact joint probability of that sequence of responses
        % and store it
        PY_MC(ss) = mean(prod(P_YgivenS_local_Resp,1));
        QY_MC(ss) = prod(mean(P_YgivenS_local_Resp,2));
        PYgivenS_MC(ss,:) = prod(P_YgivenS_local_Resp,1);
    end
    % Calculate the MC estimate of the conditional response entropy to the
    % stimulus
    %HYgivenS = - sum(sum(PYgivenS_MC./repmat(QY_MC,1,length(Stim_local)).*log2(PYgivenS_MC),1))/(N_MC*length(Stim_local));
    HYgivenS_Weight = PYgivenS_MC./repmat(QY_MC,1,length(Stim_local));
    HYgivenS_perMCSperStim = -HYgivenS_Weight.*log2(PYgivenS_MC);% This is the MC estimate of the entropy for each sample
    HYgivenS_perStim = sum(HYgivenS_perMCSperStim,1)/N_MC;
    HYgivenS = sum(HYgivenS_perStim)/length(Stim_local);
    % Calculate the MC estimate of the response entropy
    % Weights
    %HY = sum(-PY_MC./QY_MC.*log2(PY_MC))/N_MC;
    HY_Weight = PY_MC./QY_MC;% This is the MC estimate of the entropy for each sample
    HY = sum(-HY_Weight.*log2(PY_MC))/N_MC;
    % Estimate the Error on MC estimates (based on the asymptotics method
    % of Koehler et al. 2009 Am. Stat.
    %HYgivenS_MCE = ((sum((sum(HYgivenS_perMCSperStim,2)./length(Stim_local) - HYgivenS).^2))^0.5)/N_MC;
    HYgivenS_MCE = sum(sum(HYgivenS_Weight.*((-log2(PYgivenS_MC) - repmat(HYgivenS_perStim,N_MC,1)).^2).^0.5))/(N_MC*length(Stim_local));
    HY_MCE = sum(HY_Weight.*((-log2(PY_MC) - HY).^2).^0.5)/N_MC;
    %HY_MCE = sum(((-HY_Weight.*log2(PY_MC) - HY).^2).^0.5)/N_MC;
    Icum_MCE =HYgivenS_MCE +  HY_MCE;
    if DebugFig
        figure()
        plot(1:100, QY_MC(1:100), 1:100,PY_MC(1:100))
        legend('Product of probas','exact joint probas')
        xlabel('Monte Carlo samples')
        ylabel('probability')
    end
    %HY = sum(-PY_MC./QY_MC.*log2(PY_MC))./sum(PY_MC);
    %HY = sum(-PY_MC./QY_MC.*log2(PY_MC));
    fprintf(1, 'Exact entropy (entropy of the neural response)\nusing the Monte Carlo estimation with %d samples: %f\n',N_MC, HY);
    Icum=nan(2,1);
    Icum(2) = Icum_MCE;
    HYgivenS(2) = HYgivenS_MCE;
    HY(2) = HY_MCE;
end

%% Information
Icum(1) = (HY(1) - HYgivenS(1));
fprintf(1,'Cumulative information = %f bits\n',Icum);
ElapsedTime = toc;
fprintf(1,'info_cumulative_model_Calculus run for %d seconds\n',ElapsedTime);
%fprintf('info_cumulative_model_Calculus run for %d seconds or %d seconds per path for %d paths\n',ElapsedTime, ElapsedTime/NbP, NbP);
end