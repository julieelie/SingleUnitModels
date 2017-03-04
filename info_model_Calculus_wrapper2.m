function [Info]=info_model_Calculus_wrapper2(PSTH, FirstTimePoint, LastTimePoint, CatList, Response_samprate, Bootstrap_switch)

if nargin<6
    Bootstrap_switch =0;
end

NbStim = length(PSTH);
% Change the grouping parameter to a vector if not a vector
if nargin<4
    fprintf(1,'No input for categories\nNo information about category is calculated\n');
    CatList_in = [];
    NbCat = 1;
elseif iscell(CatList) && ischar(CatList{1})
    CatList_in = nan(1,NbStim);
    % Number of categories of stims
    IdCats = unique(CatList);
    NbCat = length(IdCats);
    for ss =1:NbCat
        CatList_in(find(strcmp(CatList, IdCats(ss))))=ss;
    end
elseif iscell(CatList) && isnumeric(CatList{1})
    CatList_in = cell2mat(CatList);
     % Number of categories of stims
    IdCats = unique(CatList_in);
    NbCat = length(IdCats);
elseif isnumeric(CatList)
    CatList_in = CatList;
     % Number of categories of stims
    IdCats = unique(CatList_in);
    NbCat = length(IdCats);
end

% Set the maximum value of Y investigated for the calculation of
% information
MaxY = 2*(LastTimePoint - FirstTimePoint +1)*1000/Response_samprate; % response sampling  rate should be in hertz


% Format the input
Info.InputdataStim = nan(1,NbStim);
if ~Bootstrap_switch
    if size(PSTH,1)==NbStim
        PSTH_Local = cell2mat(PSTH);
    else
        PSTH_Local = cell2mat(PSTH');
    end
    Info.InputdataStim = sum(PSTH_Local(:,FirstTimePoint:LastTimePoint),2);
else
    for ss = 1:NbStim
        PSTH_Local = cell2mat(PSTH{ss});
        NJK = size(PSTH_Local,1);
        Info.InputdataStim(ss) = sum(PSTH_Local(randperm(NJK,1),FirstTimePoint:LastTimePoint),2);
    end
end

% Entropy of the stimulus dataset
Info.stim_entropy = log2(NbStim);
    
% Calculate information about stimuli   
[Info.stim_value,Info.P_YgivenS,~] = info_model_Calculus(Info.InputdataStim, MaxY);

if ~isempty(CatList_in)
    % Derive P_YgivenC and calculate information about categories
    Info.P_YgivenC = nan(MaxY+1,NbCat);
    P_YgivenS_scaled = Info.P_YgivenS./sum(Info.P_YgivenS,1); % Make sure that all distributions of probabilities sum to 1 for each stimulus
    for Cat=1:NbCat
        Info.P_YgivenC(:,Cat) = mean(P_YgivenS_scaled(:,find(CatList_in == Cat)),2);
    end
    % Conditional entropy
    H_ycat = sum(sum(-Info.P_YgivenC.*log2(Info.P_YgivenC+(Info.P_YgivenC==0))))/NbCat;
    % entropy of the response
    P_y_i = mean(Info.P_YgivenC,2);
    % Calculate log probability
    H_y = sum(-P_y_i.*log2(P_y_i + (P_y_i==0)));% entropy of the response

    %% Calculate the model information
    Info.cat_value = (H_y - H_ycat);
    % Entropy of the categories in the dataset
    Info.cat_entropy = log2(NbCat);
end
end