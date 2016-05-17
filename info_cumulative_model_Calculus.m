function [Icum]=info_cumulative_model_Calculus(P_YgivenS, win,StimIndices_AllWin, Stim_local)
tic
MinProb = 1/(8*365*24*60*60*5); %1/number of 20ms bins in a life of a zebra finch in the lab
%% Construct the data set of individual p
P_YgivenS_local = cell(win,1);
P_YgivenS_local_rescaled1 = nan(size(P_YgivenS{1},1),length(Stim_local),win);
YIndMax = nan(win,1);
for ww=1:(win-1)
    P_YgivenS_local_w = nan(size(P_YgivenS{1},1),length(Stim_local));
    for ss=1:length(Stim_local)
        st = Stim_local(ss);
        st_local = find(StimIndices_AllWin{ww}==st);
        P_YgivenS_local_w(:,ss) = P_YgivenS{ww}(:,st_local);
        P_YgivenS_local_rescaled1(:,ss,ww)= P_YgivenS{ww}(:,st_local) ./sum(P_YgivenS{ww}(:,st_local));
    end
    Badyy = find(sum((P_YgivenS_local_w < MinProb),2)==length(Stim_local));
    Goodyy=setdiff(1:size(P_YgivenS{1},1), Badyy);
    YIndMax(ww) = length(Goodyy);
    P_YgivenS_local{ww} = P_YgivenS_local_w(Goodyy,:);
        
end
Badyy = find(sum((P_YgivenS{win} < MinProb),2)==length(Stim_local));
Goodyy=setdiff(1:size(P_YgivenS{1},1), Badyy);
YIndMax(win) = length(Goodyy);
P_YgivenS_local{win}=P_YgivenS{win}(Goodyy,:);
P_YgivenS_local_rescaled1(:,:,win)=P_YgivenS{win}(:,:)./repmat(sum(P_YgivenS{win}(:,:),1),size(P_YgivenS{win},1),1);

%% Calculating conditional entropy
HYgivenS = sum(sum(sum(-P_YgivenS_local_rescaled1.*log2(P_YgivenS_local_rescaled1))))/length(Stim_local);

% %% First strategy for coding the exact entropy
% % Construct the dataset of p products (all the different possible combinations of possible responses for each stim)
% P_YgivenS_perstim_wins = insideMultiMat(P_YgivenS_local);%The dimension of this should be size(P_YgivenS{1},1),length(Stim_local)
% 
% % Calculating exact entropy
% P_Y_wins = sum(P_YgivenS_perstim_wins,2)/length(Stim_local);%The dimension of this should be size(P_YgivenS{1},1)^win,1
% % rescale p-values so they sum to 1
% P_Y_wins_rescaled = P_Y_wins./sum(P_Y_wins);
% % entropy
% HY = sum(-P_Y_wins_rescaled.*log2(P_Y_wins_rescaled));

%% Second strategy for coding exact entropy
% Summing up the probabilities on the fly
% Define a minimum probability under which the proba should be considered 0

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
        if mean(P_path)<MinProb
            Aborted_path = 1;
            break % No need to pursue this path it's too improbable
        else
            Aborted_path = 0;
        end
    end
    if ~Aborted_path
        % cumulate entropy
        Mean_P_path = mean(P_path);
        HY = HY - Mean_P_path*log2(Mean_P_path);
        % update path
        [Path, NotLastPath] = updatePath(Path,YIndMax);
    else
        % No need to cumulate entropy
        % update path at a given window
        [Path, NotLastPath] = updatePath(Path,YIndMax,ww);
    end
end

    
    
%% Information
Icum = (HY - HYgivenS);
ElapsedTime = toc;
fprintf('info_cumulative_model_Calculus run for %d seconds or %d seconds per path for %d paths\n',ElapsedTime, ElapsedTime/NbP, NbP);
end
