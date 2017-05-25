function slurm_script_AN_CumInfo_patchcorrection(Cell)
addpath(genpath('/global/home/users/jelie/CODE/SingleUnitModels'));
addpath(genpath('/global/home/users/jelie/CODE/GeneralCode'));
addpath(genpath('/global/home/users/jelie/CODE/tlab/src'));
rmpath(genpath('/global/home/users/jelie/CODE/tlab/src/hedi'));
Storage_path = '/global/scratch/jelie/MatFiles/';

fprintf('------------------------------------------------------\n')
fprintf('---------- %s cell Dataset-----------\n', Cell)
fprintf('------------------------------------------------------\n')
fprintf('**************** Calculating cumulative information *****************\n')

%% Calculate theoretical values if the spike rate was exactly known
fprintf('**** Theoretical values: Monte Carlo 10^6, Exact calculation and Markov with 5 bins memory *****\n');
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'P_YgivenS_Theo', 'Nb_Win','Cum_boot');
Nb_Win = length(P_YgivenS_Theo);
Icum_EstMonteCarlo6_Theo = nan(1,Nb_Win);
Icum_ExactMem0_5_Theo = nan(1,Nb_Win);
Icum_EstMarkov5_Theo = nan(1,Nb_Win);
HY_Markov5 = nan(1,Nb_Win);
for tt=2:Nb_Win
    tstart = tic;
    fprintf('Time point %d/%d\n', tt, Nb_Win);
    P_YgivenS_local = P_YgivenS_Theo(1:tt);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo6_local,~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MonteCarlo', 'MCParameter',10^6);
    Icum_EstMonteCarlo6_Theo(tt) = Icum_EstMonteCarlo6_local(1);
    
    % Exact calculation with 5 bins memory
    [Icum_ExactMem0_5_Theo(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','Exact_Mem', 'Exact_history',5);
    if tt==2
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5_Theo(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
    else
        % Markov chain estimation 5 bins memory
        [Icum_EstMarkov5_Theo(tt),HY_Markov5(tt),~]=info_cumulative_model_Calculus(P_YgivenS_local,'CalMode','MarkovChain', 'MarkovParameters',[5,1], 'HY_old', HY_Markov5(tt-1));
    end
    telapsed = toc(tstart);
    fprintf('Markov + Exact calculation: total elapsed time: %d s\n', telapsed)
end
load(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Info_bcorr');
Icum_EstMonteCarlo6_Theo(1) = Info_bcorr(1); % Initializing the first value of cumulative info
Icum_ExactMem0_5_Theo(1) =  Info_bcorr(1);
Icum_EstMarkov5_Theo(1) = Info_bcorr(1);
save(sprintf('%sInfoCumInfoSpikeCount_AN_JK_KDE_%s.mat',Storage_path,Cell),'Icum_EstMonteCarlo6_Theo','Icum_ExactMem0_5_Theo','Icum_EstMarkov5_Theo', '-append');
clear P_Y* Icum* HY*
end