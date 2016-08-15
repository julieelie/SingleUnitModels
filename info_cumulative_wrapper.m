function [Model] = info_cumulative_wrapper(ParamModel,SWITCH,Model,mm,x_stim_indices_wholeset, Stim_local)
fprintf('Pre-process data for cumulative information calculations\n')
X_stim_indices_wholeset = x_stim_indices_wholeset(1:mm);
mm_local=0;
%Pre-process data for a parfor loop
P_YgivenS_allModel = cell(10,1);
if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
    % ACoustic Model
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.Acoustic.P_YgivenS_all1(1:mm,1);
end
if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
    % Semantic Model
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.Semantic.P_YgivenS_all1(1:mm,1);
end
if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
    % AcSemAc
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AcSemAc.P_YgivenS_all1(1:mm,1);
end
if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
    % AcSemSem
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AcSemSem.P_YgivenS_all1(1:mm,1);
end
% Floor
mm_local=mm_local+1;
P_YgivenS_allModel{mm_local}=Model.Floor.P_YgivenS_all1(1:mm,1);
%Ceiling
mm_local=mm_local+1;
P_YgivenS_allModel{mm_local}=Model.Ceiling.P_YgivenS_all1(1:mm,1);

% AR
if SWITCH.AR
    mm_local=mm_local+1;
    P_YgivenS_allModel{mm_local}=Model.AR.P_YgivenS_all1(1:mm,1);
end

Icum_ExactMem0_5 = nan(mm_local,1);
Icum_EstMonteCarlo10_7 = nan(3,mm_local);
Icum_EstMarkov2 = nan(3,mm_local);
Icum_EstMarkov3 = nan(3,mm_local);
Icum_EstMarkov4 = nan(3,mm_local);
Icum_EstMarkov5 = nan(3,mm_local);

fprintf('Calculate Cumulative information for all models from win %d\n',1)
NS = 10000000;
parfor modelrun=1:mm_local
    fprintf('Cumulative info %d/%d\n', modelrun, mm_local);
    % Monte Carlo estimation with full memory
    [Icum_EstMonteCarlo10_7(:,modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','MonteCarlo', 'MCParameter',NS);
    % Exact calculation with 50 ms memory
    [Icum_ExactMem0_5(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','Exact_Mem', 'Exact_history',5);
    % Markov chain estimation 20 ms
    [Icum_EstMarkov2(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[2,1]);
    % Markov chain estimation 30 ms
    [Icum_EstMarkov3(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[3,1]);
    % Markov chain estimation 40 ms
    [Icum_EstMarkov4(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[4,1]);
    % Markov chain estimation 50 ms
    [Icum_EstMarkov5(modelrun),~]=info_cumulative_model_Calculus(P_YgivenS_allModel{modelrun},'StimIndicesAll',X_stim_indices_wholeset,'StimIndicesLast',Stim_local,'Model#',modelrun,'CalMode','MarkovChain', 'MarkovParameters',[5,1]);
end

mm_local=0;
%Post-process data for a parfor loop
if ParamModel.ModelChoice(1) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on Acoustic**\n')
    % ACoustic Model
    mm_local=mm_local+1;
    Model.Acoustic.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.Acoustic.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
    Model.Acoustic.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.Acoustic.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.Acoustic.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.Acoustic.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(2) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on Semantic**\n')
    % Semantic Model
    mm_local=mm_local+1;
    Model.Semantic.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.Semantic.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
    Model.Semantic.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.Semantic.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.Semantic.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.Semantic.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(4) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on AcSemAc**\n')
    % AcSemAc
    mm_local=mm_local+1;
    Model.AcSemAc.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AcSemAc.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
    Model.AcSemAc.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AcSemAc.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AcSemAc.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AcSemAc.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
if ParamModel.ModelChoice(5) && ~SWITCH.AllAlpha
    fprintf('**CumInfo on AcSemSem**\n')
    % AcSemSem
    mm_local=mm_local+1;
    Model.AcSemSem.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AcSemSem.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
    Model.AcSemSem.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AcSemSem.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AcSemSem.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AcSemSem.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end
% Floor
fprintf('**CumInfo on Floor**\n')
mm_local=mm_local+1;
Model.Floor.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
Model.Floor.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
Model.Floor.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
Model.Floor.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
Model.Floor.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
Model.Floor.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
%Ceiling
fprintf('**CumInfo on Ceiling**\n')
mm_local=mm_local+1;
Model.Ceiling.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
Model.Ceiling.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
Model.Ceiling.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
Model.Ceiling.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
Model.Ceiling.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
Model.Ceiling.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
% AR
if SWITCH.AR
    fprintf('**CumInfo on AR**\n')
    mm_local=mm_local+1;
    Model.AR.cum_info_ExactHist5(mm,1)=Icum_ExactMem0_5(mm_local);
    Model.AR.cum_info_EstMonteCarlo10_7(mm,1:3)=Icum_EstMonteCarlo10_7(mm_local,:);
    Model.AR.cum_info_EstMarkov2(mm,1)=Icum_EstMarkov2(mm_local);
    Model.AR.cum_info_EstMarkov3(mm,1)=Icum_EstMarkov3(mm_local);
    Model.AR.cum_info_EstMarkov4(mm,1)=Icum_EstMarkov4(mm_local);
    Model.AR.cum_info_EstMarkov5(mm,1)=Icum_EstMarkov5(mm_local);
end