function [calfilename, GlobalR2A] = SpectroSemantic_Neuro_model(MatfilePath, MinWin, MaxWin, Increment, ResDelay)
FIG=1; %set to 1 for debugging figs
if nargin<5
    ResDelay = 10; %predict the neural response with a 10ms delay after the end of the stimulus
end
if nargin<4
    Increment = 5; %increase the size of the spectro window with a 5ms pace
end
if nargin<3
    MaxWin = 600; %maximum values the window of analysis can reach
end
if nargin<2
    MinWin = 40; %minimum size of the window of analysis from the begining
end

%% Load the unit matfile
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/GreBlu9508M/ZS_Site2_L1100R1450_21.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L1100R1450_e21_s0_ss1.mat')
%Res=load('/Users/elie/Documents/MATLAB/data/matfile/WholeVocMat/WholeVoc_Site2_L2000R1600_e27_s1_ss1.mat')
Res = load(MatfilePath);

%% Get the data ready
% Select first sections
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

%% Extract Emitter ID
Ename = cell(length(DataSel),1);
Esex = cell(length(DataSel),1);
Eage = cell(length(DataSel),1);
Erelated = cell(length(DataSel),1);
for ii=1:length(DataSel)
    dd=DataSel(ii);
    [Path,File,Ext] = fileparts(Res.Original_wavfiles{dd});
    Ename{ii} = File(1:10);
    Esex{ii} = File(12);
    Eage{ii} = File(13);
    Erelated{ii} = File(14);
end
Emitter.Ename = Ename;
Emitter.Esex = Esex;
Emitter.Eage = Eage;
Emitter.Erelated = Erelated;
%% Inject the data in the models
Spectro.spec = Res.Spectro(DataSel);
Spectro.to = Res.Spectroto(DataSel);
Spectro.fo = Res.Spectrofo(DataSel);
%[R2A, SSres, SSexp, SStot, ModelPredict, LL, NEC, PValLRatio, HLRatio, NeuroRes, voc, Best_nbPC, Pvalue,Wins, NeuralResponse, STRF_time, STRF_to, STRF_fo, ModSem] = GrowingModels(Spectro, Res.VocType(DataSel), Res.PSTH(DataSel), MinWin, MaxWin, Increment, ResDelay);
[R2A, SSres, SSexp, SStot, ModelPredict, LL, NEC, PValLRatio, HLRatio, NeuroRes, voc, Best_nbPC, Pvalue,Wins, NeuralResponse, STRF_time, STRF_to, STRF_fo, ModSem] = GrowingModelsRidge(Spectro, Res.VocType(DataSel), Res.PSTH(DataSel), Emitter, MinWin, MaxWin, Increment, ResDelay);

%This calculation of Global R2 needs to be revised and done as calculated in FindSemanticNeuronsFinal_encodingM 
GlobalR2A.Acoustic = 1 - nansum(SSres.Acoustic)/nansum(SStot);
GlobalR2A.Semantic = 1 - nansum(SSres.Semantic)/nansum(SStot);
GlobalR2A.AcSem = 1 - nansum(SSres.AcSem)/nansum(SStot);

if FIG==1 %these figures give the plots of predicted vs actual spike rates for each model
    figure()
    [path, name]=fileparts(MatfilePath);
    plot(Wins, R2A.Acoustic, 'ko-', Wins, R2A.Semantic, 'bo-', Wins,R2A.AcSem, 'ro-')
    legend('Acoustic', 'Semantic', 'Acoustic + Semantic', 'Location', 'SouthWest')
    hold on
    plot(Wins(find(PValLRatio.AcAcSem<0.05)), repmat(0.9,sum(PValLRatio.AcAcSem<0.05)), '*')
    hold off
    title(sprintf('%s\nR2A as a function of the window size.\n', name));
    
    figure()
    plot(Wins, Best_nbPC)
    title('Number of PC in the Acoustic model')
    axis([0 max(Wins) 0 nvoc])
    
    NCat = nan(length(voc),1);
    for vv=1:length(voc)
        NCat(vv)=length(unique(voc{vv}));
    end
    figure()
    plot(Wins,NCat)
    title('Number of remaining categories in the dataset');
    axis([0 max(Wins) 0 10])
    pause
    
    
    win = find(Wins==80);
    Legend= cell(3,1);
    Legend{1}='Acoustic only';
    Legend{2}='Semantic only';
    Legend{3}='Acoustic + Semantic';
    R2A_local = [R2A.Acoustic(win) R2A.Semantic(win) R2A.AcSem(win)];
    VocType = voc{win};
    UV = unique(VocType);
    VT = nan(length(UV),1);
    for vv = 1:length(UV)
        VT(vv) = sum(strcmp(VocType, UV(vv)));
    end
    VTMinNb=min(VT);
    MAXY = max([max(cell2mat(ModelPredict.Acoustic(win))) max(cell2mat(ModelPredict.Semantic(win))) max(cell2mat(ModelPredict.AcSem(win)))]);
    MAXX = max(NeuralResponse{win});
    Win_ModelPredict = cell(3,1);
    Win_ModelPredict(1) = ModelPredict.Acoustic(win);
    Win_ModelPredict(2) = ModelPredict.Semantic(win);
    Win_ModelPredict(3) = ModelPredict.AcSem(win);
    IndicesSelect=zeros(9,VTMinNb);
    LocalNeuralResponses =  zeros(length(UV)*VTMinNb,1);
    LocalModelPredict = cell(3,1);
    for mm = 1:3
        LocalModelPredict{mm} = zeros(length(UV)*VTMinNb,1);
    end
    LocalVoc = cell(length(UV)*VTMinNb,1);
        [C,IA,IC] = unique(VocType);
        for ii = 1:length(UV)
            LocalDataSel=find(IC==ii);
            ID = randperm(length(LocalDataSel));
            Vocrand=LocalDataSel(ID);
            IndicesSelect(ii,1:VTMinNb) = Vocrand(1:VTMinNb);
            LocalNeuralResponses(((ii-1)*VTMinNb+1):ii*VTMinNb)=NeuralResponse{win}(IndicesSelect(ii,:)); 
            for mm = 1:3
                LocalModelPredict{mm}(((ii-1)*VTMinNb+1):ii*VTMinNb) =  Win_ModelPredict{mm}(IndicesSelect(ii,:));
            end
            LocalVoc(((ii-1)*VTMinNb+1):ii*VTMinNb) = VocType(IndicesSelect(ii,:));
        end
       
        for jj=1:3
            figure()
            %subplot(1,3,jj);
            h=gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
            %gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'rkmgcbryy', '......dd.',[20 20 20 20 20 20 10 10 20]);
            title(sprintf('%s Adjusted R^2 = %.2f', Legend{jj}, R2A_local(jj)));
            MAX=max(MAXX,MAXY);
            axis([0 MAX+MAX/4 0 MAX]);
            hold on
            plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
            hold off
        end
end

%% Boostrap R2A difference global values Not doing boostrap anymore as of 06(june).08.2014
% Nmod = length(Wins);
% Bootstrap =10000;
% R2A_Shuff.Acoustic = nan(Bootstrap,1);
% R2A_Shuff.AcSem = R2A_Shuff.Acoustic;
% R2A_Shuff.Semamtic = R2A_Shuff.Acoustic;
% R2A_Shuff.SemAc = R2A_Shuff.Acoustic;
% 
% for bb=1:Bootstrap
%     Indices = [repmat(0, floor(Nmod/2),1); repmat(1,Nmod-floor(Nmod/2),1)];
%     ShuffIndices = Indices(randperm(Nmod));
%     SSMat = [SSres.Acoustic SSres.AcSem];
%     SSMat2 = [SSres.Semantic SSres.AcSem];
%     SSAcoustic_Shuff = 0;
%     SSAcSem_Shuff = 0;
%     SSSemantic_Shuff = 0;
%     SSSemAc_Shuff = 0;
%     % find where the time prediction of PSTH stopped due to not enough
%     % stims
%     NonCalWin = find(isnan(SSres.Acoustic));
%     if ~isempty(NonCalWin)
%         LastWin = NonCalWin(1) - 1;
%     else
%         LastWin = length(Wins);
%     end
%     for ii=1:LastWin
% %         SSAcoustic_Shuff = SSAcoustic_Shuff + SSMat(ii,(ShuffIndices(ii)+1));
% %         SSAcSem_Shuff = SSAcSem_Shuff + SSMat(ii,(~ShuffIndices(ii)+1));
% %         SSSemantic_Shuff = SSSemantic_Shuff + SSMat2(ii,(ShuffIndices(ii)+1));
% %         SSSemAc_Shuff = SSSemAc_Shuff + SSMat2(ii,(~ShuffIndices(ii)+1));
%         
%         ind_shuff1 = randi(2,1);
%         ind_shuff2 = ind_shuff1+1;
%         if (ind_shuff2 == 3) 
%             ind_shuff2 = 1;
%         end
%         
%         SSAcoustic_Shuff = SSAcoustic_Shuff + SSMat(ii,ind_shuff1);
%         SSAcSem_Shuff = SSAcSem_Shuff + SSMat(ii,ind_shuff2);
%         SSSemantic_Shuff = SSSemantic_Shuff + SSMat2(ii,ind_shuff1);
%         SSSemAc_Shuff = SSSemAc_Shuff + SSMat2(ii,ind_shuff2);
%     end
%     if ((SSAcoustic_Shuff + SSAcSem_Shuff) < (nansum(SSres.Acoustic) + nansum(SSres.AcSem) - 0.0001)) || ((SSAcoustic_Shuff + SSAcSem_Shuff) > (nansum(SSres.Acoustic) + nansum(SSres.AcSem) + 0.0001))
%         sprintf('WARNING!!! PB in Bootstrap calculations')
%     end
%     if ((SSSemantic_Shuff + SSSemAc_Shuff) < (nansum(SSres.Semantic) + nansum(SSres.AcSem) - 0.0001)) || ((SSSemantic_Shuff + SSSemAc_Shuff) > (nansum(SSres.Semantic) + nansum(SSres.AcSem) + 0.0001))
%         sprintf('WARNING!!! PB in Bootstrap calculations')
%     end
% 
%     R2A_Shuff.Acoustic(bb) = 1 - SSAcoustic_Shuff/nansum(SStot);
%     R2A_Shuff.AcSem(bb) = 1 - SSAcSem_Shuff/nansum(SStot);
%     R2A_Shuff.Semamtic(bb) = 1 - SSSemantic_Shuff/nansum(SStot);
%     R2A_Shuff.SemAc(bb) = 1 - SSSemAc_Shuff/nansum(SStot);
% end
% 
% DiffR2A_AcSem_Shuff = R2A_Shuff.AcSem - R2A_Shuff.Acoustic;
% DiffR2A_SemAc_Shuff = R2A_Shuff.SemAc - R2A_Shuff.Semamtic;
% DiffR2A_AcSem_Obs = GlobalR2A.AcSem - GlobalR2A.Acoustic;
% DiffR2A_SemAc_Obs = GlobalR2A.AcSem - GlobalR2A.Semantic;
% 
% GlobalR2A.PvalAcSem = sum(DiffR2A_AcSem_Shuff > DiffR2A_AcSem_Obs)/Bootstrap;
% GlobalR2A.PvalSemAc = sum(DiffR2A_SemAc_Shuff > DiffR2A_SemAc_Obs)/Bootstrap;
 

fprintf(1, 'Storing values\n');
ModelFilters.STRF_time = STRF_time;
ModelFilters.STRF_to = STRF_to;
ModelFilters.STRF_fo = STRF_fo;
ModelFilters.Sem = ModSem;
Cal.R2A = R2A;
Cal.GlobalR2A = GlobalR2A;
Cal.ModelPredict = ModelPredict;
Cal.LL = LL;
Cal.NEC = NEC;
Cal.PValLRatio = PValLRatio;
Cal.HLRatio = HLRatio;
Cal.NeuroRes = NeuroRes;
Cal.NeuralResponse = NeuralResponse;
Cal.ModelFilters=ModelFilters;
Cal.voc = voc;
Cal.Best_nbPC = Best_nbPC;
Cal.Pvalue = Pvalue;
Cal.Wins = Wins;
Cal.SSres=SSres;
Cal.SSexp = SSexp;
Cal.SStot = SStot;

if ismac()
        [status username] = system('who am i');
        if strcmp(strtok(username), 'frederictheunissen')
            if strncmp('/auto/fdata/solveig',stim_name, 19)
            elseif strncmp('/auto/fdata/julie',stim_name, 17)
                calfilename = fullfile('/Users','frederictheunissen','Documents','Data','Julie','matfile',Res.subject,['Models_' Res.Site '.mat']);
            end
        elseif strcmp(strtok(username), 'elie')
            calfilename = fullfile('/Users','elie','Documents','MATLAB','data','matfile',['Models_' Res.Site '.mat']);
        end
else
    calfilename=fullfile('/auto','k6','julie','matfile',Res.subject,['Models_' Res.Site '.mat']);
end

save(calfilename, '-struct', 'Cal');
end

