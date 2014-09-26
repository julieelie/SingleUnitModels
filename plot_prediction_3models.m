function [ output_args ] = plot_prediction_3models( ModelPath, MAXWIN,WIN )
%This function plot the figures of the observed spike rate vs predictions
%by the three models at a given time point win
%   Detailed explanation goes here
load(ModelPath)
Maxwin = find(Wins==MAXWIN);

nvoc = size(voc{1},1);

%% Figure of the R2A profils
figure(1)
subplot(4,1,1)
[path, name]=fileparts(ModelPath);
plot(Wins(1:Maxwin), R2A.Acoustic(1:Maxwin), 'ko-', Wins(1:Maxwin), R2A.Semantic(1:Maxwin), 'bo-', Wins(1:Maxwin),R2A.AcSem(1:Maxwin), 'ro-')
legend('Acoustic', 'Semantic', 'Acoustic + Semantic', 'Location', 'SouthWest')
hold on
plot(Wins(find(PValLRatio.AcAcSem(1:Maxwin)<0.05)), repmat(0.9,sum(PValLRatio.AcAcSem(1:Maxwin)<0.05)), '*')
hold off
axis([0 MAXWIN+50 0 1])
title(sprintf('%s\nR2A as a function of the window size.\n', name));


%figure()
subplot(4,1,2)
plot(Wins(1:Maxwin), Best_nbPC(1:Maxwin))
title('Number of PC in the Acoustic model')
axis([0 MAXWIN+50 0 nvoc])

NCat = nan(length(voc),1);
NStim = nan(length(voc),1);
for vv=1:length(voc)
    NCat(vv)=length(unique(voc{vv}));
    NStim(vv) = length(voc{vv});
end
subplot(4,1,3)
plot(Wins(1:Maxwin),NCat(1:Maxwin))
title('Number of remaining categories in the dataset');
axis([0 MAXWIN+50 0 10])

subplot(4,1,4)
plot(Wins(1:Maxwin),NStim(1:Maxwin))
title('Number of remaining vocalizations in the dataset');
xlim([0 MAXWIN+50])

if nargin<3
    All=input('Do you want to see all time points? (Y vs N)\n', 's');
    if strcmp(All, 'N')
        WIN=input('For which time point do you want the 3 models predictions and the STRF?\n');
        ww=find(Wins==WIN);
    else
        ww=1:length(Wins);
    end
else
    ww=find(Wins==WIN);
end

for win=ww
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
    figure(2)
    for jj=1:3
        %figure()
        subplot(1,3,jj);
        h=gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'mgcbrkyyr', '......d.d',[20 20 20 20 20 20 10 20 10]);
        %gscatter(LocalNeuralResponses, LocalModelPredict{jj}, LocalVoc, 'rkmgcbryy', '......dd.',[20 20 20 20 20 20 10 10 20]);
        title(sprintf('%s Adjusted R^2 = %.2f at %d ms', Legend{jj}, R2A_local(jj),Wins(ww)));
        MAX=max(MAXX,MAXY);
        axis([0 MAX+MAX/4 0 MAX]);
        hold on
        plot(0:MAX/10:MAX,0:MAX/10:MAX, 'k');
        hold off
    end
    %% Plot the STRF and model of semantic model at that particular time point
    figure(3)
    subplot(2,1,1)
    imagesc(ModelFilters.STRF_to{win}, ModelFilters.STRF_fo{win}, ModelFilters.STRF_time{win})
    axis xy
    title(sprintf('STRF at %d ms',Wins(ww)));
    subplot(2,1,2)
    plot(ModelFilters.Sem{win}, 1:length(ModelFilters.Sem{win}));
    title('Predicted SR with semantic Model');
    set(gca,'YTickLabel', unique(voc{win}));
    pause
end
end

