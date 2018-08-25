load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_ShuffPredict_PCCAbCh_Invariance.mat', 'List_matfilepath', 'List_anat', 'ZONES', 'ZONES_List', 'SUBJ', 'Spike_shape', 'MeanSR', 'StimTypeCM','SUBJECTS')
load('/Users/elie/Documents/MATLAB/data/matfile/SemanticAnalysis_Encod.mat');
NU=length(List_matfilepath);
plevel = 0.05;

%% Sort cells and select those for which we have at least one significant model
Expcells = find((GlobalR2A.PAc<plevel) + (GlobalR2A.PSem<plevel) + (GlobalR2A.PAcSem<plevel));
figure()
subplot(3,1,1)
hist(GlobalR2A.Acoustic(find(GlobalR2A.PAc<plevel)), 0:0.05:1);
subplot(3,1,2)
hist(GlobalR2A.Semantic(find(GlobalR2A.PSem<plevel)),0:0.05:1);
subplot(3,1,3)
hist(GlobalR2A.AcSem(find(GlobalR2A.PAcSem<plevel)),0:0.05:1);

NEx=length(Expcells)
IncSem = GlobalR2A.AcSem(Expcells) - GlobalR2A.Acoustic(Expcells);
IncAc = GlobalR2A.AcSem(Expcells) - GlobalR2A.Semantic(Expcells);
PAS_Ac = GlobalR2A.PAS_Ac(Expcells);
PAS_Sem = GlobalR2A.PAS_Sem(Expcells);


% Get rid of zero values of pvalues to obtain a finite log10 for the AC+Se
% vs Ac model comparison
ZeroPAS_Ac = find(GlobalR2A.PAS_Ac==0);
MinPAS_Ac = min(GlobalR2A.PAS_Ac(setdiff(1:NU,ZeroPAS_Ac)));
PAS_Ac_zcorrected = GlobalR2A.PAS_Ac;
PAS_Ac_zcorrected(ZeroPAS_Ac)=MinPAS_Ac;
LogPVAS_Ac = -log10(PAS_Ac_zcorrected(Expcells));

% Get rid of zero values of pvalues to obtain a finite log10 for the AC+Se
% vs Se model comparison
ZeroPAS_Sem = find(GlobalR2A.PAS_Sem==0);
MinPAS_Sem = min(GlobalR2A.PAS_Sem(setdiff(1:NU,ZeroPAS_Sem)));
PAS_Sem_zcorrected = GlobalR2A.PAS_Sem;
PAS_Sem_zcorrected(ZeroPAS_Sem)=MinPAS_Sem;
LogPVAS_Sem = -log10(PAS_Sem_zcorrected(Expcells));

% divide the explained cells population in 4 sets
% The synergic cells (both parameter categories bring more
% than just the overlap)
Expcells_SYN = find((PAS_Ac<plevel) .* (PAS_Sem<plevel));
% The semantic cells (semantic parameters bring more than the overlap)
Expcells_S = find((PAS_Ac<plevel) .* (PAS_Sem>=plevel));
% The acoustic cells (acoustic parameters bring more than the overlap)
Expcells_A = find((PAS_Ac>=plevel) .* (PAS_Sem<plevel));
% The acoustic and semantic cells (both parameter categories are equivalent
% and explanation by the model completely overlap
Expcells_as = find((PAS_Ac>=plevel) .* (PAS_Sem>=plevel));


%% Plots differences of Global R2A and significance of full models compare to
% acoustic only for Expxells
GRAD1=cubehelix(ceil(max(LogPVAS_Ac*1000)+1), 0.5, -0.8, 1.5, 1.7, [1,0]);
GRAD2=cubehelix(ceil(max(LogPVAS_Sem*1000)+1), 0.5, -0.8, 1.5, 1.7, [1,0]);
figure()
subplot(1,2,1)
for jj = 1:NEx
    plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor',GRAD1(round(LogPVAS_Ac(jj).*1000)+1,:));
    hold on
end
XL = xlim;
YL = ylim;
L1=line([0 0], YL);
set(L1, 'LineStyle','-.','Color',[1 0 0]);
L2=line(XL, [0 0]);
set(L2, 'LineStyle','-.','Color',[1 0 0]);
xlabel('Acoustic: Increase in R2A between AcSem and Sem')
ylabel('Semantic: Increase in R2A between AcSem and Ac')
title(sprintf('SignifPredicted Cells\nz-axis: Signif AcSem compare to Ac'));
colorbar()
colormap(GRAD1)
hold off
subplot(1,2,2)
for jj = 1:NEx
    plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor',GRAD2(round(LogPVAS_Sem(jj).*1000)+1,:));
    hold on
end
XL = xlim;
YL = ylim;
L1=line([0 0], YL);
set(L1, 'LineStyle','-.','Color',[1 0 0]);
L2=line(XL, [0 0]);
set(L2, 'LineStyle','-.','Color',[1 0 0]);
xlabel('Acoustic: Increase in R2A between AcSem and Sem')
ylabel('Semantic: Increase in R2A between AcSem and Ac')
title(sprintf('SignifPredicted Cells\nz-axis: Signif AcSem compare to Sem'));
colorbar()
colormap(GRAD2)
hold off

figure()
nn=0;
for jj = 1:NEx
    if PAS_Ac(jj)<plevel && PAS_Sem(jj)<plevel
        plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor','b', 'MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)<plevel && PAS_Sem(jj)>=plevel
        plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor',[0 0.9 0.9],'MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)>=plevel && PAS_Sem(jj)<plevel
        plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor',[1 0.75 0],'MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)>=plevel && PAS_Sem(jj)>=plevel
        plot(IncAc(jj), IncSem(jj), 'ko','MarkerFaceColor',[0.5 0 0.1],'MarkerSize',10);
        nn=nn+1;
    end
    hold on
end
XL = xlim;
YL = ylim;
L1=line([0 0], YL);
set(L1, 'LineStyle','-.','Color',[1 0 0]);
L2=line(XL, [0 0]);
set(L2, 'LineStyle','-.','Color',[1 0 0]);
xlabel('Acoustic: Increase in R2A between AcSem and Sem')
ylabel('Semantic: Increase in R2A between AcSem and Ac')
title(sprintf('SignifPredicted Cells\nz-axis: Signif AcSem compare to Sem and Ac'));
hold off

figure()
nn=0;
for jj = 1:NEx
    ii=Expcells(jj);
    if PAS_Ac(jj)<plevel && PAS_Sem(jj)<plevel
        plot3(GlobalR2A.Acoustic(ii), GlobalR2A.Semantic(ii), GlobalR2A.AcSem(ii), 'ko','MarkerFaceColor','b','MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)<plevel && PAS_Sem(jj)>=plevel
        plot3(GlobalR2A.Acoustic(ii), GlobalR2A.Semantic(ii), GlobalR2A.AcSem(ii), 'ko','MarkerFaceColor',[0 0.9 0.9],'MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)>=plevel && PAS_Sem(jj)<plevel
        plot3(GlobalR2A.Acoustic(ii), GlobalR2A.Semantic(ii), GlobalR2A.AcSem(ii), 'ko','MarkerFaceColor',[1 0.75 0],'MarkerSize',10);
        nn=nn+1;
    elseif PAS_Ac(jj)>=plevel && PAS_Sem(jj)>=plevel
        plot3(GlobalR2A.Acoustic(ii), GlobalR2A.Semantic(ii), GlobalR2A.AcSem(ii), 'ko','MarkerFaceColor',[0.5 0 0.1],'MarkerSize',10);
        nn=nn+1;
    end
    hold on
end
set(gca, 'XGrid', 'ON')
set(gca, 'YGrid', 'ON')
set(gca, 'ZGrid', 'ON')
xlabel('Global R2A Acoustic')
ylabel('Global R2A Semantic')
zlabel('Global R2A Acoustic + Semantic')
hold off


figure()
hist([GlobalR2A.Acoustic(Expcells) GlobalR2A.Semantic(Expcells) GlobalR2A.AcSem(Expcells)], (-0.3:0.05:1))
legend('A', 'S', 'A+S')
%% Find examples of Synergistic cells, Semantic, Acoustic and merged cells and plot their R2 profils
SYN=find((IncAc>0.3) .* (IncSem>0.1));
SYN=SYN(2);%SYN=SYN(3) is great too
GlobalR2A.Acoustic(Expcells(SYN))
GlobalR2A.Semantic(Expcells(SYN))
GlobalR2A.AcSem(Expcells(SYN))
IncAc(SYN)
IncSem(SYN)
Path_SYN=List_matfilepath{Expcells(SYN)}
[P,FSYN,ext]=fileparts(Path_SYN);
plot_prediction_3models(['/Users/elie/Documents/MATLAB/data/matfile/Models' FSYN(8:end) ext], 300,125)

S=find(IncSem==max(IncSem(Expcells_S)))
S=find(GlobalR2A.Semantic(Expcells)==max(GlobalR2A.Semantic(Expcells(Expcells_S))))
S=Expcells_S(5)%S(8) is neat too!! PB: Profil of Semantic model only is highly correlated to the shrinkage of stims in each category as time increases
GlobalR2A.Acoustic(Expcells(S))
GlobalR2A.Semantic(Expcells(S))
GlobalR2A.AcSem(Expcells(S))
Path_S=List_matfilepath{Expcells(S)}
[P,FS]=fileparts(Path_S)
plot_prediction_3models(['/Users/elie/Documents/MATLAB/data/matfile/Models' FS(8:end) ext], 300)

Expcellsbin_A = zeros(length(Expcells),1);
Expcellsbin_A(Expcells_A)=1;
A=find((IncAc>0.45).*(Expcellsbin_A));
A=A(2);%selected A(1) too
GlobalR2A.Acoustic(Expcells(A))
GlobalR2A.Semantic(Expcells(A))
GlobalR2A.AcSem(Expcells(A))
Path_A=List_matfilepath{Expcells(A)}
[P,FA]=fileparts(Path_A)
plot_prediction_3models(['/Users/elie/Documents/MATLAB/data/matfile/Models' FA(8:end) ext], 300,105)

Expcellsbin_as = zeros(length(Expcells),1);
Expcellsbin_as(Expcells_as)=1;
as=find((GlobalR2A.Semantic(Expcells)>0.05).*(Expcellsbin_as));
as=as(end);%as(2) is also selected
GlobalR2A.Acoustic(Expcells(as))
GlobalR2A.Semantic(Expcells(as))
GlobalR2A.AcSem(Expcells(as))
Path_as=List_matfilepath{Expcells(as)}
[P,Fas]=fileparts(Path_as)
plot_prediction_3models(['/Users/elie/Documents/MATLAB/data/matfile/Models' Fas(8:end) ext], 300,65)


%% Investigate predicted and observed response at particular time point for example cells


%% Plot the histogram of the 3 Global R2A and distinguish the 4 cell
% populations
MAX=max([length(Expcells_SYN) length(Expcells_S) length(Expcells_A) length(Expcells_as)]);
figure()
subplot(3,1,1)
Vect1 = [GlobalR2A.Acoustic(Expcells(Expcells_SYN))' nan(1,MAX-length(Expcells_SYN))];
Vect2 = [GlobalR2A.Acoustic(Expcells(Expcells_S))' nan(1,MAX-length(Expcells_S))];
Vect3 = [GlobalR2A.Acoustic(Expcells(Expcells_A))' nan(1,MAX-length(Expcells_A))];
Vect4 = [GlobalR2A.Acoustic(Expcells(Expcells_as))' nan(1,MAX-length(Expcells_as))];
hist([Vect1' Vect2' Vect3' Vect4'], -0.2:0.05:1)
hold on
axis([-0.2 1 0 100])
YL=ylim;
L1=line([nanmean(Vect1) nanmean(Vect1)], YL);
set(L1, 'Color', [0 0 1],'LineWidth',1);
L2=line([nanmean(Vect2) nanmean(Vect2)], YL);
set(L2, 'Color', [0 0.9 0.9],'LineWidth',1);
L3=line([nanmean(Vect3) nanmean(Vect3)], YL);
set(L3, 'Color', [1 0.75 0],'LineWidth',1);
L4=line([nanmean(Vect4) nanmean(Vect4)], YL);
set(L4, 'Color', [0.5 0 0.1],'LineWidth',1);
legend('A+S', 'S', 'A', 'as')
title('Global R2A Acoustic')
hold off

subplot(3,1,2)
Vect1 = [GlobalR2A.Semantic(Expcells(Expcells_SYN))' nan(1,MAX-length(Expcells_SYN))];
Vect2 = [GlobalR2A.Semantic(Expcells(Expcells_S))' nan(1,MAX-length(Expcells_S))];
Vect3 = [GlobalR2A.Semantic(Expcells(Expcells_A))' nan(1,MAX-length(Expcells_A))];
Vect4 = [GlobalR2A.Semantic(Expcells(Expcells_as))' nan(1,MAX-length(Expcells_as))];
hist([Vect1' Vect2' Vect3' Vect4'], -0.2:0.05:1)
hold on
axis([-0.2 1 0 200])
YL=ylim;
L1=line([nanmean(Vect1) nanmean(Vect1)], YL);
set(L1, 'Color', [0 0 1],'LineWidth',1);
L2=line([nanmean(Vect2) nanmean(Vect2)], YL);
set(L2, 'Color', [0 0.9 0.9],'LineWidth',1);
L3=line([nanmean(Vect3) nanmean(Vect3)], YL);
set(L3, 'Color', [1 0.75 0],'LineWidth',1);
L4=line([nanmean(Vect4) nanmean(Vect4)], YL);
set(L4, 'Color', [0.5 0 0.1],'LineWidth',1);
legend('A+S', 'S', 'A', 'as')
title('Global R2A Semantic')

subplot(3,1,3)
Vect1 = [GlobalR2A.AcSem(Expcells(Expcells_SYN))' nan(1,MAX-length(Expcells_SYN))];
Vect2 = [GlobalR2A.AcSem(Expcells(Expcells_S))' nan(1,MAX-length(Expcells_S))];
Vect3 = [GlobalR2A.AcSem(Expcells(Expcells_A))' nan(1,MAX-length(Expcells_A))];
Vect4 = [GlobalR2A.AcSem(Expcells(Expcells_as))' nan(1,MAX-length(Expcells_as))];
hist([Vect1' Vect2' Vect3' Vect4'], -0.2:0.05:1)
hold on
axis([-0.2 1 0 100])
YL=ylim;
L1=line([nanmean(Vect1) nanmean(Vect1)], YL);
set(L1, 'Color', [0 0 1],'LineWidth',1);
L2=line([nanmean(Vect2) nanmean(Vect2)], YL);
set(L2, 'Color', [0 0.9 0.9],'LineWidth',1);
L3=line([nanmean(Vect3) nanmean(Vect3)], YL);
set(L3, 'Color', [1 0.75 0],'LineWidth',1);
L4=line([nanmean(Vect4) nanmean(Vect4)], YL);
set(L4, 'Color', [0.5 0 0.1],'LineWidth',1);
legend('A+S', 'S', 'A', 'as')
title('Global R2A AcSem')


%% Barplot of the mean +- 2SE of Global R2 for each cell type
GR2A_mean = [mean(GlobalR2A.Acoustic(Expcells(Expcells_SYN))) mean(GlobalR2A.Acoustic(Expcells(Expcells_S))) mean(GlobalR2A.Acoustic(Expcells(Expcells_A))) mean(GlobalR2A.Acoustic(Expcells(Expcells_as)))
    mean(GlobalR2A.Semantic(Expcells(Expcells_SYN))) mean(GlobalR2A.Semantic(Expcells(Expcells_S))) mean(GlobalR2A.Semantic(Expcells(Expcells_A))) mean(GlobalR2A.Semantic(Expcells(Expcells_as)))
    mean(GlobalR2A.AcSem(Expcells(Expcells_SYN))) mean(GlobalR2A.AcSem(Expcells(Expcells_S))) mean(GlobalR2A.AcSem(Expcells(Expcells_A))) mean(GlobalR2A.AcSem(Expcells(Expcells_as)))];

GR2A_se = [std(GlobalR2A.Acoustic(Expcells(Expcells_SYN)))./sqrt(length(Expcells_SYN)) std(GlobalR2A.Acoustic(Expcells(Expcells_S)))./sqrt(length(Expcells_S)) std(GlobalR2A.Acoustic(Expcells(Expcells_A)))./sqrt(length(Expcells_A)) std(GlobalR2A.Acoustic(Expcells(Expcells_as)))./sqrt(length(Expcells_as))
    std(GlobalR2A.Semantic(Expcells(Expcells_SYN)))./sqrt(length(Expcells_SYN)) std(GlobalR2A.Semantic(Expcells(Expcells_S)))./sqrt(length(Expcells_S)) std(GlobalR2A.Semantic(Expcells(Expcells_A)))./sqrt(length(Expcells_A)) std(GlobalR2A.Semantic(Expcells(Expcells_as)))./sqrt(length(Expcells_as))
    std(GlobalR2A.AcSem(Expcells(Expcells_SYN)))./sqrt(length(Expcells_SYN)) std(GlobalR2A.AcSem(Expcells(Expcells_S)))./sqrt(length(Expcells_S)) std(GlobalR2A.AcSem(Expcells(Expcells_A)))./sqrt(length(Expcells_A)) std(GlobalR2A.AcSem(Expcells(Expcells_as)))./sqrt(length(Expcells_as))];
figure()
%bar(GR2A_mean) 
errorb(GR2A_mean,GR2A_se)
ylabel('Adjusted Global R2');
xlabel('Model')
set(gca,'XTickLabel', {'Acoustic' 'Semantic' 'Acoustic + Semantic'});
legend('A+S', 'S', 'A', 'as')

%% Barplot of the number of cells in each type
bar([length(Expcells_SYN) length(Expcells_S) length(Expcells_A) length(Expcells_as)])
ylabel('Number of units')
xlabel('Unit Type')
set(gca, 'XTickLabel',{'SYN', 'S', 'A', 'as'})

%% 3D plots of where are these different populations of cells

ID_TypeCells = zeros(length(Expcells),1);
ID_TypeCells(Expcells_SYN)=1;
ID_TypeCells(Expcells_S)=2;
ID_TypeCells(Expcells_A)=3;
ID_TypeCells(Expcells_as)=4;
xyzCells = find(~isnan(DV_theo(Expcells)));
PlotSphere(LR_theo(Expcells(xyzCells)), RC_theo(Expcells(xyzCells)), DV_theo(Expcells(xyzCells)),GlobalR2A.AcSem(Expcells(xyzCells)),1,1,ID_TypeCells(xyzCells),{'A+S', 'S', 'A', 'as'},1);
xyzCells = find(~isnan(DV(Expcells)));
PlotSphere(LR(Expcells(xyzCells)), RC(Expcells(xyzCells)), DV(Expcells(xyzCells)),GlobalR2A.AcSem(Expcells(xyzCells)),1,1,ID_TypeCells(xyzCells),{'A+S', 'S', 'A', 'as'},1);
PlotSphere(LR(Expcells(xyzCells)), RC(Expcells(xyzCells)), DV(Expcells(xyzCells)),GlobalR2A.AcSem(Expcells(xyzCells)),0,1,ID_TypeCells(xyzCells),{'A+S', 'S', 'A', 'as'},1);

%% %of each cell type in each region
UZ = unique(ZONES(Expcells,1));
NZ = length(UZ);
PropCell_PerZone = nan(NZ-1,4);
for zz=1:(NZ-1)%only look at real region ("0" is unknown)
    Unitinzone = find(ZONES(Expcells,1) == zz);
    PropCell_PerZone(zz,1) = length(intersect(Unitinzone, Expcells_SYN))./length(Unitinzone);
    PropCell_PerZone(zz,2) = length(intersect(Unitinzone, Expcells_S))./length(Unitinzone);
    PropCell_PerZone(zz,3) = length(intersect(Unitinzone, Expcells_A))./length(Unitinzone);
    PropCell_PerZone(zz,4) = length(intersect(Unitinzone, Expcells_as))./length(Unitinzone);
end
figure()
bar(PropCell_PerZone,'stacked')
set(gca, 'XTickLabel',ZONES_List(1:end));

%% Plot the proportion of different cell types along time
Pval_AS_A = R2A.Pval_AS_A(Expcells,:);
Pval_AS_S = R2A.Pval_AS_S(Expcells,:);
NSynCells = sum((Pval_AS_A<0.05).*(Pval_AS_S<0.05));
NACells = sum((Pval_AS_A>=0.05).*(Pval_AS_S<0.05));
NSCells = sum((Pval_AS_A<0.05).*(Pval_AS_S>=0.05));
NasCells = sum((Pval_AS_A>=0.05).*(Pval_AS_S>=0.05));
figure()
plot(R2A.Wins(1:IndWin), NSynCells(1:IndWin), 'ko-','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
hold on
plot(R2A.Wins(1:IndWin), NSCells(1:IndWin),'ko-','MarkerFaceColor',[0 0.9 0.9],'MarkerEdgeColor',[0 0.9 0.9])
hold on
plot(R2A.Wins(1:IndWin), NACells(1:IndWin),'ko-','MarkerFaceColor',[1 0.75 0],'MarkerEdgeColor',[1 0.75 0])
hold on
plot(R2A.Wins(1:IndWin), NasCells(1:IndWin),'ko-','MarkerFaceColor',[0.5 0 0.1],'MarkerEdgeColor',[0.5 0 0.1])
hold off
ylabel('Number of units')
xlabel('Time (ms)')
axis([0 350 0 600])

%% What is the mean profil of R2A from encoding models of Expcells
MeanProfil.Acoustic = nanmean(R2A.Acoustic(Expcells,:));
MeanProfil.Semantic = nanmean(R2A.Semantic(Expcells,:));
MeanProfil.AcSem = nanmean(R2A.AcSem(Expcells,:));
MeanProfil.ASS = nanmean(R2A.AcSem(Expcells,:))-nanmean(R2A.Acoustic(Expcells,:));

figure()
%subplot(2,1,1)
plot(R2A.Wins(1:IndWin), MeanProfil.Acoustic(1:IndWin), 'ko-', R2A.Wins(1:IndWin), MeanProfil.Semantic(1:IndWin), 'bo-', R2A.Wins(1:IndWin),MeanProfil.AcSem(1:IndWin), 'ro-', R2A.Wins(1:IndWin), 10*MeanProfil.ASS(1:IndWin), 'go-')
legend('Acoustic', 'Semantic', 'Acoustic + Semantic', '10*((Acoustic + Semantic) - Acoustic)', 'Location', 'SouthWest')
title('Mean R2A as a function of the window size for Semantic Cells');
axis([0 350 -0.15 0.5])    

%% Plot the mean profil of R2A per cell population
% I choose to represent R2A AcSem for SYN cells (best model)
% R2A_Semantic for S cells (best model with fewer parameters)
% R2A_Acoustic for A cells (best model with fewer parameters)
% R2A_Semantic for as cells (best model with fewer parameters)
IdSynCells = (Pval_AS_A<0.05).*(Pval_AS_S<0.05);
IdACells = (Pval_AS_A>=0.05).*(Pval_AS_S<0.05);
IdSCells = (Pval_AS_A<0.05).*(Pval_AS_S>=0.05);
IdasCells = (Pval_AS_A>=0.05).*(Pval_AS_S>=0.05);
MeanProfil.SYN = nan(length(R2A.Wins),1);
MeanProfil.A = nan(length(R2A.Wins),1);
MeanProfil.S = nan(length(R2A.Wins),1);
MeanProfil.as = nan(length(R2A.Wins),1);


for tt=1:length(R2A.Wins)
    MeanProfil.SYN(tt) = nanmean(R2A.AcSem(find(IdSynCells(:,tt)),tt));
    MeanProfil.A(tt) = nanmean(R2A.Acoustic(find(IdACells(:,tt)),tt));
    MeanProfil.S(tt) = nanmean(R2A.Semantic(find(IdSCells(:,tt)),tt));
    MeanProfil.as(tt) = nanmean(R2A.Semantic(find(IdasCells(:,tt)),tt));
    MeanProfil.seSYN(tt) = nanstd(R2A.AcSem(find(IdSynCells(:,tt)),tt))./sqrt(sum(IdSynCells(:,tt)));
    MeanProfil.seA(tt) = nanstd(R2A.Acoustic(find(IdACells(:,tt)),tt))./sqrt(sum(IdACells(:,tt)));
    MeanProfil.seS(tt) = nanstd(R2A.Semantic(find(IdSCells(:,tt)),tt))./sqrt(sum(IdSCells(:,tt)));
    MeanProfil.seas(tt) = nanstd(R2A.Semantic(find(IdasCells(:,tt)),tt))./sqrt(sum(IdasCells(:,tt)));
end

figure()
errorb(R2A.Wins(1:IndWin),MeanProfil.SYN(1:IndWin), MeanProfil.seSYN(1:IndWin))
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.SYN(1:IndWin), 'ko-','MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
hold on
errorb(R2A.Wins(1:IndWin),MeanProfil.S(1:IndWin), MeanProfil.seS(1:IndWin))
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.S(1:IndWin),'ko-','MarkerFaceColor',[0 0.9 0.9],'MarkerEdgeColor',[0 0.9 0.9])
hold on
errorb(R2A.Wins(1:IndWin),MeanProfil.A(1:IndWin), MeanProfil.seA(1:IndWin))
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.A(1:IndWin),'ko-','MarkerFaceColor',[1 0.75 0],'MarkerEdgeColor',[1 0.75 0])
hold on
errorb(R2A.Wins(1:IndWin),MeanProfil.as(1:IndWin), MeanProfil.seas(1:IndWin))
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.as(1:IndWin),'ko-','MarkerFaceColor',[0.5 0 0.1],'MarkerEdgeColor',[0.5 0 0.1])
axis([0 350 -0.05 0.35])
hold off
ylabel('R2A')
xlabel('Time (ms)')
legend('AcSem-R2A SYN Cells', 'Sem-R2A S Cells','Ac-R2A A Cells', 'Sem-R2A as Cells')

%% Plot Av nb of categories and proportion of stims since the begining
MeanProfil.cat = nanmean(R2A.RNb_Cat);
for tt=1:IndWin
    MeanProfil.secat(tt) = nanstd(R2A.RNb_Cat(:,tt))./sqrt(sum(~isnan(R2A.RNb_Cat(:,tt))));
end
MaxStimNb = repmat(R2A.RNb_Stim(:,1),[1 IndWin]);
PropStim = 10.*R2A.RNb_Stim./MaxStimNb;
MeanProfil.propstim = nanmean(PropStim);
for tt=1:IndWin
    MeanProfil.sepropstim(tt) = nanstd(PropStim(:,tt))./sqrt(sum(~isnan(PropStim(:,tt))));
end

figure()
errorb(R2A.Wins(1:IndWin),MeanProfil.cat, MeanProfil.secat)
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.cat, 'ko-','MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0])
hold on
errorb(R2A.Wins(1:IndWin),MeanProfil.propstim, MeanProfil.sepropstim)
hold on
plot(R2A.Wins(1:IndWin), MeanProfil.propstim,'ko-','MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
legend('Number of categories', 'Proportion of stimuli')
xlabel('Time (ms)')
axis([0 350 0 10])


%% Plot nb of SYN cells in each zone over time
IdSynCells
UZ = unique(ZONES(Expcells,1));
NZ = length(UZ);

GZ = cell(3,1);
GZ{1}=3:7;%indices of L regions
GZ{2}=1:2;%indices of CM regions
GZ{3}=8;%index of NCM units
PropCell_PerZone = nan(length(GZ),length(R2A.Wins));
PropCell_PerZone_CI_H = nan(length(GZ),length(R2A.Wins));
PropCell_PerZone_CI_L = nan(length(GZ),length(R2A.Wins));

for zz=1:length(GZ)
    nz = length(GZ{zz});
    Unitinzone = zeros(length(Expcells),1);
    for zi = 1:nz
        ii=GZ{zz}(zi);
        Unitinzone = Unitinzone + ZONES(Expcells,1) == ii;
    end
    SynUnitinzone = repmat(Unitinzone,[1 length(R2A.Wins)]).*IdSynCells;
    PropCell_PerZone(zz,:) = sum(SynUnitinzone)./sum(Unitinzone);
    for tt=1:length(R2A.Wins)
        [phat,phci] = binofit(sum(SynUnitinzone(:,tt)),sum(Unitinzone));
        PropCell_PerZone_CI_H(zz,tt)=phci(2);
        PropCell_PerZone_CI_L(zz,tt)=phci(1);
    end
end
%GRAD2=[0.8 0.8 0.8; 0 1 1; 0 0.8 1; 1 0.3 0; 1 0 0; 1 0 0 ; 1 0 0.6; 1 0.2 0.5; 0 1 0; 0 0 0];
GRAD2=[0.8 0.8 0.8; 0 0 1; 0 0.5 1; 0 1 1; 0 0 0; 0 1 0 ; 1 1 0; 1 0.5 0; 1 0 0; 0.5 0 0];
GRAD2 = [1 0 0 ; 0 0.6 1; 0.1 1 0];
figure()
for zz=1:length(GZ)
    %errorbar(R2A.Wins(1:IndWin), PropCell_PerZone(zz,1:IndWin),PropCell_PerZone(zz,1:IndWin)-PropCell_PerZone_CI_L(zz,1:IndWin),PropCell_PerZone_CI_H(zz,1:IndWin)-PropCell_PerZone(zz,1:IndWin));
    hold on
    plot(R2A.Wins(1:IndWin), PropCell_PerZone(zz,1:IndWin),'ko-', 'MarkerFaceColor', GRAD2(zz,:),'MarkerEdgeColor', GRAD2(zz,:))
    hold on
end
legend({'Field L' 'CM' 'NCM'})
ylabel('Proportion of Synergistic units')
xlabel('Time (ms)')
axis([0 350 0 0.5])

set(gca, 'XTickLabel',ZONES_List(1:end));

%%
figure(6)
% Replace negative values of R2A  by 0 for Cubehilx plottings
R2AAcoustic = GlobalR2A.Acoustic;
R2AAcoustic(find(R2AAcoustic<0))=0.001;
subplot(1,2,1)
GRAD5=cubehelix(ceil(max(GlobalR2A.Acoustic*1000)), 0.5, -0.8, 1.5, 1.7, [1,0]);
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), DiffSemP(jj), 'ko', 'MarkerFaceColor',GRAD5(round(R2AAcoustic(jj)*1000),:));
    end
    hold on
end

xlabel('Mutual Information of the IV matrix (bits)')
ylabel('MiCat - MiCatRandBGP')
title(sprintf('Semantic Cells\nz-axis: Acoustic linearity (R2A acoustic model)'));
colorbar()
colormap(GRAD5)
%axis([0 6 -0.2 0.5])
hold off

subplot(1,2,2)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD5(round(R2AAcoustic(jj)*1000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Acoustic linearity (R2A acoustic model)'));
colorbar()
colormap(GRAD5)
axis([0 6 0 1.4])

figure(7)
% Replace negative values of R2A  by 0 for Cubehelix plottings
R2ASemantic = GlobalR2A.Semantic;
R2ASemantic(find(R2ASemantic<0))=0.001;
subplot(2,2,1)
GRAD6=cubehelix(ceil(max(GlobalR2A.Semantic*1000)), 0.5, -0.8, 1.5, 1.7, [1,0]);
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), DiffSemP(jj), 'ko', 'MarkerFaceColor',GRAD6(round(R2ASemantic(jj)*1000),:));
    end
    hold on
end
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('MiCat - MiCatRandBGP')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A semantic model)'));
colorbar()
colormap(GRAD6)
%axis([0 6 -0.2 0.5])
hold off

subplot(2,2,2)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD6(round(R2ASemantic(jj)*1000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A semantic model)'));
colorbar()
colormap(GRAD6)
axis([0 6 0 1.4])

subplot(2,2,3)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(DiffSemP(jj), GlobalR2A.Semantic(jj), 'ko', 'MarkerFaceColor',[0.5 0.1 0.2]);
        hold on
    else
        plot(DiffSemP(jj), GlobalR2A.Semantic(jj), 'ko', 'MarkerFaceColor',[1 1 1]);
        hold on
    end
end
hold off
xlabel('Non-linear Semantic Index')
ylabel('R2A Semantic')
title(sprintf('Semantic Cells colored\nNon-semantic cells white'));

subplot(2,2,4)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(SemanticIndex.Observed(jj), GlobalR2A.Semantic(jj), 'ko', 'MarkerFaceColor',[0.5 0.1 0.2]);
        hold on
    else
        plot(SemanticIndex.Observed(jj), GlobalR2A.Semantic(jj), 'ko', 'MarkerFaceColor',[1 1 1]);
        hold on
    end
end
hold off
xlabel('Semantic Index')
ylabel('R2A Semantic')
title(sprintf('Semantic Cells colored\nNon-semantic cells white'));

figure(8)
NonLinear = GlobalR2A.AcSem - GlobalR2A.Acoustic;
% get rid of negative values for GRAD7
NonLinear_zero = NonLinear;
NonLinear_zero(find(NonLinear_zero < 0))=0.001;
subplot(2,2,1)
GRAD7=cubehelix(ceil(max(NonLinear_zero*1000)), 0.5, -1.1, 1.5, 0.8, [1,0]);
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), DiffSemP(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
    end
    hold on
end
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('MiCat - MiCatRandBGP')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A AcSem-Acoustic model)'));
colorbar()
colormap(GRAD7)
%axis([0 6 -0.2 0.5])
hold off

subplot(2,2,2)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A AcSem-Acoustic model)'));
colorbar()
colormap(GRAD7)
axis([0 6 0 1.4])

subplot(2,2,3)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(DiffSemP(jj), NonLinear(jj), 'ko', 'MarkerFaceColor',[0.5 0.1 0.2]);
        hold on
    else
        plot(DiffSemP(jj), NonLinear(jj), 'ko', 'MarkerFaceColor',[1 1 1]);
        hold on
    end
end
hold off
xlabel('Non-linear Semantic Index')
ylabel('R2AAcSem - R2AAcoustic')
title(sprintf('Semantic Cells colored\nNon-semantic cells white'));

subplot(2,2,4)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(SemanticIndex.Observed(jj), NonLinear(jj), 'ko', 'MarkerFaceColor',[0.5 0.1 0.2]);
        hold on
    else
        plot(SemanticIndex.Observed(jj), NonLinear(jj), 'ko', 'MarkerFaceColor',[1 1 1]);
        hold on
    end
end
hold off
xlabel('Semantic Index')
ylabel('R2AAcSem - R2AAcoustic')
title(sprintf('Semantic Cells colored\nNon-semantic cells white'));

figure(9)
NonLinear2 = GlobalR2A.AcSem - GlobalR2A.Semantic;
% get rid of under zero values for r2a
NonLinear2_zero = NonLinear2;
NonLinear2_zero(find(NonLinear2_zero < 0))=0.001;
subplot(1,2,1)
GRAD8=cubehelix(ceil(max(NonLinear2*1000)), 0.5, -1.1, 1.5, 0.8, [0,1]);
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), DiffSemP(jj), 'ko', 'MarkerFaceColor',GRAD8(round(NonLinear2_zero(jj)*1000),:));
    end
    hold on
end
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('MiCat - MiCatRandBGP')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A AcSem-Semantic model)'));
colorbar()
colormap(GRAD8)
%axis([0 6 -0.2 0.5])
hold off

subplot(1,2,2)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_tot(jj), SemanticIndex.Observed(jj), 'ko', 'MarkerFaceColor',GRAD8(round(NonLinear2_zero(jj)*1000),:));
        hold on
    end
end
hold off
xlabel('Mutual Information of the IV matrix (bits)')
ylabel('Semantic Index')
title(sprintf('Semantic Cells\nz-axis: Semantic non-linearity (R2A AcSem-Semantic model)'));
colorbar()
colormap(GRAD8)
axis([0 6 0 1.4])

figure(10)
GRAD2=cubehelix(ceil(max(MeanSR*1000)), 0.5, -1.1, 1.5, 0.5, [1,0]);
GRAD9 = cubehelix(ceil(max(MI_perArea.MI_diag_uni_cat*1000)), 0.5, -1.1, 1.5, 0.5, [1,0]);
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 %&& logpvRandBG(jj)>=2 && SemanticIndex.NonLinearBG(jj)>0
        subplot(1,3,2)
        plot(GlobalR2A.Acoustic(jj), GlobalR2A.AcSem(jj), 'ko', 'MarkerFaceColor',[0.6 0.8 0.9]);
        hold on
    else
        subplot(1,3,2)
        plot(GlobalR2A.Acoustic(jj), GlobalR2A.AcSem(jj), 'ko', 'MarkerFaceColor',[0.2 0.2 0.2]);
        hold on
        
    end
    subplot(1,3,1)
    plot(GlobalR2A.Acoustic(jj), GlobalR2A.AcSem(jj), 'ko', 'MarkerFaceColor',GRAD9(ceil(MI_perArea.MI_diag_uni_cat(jj)*1000),:));
    hold on
    subplot(1,3,3)
    plot(GlobalR2A.Acoustic(jj), GlobalR2A.AcSem(jj), 'ko', 'MarkerFaceColor',GRAD2(ceil(MeanSR(jj)*1000),:));
    hold on
end
subplot(1,3,2)
xlabel('R2A Acoustic')
ylabel('R2A Acoustic + Semantic')
title(sprintf('All Cells\n'));
line([-1 1], [-1 1])
axis([-0.5 1 -0.5 1])
hold off
subplot(1,3,1)
xlabel('R2A Acoustic')
ylabel('R2A Acoustic + Semantic')
title(sprintf('All Cells\nz-axis MI_cat'));
line([-1 1], [-1 1])
axis([-0.5 1 -0.5 1])
colorbar
colormap(GRAD9)
hold off
subplot(1,3,3)
xlabel('R2A Acoustic')
ylabel('R2A Acoustic + Semantic')
title(sprintf('All Cells\nz-axis Spike Rate Hz'));
line([-1 1], [-1 1])
axis([-0.5 1 -0.5 1])
colorbar()
colormap(GRAD2)
hold off

figure(11)
subplot(2,2,1)
for jj=1:NU
    plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
    hold on
end
    xlabel('R2A Acoustic')
    ylabel('R2A Acoustic + Semantic - R2AAcoustic')
    title(sprintf('All Cells\n'));
    colorbar()
    colormap(GRAD6)
    hline(0)
    hold off
    axis([-0.2 0.8 -0.3 0.2])
subplot(2,2,2)
for jj=1:NU
    if PAS_Ac(jj)<0.01
        plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
%     else
%         plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
%         hold on
    end
end
    xlabel('R2A Acoustic')
    ylabel('R2A Acoustic + Semantic - R2AAcoustic')
    title(sprintf('Significantly Semantic R2A cells\n'));
    colorbar()
    colormap(GRAD6)
    hline(0)
    hold off
    axis([-0.2 0.8 -0.3 0.2])
    
subplot(2,2,3)
for jj=1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
%     else
%         plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
%         hold on
    end
end
    xlabel('R2A Acoustic')
    ylabel('R2A Acoustic + Semantic - R2AAcoustic')
    title(sprintf('Significantly Semantic ConfMat cells\n'));
    colorbar()
    colormap(GRAD6)
    hline(0)
    hold off
    axis([-0.2 0.8 -0.3 0.2])
    
    subplot(2,2,4)
for jj=1:NU
    if MI_perArea.MI_diag_uni_cat(jj)>=0.87
        plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
%     else
%         plot(GlobalR2A.Acoustic(jj),GlobalR2A.AcSem(jj)-GlobalR2A.Acoustic(jj), 'ko','MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
%         hold on
    end
end
    xlabel('R2A Acoustic')
    ylabel('R2A Acoustic + Semantic - R2AAcoustic')
    title(sprintf('Super Semantic ConfMat cells\n'));
    colorbar()
    colormap(GRAD6)
    hline(0)
    hold off
    axis([-0.2 0.8 -0.3 0.2])
    


% subplot(2,2,3)
% for jj = 1:NU
%     if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
%         plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD11(ceil(SemIndex_zero(jj)*1000),:));
%         hold on
%     end
% end
% xlabel('Acoustic dimension (MiDiag)')
% ylabel('Non-linear Semantic (MiCat - MiDiag)')
% title(sprintf('Significantly Semantic ConfMat cells\nz-axis: Semantic dimension (MiCat - RandMiCat)'))
% colorbar()
% colormap(GRAD11)
% 
% subplot(2,2,4)
% 
% for jj = 1:NU
%     if PAS_Ac(jj)<0.01
%         plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD11(ceil(SemIndex_zero(jj)*1000),:));
%         hold on
%     end
% end
% xlabel('Acoustic dimension (MiDiag)')
% ylabel('Non-linear Semantic (MiCat - MiDiag)')
% title(sprintf('Significantly Semantic R2A cells\nz-axis: Semantic dimension (MiCat - RandMiCat)'))
% colorbar()
% colormap(GRAD11)

figure(12)
for jj = 1:NU
    if MI_perArea.MI_diag_uni_cat(jj)>=0.863
        if PAS_Ac(jj)<0.01
            plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'kp', 'MarkerFaceColor',GRAD11(ceil(SemIndex_zero(jj)*1000),:), 'MarkerSize', 8);
            hold on
        else
            plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD11(ceil(SemIndex_zero(jj)*1000),:));
            hold on
        end
    end
end
xlabel('Individual sound information')
ylabel(sprintf('Invariance information'))
title(sprintf('Semantic Units\nz-axis: Valuable Semantic information'))
 colorbar('YTickLabel', [0.5 1 1.5 2 2.5])
colormap(GRAD11) 
axis([0 3 0 1.5])
legend('Linear Semantic Units', 'Non-Linear Semantic Units', 'Location', 'SouthEast')




figure(13)
subplot(2,2,1)
for jj = 1:NU
    plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
    hold on
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('All Cells\nz-axis: Semantic dimension (R2A Semantic)'))
colorbar()
colormap(GRAD6)

subplot(2,2,2)
for jj = 1:NU
    if PAS_Ac(jj)<0.01
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic R2A cells\nz-axis: Semantic dimension (R2A Semantic)'))
colorbar()
colormap(GRAD6)

subplot(2,2,3)
for jj = 1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic ConfMat cells\nz-axis: Semantic dimension (R2A Semantic)'))
colorbar()
colormap(GRAD6)

subplot(2,2,4)
for jj = 1:NU
    if MI_perArea.MI_diag_uni_cat(jj)>=0.87
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD6(ceil(R2ASemantic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Super Semantic ConfMat cells\nz-axis: Semantic dimension (R2A Semantic)'))
 colorbar()
colormap(GRAD6) 
axis([0 3 0 1.4])

figure(14)
subplot(2,2,1)
for jj = 1:NU
    plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD5(ceil(R2AAcoustic(jj)*1000),:));
    hold on
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('All Cells\nz-axis: Acoustic dimension (R2A Acoustic)'))
colorbar()
colormap(GRAD5)

subplot(2,2,2)
for jj = 1:NU
    if PAS_Ac(jj)<0.01
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD5(ceil(R2AAcoustic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic R2A cells\nz-axis: Acoustic dimension (R2A Acoustic)'))
colorbar()
colormap(GRAD5)

subplot(2,2,3)
for jj = 1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD5(ceil(R2AAcoustic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic ConfMat cells\nz-axis: Acoustic dimension (R2A Acoustic)'))
colorbar()
colormap(GRAD5)

subplot(2,2,4)
for jj = 1:NU
    if MI_perArea.MI_diag_uni_cat(jj)>=0.87
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD5(ceil(R2AAcoustic(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Super Semantic ConfMat cells\nz-axis: Acoustic dimension (R2A Acoustic)'))
 colorbar()
colormap(GRAD5) 
axis([0 3 0 1.4])

figure(15)
subplot(2,2,1)
for jj = 1:NU
    plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
    hold on
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('All Cells\nz-axis: NonLinear Semantic dimension (R2A AcSem - R2A Acoustic)'))
colorbar()
colormap(GRAD7)

subplot(2,2,2)
for jj = 1:NU
    if PAS_Ac(jj)<0.01
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic R2A cells\nz-axis: NonLinear Semantic dimension (R2A AcSem - R2A Acoustic)'))
colorbar()
colormap(GRAD7)

subplot(2,2,3)
for jj = 1:NU
    if logpvRandP(jj)>=2 && SemanticIndex.NonLinearP(jj)>0 && logpvRandBGP(jj)>=2 && SemanticIndex.NonLinearBGP(jj)>0
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Significantly Semantic ConfMat cells\nz-axis: NonLinear Semantic dimension (R2A AcSem - R2A Acoustic)'))
colorbar()
colormap(GRAD7)

subplot(2,2,4)
for jj = 1:NU
    if MI_perArea.MI_diag_uni_cat(jj)>=0.87
        plot(MI_perArea.MI_diag_uni(jj), MI_perArea.MI_diag_uni_cat(jj)-MI_perArea.MI_diag_uni(jj), 'ko', 'MarkerFaceColor',GRAD7(ceil(NonLinear_zero(jj)*1000),:));
        hold on
    end
end
xlabel('Acoustic dimension (MiDiag)')
ylabel('Non-linear Semantic (MiCat - MiDiag)')
title(sprintf('Super Semantic ConfMat cells\nz-axis: NonLinear Semantic dimension (R2A AcSem - R2A Acoustic)'))
 colorbar()
colormap(GRAD7) 
axis([0 3 0 1.4])

%% Investigating semantic cells and encoding models
% semantic cells are cells that present MICat>=0.863
SemNLCell = intersect(SemCell, find(PAS_Ac<=0.01));
SemLCell = intersect(SemCell, find(PAS_Ac>0.01));

mean(GlobalR2A.Acoustic(SemCell))
std(GlobalR2A.Acoustic(SemCell))
mean(GlobalR2A.Acoustic(SemNLCell))
std(GlobalR2A.Acoustic(SemNLCell))
mean(GlobalR2A.Acoustic(SemLCell))
std(GlobalR2A.Acoustic(SemLCell))
mean(GlobalR2A.Acoustic(NonSemCell))
std(GlobalR2A.Acoustic(NonSemCell))
sum(PAc(SemCell)<=0.01)/length(SemCell)
sum(PAc(NonSemCell)<=0.01)/length(NonSemCell)

mean(GlobalR2A.Semantic(SemCell))
std(GlobalR2A.Semantic(SemCell))
mean(GlobalR2A.Semantic(SemNLCell))
std(GlobalR2A.Semantic(SemNLCell))
mean(GlobalR2A.Semantic(SemLCell))
std(GlobalR2A.Semantic(SemLCell))
mean(GlobalR2A.Semantic(NonSemCell))
std(GlobalR2A.Semantic(NonSemCell))
sum(PSem(SemCell)<=0.01)/length(SemCell)
sum(PSem(NonSemCell)<=0.01)/length(NonSemCell)

mean(GlobalR2A.AcSem(SemCell))
std(GlobalR2A.AcSem(SemCell))
mean(GlobalR2A.AcSem(SemNLCell))
std(GlobalR2A.AcSem(SemNLCell))
mean(GlobalR2A.AcSem(SemLCell))
std(GlobalR2A.AcSem(SemLCell))
mean(GlobalR2A.AcSem(NonSemCell))
std(GlobalR2A.AcSem(NonSemCell))
sum(PAcSem(SemCell)<=0.01)/length(SemCell)
sum(PAcSem(NonSemCell)<=0.01)/length(NonSemCell)

mean(GlobalR2A.AcSem(SemCell)-GlobalR2A.Acoustic(SemCell))
std(GlobalR2A.AcSem(SemCell)-GlobalR2A.Acoustic(SemCell))
mean(GlobalR2A.AcSem(SemNLCell)-GlobalR2A.Acoustic(SemNLCell))
std(GlobalR2A.AcSem(SemNLCell)-GlobalR2A.Acoustic(SemNLCell))
mean(GlobalR2A.AcSem(SemLCell)-GlobalR2A.Acoustic(SemLCell))
std(GlobalR2A.AcSem(SemLCell)-GlobalR2A.Acoustic(SemLCell))
mean(GlobalR2A.AcSem(NonSemCell)-GlobalR2A.Acoustic(NonSemCell))
std(GlobalR2A.AcSem(NonSemCell)-GlobalR2A.Acoustic(NonSemCell))

% proportion of nonlinear semantic cells
sum(PAS_Ac(SemCell)<=0.01)/length(SemCell)

% Cell with the maximum boost with Semantic info add to the acoustic model
BigBoost = SemNLCell(find(NonLinear(SemNLCell)==max(NonLinear(SemNLCell))));
List_matfilepath{BigBoost}
MI_perArea.MI_diag_uni(BigBoost)
MI_perArea.MI_diag_uni_cat(BigBoost) - MI_perArea.MI_diag_uni(BigBoost)

% Cell with the max info in diagonal
SuperLinear = SemLCell(find(MI_perArea.MI_diag_uni(SemLCell)==max(MI_perArea.MI_diag_uni(SemLCell))));
List_matfilepath{SuperLinear}

% Linear Cell with the highest linear R2A ACoustic
MaxR2AA=SemLCell(find(GlobalR2A.Acoustic(SemLCell)==max(GlobalR2A.Acoustic(SemLCell))))
List_matfilepath{MaxR2AA}
MI_perArea.MI_diag_uni(MaxR2AA)
MI_perArea.MI_diag_uni_cat(MaxR2AA) - MI_perArea.MI_diag_uni(MaxR2AA)

% Cell with the max invariance
Inv=MI_perArea.MI_diag_uni_cat - MI_perArea.MI_diag_uni;
MaxInv = find(Inv==max(Inv));
List_matfilepath{MaxInv}
MI_perArea.MI_diag_uni(MaxInv)
MI_perArea.MI_diag_uni_cat(MaxInv) - MI_perArea.MI_diag_uni(MaxInv)

% Cell with R2A Sem>R2A Acoustic
NL = GlobalR2A.Semantic-GlobalR2A.Acoustic;
MaxNL=SemCell(find(NL(SemCell)==max(NL(SemCell))))
List_matfilepath{MaxNL}
MI_perArea.MI_diag_uni(MaxNL)
MI_perArea.MI_diag_uni_cat(MaxNL) - MI_perArea.MI_diag_uni(MaxNL)

