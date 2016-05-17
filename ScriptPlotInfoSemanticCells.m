load('Models_GLMPoisson_Site4_L1500R1900_e23_s0_ss2.mat')

% Find the number of windows that run
WinMax=sum(~isnan(Model.Ceiling.info));

% Calculate the mean spike rate
Mean_SR = nan(WinMax,1);
for ww=1:WinMax
    Mean_SR(ww)=mean(Data.y_wholeset{ww});
end

figure()
% Information
subplot(1,3,1)
Info_MatrixPlot_y_left = [Model.Ceiling.info(1:WinMax), Model.Semantic.info(1:WinMax), Model.Floor.info(1:WinMax)];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR;
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 0.7])

% Cumulative Information
subplot(1,3,2)
Info_MatrixPlot_y_left = [Model.Ceiling.cum_info(1:WinMax), Model.Semantic.cum_info(1:WinMax), Model.Floor.cum_info(1:WinMax)];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR;
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Cumulative Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 2])

% Cumulative sum of information
subplot(1,3,3)
Info_MatrixPlot_y_left = [cumsum(Model.Ceiling.info(1:WinMax)), cumsum(Model.Semantic.info(1:WinMax)), cumsum(Model.Floor.info(1:WinMax))];
Info_MatrixPlot_x_left = repmat((1:WinMax)',1,3);
SPA_Plot_y_right = Mean_SR;
SPA_Plot_x_right = 1:WinMax;
LegendInfo.CellType = 'One Cell';
LegendInfo.YleftAxis = 'Cumulative Sum of Information (bits)';
myplotyyInfoSPAsemOnly(Info_MatrixPlot_x_left,Info_MatrixPlot_y_left,SPA_Plot_x_right,SPA_Plot_y_right,LegendInfo,Wins,[-0.1 2])


figure(3)
plot(1:WinMax,Model.AR.info(1:WinMax),'c', 1:WinMax, Model.Ceiling.info(1:WinMax), 'g', 1:WinMax, Model.Floor.info(1:WinMax), 'k', 1:WinMax, Model.Semantic.info(1:WinMax), 'r')
legend('Auto-regressive', 'Ceiling', 'Floor', 'Semantic', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Information (bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)
figure(4)
plot(1:WinMax,Model.AR.cum_info(1:WinMax), 'c', 1:WinMax, Model.Ceiling.cum_info(1:WinMax), 'g', 1:WinMax, Model.Floor.cum_info(1:WinMax), 'k', 1:WinMax, Model.Semantic.cum_info(1:WinMax), 'r')
hold on
plot(1:WinMax,[0; cumsum(Model.AR.info(2:WinMax))],'c--', 1:WinMax, cumsum(Model.Ceiling.info(1:WinMax)), 'g--', 1:WinMax, cumsum(Model.Floor.info(1:WinMax)), 'k--', 1:WinMax, cumsum(Model.Semantic.info(1:WinMax)), 'r--')
legend('CI Auto-regressive', 'CI Ceiling', 'CI Floor', 'CI Semantic','CSI Auto-regressive', 'CSI Ceiling', 'CSI Floor', 'CSI Semantic', 'Location','NorthWest')
xlabel('Time (ms)')
ylabel('Cumulative Information (CI, bits) and Cumulative Sum of Information (CSI, bits)')
IndLab=get(gca,'XTickLabel');
XTickLab = nan(length(IndLab),1);
for ii=1:length(XTickLab)
    XTickLab(ii) = (str2num(IndLab{ii})+1)*10;
end
set(gca,'XTickLabel', XTickLab)

SpeedInfo.Semantic=[0.0604 0.0969 0.3815 2.195 11.74 49.56 181.21 552.04 4181.1 10927.3 1520.2 3193.6 18046]
SpeedInfo.Floor=[0.0616 0.0915 0.223 0.974 3.855 16.2 48.32 2500.2 1453.1 9849.6 5628.9 9744.1 6261.2]
SpeedInfo.Ceiling=[0.06466 0.1237 0.736 6.323 35.59 220.61 828.74 1938 5108.5 734.9 20707.5 34805.8 67347.3]
SpeedInfo.AR=[0.06634 0.1632 1.269 5.981 27.93 228.68 783.04 144.09  318.5 2734.4 17373.5 41155.3 0]
WinMax_local = length(SpeedInfo.Semantic);
figure()
plot(1:WinMax_local,SpeedInfo.Semantic, 1:WinMax_local,SpeedInfo.Floor,1:WinMax_local,SpeedInfo.Ceiling,1:WinMax_local,SpeedInfo.AR)
hold on
plot(1:WinMax_local, max([SpeedInfo.Semantic; SpeedInfo.Floor; SpeedInfo.Ceiling;SpeedInfo.AR],[],1), '--r')
legend('Semantic','Floor','Ceiling','AR','Max','Location','NorthWest')