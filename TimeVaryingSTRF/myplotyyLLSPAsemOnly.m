function myplotyyLLSPAsemOnly(MatrixPlot_x_left,MatrixPlot_y_left,Plot_x_right,Plot_y_right,CellType,XLabels_in, Transp)
% plot function used for ICN2016 poster
fprintf('enlarge figure window \n')
pause()
if nargin<7
    Transp=1;
end

if length(size(MatrixPlot_y_left))==3
    ElementNum = size(MatrixPlot_y_left,3);
    MinLeft = nan(ElementNum,1);
    MaxRight = nan(ElementNum,1);
else
    ElementNum=1;
    MinLeft = nan(ElementNum,1);
    MaxRight = nan(ElementNum,1);
end

for ee=1:ElementNum
        hold on
        [ax,h1,h2]=plotyy(MatrixPlot_x_left, MatrixPlot_y_left(:,:,ee),Plot_x_right,Plot_y_right(:,ee));
        MinLeft(ee) = floor((min(median(MatrixPlot_y_left(:,:,ee),2))-0.5*abs(min(median(MatrixPlot_y_left(:,:,ee),2))))*10)/10;
        MaxRight(ee) = ceil((max(Plot_y_right(:,ee)) + 0.2 * max(Plot_y_right(:,ee)))*10)/10;
        h1(1).Color = [0 0.8 0.8 Transp];
        h1(2).Color = [0 0.8 0 Transp];
        h1(3).Color = [0.4 0.2 0 Transp];
        h1(4).Color = [0.5 0.5 0.5 Transp];
        h2.Color = [0 0 0 Transp];
        h2.LineStyle = '--';
end

ax(1).YLim = [min(MinLeft) 0];
ax(1).YTick = min(MinLeft):0.5:0;
ax(1).YTickLabel = min(MinLeft):0.5:0;
xlabel('windows (ms)')
ax(2).YLim = [0 max(MaxRight)];
ax(2).YColor = 'k';
ax(2).YTick = 0:10:max(MaxRight);
ax(2).YTickLabel = 0:10:max(MaxRight);
ylabel(ax(2),sprintf('Average Spike Rate (s) %s',CellType))
Xtickposition=get(ax(1),'XTick');
set(ax(1),'XTickLabel', [0 XLabels_in(Xtickposition(2:end))])
if length(h1)>3
    legend('Ceiling', 'Semantic','Null-Model','Saturated', 'Location','SouthEast')
    ylabel(ax(1),sprintf('Average LogLikelihood %s',CellType))
else
    legend('Ceiling','Semantic','Null-Model', 'Location','SouthEast')
    ylabel(ax(1),sprintf('Average diff LogLikelihood with Saturated Model %s',CellType))
end
hold off

end
