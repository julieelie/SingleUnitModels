function myplotyyInfoSPA(MatrixPlot_x_left,MatrixPlot_y_left,Plot_x_right,Plot_y_right,CellType,XLabels_in, Transp)
fprintf('enlarge figure window \n')
pause()
if nargin<7
    Transp=1;
end

if length(size(MatrixPlot_y_left))==3
    ElementNum = size(MatrixPlot_y_left,3);
    MaxLeft = nan(ElementNum,1);
    MaxRight = nan(ElementNum,1);
else
    ElementNum=1;
    MaxLeft = nan(ElementNum,1);
    MaxRight = nan(ElementNum,1);
end

for ee=1:ElementNum
        hold on
        [ax,h1,h2]=plotyy(MatrixPlot_x_left, MatrixPlot_y_left(:,:,ee),Plot_x_right,Plot_y_right(:,ee));
        MaxLeft(ee) = ceil((max(max(MatrixPlot_y_left(:,:,ee)))+0.2*abs(max(max(MatrixPlot_y_left(:,:,ee)))))*10)/10;
        MaxRight(ee) = ceil((max(Plot_y_right(:,ee)) + 0.1 * max(Plot_y_right(:,ee)))*10)/10;
        h1(1).Color = [0 0.8 0.8 Transp];
        h1(2).Color = [0 0.5 1 Transp];
        h1(3).Color = [1 0.2 1 Transp];
        h1(4).Color = [0 0.8 0 Transp];
        h1(5).Color = [0.4 0.2 0 Transp];
        h2.Color = [0 0 0 Transp];
        h2.LineStyle = '--';
        if length(h1)>5
            h1(6).Color = [0.5 0.5 0.5 Transp];
        end
end

ax(1).YLim = [0 max(MaxLeft)];
ax(1).YTick = 0:0.1:max(MaxLeft);
ax(1).YTickLabel = 0:0.1:max(MaxLeft);
xlabel('windows (ms)')
ax(2).YLim = [0 max(MaxRight)];
ax(2).YColor = 'k';
ax(2).YTick = 0:10:max(MaxRight);
ax(2).YTickLabel = 0:10:max(MaxRight);
ylabel(ax(2),sprintf('Average Spike Rate (s) %s',CellType))
Xtickposition=get(ax(1),'XTick');
set(ax(1),'XTickLabel', [0 XLabels_in(Xtickposition(2:end))])
if length(h1)>5
    legend('Ceiling', 'Acoustic', 'AcSem', 'Semantic','Null-Model','Saturated', 'Location','NorthEast')
    ylabel(ax(1),sprintf('Average Information (bits) %s',CellType))
else
    legend('Ceiling', 'Acoustic', 'AcSem', 'Semantic','Null-Model', 'Location','NorthEast')
    ylabel(ax(1),sprintf('Average Information (bits) %s',CellType))
end
hold off

end
