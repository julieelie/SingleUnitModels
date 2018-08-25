function myplotyyInfoSPAsemOnly(MatrixPlot_x_left,MatrixPlot_y_left,Plot_x_right,Plot_y_right,LegendInfo, YLim,Transp)
% fprintf('enlarge figure window \n')
% pause()
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
        h1(2).Color = [0 0.8 0 Transp];
        h1(3).Color = [0.4 0.2 0 Transp];
        h2.Color = [0 0 0 Transp];
        h2.LineStyle = '--';
        if length(h1)>3
            h1(4).Color = [0.5 0.5 0.5 Transp];
        end
end
if nargin<6
    YLim=[0 max(MaxLeft)];
end

ax(1).YLim = YLim;
ax(1).YTick = 0:0.1:YLim(2);
ax(1).YTickLabel = 0:0.1:YLim(2);
xlabel('windows (ms)')
ax(2).YLim = [0 max(MaxRight)];
ax(2).YColor = 'k';
ax(2).YTick = 0:10:max(MaxRight);
ax(2).YTickLabel = 0:10:max(MaxRight);
ylabel(ax(2),sprintf('Average Spike Rate (s) %s',LegendInfo.CellType))
Xtickposition=get(ax(1),'XTick');
set(ax(1),'XTickLabel', Xtickposition*10)
if length(h1)>5
    legend('Ceiling', 'Semantic','Null-Model','Saturated', 'Location','NorthEast')
    ylabel(ax(1),sprintf('%s',LegendInfo.YleftAxis))
else
    legend('Ceiling', 'Semantic','Null-Model', 'Location','NorthEast')
    ylabel(ax(1),sprintf('%s',LegendInfo.YleftAxis))
end
title(LegendInfo.CellType);
hold off

end
