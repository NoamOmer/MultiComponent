
function [xdata,ydata,zdata]=getDataFromFigure(figureHandle)

xdata=[];ydata=[];zdata=[];
if  isempty(figureHandle) figureHandle = gcf; end

xdata=figureHandle.XData;
ydata=figureHandle.YData;
zdata=figureHandle.ZData;

axesObjs = get(figureHandle, 'Children');
dataObjs = get(axesObjs, 'Children');
% for i=2:length(dataObjs)
%     objTypes = get(dataObjs{i}, 'Type');
%     xdata = get(dataObjs{i}, 'XData');
%     ydata = get(dataObjs{i}, 'YData');
%     zdata = get(dataObjs{i}, 'ZData');
% end

end


