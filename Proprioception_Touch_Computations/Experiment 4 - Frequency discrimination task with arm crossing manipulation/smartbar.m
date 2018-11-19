function X = smartbar(data, barWidth, gapBetweenGroups, barGroupNames, color )
% X = smartbar(data, barWidth, interGroupWidth, color, barGroupNames )
%
% Each column of data represents one group
% barWidth = width of each bar, 1: touching; >1: overlap
% gapBetweenGroups = gap between two adjacent groups; 0: touching
% color = specify RGB values for each color; total number of colors must be
% greater or equal to the total number of bars in a single group.
% barGroupNames = specify xticklabels under each group
% Return a vector X with the location of each bar. Can be used to draw
% errorbars. Set sem for each data point, and run the command:
% errorbar(X(:),data(:),sem(:),'o')
%

[r,c] = size(data);
if nargin < 5
    color = setplotenv(2);
    color = repmat(color,ceil(r/7),1);
end
if nargin < 4
    barGroupNames = 1:c;
end
if nargin < 3
    gapBetweenGroups = 1;
end
if nargin < 2
    barWidth = 1;
end

X = [];
for i = 1:c
    x = (0 : r-1) + (i-1)*(r+gapBetweenGroups);
    X = [X; x];
    for j = 1:r
        if strcmp(color,'none')
            bar(x(j),data(j,i),barWidth,'FaceColor',color,'EdgeColor','k'), hold on
        else
            bar(x(j),data(j,i),barWidth,'FaceColor',color(j,:),'EdgeColor','none'), hold on
        end
    end
end
xlim([-barWidth/2 x(end)+barWidth/2])
set(gca,'XTick',mean(X,2),'XTickLabel',barGroupNames)
X = X';
% sem = data; 
% errorbar(X(:)',data(:),sem(:),'o')


