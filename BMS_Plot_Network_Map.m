function [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl)

color=cat(1,[0.85 0.85 0.85],jet(8));
colormap(color)
caxis([0 90])
% cb=colorbar('YTick', linspace(0,90,10),...
%     'YTickLabel',{'0','1','2','4','8','16','32','64','128','200'},...
%     'Position',[0.5 0.2 0.01 0.6]);
%edge=edge./10000;

cb=colorbar('YTick', linspace(0,90,10),...
    'YTickLabel',{num2str(edge(1)),num2str(edge(2)),num2str(edge(3)),...
    num2str(edge(4)),num2str(edge(5)),num2str(edge(6)),num2str(edge(7)),...
    num2str(edge(8)),num2str(edge(9)),num2str(edge(10))},...
    'Position',cbpos);

%convert m to km
% cb=colorbar('YTick', linspace(0,90,10),...
%     'YTickLabel',{num2str(round(edge(1)./1000)),num2str((edge(2)./1000)),num2str(round(edge(3)./1000)),...
%     num2str(round(edge(4)./1000)),num2str(round(edge(5)./1000)),num2str(round(edge(6)./1000)),num2str(round(edge(7)./1000)),...
%     num2str(round(edge(8)./1000)),num2str(round(edge(9)./1000)),num2str(round(edge(10)./1000))},...
%     'Position',cbpos);

% cb=colorbar('Location','North','XTickMode','manual','XTick', linspace(0,90,10),...
%     'XTickLabel',{'0','1','2','4','8','16','32','64','128','200'},'Position',[0.1 0.9 0.3 0.01]);
%edge=edge.*10000;

netspec = makesymbolspec('Line',...
    {'Default', 'Color','blue'}, ...
    {'sed',[edge(1) edge(2)], 'Color',color(1,:), 'LineWidth', 0.5},...
    {'sed',[edge(2) edge(3)], 'Color',color(2,:), 'LineWidth', 1.5},...
    {'sed',[edge(3) edge(4)], 'Color',color(3,:), 'LineWidth', 1.5},...
    {'sed',[edge(4) edge(5)], 'Color',color(4,:), 'LineWidth', 1.5},...
    {'sed',[edge(5) edge(6)], 'Color',color(5,:), 'LineWidth', 1.5},...
    {'sed',[edge(6) edge(7)], 'Color',color(6,:), 'LineWidth', 1.5},...
    {'sed',[edge(7) edge(8)], 'Color',color(7,:), 'LineWidth', 1.5},...
    {'sed',[edge(8) edge(9)], 'Color',color(8,:), 'LineWidth', 1.5},...
    {'sed',[edge(9) edge(10)], 'Color',color(9,:), 'LineWidth', 1.5});

a1b = mapshow(a1,boundary, 'FaceColor', 'none', 'EdgeColor', [0.4 0.4 0.4]);

xlabel(a1,'Easting in meters','FontSize',12);
ylabel(a1,'Northing in meters','FontSize',12);
%xlim(a1,[1.5e5 5e5]);
%ylim(a1,[4.75e6 5.15e6]);
hold(a1);

a1m = mapshow(a1,network,'SymbolSpec',netspec);

title(a1,ttl,'FontSize',14);

end