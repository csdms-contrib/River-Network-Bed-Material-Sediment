%% RivNet_PlotParcels

%% Plot parcels at single time 2-D
t=10;

%create figure
f1 = figure;
set(f1,'Color','w')
set(f1,'Position',[213 50 938 632]);
a1 = axes;
%set(a1,'Color','k')
set(a1,'Color','w','XColor','w','YColor','w')
box(a1);
axis equal
hold on
%plot network
for i=1:LinkNum
    %if ~Lake(i)
        n1=plot(a1,network(i).X,network(i).Y,...
            'Color',[0.7 0.7 0.7],'LineWidth',1);
        %n1=plot(a1,network(i).X,network(i).Y,'-b');
    %end
end
%plot lakes
mapshow(lakes, 'FaceColor', 'blue', 'EdgeColor', 'blue');
% for i=1:length(lakes)
%     patch(lakes(i).X(1:end-1),lakes(i).Y(1:end-1),'b');
% end
w1=plot(a1,boundary.X,boundary.Y,...
    'Color',[0.5 0.5 0.5],'LineWidth',1);
% w2=plot(boundaryW.X,boundaryW.Y,...
%     'Color',[0.5 0.5 0.5],'LineWidth',1);
% w3=plot(boundaryL.X,boundaryL.Y,...
%     'Color',[0.5 0.5 0.5],'LineWidth',1);
xlabel(a1,'Easting in meters','FontSize',12);
ylabel(a1,'Northing in meters','FontSize',12);
% ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
%     num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)]};
% ttl={['Sand transport'];['Time since initial input = ',num2str(fix((t-1)/365/0.175*4)),' years ',...
%     num2str(round(mod((t-1)/0.175*4,365))),' days']};
ttl={['Sand transport'];['Time = ',num2str(fix(time(t))),' years ',...
   num2str(round(mod(time(t)*365,365))),' days']};
% ttl={['N transport'];['Time = ',num2str(fix(time(t))),' days ',...
%     num2str(round(mod(time(t)*24,24))),' hours']};
%title(a1,ttl,'FontSize',14,'Color',[0.5 0.5 0.5]);

%plot parcels on network at current time step    
BMS_Plot_Network_Parcels_Size2

o1=plot(a1,network(OutletLinkID).X(end-1),network(OutletLinkID).Y(end-1),...
    'LineStyle','none','Marker','^','MarkerEdgeColor',[0.8 0.8 0.8],...
    'MarkerFaceColor','r','MarkerSize',10);
t1=text(4.004207964711753e+05,4.892320388236555e+06,'Basin outlet','Color',[0.5 0.5 0.5]);
t2=text(3.211261859593419e+05,4.897663690831686e+06,ttl,'Color',[0.5 0.5 0.5]);

lp1=plot(a1,3.4e+05,4.795e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[0.8 0.8 0.8],...
    'MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',5);
%tp1=text(3.42e+05,4.795e+06,'0 - 2','Color',[0.5 0.5 0.5]);
tp1=text(3.42e+05,4.795e+06,cat(2,num2str(edge(1),1),' - ',num2str(edge(2),1)),'Color',[0.5 0.5 0.5]);

lp2=plot(a1,3.4e+05,4.8e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[146/255 39/255 143/255],...
    'MarkerFaceColor',[146/255 39/255 143/255],'MarkerSize',5);
%tp2=text(3.42e+05,4.8e+06,'2 - 4','Color',[0.5 0.5 0.5]);
tp2=text(3.42e+05,4.8e+06,cat(2,num2str(edge(2),1),' - ',num2str(edge(3),1)),'Color',[0.5 0.5 0.5]);

lp3=plot(a1,3.4e+05,4.805e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[0/255 174/255 239/255],...
    'MarkerFaceColor',[0/255 174/255 239/255],'MarkerSize',5);
%tp3=text(3.42e+05,4.805e+06,'4 - 10','Color',[0.5 0.5 0.5]);
tp3=text(3.42e+05,4.805e+06,cat(2,num2str(edge(3),1),' - ',num2str(edge(4),1)),'Color',[0.5 0.5 0.5]);

lp4=plot(a1,3.4e+05,4.81e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[141/255 198/255 63/255],...
    'MarkerFaceColor',[141/255 198/255 63/255],'MarkerSize',5);
%tp4=text(3.42e+05,4.81e+06,'10 - 20','Color',[0.5 0.5 0.5]);
tp4=text(3.42e+05,4.81e+06,cat(2,num2str(edge(4),1),' - ',num2str(edge(5),1)),'Color',[0.5 0.5 0.5]);

lp5=plot(a1,3.4e+05,4.815e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[251/255 176/255 64/255],...
    'MarkerFaceColor',[251/255 176/255 64/255],'MarkerSize',5);
%tp5=text(3.42e+05,4.815e+06,'20 - 45','Color',[0.5 0.5 0.5]);
tp5=text(3.42e+05,4.815e+06,cat(2,num2str(edge(5),1),' - ',num2str(edge(6),1)),'Color',[0.5 0.5 0.5]);

lp6=plot(a1,3.4e+05,4.82e+06,'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[237/255 28/255 36/255],...
    'MarkerFaceColor',[237/255 28/255 36/255],'MarkerSize',5);
%tp6=text(3.42e+05,4.82e+06,'45 - 100','Color',[0.5 0.5 0.5]);
tp6=text(3.42e+05,4.82e+06,cat(2,num2str(edge(6),1),' - ',num2str(edge(7),2)),'Color',[0.5 0.5 0.5]);

%ltl={['Number of hillslope contributions'];['in same cluster']};
%ltl={['Concentration of parcels in'];['each link, #/m^2 \times 1000']};
%ltl='Bed-material thickness, m';
%ltl='N conc. in link, mg/L';
%tp7=text(3.306e+05,4.8259e+06,ltl,'Color',[0.5 0.5 0.5]);

r1=rectangle('Position',[4.4e+05 4.8e+06 5000 1000],'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);
rt1=text(4.395e+05,4.797e+06,'0','Color',[0.5 0.5 0.5]);
r2=rectangle('Position',[4.45e+05 4.8e+06 5000 1000],'EdgeColor',[0.5 0.5 0.5]);
%rt2=text(4.445e+05,4.797e+06,'5','Color',[0.5 0.5 0.5]);
r3=rectangle('Position',[4.5e+05 4.8e+06 5000 1000],'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);
rt3=text(4.485e+05,4.797e+06,'10','Color',[0.5 0.5 0.5]);
r4=rectangle('Position',[4.55e+05 4.8e+06 5000 1000],'EdgeColor',[0.5 0.5 0.5]);
%rt4=text(4.535e+05,4.797e+06,'15','Color',[0.5 0.5 0.5]);
rt5=text(4.585e+05,4.797e+06,'20 km','Color',[0.5 0.5 0.5]);
%nt1=text(4.39e+05,4.783e+06,'Czuba and Foufoula-Georgiou, 2015','Color',[0.5 0.5 0.5]);

