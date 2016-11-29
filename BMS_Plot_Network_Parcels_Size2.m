%determine parcel location at current time step and plot on network
clear x y z n
n=1;
x(1:sum(numpar(t,:),1))=NaN;
y(1:sum(numpar(t,:),1))=NaN;
z(1:sum(numpar(t,:),1))=NaN;

%zval = lnkconc;
zval = numpar;
%zval = seddepth;

for i=1:LinkNum
    if isempty(P_loc{t,i})
        continue
    end
    clear J
    J=length(P_loc{t,i});
    for j=1:J
        clear idx
        idx=find(network(i).L<=P_loc{t,i}(j),1,'last');
        if P_loc{t,i}(j)==1
            x(n,1)=network(i).X(idx);
            y(n,1)=network(i).Y(idx);
        else
            x(n,1)=interp1(network(i).L(idx:idx+1),network(i).X(idx:idx+1),P_loc{t,i}(j));
            y(n,1)=interp1(network(i).L(idx:idx+1),network(i).Y(idx:idx+1),P_loc{t,i}(j));
        end
        %z(n,1)=numparconc(t,i);
        z(n,1)=zval(t,i);
        n=n+1;
    end
end
%plot parcels at next time step
%white background
% p1=plot(a1,x,y,'LineStyle','none',...
%     'Marker','o','MarkerEdgeColor','k',...
%     'MarkerFaceColor','k','MarkerSize',5);
%black background
% p1=plot(a1,x,y,'LineStyle','none',...
%     'Marker','o','MarkerEdgeColor',[0.9 0.9 0.9],...
%     'MarkerFaceColor',[0.9 0.9 0.9],'MarkerSize',5);

mn=nanmin(nanmin(zval(zval>0)));
edge = exp(linspace(log(mn),log((max(max(zval)))),7));
%edge = (linspace((mn),((max(max(zval)))),7));
%edge = [mn 0.001 0.01 0.05 0.1 0.5 1.5];
%edge = [mn 1 10 30 70 100 1000];
edge(end)=edge(end)+0.01*edge(end);

b1=find(z>=edge(1)&z<edge(2));
b2=find(z>=edge(2)&z<edge(3));
b3=find(z>=edge(3)&z<edge(4));
b4=find(z>=edge(4)&z<edge(5));
b5=find(z>=edge(5)&z<edge(6));
b6=find(z>=edge(6)&z<edge(7));

p1=plot(a1,x(b1),y(b1),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[0.8 0.8 0.8],...
    'MarkerFaceColor',[0.8 0.8 0.8],'MarkerSize',5);
p2=plot(a1,x(b2),y(b2),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[146/255 39/255 143/255],...
    'MarkerFaceColor',[146/255 39/255 143/255],'MarkerSize',5);
p3=plot(a1,x(b3),y(b3),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[0/255 174/255 239/255],...
    'MarkerFaceColor',[0/255 174/255 239/255],'MarkerSize',5);
p4=plot(a1,x(b4),y(b4),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[141/255 198/255 63/255],...
    'MarkerFaceColor',[141/255 198/255 63/255],'MarkerSize',5);
p5=plot(a1,x(b5),y(b5),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[251/255 176/255 64/255],...
    'MarkerFaceColor',[251/255 176/255 64/255],'MarkerSize',5);
p6=plot(a1,x(b6),y(b6),'LineStyle','none',...
    'Marker','o','MarkerEdgeColor',[237/255 28/255 36/255],...
    'MarkerFaceColor',[237/255 28/255 36/255],'MarkerSize',5);

clear mn b1 b2 b3 b4 b5 b6 zval