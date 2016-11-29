function [xL,hl3] = BMS_Plot_Width_Function2(Y,n,XAxisLabel,color,num,t)

[N,xL]=hist(Y,n);

NN = max(size(Y));
%NN = 2427;
x1 = xL./max(xL);
y1 = N./NN;
x2 = xL;
y2 = N;
%inc=xL(73)-xL(25);
my = 0.04;

if num==1
%     hl1 = line(x1,y1,'Color',color);
%     ax1 = gca;
%     %set(ax1,'XColor','r','YColor','r')
%     axis(ax1, [0 1 0 my])
%     %axis(ax1, [0 max(xL)+inc 0 0.04])
%     xlabel(ax1,'Normalized time','FontSize',14)
%     %xlabel(ax1,'Normalized distance','FontSize',14)
%     %xlabel(ax1,XAxisLabel,'FontSize',14)
%     ylabel(ax1,'Fraction of links','FontSize',14)
    
%     ax2 = axes('Position',get(ax1,'Position'),...
%         'XAxisLocation','top',...
%         'YAxisLocation','right',...
%         'Color','none',...
%         'XColor','k','YColor','k');
    %set('CurrentAxes',a1)
    %hl2 = line(x2+inc,y2,'Color','r','Parent',ax2);
    %hl2 = line(x2,y2,'Color','k','Parent',ax2);
    hl2 = line(x2,y2,'Color','k');
    %axis(ax2, [0 max(xL) 0 my*NN])
    axis(gca, [0 max(xL) 0 100])
    %axis(ax2, [0 max(xL)+inc 0 0.04*2427])
%     xlabel(ax2,XAxisLabel,'FontSize',14)
%     ylabel(ax2,'Number of links','FontSize',14)
    xlabel(gca,XAxisLabel,'FontSize',12)
    ylabel(gca,'Number of parcels','FontSize',12)
    
    hl3 = line([t/365/0.175 t/365/0.175],[0 60],...
        'Color','k','LineStyle',':','Parent',gca);
    title(gca,{'Outlet response'},'FontSize',14);

else if num==2
        hl2 = line(x2,y2,'Color',color);
        %     else if num==3 %amplification
        %             hl2 = line(x2,y2,'Color','k');
        %         end
    end
    
end