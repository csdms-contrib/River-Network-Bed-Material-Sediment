function Plot_Area_Function(Y,xL,XAxisLabel,color,A,num)

hfinc=(xL(2)-xL(1))/2;
Awfc=zeros(max(size(xL)),1);
Awf=zeros(max(size(xL)),1);
for i=1:max(size(xL))
    Awfc(i,1)=sum(A((Y)<=(hfinc+xL(i))));
end
Awf(1,1)=Awfc(1,1);
Awf(2:end,1)=diff(Awfc);
%plot(xL,Awf)

NN = sum(A);%sum(Awf);
x1 = xL./max(xL);
y1 = Awf./NN;
x2 = xL;
y2 = Awf;
%inc=xL(73)-xL(25);
my = 0.07;

if num==1
hl1 = line(x1,y1,'Color',color);
ax1 = gca;
%set(ax1,'XColor','r','YColor','r')
axis(ax1, [0 1 0 my])
%axis(ax1, [0 max(xL)+inc 0 0.04])
xlabel(ax1,'Normalized time','FontSize',14)
%xlabel(ax1,'Normalized distance','FontSize',14)
%xlabel(ax1,XAxisLabel,'FontSize',14)
ylabel(ax1,'Fraction of area','FontSize',14)

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
       
%hl2 = line(x2+inc,y2,'Color','r','Parent',ax2);
hl2 = line(x2,y2,'Color','k','Parent',ax2);
axis(ax2, [0 max(xL) 0 my*NN])
%axis(ax2, [0 max(xL)+inc 0 0.04*2427])
xlabel(ax2,XAxisLabel,'FontSize',14)
ylabel(ax2,'Area, m^2','FontSize',14)

else if num==2
        hl2 = line(x2,y2,'Color',color);
    end
    
end
end