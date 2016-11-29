%% BMS_ClusterAnalysisPlotting
% Analyze and plot some of the data created by BMS_Pardist_Cluster
% Some code is legacy and may not work.

%%
pardist(pardist==0)=NaN;
clear nn n xout
bins=25;
nn(1:steps,1:bins)=NaN;
%pdmx=nanmax(nanmax(pardist));
pdmxx=nanmax(nanmax(pardist));
pdmxi=nanmin(nanmin(pardist));
%pdmxx=nanmax(nanmax(Clst_spdist));
%pdmxi=nanmin(nanmin(Clst_spdist));
for i=1:steps
   [n, xout]=hist(log(pardist(i,:)),linspace(log(pdmxi),log(pdmxx),bins));
   %[n, xout]=hist((Clst_spdist(i,:)),linspace((pdmxi),(pdmxx),100));
   nn(i,:)=n./nanmax(n);
   %nn(i,:)=n./nansum(n);
end
%n(isnan(n))=0;
%nn(nn==0)=NaN;
%%
nn=nn./repmat(nansum(nn,2),1,bins);
nn(isnan(nn))=0;
nn=cumsum(nn,2);
%%
[a, Mde]=nanmax(nn,[],2);
windowSize = 100;
MdeSM=filtfilt(ones(1,windowSize)/windowSize,1,Mde(~isnan(Mde(2:end,1))));
%%
figure
hold on
box on
plot(exp(xout),nn(1,:),'k');
%hist(Length,linspace(1,pdmx,100))
xlabel('Inter-parcel distance, km')
%ylabel('Count')
%title('Time step = 30')
%xlim([log(pdmxi),log(pdmxx)])
xlim([5, 210000])
set(gca,'XScale','log')
%%
pdmxx=nanmax(nanmax(nn));
pdmxi=nanmin(nanmin(nn));
figure
%hold on
%box on
%contourf(linspace(log(pdmxi),log(pdmxx),100),1:6400,(nn),'LineStyle','none')
%contourf(linspace(log(pdmxi),log(pdmxx),100),1:1000,(nn))
%contourf(linspace(log(pdmxi),log(pdmxx),100),1:1000,log(nn))
contourf(time,exp(xout)./1000,(nn)','LineStyle','none')
%set(gca,'YScale','log')
set(gca,'YScale','log','XScale','log')
ylim([5./1000, max(exp(xout))./1000])
%xlim([time(2), time(5753)])
%colorbar
%caxis([log(pdmxi), log(0.1)])
%plot((xout(MdeSM)),2:6400,'k','LineWidth',1)
%%
% pardist(pardist==0)=NaN;
% % pdmx1=nanmax(pardist(30,:));
% % pdmx2=nanmax(pardist(100,:));
% % pdmx3=nanmax(pardist(600,:));
% % pdmxx=pdmx2;
% % pdmx1=nanmin(pardist(30,:));
% % pdmx2=nanmin(pardist(100,:));
% % pdmx3=nanmin(pardist(600,:));
% % pdmxi=pdmx2;
% figure
% hold on
% box on
% %hist(log(pardist(400,:)),linspace(0,12,100));
% hist(log(pardist(100,:)),linspace(log(pdmxi),log(pdmxx),100));
% %hist(Length,linspace(1,pdmx,100))
% xlabel('Distance between parcels, log(meters)')
% ylabel('Count')
% title('Time step = 100')
% xlim([log(pdmxi),log(pdmxx)])
%%
% maxsize=max(Cluster,[],2);
% figure
% hold on
% box on
% plot(time, maxsize,'LineWidth',1)
% ylabel('Largest cluster size, links')
% xlabel('Time, years')
% windowSize = 100;
% maxsizeavg=filtfilt(ones(1,windowSize)/windowSize,1,maxsize(~isnan(maxsize)));
% plot(maxsizeavg,'r','LineWidth',1)
% [maxcluster, maxclustime] = max(maxsizeavg);
% [absmax, absmaxtime] = max(maxsize);
% %ylim([0 120])

%%
for i=1:steps
medclst(i,1)=prctile(Clst_spdist(i,Clst_spdist(i,:)>0),50);
sfclst(i,1)=prctile(Clst_spdist(i,Clst_spdist(i,:)>0),75);
tfclst(i,1)=prctile(Clst_spdist(i,Clst_spdist(i,:)>0),25);
ntclst(i,1)=prctile(Clst_spdist(i,Clst_spdist(i,:)>0),90);
stdclst(i,1)=std(Clst_spdist(i,Clst_spdist(i,:)>0),[],2);
varclst(i,1)=var(Clst_spdist(i,Clst_spdist(i,:)>0),0,2);
meanclst(i,1)=mean(Clst_spdist(i,Clst_spdist(i,:)>0),2);
end
%medclst2=median(Cluster(9,(Cluster(9,:)>0)),2);
cvclst=stdclst./meanclst;
%% Mean Cluster Size
windowSize = 10;
maxsize=meanclst;
meansizeavg=maxsize;
meansizeavg(~isnan(maxsize))=filtfilt(ones(1,windowSize)/windowSize,1,maxsize(~isnan(maxsize)));

figure
hold on
box on
plot(time,meanclst./1000,'k','LineWidth',1)
plot(time,meansizeavg./1000,'r')
%ylabel('Mean cluster size (Spanning Distance)')
%xlabel('Time step')
%xlim([time(2), time(5753)])
%set(gca,'XScale','log')
%%
Clst_spdist(Clst_spdist==0)=NaN;
Clst_sort=sort(Clst_spdist,2);
Clst_sort(isnan(Clst_sort))=0;
%Clst_pdf=Clst_sort./repmat(sum(Clst_sort,2),1,2427);
%Clst_cdf=cumsum(Clst_pdf,2);
%%
clear n xout histclst
bins=15;
histclst(1:steps,1:bins)=NaN;
pdmxx=nanmax(nanmax(Clst_spdist));
pdmxi=nanmin(nanmin(Clst_spdist));
for i=1:steps
    %histclst(i,:)=hist(Cluster(i,Cluster(i,:)>0),1:100);
    %[n, xout]=hist(log(pardist(i,:)),linspace(log(pdmxi),log(pdmxx),bins));
    [n, xout]=hist(log(Clst_spdist(i,:)),linspace(log(pdmxi),log(pdmxx),bins));
    histclst(i,:)=n;
end
histclst=histclst./repmat(nansum(histclst,2),1,bins);
histclst(isnan(histclst))=0;
histclst=cumsum(histclst,2);
%nn(isnan(nn))=0;
%histclst(histclst==0)=NaN;
figure
contourf(time(2:end),exp(xout)./1000,histclst(2:end,:)','LineStyle','none')
set(gca,'YScale','log','XScale','log')
%contourf(histclst,'LineStyle','none')
%pcolor(log(histclst))
xlabel('Time, years')
ylabel('Cluster size, km')
%%
% for i=1:steps
% numclst(i,1)=sum(Cluster(i,:)>0);    
% end
% figure
% hold on
% box on
% plot(time,numclst,'LineWidth',1)
% ylabel('Number of clusters')
% xlabel('Time, years')
%% Max Cluster Size
figure
hold on
box on
plot(time,maxclst_spdist./1000,'b')
%plot(time,maxsizeavg./1000,'k')
%ylabel('Max cluster size (Spanning Distance)')
%xlabel('Time step')
%xlim([time(2), time(5753)])
%ylim([1, 175])
%set(gca,'XScale','log')
%[a b]=nanmax(maxsizeavg);
%plot(time(b),a./1000,'ok')
%% Time
% time=(1:1:steps)./365./0.175.*4;
% time=(1:1:steps).*144./60./60./24;
%% Number of Clusters
for i=1:steps
numclst(i,1)=sum(Clst_spdist(i,:)>0);    
end
figure
hold on
box on
plot(time,numclst,'LineWidth',1)
%ylabel('Number of clusters (Spanning Distnace)')
%xlabel('Time, years')
set(gca,'YScale','log','XScale','log')
%xlim([time(2), time(5753)])

%%
% D pdf
% T=cat(1,Time(:,1),Time(:,2),Time(:,3),Time(:,4),Time(:,5));
% figure
% hold on
% [N(:,5),xT(:,5)]=Plot_Width_Function(Time(:,5),xT(:,1),'Time to the outlet, yrs','r',1);
% scale=[];
% Nsca=N.*repmat(scale,400,1);
% figure
% plot(xT(:,1),sum(Nsca,2)./sum(sum(Nsca)))
% hold on
% plot(xT(:,1),Nsca(:,3)./sum(sum(Nsca)),'k')
% %xlim([0 400])
% plot(xT(:,1),Nsca(:,1)./sum(sum(Nsca)),'c')
% plot(xT(:,1),Nsca(:,2)./sum(sum(Nsca)),'m')
% plot(xT(:,1),Nsca(:,4)./sum(sum(Nsca)),'g')
% %xlim([0 400])
% plot(xT(:,1),Nsca(:,5)./sum(sum(Nsca)),'r')
% %%
% plot([32 32],[1 10000],':k','LineWidth',2)
% plot([400 400],[1 10000],':k','LineWidth',2)
% plot([1200 1200],[1 10000],':k','LineWidth',2)

%%
Clst_spdist(Clst_spdist==0)=NaN;
Clst_met=Clst_parconc./Clst_spdist;
maxclst_met=nanmax(Clst_met,[],2);
%%
figure
hold on
box on
plot(maxclst_spdist)
xlabel('Time step')
ylabel('Max Cluster (sp dist)')
%ylabel('Max number of parcels in a cluster')
%ylabel('Max of cluster spanning distance \times number of parcels in cluster')
%ylabel('Max of number of parcels in cluster / cluster spanning distance')

%%
windowSize = 10;
maxsize=maxclst_spdist;
maxsizeavg=filtfilt(ones(1,windowSize)/windowSize,1,maxsize(~isnan(maxsize)));
plot(1:10:120,maxsizeavg(1:10:120),'r','LineWidth',2)
%%
% H=(0.0029).*(usarea').^(0.2943);
% U=(0.1976).*(usarea').^(0.0679);
% VelocityS=0.05./sqrt(9.81)./1.65./1.65./0.0001.*U.^2.*H.^(1/2).*Slope'.^(3/2)./.1;
% numpar(numpar>1)=1;
% parvel=numpar.*repmat(VelocityS,6400,1);
% parvel(parvel==0)=NaN;
%%
% figure
% hold on
% box on
% % for i=1:6400
% %     for j=1:2427
% %         plot(i,parvel(i,j),'.b')
% %     end
% % end
% plot(1:1:6400,nanmean(parvel,2),'r')
%% first appear
% ttime(1:steps,1)=NaN;
% tlnk=13;%travel through tlnk links, look for appear in tlnk+1
% for i=1:2427
%    if input(i)==0
%        continue
%    else
%        %move through first 3 links, appearance in 4
%        for j=1:steps
%            if isnan(Connect(i,tlnk+1))
%                continue
%            elseif isempty(Source{j,Connect(i,tlnk+1)})
%                continue            
%            elseif find(Source{j,Connect(i,tlnk+1)}(1,:)==i,1,'first')
%                ttime(i)=j;
%                break
%            end  
%        end
%    end  
% end
% figure
% hist(ttime,50)
% xlabel('Time step to travel through 13 links')
% nanmean(ttime)
%% Outlet histogram
% for i = 1:inputs
%     out(i,1) = min(find(State(:,i)==0));
% end
% figure
% hold on
% [xT]=Plot_Width_Function(outtime/365,100,'Time to the outlet, yrs','b',1);

%%
% HW=Sub(:,j);
% for i=1:length(extra)
% HW(find(GridID==extra(i)))=1;
% end
% %%
% active(1:steps,1:LinkNum)=1;
% active(cellfun(@isempty,Active))=0;
% cactive=cumsum(active,1);
% cactive(cactive>1)=1;
%% Plot Spatial Data
%set timestep to plot
t = 1;

%set n as # of parcels in each state
%n = hist(State(t,:),0:1:LinkNum);
%nttl = 'Number of parcels';
%nttl = 'Cluster Size (Spanning Distance)';
nttl = 'CPI';

%set n as # of inactive states
% n = hist(State(t,logical(Inactive(t,:))),0:1:LinkNum);
% nttl = 'Number of inactive parcels';

CPI=nansum(LCS_spdist,1);
% LCS=LCS_spdist;
% mn=nanmin(nanmin(LCS));
% %LCS(isnan(LCS_spdist))=0;
% LCS(LCS~=repmat(nanmax(LCS_spdist,[],2),1,LinkNum))=0;

% mn=nanmin(nanmin(numpar(numpar>0)));

% n_zero = n(1);
% n = n(2:end);

for i = 1:LinkNum
    %[network(i).sed] = n(i);
    %[network(i).sed] = removed(i,1);
    [network(i).sed] = CPI(t,i);
    %[network(i).sed] = CsizeComp(t,i);
    %[network(i).sed] = LCS(t,i);
    %[network(i).sed] = TimeS(i)/60/60/24/365;
    %[network(i).sed] = limit(i);
    %[network(i).sed] = Dgm(i);
    %[network(i).sed] = max(numpar(:,i),[],1);%max in link
    %[network(i).sed] = tleng(i);
    %[network(i).sed] = VelocityS(i).*10000;
    %[network(i).sed] = Sub(i,j);
    %[network(i).sed] = HW(i);
    %[network(i).sed] = Length(i);
    %[network(i).sed] = input(i,1);
end
%
%
f1 = figure;
set(f1,'Position',[213 50 938 632]);
%a1 = axes('Position',[0.08 0.1 0.4 0.8]);
a1 = axes;
cbpos = [0.78 0.22 0.01 0.6];
box(a1);

%edge = [0 1 2 4 8 16 32 64 128 200];
%edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
%edge = [0 1 2 4 6 8 10 12 16 20];
%edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(LCS))),9))));
%edge = cat(2,0,round(linspace(mn,max(max(LCS)),9)));
edge = cat(2,0,round((linspace((min([network.sed])),(max([network.sed])),9))));
%edge = cat(2,0,round((linspace((1),(max(max(numpar))),9))));

% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
 ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
     num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)];nttl};
 
[a1b, a1m] = BMS_Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);

%% Create Movie
% 
% %set timestep to plot
% t = 1;
% %nttl = 'Number of parcels';
% %nttl = 'Comp. Cluster Size';
% nttl = 'Max Cluster Size';
% 
% LCS=LCS_spdist;
% mn=nanmin(nanmin(LCS));
% %LCS(isnan(LCS_spdist))=0;
% LCS(LCS~=repmat(nanmax(LCS_spdist,[],2),1,LinkNum))=0;
% 
% for i = 1:LinkNum
%     %[network(i).sed] = numpar(t,i);
%     [network(i).sed] = LCS(t,i);
% end
% %
% %
% f1 = figure;
% set(f1,'Position',[213 50 938 632]);
% %a1 = axes('Position',[0.08 0.1 0.4 0.8]);
% a1 = axes;
% cbpos = [0.78 0.22 0.01 0.6];
% box(a1);
% 
% %edge = [0 1 2 4 8 16 32 64 128 200];
% %edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
% %edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
% %edge = cat(2,0,round(exp(linspace(log(1),log(max(max(CsizeComp))),9))));
% %edge = cat(2,0,round(exp(linspace(log(1),log(max(max(LCS))),9))));
% edge = cat(2,0,round(linspace(mn,max(max(LCS)),9)));
% 
% % ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %      num2str(round(mod(t/0.175,365))),' days'];nttl};
% % ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %      num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)];nttl};
% %sand 0.4 mm
% % ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
% %      num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)];nttl};
% %uniform 1 m/s
% ttl={['Time = ',num2str(fix(t.*144./60./60./24)),' days ',...
%      num2str(round(mod(t.*144./60./60,24))),' hours ',...
%      num2str(round(mod(t.*144./60,60))),' minutes'];['Time step = ',num2str(t)];nttl}; 
%  
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);
% 
% writerObj = VideoWriter('ex6','MPEG-4');
% writerObj.FrameRate = 40;
% writerObj.Quality =100;
% open(writerObj);
% 
% set(gcf,'Renderer','Painters');
% 
% frame = getframe(gcf);
% writeVideo(writerObj,frame);
% cla(a1);
% colorbar('delete');
% 
% for t=2:1:6000
%     for i = 1:LinkNum
%         %[network(i).sed] = numpar(t,i);
%         %[network(i).sed] = CsizeComp(t,i);
%         [network(i).sed] = LCS(t,i);
%     end
%     
% %     ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %      num2str(round(mod(t/0.175,365))),' days'];nttl};
% %     ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %         num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)];nttl};
% %     ttl={['Time = ',num2str(fix(t/365/0.175*4)),' years ',...
% %         num2str(round(mod(t/0.175*4,365))),' days'];['Time step = ',num2str(t)];nttl};
%     %uniform 1 m/s
%     ttl={['Time = ',num2str(fix(t.*144./60./60./24)),' days ',...
%         num2str(round(mod(t.*144./60./60,24))),' hours ',...
%         num2str(round(mod(t.*144./60,60))),' minutes'];['Time step = ',num2str(t)];nttl};
%     
%     [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);
% 
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     cla(a1);
%     colorbar('delete');
%     
% end
% 
% close(writerObj);
% %% Create Cluster and Cluster Complement 2 Pane Movie
% 
% %set timestep to plot
% t = 1;
% 
% f1=figure;
% set(f1,'Position',[160 60 1024 620]);
% a1=axes('Position',[0.065 0.1 0.4 0.8]);
% a2=axes('Position',[0.543 0.1 0.4 0.8]);
% cbpos1 = [0.474 0.2 0.01 0.6];
% cbpos2 = [0.954 0.2 0.01 0.6];
% anpos = [0.4 0.896 0.3 0.1];
% box(a1);
% box(a2);
% 
% nttl1 = 'Cluster Size';
% for i = 1:LinkNum
%     [network(i).sed] = Csize(t,i);
% end
% edge = cat(2,0,round(exp(linspace(log(1),log(max(max(Csize))),9))));
% 
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos1,boundary,network,nttl1);
% 
% nttl2 = 'Comp. Cluster Size';
% for i = 1:LinkNum
%     [network(i).sed] = CsizeComp(t,i);
% end
% edge = cat(2,0,round(exp(linspace(log(1),log(max(max(CsizeComp))),9))));
% 
% [a2b, a2m] = Plot_Network_Map(a2,edge,cbpos2,boundary,network,nttl2);
% ylabel(a2,'')
% 
% % an1 = annotation(f1,'textbox',...
% %     anpos,...
% %     'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %      num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)]},...
% %     'FitBoxToText','off',...
% %     'LineStyle','none',...
% %     'FontSize',16);
% an1 = annotation(f1,'textbox',...
%     anpos,...
%     'String',{['Time step = ',num2str(t)]},...
%     'FitBoxToText','off',...
%     'LineStyle','none',...
%     'FontSize',16);
% 
% writerObj = VideoWriter('ex7','MPEG-4');
% writerObj.FrameRate = 30;
% writerObj.Quality =100;
% open(writerObj);
% 
% set(gcf,'Renderer','Painters');
% 
% frame = getframe(gcf);
% writeVideo(writerObj,frame);
% cla(a1);
% cla(a2);
% set(an1,'Visible','off')
% colorbar('delete');
% 
% for t=2:1:1000
%     nttl1 = 'Cluster Size';
%     for i = 1:LinkNum
%         [network(i).sed] = Csize(t,i);
%     end
%     edge = cat(2,0,round(exp(linspace(log(1),log(max(max(Csize))),9))));
%     
%     [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos1,boundary,network,nttl1);
%     
%     nttl2 = 'Comp. Cluster Size';
%     for i = 1:LinkNum
%         [network(i).sed] = CsizeComp(t,i);
%     end
%     edge = cat(2,0,round(exp(linspace(log(1),log(max(max(CsizeComp))),9))));
%     
%     [a2b, a2m] = Plot_Network_Map(a2,edge,cbpos2,boundary,network,nttl2);
%     ylabel(a2,'')
%     
% %     an1 = annotation(f1,'textbox',...
% %         anpos,...
% %         'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
% %         num2str(round(mod(t/0.175,365))),' days'];['Time step = ',num2str(t)]},...
% %         'FitBoxToText','off',...
% %         'LineStyle','none',...
% %         'FontSize',16);
%     an1 = annotation(f1,'textbox',...
%         anpos,...
%         'String',{['Time step = ',num2str(t)]},...
%         'FitBoxToText','off',...
%         'LineStyle','none',...
%         'FontSize',16);
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     cla(a1);
%     cla(a2);
%     set(an1,'Visible','off')
%     colorbar('delete');
%     
% end
% 
% close(writerObj);
% %% old movie
% %aviobj = avifile('example1.avi','compression','None','fps',5);
% writerObj = VideoWriter('ex2.avi');
% writerObj.FrameRate = 20;
% open(writerObj);
% 
% figure
% set(gcf,'Position',[213 50 938 632])
% set(gcf,'Renderer','Painters');
% box on
% 
% %edge=[0 1 10 25 45 70 90];
% edge=[0 1 2 4 8 16 32 64 128 200];
% 
% %color=colormap(jet(5));
% color=cat(1,[0.85 0.85 0.85],jet(8));
% colormap(color)
% caxis([0 90])
% colorbar('YTick', linspace(0,90,10),...
%     'YTickLabel',{'0','1','2','4','8','16','32','64','128','200'});
% 
% %rainbowscale
% % color = [0 0 210;...%blue
% %     0 210 210;...%light blue
% %     0 200 0;...%green
% %     255 215 0;...%yellow
% %     220 0 0];%red
% % color=color./255;
% %bluescale
% % color = [0.7 0.7 1.0
% %     0.3 0.3 1.0;...
% %     0.0 0.0 1.0;...
% %     0.0 0.0 0.6;...
% %     0.0 0.0 0.2];
% %grayscale
% % color = [0.8 0.8 0.8
% %     0.6 0.6 0.6;...
% %     0.4 0.4 0.4;...
% %     0.2 0.2 0.2;...
% %     0.0 0.0 0.0];
% 
% netspec = makesymbolspec('Line',...
%     {'Default', 'Color','blue'}, ...
%     {'sed',[edge(1) edge(2)], 'Color',color(1,:), 'LineWidth', 1},...
%     {'sed',[edge(2) edge(3)], 'Color',color(2,:), 'LineWidth', 2},...
%     {'sed',[edge(3) edge(4)], 'Color',color(3,:), 'LineWidth', 2},...
%     {'sed',[edge(4) edge(5)], 'Color',color(4,:), 'LineWidth', 2},...
%     {'sed',[edge(5) edge(6)], 'Color',color(5,:), 'LineWidth', 2},...
%     {'sed',[edge(6) edge(7)], 'Color',color(6,:), 'LineWidth', 2},...
%     {'sed',[edge(7) edge(8)], 'Color',color(7,:), 'LineWidth', 2},...
%     {'sed',[edge(8) edge(9)], 'Color',color(8,:), 'LineWidth', 2},...
%     {'sed',[edge(9) edge(10)], 'Color',color(9,:), 'LineWidth', 2});
% 
% hb = mapshow(boundary, 'FaceColor', 'none', 'EdgeColor', [0.4 0.4 0.4]);
% xlabel('Easting in meters')
% ylabel('Northing in meters')
% xlim([1.5e5 5e5])
% ylim([4.75e6 5.15e6])
% hold on
% % plot([2e5 2e5],[4.8e6 4.8e6],'Color',[0.85 0.85 0.85], 'LineWidth', 1)
% % plot([2e5 2e5],[4.8e6 4.8e6], 'Color',color(1,:), 'LineWidth', 2)
% % plot([2e5 2e5],[4.8e6 4.8e6], 'Color',color(2,:), 'LineWidth', 2)
% % plot([2e5 2e5],[4.8e6 4.8e6], 'Color',color(3,:), 'LineWidth', 2)
% % plot([2e5 2e5],[4.8e6 4.8e6], 'Color',color(4,:), 'LineWidth', 2)
% % plot([2e5 2e5],[4.8e6 4.8e6], 'Color',color(5,:), 'LineWidth', 2)
% % legend([num2str(edge(1)),'-',num2str(edge(2))],...
% %     [num2str(edge(2)),'-',num2str(edge(3))],...
% %     [num2str(edge(3)),'-',num2str(edge(4))],...
% %     [num2str(edge(4)),'-',num2str(edge(5))],...
% %     [num2str(edge(5)),'-',num2str(edge(6))],...
% %     [num2str(edge(6)),'-',num2str(edge(7))])
% 
% for t=1:10:6400
%     display = t;
%     n = hist(State(display,:),0:1:LinkNum);
%     n_zero = n(1);
%     n = n(2:end);
%     
%     for i = 1:LinkNum
%         [network(i).sed] = n(i);
%     end
%     
%     % figure
%     % box on
%     hn = mapshow(network,'SymbolSpec',netspec);
%     
%     title({'Number of parcels';...
%         ['t = ',num2str(fix(display/365/0.175)),' years ',...
%         num2str(round(mod(display/0.175,365))),' days']})
%     
%     %     frame = getframe(gcf);
%     %     aviobj = addframe(aviobj,frame);
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     delete(hn);
%     %delete(hb);
% end
% 
% close(writerObj);
% %aviobj = close(aviobj);
% 
% %% Multi Axis Figure 3 part
% %set timestep to plot
% t=500;
% 
% f1=figure;
% set(f1,'Position',[160 60 1024 620]);
% a1=axes('Position',[0.08 0.1 0.4 0.8]);
% a2=axes('Position',[0.62 0.07 0.3 0.5]);
% a3=axes('Position',[0.6 0.685 0.35 0.25]);
% cbpos = [0.5 0.2 0.01 0.6];
% anpos = [0.31 0.886 0.3 0.1];
% box(a1);
% box(a2);
% box(a3);
% 
% edge=[0 1 2 4 8 16 32 64 128 200];
% %edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
% 
% n = hist(State(t,:),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,'Number of parcels');
% 
% n = hist(State(t,logical(Inactive(t,:))),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% [a2b, a2m] = Plot_Network_Map(a2,edge,cbpos,boundary,network,'      Number of inactive parcels');
% 
% % Outlet histogram
% for i = 1:inputs
%     out(i,1) = min(find(State(:,i)==0));
% end
% 
% [xT]=Plot_Width_Function2(out/365/0.175,100,'Time to the outlet, yrs','b',1,t);
% 
% an1 = annotation(f1,'textbox',...
%     anpos,...
%     'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days']},...
%     'FitBoxToText','off',...
%     'LineStyle','none',...
%     'FontSize',16);
% 
% %% Multi Axis Figure 2 part
% %set timestep to plot
% t = 500;
% 
% f1=figure;
% set(f1,'Position',[160 60 1024 620]);
% a1=axes('Position',[0.065 0.1 0.4 0.8]);
% a2=axes('Position',[0.58 0.1 0.4 0.8]);
% cbpos = [0.485 0.2 0.01 0.6];
% anpos = [0.4 0.87 0.3 0.1];
% box(a1);
% box(a2);
% 
% edge=[0 1 2 4 8 16 32 64 128 200];
% %edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
% 
% n = hist(State(t,:),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,'Number of parcels');
% 
% n = hist(State(t,logical(Inactive(t,:))),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% [a2b, a2m] = Plot_Network_Map(a2,edge,cbpos,boundary,network,'Number of inactive parcels');
% 
% an1 = annotation(f1,'textbox',...
%     anpos,...
%     'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days']},...
%     'FitBoxToText','off',...
%     'LineStyle','none',...
%     'FontSize',16);
% 
% %% Multi Axis Figure 2 part outlet
% %set timestep to plot
% t = 100;
% 
% %set n as # of parcels in each state
% n = hist(State(t,:),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% % Outlet histogram
% for i = 1:inputs
%     out(i,1) = min(find(State(:,i)==0));
% end
% 
% f1 = figure;
% set(f1,'Position',[160 60 1024 620]);
% a1 = axes('Position',[0.08 0.1 0.4 0.8]);
% a3 = axes('Position',[0.6 0.3 0.35 0.4]);
% cbpos = [0.5 0.2 0.01 0.6];
% anpos = [0.4 0.87 0.3 0.1];
% box(a1);
% box(a3);
% 
% edge = [0 1 2 4 8 16 32 64 128 200];
% %edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
% 
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,'Number of parcels');
% 
% [xT,h2l] = Plot_Width_Function2(out/365/0.175,100,'Time, years','b',1,t);
% 
% an1 = annotation(f1,'textbox',...
%     anpos,...
%     'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days']},...
%     'FitBoxToText','off',...
%     'LineStyle','none',...
%     'FontSize',16);
%%
% figure;
% axes('XScale','log','YScale','log')
% hold on
% box on
% plot(maxpar,'k')
% plot(numclst,'b')
% legend('Largest number of parcels in a cluster','Number of clusters')
% %%
% %set timestep to plot
% t=1;
% 
% f1=figure;
% %set(f1,'Position',[160 60 1024 620]);
% set(f1,'Units','inches','Position',[1 1 8.5 6]);
% 
% a1=axes('Position',[0.08 0.1 0.4 0.8]);
% a2=axes('Position',[0.63 0.1 0.34 0.45]);
% a3=axes('Position',[0.63 0.65 0.34 0.25]);
% cbpos = [0.5 0.2 0.01 0.6];
% anpos = [0.3 0.882 0.4 0.1];
% box(a1);
% box(a2);
% box(a3);
% 
% edge=[0 1 2 4 8 16 32 64 128 200];
% %edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
% 
% n = hist(State(t,:),0:1:LinkNum);
% n_zero = n(1);
% n = n(2:end);
% for i = 1:LinkNum
%     [network(i).sed] = n(i);
% end
% 
% [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,'Number of parcels');
% 
% % n = hist(State(t,logical(Inactive(t,:))),0:1:LinkNum);
% % n_zero = n(1);
% % n = n(2:end);
% % for i = 1:LinkNum
% %     [network(i).sed] = n(i);
% % end
% % 
% % [a2b, a2m] = Plot_Network_Map(a2,edge,cbpos,boundary,network,'      Number of inactive parcels');
% 
% axes(a2)
% set(a2,'XScale','log','YScale','log');
% hold(a2);
% plot((1:steps)./365./0.175, maxpar,'k')
% %hold(a2);
% plot(a2,(1:steps)./365./0.175, numclst,'b')
% h2l = line([t/365/0.175 t/365/0.175],[1 10000],...
%         'Color','k','LineStyle',':');
% ylabel({'Largest number of parcels in a cluster (black)';'Number of clusters (blue)'},'FontSize',12)
% xlabel('Time, yrs','FontSize',12)
% hold off
% 
% % Outlet histogram
% for i = 1:inputs
%     out(i,1) = min(find(State(:,i)==0));
% end
% 
% axes(a3)
% [xT,h3l]=Plot_Width_Function2(out/365/0.175,100,'Time to the outlet, yrs','b',1,t);
% 
% an1 = annotation(f1,'textbox',...
%     anpos,...
%     'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days']},...
%     'FitBoxToText','off',...
%     'LineStyle','none',...
%     'FontSize',16);
% %
% writerObj = VideoWriter('ex5','MPEG-4');
% writerObj.FrameRate = 20;
% writerObj.Quality =100;
% open(writerObj);
% 
% set(gcf,'Renderer','Painters');
% 
% frame = getframe(gcf);
% writeVideo(writerObj,frame);
% delete(h2l,h3l);
% cla(a1);
% colorbar('delete');
% 
% for t=2:1:1600
%     n = hist(State(t,:),0:1:LinkNum);
%     n_zero = n(1);
%     n = n(2:end);
%     for i = 1:LinkNum
%         [network(i).sed] = n(i);
%     end
%     
%     [a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,'Number of parcels');
%     
%     axes(a2);
%     hold(a2);
%     h2l = line([t/365/0.175 t/365/0.175],[1 10000],...
%         'Color','k','LineStyle',':');
%     hold off
%     
%     axes(a3);
%     hold(a3);
%     h3l = line([t/365/0.175 t/365/0.175],[0 120],...
%         'Color','k','LineStyle',':');
%     hold off
%     
%     set(an1,'String',{['Time = ',num2str(fix(t/365/0.175)),' years ',...
%         num2str(round(mod(t/0.175,365))),' days']});
%  
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     delete(h2l,h3l);
%     cla(a1);
%     colorbar('delete');
%     
% end
% 
% close(writerObj);
% %%
% set(gcf,'PaperPositionMode','auto');
% set(gcf,'PaperSize',[8.5 6]);
% print('-dpdf',sprintf('image%05.5d.pdf', 5));
% 
% %% Determine Node Locations
% % X1(1:LinkNum,1)=NaN;
% % Y1(1:LinkNum,1)=NaN;
% % Node1(1:LinkNum,1)=NaN;
% % Elev1(1:LinkNum,1)=NaN;
% % e1(1:LinkNum,1)=NaN;
% % X2(1:LinkNum,1)=NaN;
% % Y2(1:LinkNum,1)=NaN;
% % Node2(1:LinkNum,1)=NaN;
% % Elev2(1:LinkNum,1)=NaN;
% % e2(1:LinkNum,1)=NaN;
% % 
% % X(1:LinkNum,1)=NaN;
% % Y(1:LinkNum,1)=NaN;
% % Node(1:LinkNum,1)=NaN;
% % Elev(1:LinkNum,1)=NaN;
% % for i=1:2427
% %    X1(i,1)=network(i,1).X(1,1);
% %    Y1(i,1)=network(i,1).Y(1,1);
% %    Node1(i,1)=network(i,1).from_node(1,1);
% %    X2(i,1)=network(i,1).X(1,end-1);
% %    Y2(i,1)=network(i,1).Y(1,end-1);
% %    Node2(i,1)=network(i,1).to_node(1,1); 
% %    
% %    e1(i,1)=network(i,1).DL10_C10_1;
% %    e2(i,1)=network(i,1).DL10_C10_2;
% %    
% %    Elev1(i,1)=max(e1(i,1),e2(i,1));
% %    Elev2(i,1)=min(e1(i,1),e2(i,1));
% % end
% % Node=cat(1,Node1,Node2(903,1));
% % X=cat(1,X1,X2(903,1));
% % Y=cat(1,Y1,Y2(903,1));
% % Elev=cat(1,Elev1,Elev2(903,1));