%% Interarrival
% This script computes the statistics on the interarrival times for
% parcel first entering a link.

% Jon Czuba
% February 25, 2015

%% Compute width function for every link
% for every row of Connect i (every link in the basin),
% find which links connect through that point and mark these
% links j in the column of Sub
% diagonal = 1

% Sub(i,j)
% *each row i indicates which indicies that node connects
% through (including itself) to reach the outlet
% *each column j corresponds to all indicies upstream of
% and including j
Sub(1:LinkNum,1:LinkNum)=0;
for j=1:LinkNum
    for i=1:LinkNum
        if ~isempty(find(Connect(i,:)==j,1))
            Sub(i,j)=1;
        end
    end
end
%% Remove Lakes from Sub
% remove from Sub the contributions upstream of lakes for locations
% downstream of lakes
for i=1:LinkNum
    if Lake(i)~=0 %for lake links
        for j=1:LinkNum
            if Sub(i,j) 
                if j~=i
                Sub(logical(Sub(:,i)),j)=0;
                end
            end
        end
    end
end
%% compute PSWF
%same as PSWFi{OutletLinkID,1}
PSWFout(1:LinkNum,1)=NaN;
WFout(1:LinkNum,1)=NaN;
%Remout(1:LinkNum,1)=NaN;
for i=1:LinkNum
    lastConn=find(Connect(i,:)==OutletLinkID);
    PSWFout(i,1)=sum(STTime(Connect(i,1:lastConn),1))./60./60./24./365;%years to outlet 
    WFout(i,1)=sum(Length(Connect(i,1:lastConn),1))./1000;%km to outlet
    %Remout(i,1)=sum(Rem(Connect(i,1:lastConn),1));%kg removed from input to outlet
end
clear lastConn

%% compute max PSWF for each link
PSWFi=cell(LinkNum,1);%years, PSWF for each link
%WFi=cell(LinkNum,1);%km, WF for each link
%WFidx=cell(LinkNum,1);%km, WF for each link
for i=1:LinkNum
    uslnks=find(Sub(:,i)==1);
    for ii=1:length(uslnks)    
        lastConn=find(Connect(uslnks(ii),:)==i);
        trt=sum(STTime(Connect(uslnks(ii),1:lastConn),1))./60./60./24./365;%years to outlet
        PSWFi{i,1}=cat(2,PSWFi{i,1},trt);%years
        %dst=sum(Length(Connect(uslnks(ii),1:lastConn),1))./1000;%km to outlet
        %WFi{i,1}=cat(2,WFi{i,1},dst);%km
        %WFidx{i,1}=cat(2,WFidx{i,1},uslnks(ii));%index corresponding to distance
    end 
end

PSWFmx=cellfun(@max,PSWFi);
clear uslnks ii lastConn trt dst

%hist(PSWFi{OutletLinkID,1},50)

%% Compute Statistics of Interarrival Times (IT)
IT_n(1:LinkNum,1)=NaN; %number of interarrival times (IT)
IT_mu(1:LinkNum,1)=NaN; %mean of IT
IT_std(1:LinkNum,1)=NaN; %std of IT
IT_cov(1:LinkNum,1)=NaN; %cov of IT

for i = 1:LinkNum
    art=L_arrival{i,1}; %yrs, arrival time
    art=art(art>(PSWFmx(i,1)+mu));%use data greater than PSWFi+mu
    sart=sort(art); %sorted arrival time
    iart=diff(sart); %interarrival time
    
    IT_n(i,1)=length(iart); %number of interarrival times (IT)
    IT_mu(i,1)=mean(iart); %mean of IT
    IT_std(i,1)=std(iart); %std of IT
    IT_cov(i,1)=IT_std(i,1)/IT_mu(i,1); %cov of IT
    
    clear art sart iart
end

%% Number of upstream links
US_links=sum(Sub,1)';%number of upstream links including itself

%% Plot IT statistics
figure; hold on; box on
%plot(US_links,IT_mu,'^r')
plot(VpLi_tni,IT_mu,'^r')
plot([1 10^3],[1 10^-3],'-k')
set(gca,'YScale','log','XScale','log')
ylabel('Mean of interarrival times')
xlabel('Number of contributing links')

figure; hold on; box on
%plot(US_links,IT_std,'^r')
plot(VpLi_tni,IT_std,'^r')
plot([1 10^3],[10 10^-2],'-k')
set(gca,'YScale','log','XScale','log')
ylabel('Standard deviation of interarrival times')
xlabel('Number of contributing links')

figure; hold on; box on
%plot(US_links,IT_cov,'ob')
plot(VpLi_tni,IT_cov,'ob')
plot([1 10^3],[1 1],'-k')
set(gca,'XScale','log')
ylim([0 2])
ylabel('COV of interarrival times')
xlabel('Number of contributing links')

%% compute arrival times at outlet for all inputs
%likely double counting because of how L_arrival assigned for inputs
% OutArrival_inputs=[];
% for i=1:LinkNum
%     OutArrival_inputs=cat(1,OutArrival_inputs,L_arrival{i,:}'+PSWF(i,1));%years
% end

%% Compute Statistics of Sediment Depth (SD)

SD_n(1:LinkNum,1)=NaN; %number of values of sed depth (SD)
SD_m(1:LinkNum,1)=NaN; %mean of SD
SD_std(1:LinkNum,1)=NaN; %std of SD
SD_cov(1:LinkNum,1)=NaN; %cov of SD

for i = 1:LinkNum
    art=seddepth(time>(PSWFmx(i,1)+1),i);%use data greater than PSWFi+mu
    %art=seddepth(time>(100),i);%use data greater than PSWFi+mu
    
    SD_n(i,1)=length(art); %number of interarrival times (SD)
    SD_m(i,1)=mean(art); %mean of SD
    SD_std(i,1)=std(art); %std of SD
    SD_cov(i,1)=SD_std(i,1)/SD_m(i,1); %cov of SD
    
    clear art
end
%% Plot Statistics of SD
figure; hold on; box on
%plot(rho.*Vp,SD_m,'.r')
%plot(rho.*Vp,hs,'ok')
%plot(rho.*Vp./B./Length,SD_std.^2,'.b')
plot(Slope,SD_m,'.b')
%plot([10^-4 10^4],[10^-4 10^4],'-k')
set(gca,'YScale','log','XScale','log')
ylabel('Mean of sediment depth')
%xlabel('Variance of sediment depth')
%xlim([10^-4 10^4])
%%
herr=(SD_m-hs)./hs.*100;
herr(logical(Lake))=NaN;
%%
figure
plot(VpLi_gen./VpLi_arr,herr,'.b')
set(gca,'XScale','log')
%%
% CapFact=US_links.*(1/10).*1./(Velocity.*60.*60.*24.*365.*B.*theta.*H);
% CapRate=Velocity.*60.*60.*24.*365.*B.*theta.*H./US_links;%m3/yr
%% Analytical Mean SD

Hs=theta.*H;%sediment depth at capacityhs(1:LinkNum,1)=NaN;
% for i=1:LinkNum
%     hs(i,1)=US_links(i)/mu/365/24/60/60/Velocity(i)/B(i);%mean bed thickness
% end

hs=VpLi_arr./Velocity./B./365./24./60./60;


RelCap=hs./Hs;%sediment depth relative to capacity
%% Numerical SD timeseries
%seddepth=numpar.*1000./repmat(Length',timesteps,1)./repmat(B',timesteps,1);%m

%% Plot numerical SD timeseries
i=630%284
figure
hold on
box on
plot(time,seddepth(:,i))
plot([PSWFmx(i,1)+mu PSWFmx(i,1)+mu], [0 0.1],'k')
%ylim([0 1])
xlabel('Time, years')
ylabel('Bed sediment thickness, m')
%% Plot comparison of numerical and analytical mean SD
figure
plot(hs,Vp.*rho./B./Length,'.b')
%plot(SD_std.^2,Vp_invar,'.b')
set(gca,'YScale','log','XScale','log')
xlabel('Mean of sediment depth')
ylabel('hs')
%% Compute Analytical PDF of SD
STTimeyr=STTime./365./24./60./60;
%rho=STTimeyr.*US_links./mu;
rho=STTimeyr.*VpLi_tni./1;

clear ph cph n
n=(0:1:10000)';
ph(1:length(n),1:LinkNum)=NaN;
for i=1:LinkNum
    %ph(:,i)=exp(-rho(i)).*rho(i).^(n)./factorial(n);
    ph(:,i)=poisspdf(n,rho(i));
    cph(:,i)=poisscdf(n,rho(i));
end

Vp=10;%m3
%hsn=repmat(n,1,LinkNum).*Vp./repmat(Length',length(n),1)./repmat(B',length(n),1);
%Vp_bar=VpLi_arr./VpLi_tni;%m3
%Vp_bar=cellfun(@geomean,VpLi_inp);
%hsn=repmat(n,1,LinkNum).*repmat(Vp_bar',length(n),1)./repmat(Length',length(n),1)./repmat(B',length(n),1);
hsn=repmat(n,1,LinkNum).*Vp./repmat(Length',length(n),1)./repmat(B',length(n),1);

%% Plot histogram of bed sediment thickness
%hval=hist(numpar(:,1147),0:1:50);
i=630%159%284
[hval, cent]=hist(seddepth(time>200,i),15);
%[hval, cent]=hist(seddepth(time>(PSWFmx(i,1)+mu),i),15);
%[hval, cent]=hist(lnkvol(time>(PSWFmx(i,1)+mu),i)./Vp,15);
%[hval, cent]=hist(seddepth(time>100,i),15);

abw=diff(hsn(1:2,i));
nbw=diff(cent(1:2));
bwsf=nbw/abw;

figure
hold on
box on
bar(cent,hval./sum(hval))
plot(hsn(:,i)+Elevch(end-1,i),ph(:,i).*bwsf,'k','LineWidth',3)
%plot(n,ph(:,i).*bwsf,'k','LineWidth',3)
%xlim([0 1])
xlabel('Bed sediment thickness, m')
ylabel('Probability')
%%
Hs=theta.*H;%sediment depth at capacity

%% Plot analytical SD pdf
i=630;%513;
figure
plot(hsn(:,i),ph(:,i))
Hs(i)
%%
figure
plot(hsn(:,i),1-cph(:,i))
%xlim([0 0.3])
%%
% %looking at exp pdf for erosion rate development
% mu=1.46380;
% X=0.1:0.1:5;
% Y=pdf('exp',X.*mu,mu);
% 
% figure
% plot(X.*mu,Y)

%% Determine fraction of time capacity is exceeded
pHs(1:LinkNum,1)=NaN;
for i=1:LinkNum
    %if isnan(Vp_bar(i,1))
    %    pHs(i,1)=NaN;
    %else
        pHs(i,1)=interp1(hsn(:,i),cph(:,i),Hs(i));
    %end
end
exceedHs=1-pHs;

%% Moved to Analytical_bed_sed_pdf.m
% %% VpLi generated, arrival, capacity
% VpLi_gen(1:LinkNum,1)=0;%vol rate generated in each link
% VpLi_arr(1:LinkNum,1)=0;%vol rate arriving to and generated in each link
% VpLi_cap(1:LinkNum,1)=0;%vol rate capacity of each link
% VpLi_ni(1:LinkNum,1)=0;%number of effective inputs in each link
% VpLi_tni(1:LinkNum,1)=0;%number of effective inputs arriving and generated in each link
% %VpLi_inp=cell(LinkNum,1);
% 
% for i=1:LinkNum
% %     VpLi_gen(i,1)=sum(Rav_In_Vol(find(Rav_Link==i),1))+...
% %         sum(Up_In_Vol(find(Up_Link==i),1))+...
% %         sum(Blf_In_Vol(find(Blf_Link==i),1));
% %     VpLi_inp{i,1}=cat(2,(Rav_In_Vol(find(Rav_Link==i),1))',...
% %         (Up_In_Vol(find(Up_Link==i),1))',...
% %         (Blf_In_Vol(find(Blf_Link==i),1))');
% %     VpLi_ni(i,1)=sum(Rav_Link==i)+...
% %         sum(Up_Link==i)+...
% %         sum(Blf_Link==i);
%     VpLi_gen(i,1)=sum(In_Vol(find(In_Link==i),1));
%     VpLi_ni(i,1)=sum(In_Link==i);
% end
% 
% VpLi_cap=Velocity.*60.*60.*24.*365.*B.*theta.*H;%m3/yr
% 
% VpLi_arr=sum(Sub.*repmat(VpLi_gen,1,LinkNum),1)';%volume generated in upstream links including itself
% 
% VpLi_tni=sum(Sub.*repmat(VpLi_ni,1,LinkNum),1)';%number of discrete input locations in upstream links including itself
% 
% figure
% plot(VpLi_arr./VpLi_cap)
% 
% %%
% Elevch=Elev-repmat(Elev(1,:),timesteps,1);
% %%
% figure
% plot(time,Elevch(:,230))
% %%
% %compute/update slope
% sslope(1:timesteps,1:LinkNum)=NaN;
% for t=1:timesteps
%     for i=1:LinkNum
%         if i==OutletLinkID
%             sslope(t,i)=(Elev(t,i)-mnelev(i,1))./Length(i,1);
%         else
%             sslope(t,i)=(Elev(t,i)-Elev(t,Connect(i,2)))./Length(i,1);
%         end
%     end
% end
% Slope(Slope<1e-4)=1e-4;
% %%
% sslopech=sslope-repmat(sslope(1,:),timesteps,1);
% %% slope at capacity
% Sc=((Vp.*VpLi_tni./mu./365./24./60./60.*sqrt(9.81).*1.65.*1.65.*D)./...
%     (0.05.*0.175.*B.*H.^(3/2).*U.^(2))).^(2/3);
% 
% %%
% i=1296
% figure
% hold on
% plot(time,sslope(:,i))
% plot([0 tmax],[Sc(i) Sc(i)],'k')
% %%
% figure
% hold on
% plot((Vp.*VpLi_tni./mu)./(capacity./STTimeyr))
% %%
% vvelocity=0.05./sqrt(9.81)./1.65./1.65./D.*repmat(U',timesteps,1).^2.*...
%     repmat(H',timesteps,1).^(1/2).*repmat(Sc',timesteps,1).^(3/2)./theta.*0.175;%m/s, realtime
% %vvelocity=0.05./sqrt(9.81)./1.65./1.65./D.*repmat(U',timesteps,1).^2.*...
% %    repmat(H',timesteps,1).^(1/2).*sslope.^(3/2)./theta.*0.175;%m/s, realtime
%         
% ccapacityrate=vvelocity.*60.*60.*24.*365.*repmat(B',timesteps,1).*...
%     theta.*repmat(H',timesteps,1);%m3/yr
% %%
% figure
% plot(VpLi_arr./ccapacityrate(end-1,:)')
% %%
% i=717;
% figure
% hold on
% plot(time,ccapacityrate(:,i))
% plot([0 tmax],[VpLi_arr(i) VpLi_arr(i)],'k')
% %%
% STTime=Length./vvelocity(end-1,:)';%seconds, travel time through each link
% STTimeyr=STTime./365./24./60./60;

%% Display Results
%% Correct Elev for initial trial runs
% Elevcorr=Elev;
% 
% for i=1:LinkNum
%     idxcorr=find(Elev(:,i)~=mxelevmod(i,1));
%     if ~isempty(idxcorr)
%         Elevcorr(idxcorr,i)=Elev(idxcorr,i)-mxelev(i,1)+mxelevmod(i,1);
%     end
% end
% clear idxcorr
%% Plot long profile
i=446;
t=1;
figure
%hold on
plot(WFout(Connect(i,~isnan(Connect(i,:)),1),1),...
    capacity(Connect(i,~isnan(Connect(i,:))),1),'k','LineWidth',2)

%%
figure
plot(Elev(5999,:)-Elev(1,:));%,Length,'.b')

%%
seddepthd=seddepth-seddepth_base;
%% Plot x-axis time
i=446;%446;
ct=11;%2,6,7,11
clear y ys
y=seddepthd(:,Connect(i,ct));
y(y<0)=0;
ys=y;

windowSize=50;
ys=filter(ones(1,windowSize)/windowSize,1,y);

figure
hold on
box on
plot(time, 0.5*Ppulsep(:,Connect(i,ct)),'r')
plot(time, ys, 'b')
plot(time, overcap(:,Connect(i,ct)),'k','LineWidth',2)
%plot(time, ys, 'k')
%plot(time, repmat(Hs(Connect(i,ct),1),timesteps,1),'k')
xlabel('Time, years')
ylabel('Depth of bed sediment, meters')
xlim([0 tmax])
ylim([0 0.5])
title(num2str(Connect(i,ct)))
%%
figure
hold on
plot(time, seddepth(:,Connect(i,ct))./repmat(Hs(Connect(i,ct),1),timesteps,1))
%%
trate=capacity./STTime*60*60*24*365./B;

%%
Hs=theta.*H;%sediment depth at capacity
relcap(timesteps,LinkNum)=NaN;
for i=1:LinkNum
relcap(time>(PSWFmx(i,1)+mu),i)=seddepth(time>(PSWFmx(i,1)+mu),i)./Hs(i);
end
%relcap=seddepth./repmat(Hs',timesteps,1);
overcap=relcap;
overcap(overcap<1)=NaN;
overcap(overcap>=1)=0;
%%
% Ppulse(1:timesteps,1:LinkNum)=NaN;
% Ppvol(1:timesteps,1:LinkNum)=NaN;
% %ppvol=P_vol{2001,446}(1,1);
% for t=1:timesteps
%     t
%     for i=1:LinkNum        
%         Ppulse(t,i)=sum(P_idx{t,i}>=ppidxstart);  
%         Ppvol(t,i)=sum(P_vol{t,i}(P_idx{t,i}>=ppidxstart));
%     end
% end
% 
% save('BE_NHD_MartinLakes4_BRU_pulseall_300yr_v2_Ppulse.mat',...
%     'Ppulse','Ppvol');

%%
% Pidx=P_idx(6001,:);
% Ploc=P_loc(6001,:);
% Pstor=P_storage(6001,:);
% Pvol=P_vol(6001,:);
% elev=Elev(6001,:);
% 
% save('BE_NHD_MartinLakes4_BRU_pulseall_300yr_v2_seed.mat',...
%     'Pidx','Ploc','Pstor','Pvol','elev','Slope');

%%
Ppulsep=Ppulse./max(max(Ppulse));
%%
Pout=diff(1-sum(Ppulsep,2));
%%
Ppulsed=Ppulse.*ppvol./repmat((Length.*B)',timesteps,1);
%%
Ppulsep=Ppvol./sum(Ppvol(2001,:));
%% Plot x-axis distance
% writerObj = VideoWriter('ex1','MPEG-4');
% writerObj.FrameRate = 15;
% writerObj.Quality =100;
% open(writerObj);

f1=figure;
hold on
box on
i=446;
t=3101;

% set(gcf,'Renderer','Painters');
% for t=2000:12001

clear x y xs is y2 y3
x=WFout(Connect(i,~isnan(Connect(i,:)),1),1);
x=cat(1,x,x(2:end),0);
y=seddepthd(t,Connect(i,~isnan(Connect(i,:))))';
y(y<0)=0;
%y=trate(Connect(i,~isnan(Connect(i,:))),1);
y=cat(1,y,y);
y2=overcap(t,Connect(i,~isnan(Connect(i,:))))';
y2=cat(1,y2,y2);
y3=0.5*Ppulsep(t,Connect(i,~isnan(Connect(i,:))))';
y3=cat(1,y3,y3);
[xs is]=sort(x);

plot(xs,y3(is),'r')
hold on
plot(xs,y(is),'b')
plot(xs,y2(is),'k','LineWidth',5)
%ylim([-0.5 0.5])
ylim([0 0.5])
xlim([0 max(x)])
title(cat(2,num2str(time(t)),' years'))
set(gca,'XDir','reverse')
% 
% frame = getframe(gcf);
% writeVideo(writerObj,frame);
% delete(gca)
% end
% close(writerObj);

%% 
overcap(overcap==0)=1;
overcap(isnan(overcap))=0;
%%
fracocap=sum(overcap(1:12000,:),1)./12000;

%% Plot elevation profile
figure
hold on
box on
plot(cat(1,WFout(Connect(i,~isnan(Connect(i,:)),1),1),0),...
    cat(1,Elev(t,Connect(i,~isnan(Connect(i,:))),1)',mnelev(OutletLinkID)),'b')
xlim([0 max(x)])
set(gca,'XDir','reverse')
%plot(WFout(Connect(i,~isnan(Connect(i,:)),1),1),...
%    trate(Connect(i,~isnan(Connect(i,:))),1),'k')
%plot(WFout(Connect(i,~isnan(Connect(i,:)),1),1),...
%    seddepthd(t,Connect(i,~isnan(Connect(i,:)))),'k')
%%
figure
contourf(seddepthd(:,Connect(i,~isnan(Connect(i,:)))))
%%
OutVold=OutVol;%-OutVol_base;

OutVolyr(1:ceil(tmax)-1,1)=NaN;
for i=1:ceil(tmax)-1
OutVolyr(i,1)=sum(OutVold(find(ceil(time)==i)));
end
%%
Pout(Pout==-1)=0;
Pout=cat(1,Pout,0);
Poutyr(1:ceil(tmax)-1,1)=NaN;
for i=1:ceil(tmax)-1
Poutyr(i,1)=sum(Pout(find(ceil(time)==i)));
end
%%
clear y ys
y=OutVolyr;
%y2=Poutyr;
%ys=y;
%windowSize=5;
%ys=filter(ones(1,windowSize)/windowSize,1,y);
figure
hold on
box on
plot(1:ceil(tmax)-1,y.*2.65,'b')
%plot([PSWFmx(OutletLinkID,1) PSWFmx(OutletLinkID,1)], [0 500000],'k')
plot([100 100], [0 500000],'k')
%plot(1:ceil(tmax)-1,y2.*-50000,'r')
%plot(1:ceil(tmax)-1,ys,'k','LineWidth',2)

%%
clear xbound wfx wfy
xbound=0:10:250;
wfx=xbound+xbound(2)/2;
wfx(end)=[];
for i=1:length(xbound)-1
    wfy(i,1)=sum(and(WFout>xbound(i),WFout<xbound(i+1)));
end
figure
plot(wfx,wfy./sum(wfy))
%%
clear pwfy
pwfy(timesteps,length(xbound)-1)=NaN;
for ti=2001:12001
    for i=1:length(xbound)-1
        pwfy(ti,i)=sum(Ppulsep(ti,and(WFout>xbound(i),WFout<xbound(i+1))));
    end
end
%%
ti=3001;
figure
plot(wfx,pwfy(ti,:))
title(num2str(time(ti)))
ylim([0 0.1])
xlim([0 250])
%%
[N, xL]=hist(WFout,50);

%NN = max(size(Y));
%NN = 2427;
%x1 = xL./max(xL);
%y1 = N./NN;
x2 = xL;
y2 = N;
figure
plot(x2,y2)


%% Input shredding?
%%
in=max(PSWFout)./100;
PSWF=PSWFout;
for i=1:100
    PSWF=cat(1,PSWF,in.*i+PSWFout);
end
%%
figure
hold on
%[y x]=hist(WF,1000);
[y x]=hist(WFi{1034,1},20);
plot(x,y,'k')
%[y x]=hist(WFout(logical(Lake)),50);
%plot(x,y,'b')
%%
%dt=diff(x(1:2));
y=art;
Fs=1/dt;
N = length(y);
% FFT compute
xdft = fft(y);
xdft = xdft(1:floor(N/2)+1);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(y):Fs/2;
figure
hold on
box on
set(gca,'YScale','log','XScale','log')
% plot(freq,psdx,'k')%'LineWidth',2)
% xlabel('Frequency (Hz)'); ylabel('Power');
%plot(1./freq,psdx,'k')%'LineWidth',2)
plot(1./freq./60./60./24./365,psdx,'k')%'LineWidth',2)
xlabel('Period (yr)'); ylabel('Power');
%
%% WT
clear coefs scales
%effective support mexh=5, morl=4
effsup=4;
sig.val=y;
sig.period=dt;
cwtstruct = cwtft(sig);%default morlet wo=6
wfreq=1./(4.*pi.*cwtstruct.scales./(6+sqrt(2+6^2)));%wo=6, morl
clear wpsd
%plot
mfreq=1./(4.*pi.*(N/effsup/2)./(6+sqrt(2+6^2)));
%vscal=find(wfreq<mfreq,1,'first')-1;
vscal=wfreq;
wpsd=(0.776/(N))*nansum(abs(cwtstruct.cfs).^2,2);
 figure
% hold on
% box on
%plot(wfreq(1:vscal),wpsd(1:vscal),'LineWidth',2)
plot(1./wfreq./60./60./24./365,wpsd,'LineWidth',2)
set(gca,'YScale','log','XScale','log')
xlabel('Frequency (Hz)'); ylabel('Power');

%% Basis for theta via Rouse profile
D=0.0004; %m
R=1.65; %submerged specific gravity
g=9.81; %m/s2, gravitational acceleration
nu=10^-6; %m2/s, kinematic viscosity
kappa=0.41; %von Karman constant
ustar=sqrt(g.*H.*Slope);%m/s
Rep=(sqrt(g.*R.*D).*D)./nu; % particle Reynolds number (dimensionless)

% Dietrich (1982) relation for fall velocity
% coefficients
b1=2.891394;
b2=0.95296;
b3=0.056835;
b4=0.002892;
b5=0.000245;
% functional relation
Rf=exp(-b1 + b2.*log(Rep) - b3.*(log(Rep)).^2 - ...
    b4.*(log(Rep)).^3 + b5.*(log(Rep)).^4);

vs=Rf.*sqrt(g.*R.*D);%m/s, fall velocity

ZR=vs./kappa./ustar;% Rouse Number (dimensionless, m/s/m/s)
%%
f10ccb(1:LinkNum,1)=NaN;
f10ccbuus(1:LinkNum,1)=NaN;
zHq(1:LinkNum,1)=NaN;
for i=1:LinkNum
    clear zH ccb zH ccb uus Fccb Fccbuus
%ZR=2;%0.5;
zH=(0.05:0.01:0.99)';
ccb=(((1./zH)-1)./((1./0.05)-1)).^ZR(i);% Rousean distribution for suspended sediment
zH=cat(1,(0.01:0.01:0.04)',zH);
ccb=cat(1,repmat(1,4,1),ccb);
%
ks=2*D;
uus=kappa.*log(30.*zH.*H(i)./ks);
Fccb=cumsum((ccb).*0.01);
f10ccb(i)=Fccb(10)./max(Fccb);
Fccbuus=cumsum((ccb.*uus).*0.01);
f10ccbuus(i)=Fccbuus(10)./max(Fccbuus);
Fccbuusnorm=Fccbuus./max(Fccbuus);

sth=0.9;
Fccbuusnorm=cat(1,0,Fccbuusnorm);
zH=cat(1,0,zH);
for j=1:length(zH)-1
    if Fccbuusnorm(j,1)<sth && Fccbuusnorm(j+1,1)>=sth
        zHq(i,1)=interp1(Fccbuusnorm(j:j+1,1),zH(j:j+1,1),sth);
    end
end
%zHq(i) = interp1(Fccbuusnorm,zH,0.1);
end
%%
figure; hold on; box on
xlim([0 1])
ylim([0 1])
plot(uus./max(uus),zH,'-b')
%
figure; hold on; box on
%set(gca,'XScale','log','YDir','reverse')
set(gca,'XScale','log')
%xlim([0.001 1])
ylim([0 1])
plot(ccb,zH,'-b')
%plot(ccbest,zHest,':r')
xlabel('Relative concentration, c/c_b')
ylabel('Relative depth, z/H')
%title([strcat('Z_R = ',num2str(round(ZR.*100)./100)),...
%    sprintf('\n'),...
%    strcat('particle diameter = ',num2str(dsf*1000000),' microns')])
%
figure; hold on; box on
xlim([0 1])
ylim([0 1])
plot(Fccb./max(Fccb),zH,'-b')
plot([0 1],[0.1 0.1],'k')
%
figure; hold on; box on
set(gca,'XScale','log')
%xlim([0 1])
ylim([0 1])
plot(ccb.*uus,zH,'b')
%
figure; hold on; box on
xlim([0 1])
ylim([0 1])
plot(Fccbuus./max(Fccbuus),zH,'-b')
plot([0 1],[0.1 0.1],'k')

%%


%%
% % estimate the slope of the distribution between 0.2 and 0.4 depth
% whichstats = {'tstat','yhat','rsquare'};
% clear zHest ccbest
% %for i=1:size(ccb,2)
%     clear Ystat Xstat stats1 m b idx
%     idx=(and(zH>=0.6,zH<=0.8));
%     Xstat = zH(idx,1);
%     Ystat = log(ccb(idx,1));
%     stats1 = regstats(Ystat,Xstat,'linear',whichstats);
%     m = stats1.tstat.beta(2);
%     b = stats1.tstat.beta(1);
%     % estimate the distribution near the surface to 0.2 depth
%     zHest=(0.8:0.01:0.99);
%     ccbest=exp(m.*zHest+b);
% %end
