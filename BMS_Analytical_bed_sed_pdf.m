%% Analytical Bed Sediment PDF
% This code computes an analytical bed-sediment depth pdf when
% in-channel storage effects slope and thus transport.

% January 21, 2016
% Jon Czuba

%% VpLi generated, arrival, capacity
VpLi_gen(1:LinkNum,1)=0;%vol rate generated in each link
VpLi_arr(1:LinkNum,1)=0;%vol rate arriving to and generated in each link
VpLi_cap(1:LinkNum,1)=0;%vol rate capacity of each link
VpLi_ni(1:LinkNum,1)=0;%number of effective inputs in each link
VpLi_tni(1:LinkNum,1)=0;%number of effective inputs arriving and generated in each link

for i=1:LinkNum
    VpLi_gen(i,1)=sum(In_Vol(find(In_Link==i),1));
    VpLi_ni(i,1)=sum(In_Link==i);
end
%VpLi_cap=Velocity.*60.*60.*24.*365.*B.*theta.*H;%m3/yr, initial capacity, no limiter
VpLi_arr=sum(Sub.*repmat(VpLi_gen,1,LinkNum),1)';%volume generated in upstream links including itself
VpLi_tni=sum(Sub.*repmat(VpLi_ni,1,LinkNum),1)';%number of discrete input locations in upstream links including itself

%figure
%plot(VpLi_arr./VpLi_cap)

%%
%load('Poisson\BE_NHD_ML4_v4_setupfor_600Vpconst_loc0.mat')
%load('Poisson\BE_NHD_ML4_v4_600Vpconst_loc0_sedbudinchstor2.mat', 'Elev', 'OutVol', 'seddepth')
%load('Poisson\BE_NHD_ML4_v4_600Vpconst_loc0_sedbudnostor.mat', 'Elev', 'OutVol', 'seddepth')
%% Simulated Elevation Change
Elevch=Elev-repmat(Elev(1,:),timesteps,1);
%figure
%plot(time,Elevch(:,230))
%% Simulated Slope Change
%compute/update slope
sslope(1:timesteps,1:LinkNum)=NaN;
for t=1:timesteps
    for i=1:LinkNum
        if i==OutletLinkID
            sslope(t,i)=(Elev(t,i)-mnelev(i,1))./Length(i,1);
        else
            sslope(t,i)=(Elev(t,i)-Elev(t,Connect(i,2)))./Length(i,1);
        end
    end
end
Slope(Slope<1e-4)=1e-4;
%
sslopech=sslope-repmat(sslope(1,:),timesteps,1);

%% slope at capacity
mu=1;
Sc=((Vp.*VpLi_tni./mu./365./24./60./60.*sqrt(9.81).*1.65.*1.65.*D)./...
    (0.05.*0.175.*B.*H.^(3/2).*U.^(2))).^(2/3);

%% Compute capacity slopes in network
%initialize
SlopeNew=slope;%lowercase slope is original
ElevNew=Elev(1,:)';
hstor(1:LinkNum,1)=0;
idlac=1;
itcntr=0;
%iterate to find capacity slopes
while ~isempty(idlac)
    %compute capacity
    Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*SlopeNew.^(3/2)./theta.*0.175;%m/s
    Velocity(SlopeNew<0)=0;%set to zero if negative slope
    STTime=Length./Velocity;%seconds, travel time through each link
    STTimeyr=STTime./365./24./60./60;%years
    VpLi_cap=capacity./STTimeyr;%m3/yr, initial capacity, including limiter
    % compute relative capacity
    RelCap=VpLi_arr./VpLi_cap;
    figure
    plot(RelCap)
    % identify links above capacity
    idlac=find(RelCap>1.001);
    % adjust elevations based on new slopes
    en=Length(idlac).*(Sc(idlac)-SlopeNew(idlac));
    %ElevNew(idlac)=ElevNew(idlac)+Length(idlac).*(Sc(idlac)-SlopeNew(idlac)); 
    ElevNew(idlac)=ElevNew(idlac)+en;
    % compute depth equivalent volumes in storage
    for i=1:length(idlac)
        usid=find(Connect(:,2)==idlac(i));%determine US links           
        nlusid=cat(1,idlac(i),usid(Lake(usid)==0));%include current link and upstream links that are not lakes
        hstor(idlac(i))=hstor(idlac(i))+en(i)./2.*(1-Lp).*...
            sum(B(nlusid,1).*Length(nlusid,1))./B(idlac(i))./Length(idlac(i));
    end
    % comupte new slopes for all links
    for i=1:LinkNum
        if i==OutletLinkID
            SlopeNew(i,1)=(ElevNew(i,1)-mnelev(i,1))./Length(i,1);
        else
            SlopeNew(i,1)=(ElevNew(i,1)-ElevNew(Connect(i,2),1))./Length(i,1);
        end
    end
    itcntr=itcntr+1;%increment counter
end

SlopeNew(SlopeNew<1e-4)=1e-4;
Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*SlopeNew.^(3/2)./theta.*0.175;%m/s
STTime=Length./Velocity;%seconds, travel time through each link
STTimeyr=STTime./365./24./60./60;%years
VpLi_cap=capacity./STTimeyr;%m3/yr, initial capacity, including limiter
RelCap=VpLi_arr./VpLi_cap;

%% Compute simulated capacity
%vvelocity=0.05./sqrt(9.81)./1.65./1.65./D.*repmat(U',timesteps,1).^2.*...
%    repmat(H',timesteps,1).^(1/2).*repmat(Sc',timesteps,1).^(3/2)./theta.*0.175;%m/s, realtime
vvelocity=0.05./sqrt(9.81)./1.65./1.65./D.*repmat(U',timesteps,1).^2.*...
    repmat(H',timesteps,1).^(1/2).*sslope.^(3/2)./theta.*0.175;%m/s, realtime
ccapacityrate=vvelocity.*60.*60.*24.*365.*repmat(B',timesteps,1).*...
    theta.*repmat(H',timesteps,1);%m3/yr

%%
i=17;
figure
hold on
plot(time,ccapacityrate(:,i),'b')
plot([0 tmax],[VpLi_arr(i) VpLi_arr(i)],'k')
ylabel('Volume flux, m^3/yr')
legend('simulated transport capacity','arrival')

figure
hold on
plot(time,sslope(:,i),'b')
plot([0 tmax],[Sc(i) Sc(i)],'k')
plot([0 tmax],[SlopeNew(i) SlopeNew(i)],'g')
ylabel('Slope')
legend('simulated','at-capacity','target')

% figure
% hold on
% plot(time,Elevch(:,i),'b')
% ylabel('Change in elevation, m')

%%
%% 

%% Analytical sediment depth at capacity
Hs_theo=theta.*H;%m, sediment depth at theoretical capacity
Hs_adj=capacity./Length./B;%m, sediment depth at adjusted capacity
%Hs_sim

%% Compute Statistics of Sediment Depth (SD)
SD_n(1:LinkNum,1)=NaN; %number of values of sed depth (SD)
SD_m(1:LinkNum,1)=NaN; %mean of SD
SD_std(1:LinkNum,1)=NaN; %std of SD
SD_var(1:LinkNum,1)=NaN; %var of SD
SD_mx(1:LinkNum,1)=NaN; %max of SD
%SD_fexHs(1:LinkNum,1)=NaN; %fraction of time of SD is above capacity
SD_iqr(1:LinkNum,1)=NaN; %iqr of SD

for i = 1:LinkNum
    %art=seddepth(time>(PSWFmx(i,1)+mu),i);%use data greater than PSWFi+mu
    art=seddepth(time>200,i)./(1-Lp);%use data greater than PSWFi+mu
    SD_n(i,1)=length(art); %number of interarrival times (SD)
    SD_m(i,1)=mean(art); %mean of SD
    SD_std(i,1)=std(art); %std of SD
    SD_var(i,1)=var(art); %var of SD
    SD_mx(i,1)=max(art); %max of SD
    %SD_fexHs(i,1)=sum(art>Hs_adj(i))./length(art);
    SD_iqr(i,1)=diff(quantile(art,[0.25 0.75])); %iqr of SD
    clear art
end

%% Analytical sediment depth
% analytical mean sediment depth, accounts for bed storage
hs_theo_avg=(VpLi_arr./Velocity./B./365./24./60./60+hstor)./(1-Lp);
%hs_theo_pdf
%hs_sim=seddepth;

%% Compute Analytical PDF of SD
rho=STTimeyr.*VpLi_tni./1;

clear ph cph n
n=(0:1:10000)';
ph(1:length(n),1:LinkNum)=NaN;
for i=1:LinkNum
    ph(:,i)=poisspdf(n,rho(i));
    cph(:,i)=poisscdf(n,rho(i));
end
%Vp=10;%m3
hs_theo_pdf=repmat(n,1,LinkNum).*Vp./repmat(Length',length(n),1)./repmat(B',length(n),1);

%% Plot histogram of bed sediment thickness
%hval=hist(numpar(:,1147),0:1:50);
i=406%159%284
[hval, cent]=hist(seddepth(time>200,i),15);
%[hval, cent]=hist(seddepth(time>(PSWFmx(i,1)+mu),i),15);
%[hval, cent]=hist(lnkvol(time>(PSWFmx(i,1)+mu),i)./Vp,15);
%[hval, cent]=hist(seddepth(time>100,i),15);

abw=diff(hs_theo_pdf(1:2,i));
nbw=diff(cent(1:2));
bwsf=nbw/abw;

clear hsx hsy
hsx(1:length(cent)*2,1)=NaN;
hsy(1:length(cent)*2,1)=NaN;
for ii=1:length(cent)
   hsx(ii*2-1,1)=cent(1,ii)-nbw/2;
   hsx(ii*2,1)=cent(1,ii)+nbw/2;
   hsy(ii*2-1:ii*2,1)=hval(1,ii);
end
hsx=cat(1,hsx(1),hsx,hsx(end));
hsy=cat(1,0,hsy,0);

figure
hold on
box on
%bar(cent./(1-Lp),hval./sum(hval))
plot(hsx./(1-Lp),hsy./sum(hval),'-b')
%Elevch is not the correct depth adjustment
plot((hs_theo_pdf(:,i)+hstor(i))./(1-Lp),ph(:,i).*bwsf,'k','LineWidth',1)
%plot((hs_theo_pdf(:,i))./(1-Lp),ph(:,i).*bwsf,'k','LineWidth',1)
%plot(n,ph(:,i).*bwsf,'k','LineWidth',3)
%plot([(Hs_adj(i)+hstor(i))./(1-Lp) (Hs_adj(i)+hstor(i))./(1-Lp)],[0 max(ph(:,i).*bwsf)],'r')
%plot([hs_theo_avg(i)./(1-Lp) hs_theo_avg(i)./(1-Lp)],[0 max(ph(:,i).*bwsf)],'r')
plot([hs_theo_avg(i) hs_theo_avg(i)],[0 max(ph(:,i).*bwsf)],'r')
%xlim([0 1])
xlabel('Bed sediment thickness, m')
ylabel('Probability')
title(num2str(i))

% Determine fraction of time capacity is exceeded
% pHs(1:LinkNum,1)=NaN;
% for i=1:LinkNum
%     pHs(i,1)=interp1(hs_theo_pdf(:,i)+hstor(i),cph(:,i),Hs_adj(i)+hstor(i));
% end
% fexHs=1-pHs;

%% Plot numerical SD timeseries
%i=376%17%284
%prd2(1:LinkNum,1)=0;
%dt=diff(x(1:2));
%for j=1:48
%close all
%j=j+1
%i=lnks2(j)
i=406
figure; hold on; box on
%set(gcf,'Position',p1);
plot(time,seddepth(:,i)./(1-Lp))
%plot(time,sdsave(i,:))
plot([PSWFmx(i,1)+mu PSWFmx(i,1)+mu], [0 max(seddepth(:,i)./(1-Lp))],'k')
 %plot([0 tmax],[(Hs_adj(i)+hstor(i))./(1-Lp) (Hs_adj(i)+hstor(i))./(1-Lp)],'k')
 %plot([0 tmax],[hs_theo_avg(i)./(1-Lp) hs_theo_avg(i)./(1-Lp)],'k')
plot([0 tmax],[hs_theo_avg(i) hs_theo_avg(i)],'k')
 %plot([0 tmax],[hstor(i) hstor(i)],'k')
 %plot([0 tmax],[Hs_adj(i) Hs_adj(i)],'k')
 %ylim([0 1])
title(num2str(i))
xlabel('Time, years'); ylabel('Bed sediment thickness, m')
%%
% FFT of sediment depth
y=seddepth(4001:end,i);
%y=sdsave(i,4001:end);
Fs=1/(daystp./365);%1/dt;
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
set(gcf,'Position',p2);
% plot(freq,psdx,'k')%'LineWidth',2)
% xlabel('Frequency (Hz)'); ylabel('Power');
%plot(1./freq,psdx,'k')%'LineWidth',2)
plot(1./freq,psdx,'k')%'LineWidth',2)
%plot(1./freq./60./60./24./365,psdx,'k')%'LineWidth',2)
xlabel('Period (yr)'); ylabel('Power');
title(i)

[dummy, I]=max(psdx(14:601));
plot(1./freq(14+I-1),psdx(14+I-1),'or')
%prd2(i,1)=1./freq(14+I-1);
%waitforbuttonpress
%close all
%end
%%
prd1(prd1==0)=NaN;
prd2(prd2==0)=NaN;
%%
prd1(isnan(prd1))=0;
prd2(isnan(prd2))=0;
%% Periodicity
%val=usarea;
%val=Length.*0.1.*H.*B./STTimeyr;
val=VpLi_arr;
val(val==0)=NaN;
%figure; hold on; box on
plot(val,prd1,'.r')
%set(gca,'YScale','log','XScale','log')

stats=regstats(log10(prd1),log10(val),'linear','all');
mnmx=cat(1,min(val),max(val));
plot(mnmx,10^stats.beta(1).*(mnmx).^stats.beta(2),'-r')
title(['y = ',num2str(10^stats.beta(1)),' x ',num2str(stats.beta(2)),...
    '; R2 = ',num2str(stats.rsquare)])
  
%% Plot comparison of numerical and analytical mean SD
figure; hold on; box on
%plot(SD_m./(1-Lp),hs_theo_avg./(1-Lp),'.b')
plot(SD_m,hs_theo_avg,'.b')
%plot(1000*(hs_theo_avg),(SD_std*1000).^2,'.b')
set(gca,'YScale','log','XScale','log')
xlabel('Mean of simulated sediment depth')
ylabel('Theoretical mean sediment depth')
plot([10^-5 10^1],[10^-5 10^1],'-k')
%xlim([10^-4 10^4])
plot(Param(:,8),Param_hs_theo_avg,'.r')
%text(Param(:,8),hs_theo_avg,num2str(Param(:,1)))

%% Compute capacity slopes in single reach storage dynamic
%initialize
Param_SlopeNew=slope;%lowercase slope is original
Param_hstor(1:LinkNum,1)=0;
%compute capacity
Param_Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Param_SlopeNew.^(3/2)./theta.*0.175;%m/s
Param_STTime=Length./Param_Velocity;%seconds, travel time through each link
Param_STTimeyr=Param_STTime./365./24./60./60;%years
Param_VpLi_cap=capacity./Param_STTimeyr;%m3/yr, initial capacity, including limiter
% compute relative capacity
Param_RelCap=VpLi_arr./Param_VpLi_cap;
figure
plot(Param_RelCap)
%%
% identify links above capacity
idlac=find(Param_RelCap>1.00001);
en=Length(idlac).*(Sc(idlac)-Param_SlopeNew(idlac));
for i=1:length(idlac)
    %Param_hstor(idlac(i))=en(i)./2.*(1-Lp); %height will need /(1-Lp)
    Param_hstor(idlac(i))=en(i)./2.*(1-Lp).*...
            Param(idlac(i),6)./B(idlac(i))./Length(idlac(i));
end
Param_SlopeNew(idlac)=Sc(idlac);
Param_SlopeNew(Param_SlopeNew<1e-4)=1e-4;
Param_Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Param_SlopeNew.^(3/2)./theta.*0.175;%m/s
Param_STTime=Length./Param_Velocity;%seconds, travel time through each link
Param_STTimeyr=Param_STTime./365./24./60./60;%years
Param_VpLi_cap=capacity./Param_STTimeyr;%m3/yr, initial capacity, including limiter
% compute relative capacity
Param_RelCap=VpLi_arr./Param_VpLi_cap;
figure
plot(Param_RelCap)
%%
Param_hs_theo_avg=(VpLi_arr./Param_Velocity./B./365./24./60./60+Param_hstor)./(1-Lp);
Param_rho=Param_STTimeyr.*VpLi_tni./1;

%%
%BL=B.*Length;
figure; hold on; box on
plot(rho,SD_std./(SD_m-(hstor)./(1-Lp)),'.b')
%plot(SD_m,SD_var,'.b')
%plot(SD_m(BL>22000),SD_var(BL>22000),'.r')
%plot(1000*(hs_theo_avg),(SD_std*1000).^2,'.b')
set(gca,'YScale','log','XScale','log')
xlabel('rho=t*n*lamdaprime')
ylabel('COV of simulated sediment depth, std/mean')
%plot([10^-3 10^5],[10^5 10^-3],'-k')
plot([10^-3 10^5],[(10^-3)^-0.5 (10^5)^-0.5],'k')
xlim([10^-3 10^4])
ylim([10^-2 10^2])
%%
plot(Param_rho,Param(:,9)./(Param(:,8)-(Param_hstor)./(1-Lp)),'.r')
%plot(rho(Param(:,1)),Param(:,9)./(Param(:,8)),'or')
%text(rho(Param(:,1)),Param(:,9)./(Param(:,8)-(hstor(Param(:,1)))./(1-Lp)),num2str(Param(:,1)))
%%
p1=get(gcf,'Position');
%%
set(gcf,'Position',p1);
%%
% lnks=[717 507 213 17];
% plot(rho(lnks),SD_std(lnks)./(SD_m(lnks)-(hstor(lnks))./(1-Lp)),'.r')
% for ii=1:length(lnks)
% text(rho(lnks(ii)),SD_std(lnks(ii))./(SD_m(lnks(ii))-(hstor(lnks(ii)))./(1-Lp)),num2str(lnks(ii)))
% end
%%
dev1=((SD_std./(SD_m-(hstor)./(1-Lp)))-rho.^-0.5)./rho.^-0.5;
dev2=((Param(:,9)./(Param(:,8)-(Param_hstor)./(1-Lp)))-Param_rho.^-0.5)./Param_rho.^-0.5;
%dev2=(Param(:,9)./(Param(:,8))-rho(Param(:,1)).^-0.5)./rho(Param(:,1)).^-0.5;
figure; hold on; box on
set(gca,'XScale','log')
plot(rho,dev1,'.b')
plot(Param_rho,dev2,'.r')
%%
plot([10^-3 10^4],[0.2 0.2],'-k')
%%
%plot(rho(hstor>0),dev(hstor>0),'.r')

% plot(rho(lnks),dev(lnks),'.r')
% for ii=1:length(lnks)
% text(rho(lnks(ii)),dev(lnks(ii)),num2str(lnks(ii)))
% end
%clear lnks2
lnks1=find(dev1>0.2);
lnks2=find(dev2>0.2);
%plot(rho(lnks1),dev(lnks1),'.r')
%plot(Param_rho(lnks2),dev2(lnks2),'.k')
% for ii=1:length(lnks2)
% text(rho(lnks2(ii)),dev(lnks2(ii)),num2str(lnks2(ii)))
% end
% xlabel('rho=t*n*lamdaprime')
% ylabel('Relative COV deviation from Poisson')
%%
figure; hold on; box on
plot(Param_RelCap,dev2,'.r')
%plot(Param_rho(lnks2),Param_RelCap(lnks2),'.k')
%set(gca,'YScale','linear','XScale','log')
%%
plot([0.995 0.995],[0 20],'-k')
%%
val=BiLi;
%val=B.*Length;
%val=Param(:,6).*Length;
%val=Param_rho;
%BL=BiLi.*Slope;
lnks3=find(Param_RelCap>0.995);
figure; hold on; box on
plot(val(lnks3),dev2(lnks3),'.r')
%plot(val(lnks2),dev2(lnks2),'.k')
% for ii=1:length(lnks3)
% text(BiLi(lnks3(ii)),dev2(lnks3(ii)),num2str(lnks3(ii)))
% end
set(gca,'YScale','log','XScale','log')
%plot(B(Param(lnks2,1)).*Length(Param(lnks2,1)),dev2(lnks2),'.k')

stats=regstats(log10(dev2(lnks3)),log10(val(lnks3)),'linear','all');
mnmx=cat(1,min(val(lnks3)),max(val(lnks3)));
plot(mnmx,10^stats.beta(1).*(mnmx).^stats.beta(2),'-k')
title(['y = ',num2str(10^stats.beta(1)),' x ',num2str(stats.beta(2)),...
    '; R2 = ',num2str(stats.rsquare)])
    
%%
var1(1:LinkNum,1)=0;
var2(1:LinkNum,1)=0;
var1(lnks1)=1;
var2(lnks2)=1;
class(1:LinkNum,1)=NaN;
class(and(var1==1,var2==1))=1;%gen
class(and(var1==0,var2==0))=0;%none
class(and(var1==1,var2==0))=2;%prop
class(and(var1==0,var2==1))=3;%damp
% 
% gen(1:LinkNum,1)=0;
% non(1:LinkNum,1)=0;
% pro(1:LinkNum,1)=0;
% damp(1:LinkNum,1)=0;
% gen(Param(lnks2))=1;%generate variability
% non(lnks1)=1;%no variability
% pro(~(gen+non))=1;%propagate or amplify variability
% gen=logical(gen);
% non=logical(non);
% pro=logical(pro);
%%
non1=and(non,RelCap>0.95);
val=rho;
figure; hold on; box on
plot(val(lnks1),dev(lnks1),'.b')
plot(val(non1),dev(non1),'.k')
%plot(BiLi(RelCap>0.9),dev(RelCap>0.9),'.r')
plot(val(Param(lnks2,1)),dev2(lnks2),'.r')
set(gca,'XScale','log')
text(val(non1),dev(non1),num2str(GridID(non1)))
%%
figure; hold on
plot(BiLi(lnks1),RelCap(lnks1),'.b')
plot(BiLi(Param(lnks2,1)),RelCap(Param(lnks2,1)),'.k')
set(gca,'XScale','log')
%%
x=RelCap;
y=BiLi;
figure; hold on; box on
plot(x(non),y(non),'.b')
plot(x(pro),y(pro),'.k')
plot(x(gen),y(gen),'.r')
set(gca,'XScale','log','YScale','log')
%%
BiLi(1:LinkNum,1)=NaN;
for ii=1:LinkNum
    usid=find(Connect(:,2)==ii);%determine US links
    nlusid=cat(1,ii,usid(Lake(usid)==0));%include current link and upstream links that are not lakes
    %BiLi(ii,1)=2./(1-Lp)./sum(B(nlusid,1).*Length(nlusid,1))./Length(ii);
    BiLi(ii,1)=sum(B(nlusid,1).*Length(nlusid,1)).*Length(ii);
    %BiLi(ii,1)=sum(B(nlusid,1).*Length(nlusid,1));
    clear usid nlusid
end
%%
%find(dev>0.2)
ans=GridID;
Param(1:length(ans),1:9)=NaN;
Param(:,1)=ans;
Param(:,2)=usarea(ans);
Param(:,3)=capacity(ans);
Param(:,4)=slope(ans);%was SlopeNew in early trials
Param(:,5)=Length(ans);
Param(:,6)=BiLi(ans);
Param(:,7)=VpLi_tni(ans);

sdsave(1:length(ans),1:12001)=NaN;

%%
%%
%% Plot Spatial Data
%set timestep to plot
t = timesteps-1;

for i = 1:LinkNum
      %[network(i).sed] = hs_theo_avg(i).*1000./(1-Lp);
      %[network(i).sed] = hstor(i).*1000./(1-Lp);
      %[network(i).sed] = uptail(i).*1000./(1-Lp);
      %[network(i).sed] = SD_fexHs(i).*1000;
      %[network(i).sed] = dev(i);
      [network(i).sed] = prd1(i);
end

f1 = figure;
%set(f1,'Position',[213 50 938 632]);
%a1 = axes('Position',[0.08 0.1 0.4 0.8]);
a1 = axes;
%cbpos = [0.78 0.22 0.01 0.6];%mrb,hw,c
cbpos = [0.90 0.22 0.01 0.6];%be
%cbpos = [0.93 0.22 0.01 0.6];%w
box(a1);

%edge = [0 1 2 4 8 16 32 64 128 200];
%edge=[0 0.0004 0.0007 0.001 0.004 0.007 0.01 0.04 0.07 0.1];
%edge = [0 1 2 4 6 8 10 12 16 20];
edge = [0 1 2 3 4 5 10 20 50 100];
%edge = [0 1 10 50 100 500 1000 5000 10000 100000];
%edge = [0 1 10 50 70 100 200 300 500 1000];
%edge = [0 .5 .75 1 1.25 1.5 1.75 2 2.25 2.5];
%edge = cat(2,0,logspace(-2,0,9));
%edge = cat(2,0,round(exp(linspace(log(1),log(max([network.sed])),9))));
%edge = cat(2,0,(exp(linspace(log(min([network.sed])),log(max([network.sed])),9))));
%edge = cat(2,0,(exp(linspace(log(1),log(max(max(LCS))),9))));
%edge = cat(2,0,round(exp(linspace(log(1),log(max(max(numpar))),9))));
%edge = cat(2,0,round(linspace(mn,ceil(max(max(LCS))),9)));
%edge = cat(2,0,round(exp(linspace(log(mn),log(ceil(max(max(LCS)))),9))));

%edge = cat(2,0,round((linspace((min([network.Nconc])),(max([network.Nconc])),9))));
%edge = cat(2,0,round((linspace((1),(max(max(numpar))),9))));
%edge = cat(2,0,round(quantile([network.sed],0.1:0.1:0.8)),round(max([network.sed])));

% ttl={['Time = ',num2str(fix(t/365/0.175)),' years ',...
%      num2str(round(mod(t/0.175,365))),' days'];nttl};
%ttl={['Time = ',num2str(fix(time(t))),' days ',...
%    num2str(round(mod(time(t)*24,24))),' hours'];['Time step = ',num2str(t)];nttl};
%ttl=''; 
%[a1b, a1m] = Plot_Network_Map(a1,edge,cbpos,boundary,network,ttl);

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
    {'Default', 'Color','black'}, ...
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
a1l = mapshow(a1,lakes,'FaceColor','blue');

%title(a1,ttl,'FontSize',14);

%% Export data
%network = rmfield(network, 'L');

for i=1:LinkNum
   %[network(i).seddepth] = seddepth(6001,i);
   %[network(i).hstheoavg] = hs_theo_avg(i,1)./(1-Lp);
   %[network(i).VpLigen] = VpLi_gen(i,1).*2.65;
   %[network(i).numfracocap] = numfracocap(i); 
   %[network(i).exceedHs] = exceedHs(i);
   %[network(i).wfpeak1] = lnkpl(i);
   %[network(i).relcap1] = out(i,1);
   %[network(i).relcap2] = out(i,2);
   %[network(i).hstor1Lp] = out(i,3);
   [network(i).prd] = prd2(i);
end
% Write Shapefile
%shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att_overcap.shp');
shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_prd2.shp');

