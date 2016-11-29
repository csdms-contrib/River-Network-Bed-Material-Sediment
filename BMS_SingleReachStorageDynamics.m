%% Single Reach Storage Dynamics
% This code contains all the transport and storage dynamics of the network
% model, but for a single reach to allow for the exploration of these
% dynamics versus relevant link parameters.

% Jon Czuba
% January 25, 2016

%%
clear all
close all
clc

% load Param_sdsave4.mat
% 
% for iter=1:size(Param,1)
%     iter
%     usarea=Param(iter,2);
%     capacity=Param(iter,3);
%     Si=Param(iter,4);
%     Length=Param(iter,5);
%     sBiLi=Param(iter,6);
%     n=Param(iter,7);

%% Link 17 parameters
usarea=73636200;%m2
capacity=1.225884150335781e+03;%m3
%Si=1.027385983000258e-04;
Si=5e-5;
Length=1827;%m
sBiLi=1.772041739038366e+05;%m2, sum(B*L) for current reach and imm. US reaches
n=20;%VpLi_tni

%%
daystp=18.25; %number of days per timestep
dt = 1*24*60*60*daystp; %seconds in daystp number of days
timesteps=12001;%12001;%4001;%4000;%2115; %timesteps to simulate sand
tmax=dt*(timesteps-1)/60/60/24/365; %years, max time of simulation
time=(0:daystp/365:tmax)'; %years, time of current step

theta=0.1;
D=0.0004;%m
Lp=0.4;%porosity
Vp=10;%m3

% Determine flow and its scaling
B=(0.0238).*(usarea).^(0.3397);%m, width at Q2
H=(0.0029).*(usarea).^(0.2943);%m, water depth
U=(0.1976).*(usarea).^(0.0679);%m/s, water velocity

mu=1;
Sc=((Vp.*n./mu./365./24./60./60.*sqrt(9.81).*1.65.*1.65.*D)./...
    (0.05.*0.175.*B.*H.^(3/2).*U.^(2))).^(2/3);

InArrival=[];%years, arrival time of parcels to each link
OutArrival=[];%years, arrival time of parcels to DS link

Elev(1:timesteps,1)=0;
Slope(1:timesteps,1)=Si;

Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Si.^(3/2)./theta.*0.175;%m/s, realtime
STTime=Length./Velocity;%seconds, travel time through each link
STTimeyr=STTime./365./24./60./60;%years

% Initialize
% Parcel properties for 
P_idx=cell(timesteps,1);% index of each parcel (Source)
P_loc=cell(timesteps,1);% location of each parcel within a link (Loc)
P_storage=cell(timesteps,1);% parcel in storage (inactive)
P_vol=cell(timesteps,1);%m3, parcel volume

% Inputs
% Poisson inputs in continuous time at upstream end of all links
% constant parcel volume version 1

% change the seed to give truer random numbers
RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));

mu=1; %years, mean of interarrival time, exp. dist., =1/lambda
pidx=1; %set unique parcel index to 1 for first parcel
i=1;%one location, one input
R = exprnd(mu./n); %exponential random variables, interarrival time, years
cR=R; %arrival time, years
while cR<tmax, %loop through inputs
    % concatenate arrival times of parcels to each link
    %L_arrival{1,1}=cat(2,L_arrival{1,1},cR);
    InArrival=cat(1,InArrival,cR);
    
    % locate continuous time input at next timestep
    add = ceil(cR./(daystp./365))+1; %arrival at discrete timestep
    pdt=(time(add,1)-cR)*365*24*60*60; %parcel dt, convert years to seconds
    
    % Move parcels downstream
    sm=[];%time to move through current and downstream links
    csm=[];%cumulative time to move through current and ds links
    sm=STTime;%seconds, time to move out of current link
    csm=sm;%seconds
    if csm(end)<=pdt %if leave before dt ends
        sm=cat(1,sm,NaN);
        csm=cat(1,csm,sum(sm));
        OutArrival=cat(1,OutArrival,cR+csm(1)/60/60/24/365);
    else % update parcel location
        pi=1; %parcel link
        pl=pdt/STTime;%update location from us end
        
        P_loc{add,pi}=cat(2,P_loc{add,pi},pl);
        P_idx{add,pi}=cat(2,P_idx{add,pi},pidx);
        P_storage{add,pi}=cat(2,P_storage{add,pi},0);%activate
        P_vol{add,pi}=cat(2,P_vol{add,pi},Vp);
    end
    
    pidx=pidx+1; %increment unique parcel index
    
    % next input
    R = exprnd(mu./n);
    cR = cR + R;
end
%pidx
clear R cR add sm csm pi pl pdt
%
P_idx_IN=P_idx;
P_loc_IN=P_loc;
P_storage_IN=P_storage;
P_vol_IN=P_vol;

P_idx=cell(timesteps,1);
P_loc=cell(timesteps,1);
P_storage=cell(timesteps,1);
P_vol=cell(timesteps,1);
% Run RivNet
tic
i=1;
for t = 1:timesteps-1 %step through each timestep
    %t 
    % Generate autogenic inputs here? If desired.
    
    % Add new inputs at given timestep on top of what is already in system
    % background inputs
    P_loc{t,i}=cat(2,P_loc{t,i},P_loc_IN{t,i});
    P_idx{t,i}=cat(2,P_idx{t,i},P_idx_IN{t,i});
    P_storage{t,i}=cat(2,P_storage{t,i},P_storage_IN{t,i});
    P_vol{t,i}=cat(2,P_vol{t,i},P_vol_IN{t,i});

    % Determine capacity of links, active parcels, and adjust bed slopes
    % based on storage volume
    %only do this check capacity if parcels are in link
    if ~isempty(P_vol{t,i})
        %First In Last Out
        %compute cumulative volume in link beginning with last in
        cvol=fliplr(cumsum(fliplr(P_vol{t,i})));
        %determine which parcels that were the first to enter are above capacity
        exc=find(cvol>capacity(i,1),1,'last');
        %if parcels have been identified above capacity then
        if ~isempty(exc)
            %set their status to inactive, 1
            P_storage{t,i}(1,1:exc)=1;
            
            %UPDATE Elevations beginning here
            %compute volume inactive
            vstor=sum(P_vol{t,i}(1:exc))./(1-Lp);%m3 vol in storage, porosity Lp=0.4
            %update elevation at upstream end of link
            %volume is placed at upstream end of link and along current link
            %length and upstream link lengths to compute new elev
            Elev(t,i)=0+(2.*vstor)./sBiLi;
            %UPDATE Elevations complete
            Slope(t,i)=Si+(2.*vstor)./sBiLi./Length;
            
            clear cvol exc vstor
        end
    end
    %Slope(Slope<1e-4)=1e-4;

    Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope(t,i).^(3/2)./theta.*0.175;%m/s, realtime
    STTime=Length./Velocity;%seconds, travel time through each link

    % if current link is empty go to the next link
    if isempty(P_loc{t,i})
        continue
    end
           
    for p = 1:length(P_loc{t,i})

        % if parcel is inactive then update and go to the next parcel
        if P_storage{t,i}(1,p)
            P_loc{t+1,i}=cat(2,P_loc{t+1,i},P_loc{t,i}(1,p));%keep loc
            P_idx{t+1,i}=cat(2,P_idx{t+1,i},P_idx{t,i}(1,p));%keep idx
            P_storage{t+1,i}=cat(2,P_storage{t+1,i},0);% make active as
            %this is determined for each timestep
            P_vol{t+1,i}=cat(2,P_vol{t+1,i},P_vol{t,i}(1,p));%keep vol
            continue
        end
        
        % Move parcels downstream
        sm=[];%time to move through current and downstream links
        csm=[];%cumulative time to move through current and ds links
        sm=STTime(i,1)*(1-P_loc{t,i}(1,p));%seconds, time to move out of current link
        csm=sm;%seconds     
        if csm(end)<=dt %loop until cumulative time toward outlet exceeds timestep
            sm=cat(1,sm,NaN);
            csm=cat(1,csm,sum(sm));
            OutArrival=cat(1,OutArrival,time(t,1)+csm(1)/60/60/24/365);
        else
            pi=1;
            pl=(STTime(i,1)*P_loc{t,i}(1,p)+dt)/STTime(i,1);%update location from current ploc
            P_loc{t+1,pi}=cat(2,P_loc{t+1,pi},pl);
            P_idx{t+1,pi}=cat(2,P_idx{t+1,pi},P_idx{t,i}(1,p));
            P_storage{t+1,pi}=cat(2,P_storage{t+1,pi},0);%activate
            P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},P_vol{t,i}(1,p));   
        end
        
    end %p loop
end %t loop
clear sm csm ii pi pl ai
toc
%
lnkvol=cellfun(@sum,P_vol); %sum the volume in each link
seddepth=lnkvol./Length./B;%m
% 
% sdsave(iter,:)=seddepth;
% Param(iter,8)=mean(seddepth(time>200)./(1-Lp)); %mean of SD
% Param(iter,9)=std(seddepth(time>200)./(1-Lp)); %std of SD
% 
% save('Param_sdsave4.mat','Param', 'sdsave');
% clearvars -except Param sdsave iter
%end
%%
% 
% %%
hstor=(Length.*(Sc-Si))./2.*sBiLi.*(1-Lp)./B./Length;
Hs_adj=capacity./Length./B;%m, sediment depth at adjusted capacity
figure; hold on; box on
plot(time,seddepth./(1-Lp))
%plot([PSWFmx(i,1)+mu PSWFmx(i,1)+mu], [0 0.1],'k')
plot([0 tmax],[(Hs_adj+hstor)./(1-Lp) (Hs_adj+hstor)./(1-Lp)],'k')
xlabel('Time, years'); ylabel('Bed sediment thickness, m')
% 
% %
%% %
Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Sc.^(3/2)./theta.*0.175;%m/s, realtime
STTime=Length./Velocity;%seconds, travel time through each link
STTimeyr=STTime./365./24./60./60;%years

%
Hs_adj=capacity./Length./B;%m, sediment depth at adjusted capacity
rho=STTimeyr(1).*n./1;

clear ph cph nn
nn=(0:1:800)';
ph(1:length(nn),1)=NaN;

ph(:,1)=poisspdf(nn,rho(1));
cph(:,1)=poisscdf(nn,rho(1));

hs_theo_pdf=nn.*Vp./Length./B;

% Plot histogram of bed sediment thickness
%hval=hist(numpar(:,1147),0:1:50);

%[hval, cent]=hist(seddepth(time>100,i),20);
[hval, cent]=hist(seddepth(time>200,1));%,0.005:0.01:0.205);
%[hval, cent]=hist(seddepth(time>(PSWFmx(i,1)+mu),i),15);
%[hval, cent]=hist(lnkvol(time>(PSWFmx(i,1)+mu),i)./Vp,15);
%[hval, cent]=hist(seddepth(time>100,i),15);

abw=diff(hs_theo_pdf(1:2,1));
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
plot((hs_theo_pdf(:,1)+hstor)./(1-Lp),ph(:,1).*bwsf,'k','LineWidth',3)
%plot(n,ph(:,i).*bwsf,'k','LineWidth',3)
plot([(Hs_adj+hstor)./(1-Lp) (Hs_adj+hstor)./(1-Lp)],[0 max(ph(:,i).*bwsf)],'r')
%plot(nn(123:394).*Vp./Length./B,PI(123:394).*bwsf.*.5./sum(PI(123:394)),'b','LineWidth',3)
%xlim([0 1])
xlabel('Bed sediment thickness, m')
ylabel('Probability')
%% 
% %
% %pHs=interp1(hs_theo_pdf(:,1),cph(:,i),Hs_adj(i));
% %fexHs=1-pHs;
% art=seddepth(time>100,i);%use data greater than PSWFi+mu
% SD_n=length(art); %number of interarrival times (SD)
% SD_min=min(art); %min of SD
% SD_max=max(art); %max of SD
% SD_m=mean(art); %mean of SD
% SD_med=median(art); %median of SD
% SD_std=std(art); %std of SD
% SD_fexHs=sum(art>Hs_adj+hstor)./length(art);
% 
% output(1,1)=SD_min;
% output(2,1)=SD_med;
% output(3,1)=SD_m;
% output(4,1)=SD_max;
% output(5,1)=SD_std;
% output(6,1)=Hs_adj;
% output(7,1)=hstor;
% output(8,1)=Hs_adj+hstor;
% output(9,1)=SD_fexHs;
% output(10,1)=capacity;
% output(11,1)=STTimeyr;
% output(12,1)=capacity./STTimeyr;
% output(13,1)=sBiLi;
% %%
% Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
% STTime=Length./Velocity;%seconds, travel time through each link
% STTimeyr=STTime./365./24./60./60;%years
% %%
% figure
% hold on
% plot(art,'.-b')
% find(art>Hs_adj+hstor);
% plot(ans,art(ans),'.r')
% find(art>Hs_adj+hstor);
% plot(ans,art(ans),'.g')
% %%
% gc=(art>Hs_adj+hstor);
% 
% lastidx=find(diff(gc)~=0);
% zo=gc(lastidx+1);
% zo(end)=[];
% intar=diff(lastidx);
% 
% tgc=intar(zo);
% tlc=intar(~zo);    
% %%
% figure; hold on; box on
% plot(hs_theo_pdf(:,i),1-cph(:,i))
% plot([Hs_adj(i) Hs_adj(i)],[0 1],'r')
% %% Plot numerical SD timeseries
% figure; hold on; box on
% plot(time,seddepth)
% %plot([PSWFmx(i,1)+mu PSWFmx(i,1)+mu], [0 0.1],'k')
% plot([0 tmax],[Hs_adj(i) Hs_adj(i)],'k')
% %ylim([0 1])
% xlabel('Time, years'); ylabel('Bed sediment thickness, m')
% %%
% hstor=(Length.*(Sc-Si))./2.*sBiLi.*(1-Lp)./B./Length;
% (hstor+Hs_adj).*Length.*B
% %%
% clear Velocity STTime Vss Sn STTimeyr
% Vss(1:4000,1)=0;
% Sn(1:4000,1)=NaN;
% STTimeyr(1:4000,1)=NaN;
% Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Si.^(3/2)./theta.*0.175;%m/s, realtime
% STTime=Length./Velocity;%seconds, travel time through each link
% STTimeyr(1)=STTime./365./24./60./60;%years
% for i=1:4000
% Vss(i)=((Vp.*n./mu)-(capacity./STTimeyr(i))).*dt./60./60./24./365;
% Sn(i)=Si+(2.*sum(Vss))./sBiLi./Length./(1-Lp);
% Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Sn(i).^(3/2)./theta.*0.175;%m/s, realtime
% STTime=Length./Velocity;%seconds, travel time through each link
% STTimeyr(i+1)=STTime./365./24./60./60;%years
% end
% figure
% plot((1:4000).*dt./60./60./24./365,Sn)%cumsum(Vss))
% %%
% figure
% plot(seddepth(1:end-1),STTimeyr(1:end-1),'.b')
% %%
% IN(1:601,1)=0;
% OUT(1:601,1)=0;
% for t=0:600
% IN(t+1,1)=sum(and(InArrival>=t,InArrival<t+1));
% OUT(t+1,1)=sum(and(OutArrival>=t,OutArrival<t+1));
% end
% %%
% IN(1:timesteps,1)=0;
% OUT(1:timesteps,1)=0;
% for t=1:timesteps-1
% IN(t+1,1)=sum(and(InArrival>=time(t),InArrival<time(t+1)));
% OUT(t+1,1)=sum(and(OutArrival>=time(t),OutArrival<time(t+1)));
% end
% %%
% %%
% Hs_adj=capacity./Length./B;%m, sediment depth at adjusted capacity
% 
% nn=0:1:1000;
% Pi(1:length(nn),1:length(nn))=0;
% for i=1:length(Pi)
%     i
%     if ((i-1).*Vp)>capacity
%         vstor=(((i-1).*Vp)-capacity)./(1-Lp);%m3
%         Sp=Si+(2.*vstor)./sBiLi./Length;
%         Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Sp.^(3/2)./theta.*0.175;%m/s, realtime
%     else
%         Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Si.^(3/2)./theta.*0.175;%m/s, realtime
%     end
%     STTime=Length./Velocity;%seconds, travel time through each link
%     STTimeyr=STTime./365./24./60./60;%years
%     rho=dt./60./60./24./365.*n./1;
%        
%     %dt or STTimeyr for rho
%     %+ vs -
%     
%     ph=poisspdf(nn,rho);
%     Pi(i,i:end)=ph(1,1:end-i+1);
%     
%     (dt./60./60./24./365)/STTimeyr
%     
% end
% 
% hs_theo_pdf=nn.*Vp./Length./B;
% %%
% %%
% nn=0:1:400;
% pin(1:length(nn),1:length(nn))=0;
% pout(1:length(nn),1:length(nn))=0;
% 
% for i=1:length(nn)
%     i
%     if ((i-1).*Vp)>capacity
%         vstor=(((i-1).*Vp)-capacity)./(1-Lp);%m3
%         Sp=Si+(2.*vstor)./sBiLi./Length;
%         V=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Sp.^(3/2)./theta.*0.175;%m/s, realtime
%     else
%         V=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Si.^(3/2)./theta.*0.175;%m/s, realtime
%     end
%     STT=Length./V;%seconds, travel time through each link
%     STTyr=STT./365./24./60./60;%years
%     rho_in=dt./60./60./24./365.*n./1;
%     rho_out=dt./60./60./24./365.*capacity./Vp./STTyr;
%     
%     phi=poisspdf(nn,rho_in);
%     pho=poisspdf(nn,rho_out);
%     
%     pin(i,i:end)=phi(1:end-i+1);
%     pout(i,1:i)=fliplr(pho(1:i));
% end
% % figure; hold on
% % plot(phi,'b');plot(pho,'r')
% %%
% pin(1:length(nn),1:length(nn))=0;
% pout(1:length(nn),1:length(nn))=0;
% for i=1:11
%    pin(10+i,10+i:end)=phi(1:end-i+1);
%    pout(10+i,i:10+i)=fliplr(pho(1:end));
% end
% %%
% ps=pin*pout;
% clear A
% A=ps-eye(length(ps));
% A(:,end)=1;
% Ainv=inv(A);
% PI=Ainv(end,:);
% figure
% plot(nn.*Vp./Length./B,PI)
% %%
% figure
% plot(nn.*Vp./Length./B,cumsum(PI))
% %%
% PIp=PI(123:394);
% %%
% figure
% %plot(seddepth,capacity./STTimeyr./B./Length,'.b')
% %plot(seddepth,capacity.*Velocity./dt./B,'.b')
% %plot(seddepth.*B.*Length./Vp,capacity./STTimeyr./Vp,'.b')
% plot(seddepth.*B.*Length./Vp,capacity./Vp./STTimeyr./365./24./60./60.*dt,'.b')
% xlabel('# parcels in system');ylabel('Output rate,# parcels per dt')
% %plot(seddepth.*B.*Length,capacity./STTimeyr,'.b')
% %%
% f1=figure; a1=axes; hold on; box(a1);
% xlabel(a1,'Relative distance through link')
% ylabel(a1,'Bed-sediment depth, m')
% ylim(a1,[0 0.16])
% xlim(a1,[0 1])
% 
% writerObj = VideoWriter('lnkstrdyn2','MPEG-4');
% writerObj.FrameRate = 50;
% writerObj.Quality =100;
% open(writerObj);
% set(gcf,'Renderer','Painters');
% 
% for t=1:timesteps-1;   
%     clear yd
%     yd=1:length(P_loc{t,1});
%     yd=yd.*Vp./Length./B;
%     
%     p1=plot(a1,P_loc{t,1}(logical(P_storage{t,1})),...
%         yd(logical(P_storage{t,1})),'ok','MarkerSize',5);
%     p2=plot(a1,P_loc{t,1}(P_storage{t,1}==0),...
%         yd(P_storage{t,1}==0),'ob','MarkerSize',5);
%     t1=title(a1,{['Time = ',num2str(fix(time(t))),' years ',...
%         num2str(round(mod(time(t)*365,365))),' days']});
%     
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
%     delete(p1)
%     delete(p2)
%     delete(t1)
%     
% end
% 
% close(writerObj);
