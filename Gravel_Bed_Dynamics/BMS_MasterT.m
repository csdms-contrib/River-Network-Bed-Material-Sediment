%% BMS_Master
% This program explores bed-material sediment dynamics on river networks. 
% This file is the main file that drives the entire program.

% Jon Czuba
% February 17, 2015

%% Functions/Files Used


%% Variables Used
% t - time index
% i - space (link) index
% p - parcel index
% timesteps - number of timesteps
% LinkNum - number of links in network

%% Import Network, Initialize Variables
clear all
close all
clc

% Import a network
% Determine connectivity of links
%load('BE_NHD_MartinLakes4.mat');%load river network with lakes
%load('BE_NHD_MartinLakes2_LKE.mat');%load lake attributes

% Tushar
load Tushar\TusharPreProcData.mat
boundary = shaperead('Tushar\Data\112717_update\ClearCreekWatershedBndy.shp');
fire = shaperead('Tushar\Data\112717_update\Twitchell_Fire_Boundary.shp');
network = shaperead('Tushar\Tushar_Network.shp');
burn=[network.Burn]';
% Methow
% load Methow\MethowPreProcData.mat

% Nisqually
% load C:\Users\jczuba\Documents\Projects\Nisqually\MATLAB\Nisqually_NHDslope.mat

%%
% Initialize variables
daystp=1;%182.5;%18.25; %number of days per timestep
dt = 1*24*60*60*daystp; %seconds in daystp number of days

% hrstp=1; %number of hours per timestep
% dt = 1*60*60; %seconds in hrstp

% Tushar
timesteps=2557;%1912;%14610;%1912*6;%365*2;%12001;%12001;%4001;%4000;%2115; %timesteps to simulate sand
% Methow
%this runs 5 years
%timesteps=1826;
% now it runs for 30 days so you get an idea of how the model works
%timesteps=365;

tmax=dt*(timesteps-1)/60/60/24/365; %years, max time of simulation
time=(0:daystp/365:tmax)'; %years, time of current step
% tmax=dt*(timesteps-1)/60/60/24; %days, max time of simulation
% time=(0:hrstp/24:tmax)'; %days, time of current step


theta=0.064;%0.1;%2D90, 64mm, units in meters
g=9.81;%m/s2
R=1.65;
rho=1000;%kg/m3
%D=0.0004;%m
%D=0.02;%m

BMS_FlowT


% Determine flow and its scaling
%B=(0.0238).*(usarea).^(0.3397);%m, width at Q2
%H=(0.0029).*(usarea).^(0.2943);%m, water depth
%U=(0.1976).*(usarea).^(0.0679);%m/s, water velocity

OutArrival=[];
OutVol(1:timesteps,1)=0;

%compute capacity for all links
capacity=Btmax.*Length.*theta;%m3
% caplim=50;%m3, lower limit for transport capacity
% capacity(capacity<caplim)=caplim;%m3, set lower limit for transport capacity
% for i=1:LinkNum
%     if Length(i,1)<300
%         %set capacity to max of current capacity, directly US capacity, or 100 m3
%         capacity(i,1)=max([capacity(i,1);capacity(find(Connect(:,2)==i),1);100]);
%     end
% end
% clear i

% H=((0.0029).*(usarea).^(0.2943))./2;%Use 1/2 ~bankfull, **UPDATE**
% Velocity=0.0354.*(usarea_km).^0.3543;%Use arbitrary for Q velocity, **UPDATE**
% REM=200/60/60;%mg/m2/hr->mg/m2/sec
% kd=0.00001;%0.1*24*60*60;%1/day->1/sec;

%remove the slope effect of the dam
%1259 connects to 1278 and together span the dam
% Slope(1259,1)=(mxelev(1259)-mnelev(1278))./(Length(1259)+Length(1278));
% Slope(1278,1)=Slope(1259,1);

Slope(Slope<0.001)=0.001;
%modify elevations so they are consistent with adjusted slopes
mxelevmod(1:LinkNum,1)=NaN;
for i = 1:LinkNum
    idx=fliplr(Connect(i,~isnan(Connect(i,:))));
    for j=1:length(idx)
        if idx(j)==OutletLinkID
            mxelevmod(idx(j),1)=mnelev(OutletLinkID)+Slope(idx(j)).*Length(idx(j));
        else
            mxelevmod(idx(j),1)=mxelevmod(idx(j-1),1)+Slope(idx(j)).*Length(idx(j));
        end
    end
    clear idx
end
clear i j

Elev=repmat(mxelevmod',timesteps,1);
% %Slope - slope at current timestep
% %slope - original slope
% slope=Slope;
% %Slope=repmat(slope',timesteps,1);


%Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
%STTime=Length./Velocity;%seconds, travel time through each link


Lp=0.25;%0.21;%0.4; %porosity
%Lp=0.13+0.21/(((D*1000)+0.002)^0.21);%Wu and Wang 2006
Lake(1:LinkNum,1)=0;
%Ppulse(1:timesteps,1:LinkNum)=NaN;
%Ppvol(1:timesteps,1:LinkNum)=NaN;

lnkvol(1:timesteps,1:LinkNum)=NaN;
Dg(1:timesteps,1:LinkNum)=NaN;

%% Assign inputs
BMS_Inputs

%% Map input structure
% so inputs can always be added on top of what is in the system

P_idx_IN=P_idx;
P_loc_IN=P_loc;
P_storage_IN=P_storage;
P_vol_IN=P_vol;
P_d_IN=P_d;
P_tt_IN=P_tt;

P_idx=cell(timesteps,LinkNum);
P_loc=cell(timesteps,LinkNum);
P_storage=cell(timesteps,LinkNum);
P_vol=cell(timesteps,LinkNum);
P_d=cell(timesteps,LinkNum);
P_tt=cell(timesteps,LinkNum);

%% Run Model

tic

for t = 1:timesteps-1 %step through each timestep
    t 
    % Generate autogenic inputs here? If desired.
    
    % Add new inputs at given timestep on top of what is already in system
    % background inputs
    for i=1:LinkNum
        P_loc{t,i}=cat(2,P_loc{t,i},P_loc_IN{t,i});
        P_idx{t,i}=cat(2,P_idx{t,i},P_idx_IN{t,i});
        P_storage{t,i}=cat(2,P_storage{t,i},P_storage_IN{t,i});
        P_vol{t,i}=cat(2,P_vol{t,i},P_vol_IN{t,i});
        P_d{t,i}=cat(2,P_d{t,i},P_d_IN{t,i});
        P_tt{t,i}=cat(2,P_tt{t,i},P_tt_IN{t,i});      
    end
    
%     %input pulse here
%     if t==2001 %add pulse at 100 yrs
%         for i=1:LinkNum
%             P_loc{t,i}=cat(2,P_loc{t,i},P_loc_PIN{t,i});
%             P_idx{t,i}=cat(2,P_idx{t,i},P_idx_PIN{t,i});
%             P_storage{t,i}=cat(2,P_storage{t,i},P_storage_PIN{t,i});
%             P_vol{t,i}=cat(2,P_vol{t,i},P_vol_PIN{t,i});
%         end
%         clear P_loc_PIN P_idx_PIN P_storage_PIN P_vol_PIN
%     end
    
    % Determine capacity of links, active parcels, and adjust bed slopes
    % based on storage volume
    BMS_CapacitySlopeT

    %Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
    %STTime=Length./Velocity;%seconds, travel time through each link
    %STTime -> P_tt{t,i}(p)

    for i = 1:LinkNum %step through each link
        
        % if current link is empty go to the next link
        if isempty(P_loc{t,i})
            continue
        end
        
%         if Lake(i) % if a lake, do not move sediment
%             continue
%         end
        
        %capacity and slope here? no, needs to be above i loop because
        %slopes are updated from more than just its own link
        
        for p = 1:length(P_loc{t,i})
            
            % if parcel is inactive then update and go to the next parcel
            if P_storage{t,i}(1,p)
                P_loc{t+1,i}=cat(2,P_loc{t+1,i},P_loc{t,i}(1,p));%keep loc
%                 if P_loc{t,i}(1,p)<0
%                     1;
%                 end
                P_idx{t+1,i}=cat(2,P_idx{t+1,i},P_idx{t,i}(1,p));%keep idx
                P_storage{t+1,i}=cat(2,P_storage{t+1,i},0);% make active as 
                %this is determined for each timestep
                P_vol{t+1,i}=cat(2,P_vol{t+1,i},P_vol{t,i}(1,p));%keep vol
                P_d{t+1,i}=cat(2,P_d{t+1,i},P_d{t,i}(1,p));%keep d
                P_tt{t+1,i}=cat(2,P_tt{t+1,i},0);% make 0 as 
                %this is determined for each timestep
                continue
            end
            
            % Enter long term storage here?
            
            % Move parcels downstream
            sm=[];%time to move through current and downstream links
            csm=[];%cumulative time to move through current and ds links
            sm=P_tt{t,i}(p)*(1-P_loc{t,i}(1,p));%seconds, time to move out of current link
            csm=sm;%seconds
            ii=i;
            while csm(end)<=dt %loop until cumulative time toward outlet exceeds timestep
                ii=Connect(ii,2); %set ds link
                if isnan(ii) %if downstream is the outlet
                    %leave outlet before dt ends
                    sm=cat(1,sm,NaN);
                    csm=cat(1,csm,sum(sm));
                    OutVol(t,1)=OutVol(t,1)+P_vol{t,i}(1,p);
                    break
                end
%                 if Lake(ii) %if move into a lake
%                     sm=cat(1,sm,NaN);
%                     csm=cat(1,csm,sum(sm));                   
%                    break 
%                 end
                %sm=cat(1,sm,STTime(ii,1));%add time through next link
                %EDIT IN THE FUTURE
                sm=cat(1,sm,P_tt{t,i}(p)./Length(i).*Length(ii));%travel at same velocity but in ds link
                csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
            end
            % update parcel location
            if ~isnan(csm(end)) %check to make sure parcel is still in the system
                pi=Connect(i,length(csm)); %parcel link
                if length(csm)==1 %still in same link
                    pl=(P_tt{t,i}(p)*P_loc{t,i}(1,p)+dt)/P_tt{t,i}(p);%update location from current ploc
                else %moved to a ds link
                    pl=1-((csm(end)-dt)/P_tt{t,i}(p));%update location from beginning of link
                    if pl<1 %overshooting the next link, UPDATE IN FUTURE
                        pl=1;
                    end
                    %note csm(end)-dt computes time remaining to move
                    %through the rest of the link
                end
                              
                P_loc{t+1,pi}=cat(2,P_loc{t+1,pi},pl);
%                 if pl<0
%                     1;
%                 end
                P_idx{t+1,pi}=cat(2,P_idx{t+1,pi},P_idx{t,i}(1,p));
                P_storage{t+1,pi}=cat(2,P_storage{t+1,pi},0);%activate
                P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},P_vol{t,i}(1,p));
                P_d{t+1,pi}=cat(2,P_d{t+1,pi},P_d{t,i}(1,p));
                P_tt{t+1,pi}=cat(2,P_tt{t+1,pi},0);%set to 0
%                P_mass{t+1,pi}=cat(2,P_mass{t+1,pi},P_mass{t,i}(1,p).*exp(-kd*dt));
            end
%             %assign arrival times
%             for ai=1:length(csm)-1 %as long the parcel moved to a new link
%                 if isnan(Connect(i,ai+1)) %left the system
%                     OutArrival=cat(1,OutArrival,...
%                         time(t,1)+csm(ai,1)/60/60/24/365);
%                 else %arrival times at each link
%                     L_arrival{Connect(i,ai+1),1}=cat(2,L_arrival{Connect(i,ai+1),1},...
%                         time(t,1)+csm(ai,1)/60/60/24/365);
%                 end
%             end
            
        end %p loop
    end %i loop 
    
%     if t==2001
%         lnkvol=cellfun(@sum,P_vol); %sum the volume in each link
%         numpar=cellfun(@length,P_loc); %compute number of parcel in each link
%         seddepth=lnkvol./repmat(Length',timesteps,1)./repmat(B',timesteps,1);%m
%         save('BE_NHD_MartinLakes4_BRU_baseline_100yr_v1.mat','seddepth','lnkvol','Elev','OutVol','OutArrival','L_arrival','-v7.3');
%         clear seddepth
%     end

lnkvol(t,:)=cellfun(@sum,P_vol(t,:));
% 
% for i=1:LinkNum
%     Ppulse(t,i)=sum(P_idx{t,i}>=ppidxstart);
%     Ppvol(t,i)=sum(P_vol{t,i}(P_idx{t,i}>=ppidxstart));
% end
% 
% %clear contents to conserve space
% for i=1:LinkNum
%     P_idx{t,i}=[];
%     P_loc{t,i}=[];
%     P_storage{t,i}=[];
%     P_vol{t,i}=[];
%     
%     P_idx_IN{t,i}=[];
%     P_loc_IN{t,i}=[];
%     P_storage_IN{t,i}=[];
%     P_vol_IN{t,i}=[];
% end

end %t loop

clear sm csm ii pi pl ai

toc
%%
%
% lnkmass=cellfun(@sum,P_mass); %sum the mass in each link
% lnkconc=lnkmass./(repmat((Length.*B.*H)',timesteps,1)).*1000;%mg/L
numpar=cellfun(@length,P_loc); %compute number of parcel in each link

%
%lnkvol=cellfun(@sum,P_vol); %sum the volume in each link
%numpar=cellfun(@length,P_loc); %compute number of parcel in each link
%B=(0.0238).*(usarea).^(0.3397);
%numparconc=numpar./repmat(Length',timesteps,1)./repmat(B',timesteps,1).*1000;
lnkvol(t+1,:)=cellfun(@sum,P_vol(t+1,:));
seddepth=lnkvol./repmat(Length',timesteps,1)./repmat(Btmax',timesteps,1);%m

%save MethowOutput.mat

%sum(sum(cellfun(@sum,P_storage),2))
%
%save('TusharOutput.mat','lnkvol','Dg','P_loc','P_vol','P_d','-v7.3');

%save('TusharOutput14.mat','-v7.3');

%%
Fs(1:timesteps,1:LinkNum)=NaN;
for t=1:timesteps
    for i=1:LinkNum
        if isempty(P_loc{t,i})
            continue
        end
        actidx=P_storage{t,i}==0;
        actvol=sum(P_vol{t,i}(actidx));
        
        actsandidx=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(1));
        actsandvol=sum(P_vol{t,i}(actsandidx));
        Fs(t,i)=actsandvol./actvol;
    end
end
%%
figure; hold on; box on
plot((1:length(Qgage))./365,Qgage)
%set(gca,'YScale','log')
ylabel('Daily streamflow, m^3/s')
xlabel('Time since debris flow input, years')
xlim([0 7])

%%
% Methow: 61, 689, 444
%%
figure; hold on; box on
%plot((1:length(Qgage))./365,seddepth,'b')
%set(gca,'YScale','log')
%ylabel('Total sediment depth within a link, m')
xlabel('Time since debris flow input, years')
%xlim([0 5])
plot((1:length(Qgage))./365,seddepth(:,316),'b')
plot((1:length(Qgage))./365,seddepth(:,438),'r')
plot((1:length(Qgage))./365,seddepth(:,416),'k')
% plot((1:length(Qgage))./365,seddepthdf(:,161),'b')
% plot((1:length(Qgage))./365,seddepthdf(:,380),'r')
% plot((1:length(Qgage))./365,seddepthdf(:,424),'k')
%plot((1:length(Qgage))./365,seddepth(:,267),'k')
%%
BMS_PlotParcels
%%
figure; hold on; box on
plot((1:length(Qgage))./365,seddepthdf)
set(gca,'YScale','log')
ylabel('Sediment depth within a link from DF inputs, m')
xlabel('Time since debris flow input, years')
%xlim([0 5])
%%
figure; hold on; box on
%plot((1:length(Qgage))./365,Dg)
%set(gca,'YScale','log')
ylabel('D50 within a link, m')
xlabel('Time since debris flow input, years')
%xlim([0 5])
plot((1:length(Qgage))./365,Dg(:,316),'b')
plot((1:length(Qgage))./365,Dg(:,438),'r')
plot((1:length(Qgage))./365,Dg(:,416),'k')
%plot((1:length(Qgage))./365,Dg(:,477),'k')

%%
figure; hold on; box on
%plot((1:length(Qgage))./365,Fs)
%set(gca,'YScale','log')
ylabel('Sand fraction')
xlabel('Time since debris flow input, years')
%xlim([0 5])
plot((1:length(Qgage))./365,Fs(:,316),'b')
plot((1:length(Qgage))./365,Fs(:,438),'r')
plot((1:length(Qgage))./365,Fs(:,416),'k')
%plot((1:length(Qgage))./365,Fs(:,477),'k')

%%
Dist(1:LinkNum,1)=NaN;
for i=1:LinkNum
    lastConn=find(Connect(i,:)==OutletLinkID);
    Dist(i,1)=sum(Length(Connect(i,1:lastConn)),1);
end
%% plot long profile of given value

loc=3;%middle
%loc=76;%right
figure; hold on; box on
color=['b','r','k'];
yr=[0,1,7];
for j=1:3
t=1+365*yr(j);

index=Connect(loc,~isnan(Connect(loc,:)));

%y=seddepthdf(t,index);
%[y,yt]=max(seddepthdf(:,index),[],1);
%y=time(yt);
y=Slope(index);
%y=max(Fs(:,index),[],1);
%y=min(Dg(:,index),[],1);
%y=Dg(t,index);%
%y(isnan(y))=0;
%y(seddepth(t,index)<0.1)=NaN;

%x=cumsum(Length(index))./1000;%km
x=Dist(index)./1000;%km

xx=[];
yy=[];
% xx=cat(1,xx,0,x(1));
% yy=cat(1,yy,y(1),y(1));
% for i=2:length(x)
%     xx=cat(1,xx,x(i-1:i,1));
%     yy=cat(1,yy,y(i),y(i));
% end
for i=1:length(x)-1
    xx=cat(1,xx,x(i:i+1,1));
    yy=cat(1,yy,y(i),y(i));
end

plot(xx,yy,color(j));

end

xlabel('Distance upstream from basin outlet, km')
legend('0 years','1 years','7 years')

%%
figure; hold on; box on
plot(Dist(index)./1000,mxelev(index))
xlabel('Distance upstream from basin outlet, km')
%%
figure; hold on; box on
plot(Dist(index)./1000,Slope(index))
xlabel('Distance upstream from basin outlet, km')
%% vol by grainsize

lnkvoldfgs(1:timesteps,1:LinkNum,1:gsclass)=NaN;
for j=1:gsclass
    for i=1:LinkNum
        for tt=1:timesteps
            lnkvoldfgs(tt,i,j)=sum(P_vol{tt,i}(and(P_idx{tt,i}>=pidxdf,P_d{tt,i}==Dpsd(j))));
        end
    end
end
%seddepthdf=lnkvoldf./repmat(Length',timesteps,1)./repmat(Btmax',timesteps,1);%m

volinbasin=sum(lnkvoldfgs,2);
totalvolinbasin=sum(volinbasin,3);
%%
figure; hold on; box on
for i=1:gsclass
plot((1:length(Qgage))./365,volinbasin(:,1,i)./volinbasin(1,1,i))
end
ylabel('Fraction of DF input within watershed')
xlabel('Time since debris flow input, years')
legend(num2str(Dpsd))
plot((1:length(Qgage))./365,totalvolinbasin(:,1,1)./totalvolinbasin(1,1,1),'k')
%%
% figure
% plot(volinbasin(:,1,7))
% %%
% volinbasin(end,1,:)./volinbasin(1,1,:)
%%
lnkvoldfgs_burn=lnkvoldfgs;
lnkvoldfgs_burn(:,burn==0,:)=0;
volinburn=sum(lnkvoldfgs_burn,2);
%%
figure; hold on; box on
for i=1:gsclass
plot((1:length(Qgage))./365,volinburn(:,1,i)./volinburn(1,1,i))
end
ylabel('Fraction of DF input within burned area')
xlabel('Time since debris flow input, years')
legend(num2str(Dpsd))
%%
% volinburn(end,1,:)./volinburn(1,1,:)

%% surface grain size distribution
t=1
i=416%359%205%424
Dft=[2,4,8,32,64,128,256,512]'./1000;%m
cpsd(1:timesteps,1:gsclass,1)=NaN;
for t=1:timesteps
%for i=1:LinkNum
if isempty(P_loc{t,i})
    continue
end

actidx=P_storage{t,i}==0;
actvol=sum(P_vol{t,i}(actidx));

for j=1:gsclass
    actjidx=and(P_storage{t,i}==0,P_d{t,i}<=Dpsd(j));
    actftjvol=sum(P_vol{t,i}(actjidx));
    cpsd(t,j)=actftjvol./actvol;
end
end
%
figure; hold on; box on
jj=1:10:(timesteps-7);
ci=linspace(0.8,0.2,timesteps); 
for jjj=1:length(jj)
j=jj(jjj);
plot(Dft,cpsd(j,:),'Color',[ci(j) ci(j) ci(j)])
end
plot(Dft,cpsd(1+365*0,:),'b','LineWidth',3)
plot(Dft,cpsd(1+365*1,:),'r','LineWidth',3)
plot(Dft,cpsd(1+365*7,:),'k','LineWidth',3)
set(gca,'XScale','log')
ylim([0 1])

inisspsd=cumsum(Fsspsd);
inisfpsd=cumsum(Fsfpsd);
inidfpsd=cumsum(Fdfpsd);
plot(Dft,inisspsd,'g','LineWidth',2)
plot(Dft,inisfpsd,'g','LineWidth',2)
plot(Dft,inidfpsd,'g','LineWidth',2)
title(num2str(i))

%% contour
figure; hold on; box on
loc=76;%middle
index=Connect(loc,~isnan(Connect(loc,:)));
y=seddepthdf(:,index);%
x=cumsum(Length(index))./1000;%km
contourf(repmat(x',timesteps,1),repmat(time,1,length(x)),(y),'LineStyle','none')
ylabel('Time since debris flow input, years')
xlabel('Distance downstream, km')

%% parcel travel distances and travel times
DF_idx=[];
DF_d=[];
DF_loc0=[];

for i=1:LinkNum

    dfcellidx=find(P_idx{1,i}>=pidxdf)';
    if isempty(dfcellidx)
        continue
    end
    DF_idx=cat(1,DF_idx,P_idx{1,i}(dfcellidx)');
    DF_d=cat(1,DF_d,P_d{1,i}(dfcellidx)');
    DF_loc0=cat(1,DF_loc0,repmat(i,length(dfcellidx),1));
    
end
%%
DF_loc7(1:length(DF_idx),1)=NaN;
DF_tt(1:length(DF_idx),1)=NaN;
DF_ll(1:length(DF_idx),1)=NaN;

for dfp=1:length(DF_idx)
dfp
index=Connect(DF_loc0(dfp),~isnan(Connect(DF_loc0(dfp),:)));

t=1;
i=1;
while isnan(DF_loc7(dfp))
    if ~isempty(find(DF_idx(dfp)==P_idx{t,index(i)},1))
        t=t+1;
    else
        i=i+1;
    end

    if i>length(index)
        DF_loc7(dfp)=index(i-1);
        DF_tt(dfp)=t;      
    end
    
    if t>timesteps
        DF_loc7(dfp)=index(i);
        DF_tt(dfp)=t-1;
    end
    
    if i>length(index) && t<timesteps
        DF_loc7(dfp)=index(i-1);
        DF_tt(dfp)=t-1;
    end
end
tosum=Connect(DF_loc0(dfp),1:find(Connect(DF_loc0(dfp),:)==DF_loc7(dfp),1));
if numel(tosum)==1
    DF_ll(dfp)=Length(tosum(1)).*...
        P_loc{DF_tt(dfp),DF_loc7(dfp)}(find(P_idx{DF_tt(dfp),DF_loc7(dfp)}==DF_idx(dfp),1));
elseif isempty(find(P_idx{DF_tt(dfp),DF_loc7(dfp)-1}==DF_idx(dfp),1))
    DF_ll(dfp)=sum(Length(tosum(1:end)));
else
    DF_ll(dfp)=sum(Length(tosum(1:end-1)))+Length(tosum(end)).*...
        P_loc{DF_tt(dfp),DF_loc7(dfp)}(find(P_idx{DF_tt(dfp),DF_loc7(dfp)-1}==DF_idx(dfp),1));
end
end
%%
i=5;
figure
hist(log10(DF_ll(and(DF_d==Dpsd(i),DF_ll>0))),30)
xlim([-5 5])
title(num2str(Dpsd(i)))

%%
figure; hold on; box on
for i=1:gsclass
x=sort(DF_ll(and(DF_d==Dpsd(i),DF_ll>0)));
y=(1:length(x))'./length(x);
plot(x,y)
end
set(gca,'XScale','log')
legend(num2str(Dpsd))

%%

DF_tt(DF_tt<timesteps)
DF_d(DF_tt<timesteps)

%%
figure
plot(DF_ll,DF_tt,'.b')

