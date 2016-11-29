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
load('BE_NHD_MartinLakes4.mat');%load river network with lakes
%load('BE_NHD_MartinLakes2_LKE.mat');%load lake attributes
%%
% Initialize variables
daystp=182.5;%18.25; %number of days per timestep
dt = 1*24*60*60*daystp; %seconds in daystp number of days

% hrstp=1; %number of hours per timestep
% dt = 1*60*60; %seconds in hrstp

timesteps=401;%12001;%12001;%4001;%4000;%2115; %timesteps to simulate sand
tmax=dt*(timesteps-1)/60/60/24/365; %years, max time of simulation
time=(0:daystp/365:tmax)'; %years, time of current step
% tmax=dt*(timesteps-1)/60/60/24; %days, max time of simulation
% time=(0:hrstp/24:tmax)'; %days, time of current step

theta=0.1;
D=0.0004;%m

% Determine flow and its scaling
B=(0.0238).*(usarea).^(0.3397);%m, width at Q2
H=(0.0029).*(usarea).^(0.2943);%m, water depth
U=(0.1976).*(usarea).^(0.0679);%m/s, water velocity

OutArrival=[];
OutVol(1:timesteps,1)=0;

%compute capacity for all links
capacity=B.*Length.*theta.*H;%m3
caplim=50;%m3, lower limit for transport capacity
capacity(capacity<caplim)=caplim;%m3, set lower limit for transport capacity
for i=1:LinkNum
    if Length(i,1)<300
        %set capacity to max of current capacity, directly US capacity, or 100 m3
        capacity(i,1)=max([capacity(i,1);capacity(find(Connect(:,2)==i),1);100]);
    end
end
clear i

% H=((0.0029).*(usarea).^(0.2943))./2;%Use 1/2 ~bankfull, **UPDATE**
% Velocity=0.0354.*(usarea_km).^0.3543;%Use arbitrary for Q velocity, **UPDATE**
% REM=200/60/60;%mg/m2/hr->mg/m2/sec
% kd=0.00001;%0.1*24*60*60;%1/day->1/sec;

%remove the slope effect of the dam
%1259 connects to 1278 and together span the dam
Slope(1259,1)=(mxelev(1259)-mnelev(1278))./(Length(1259)+Length(1278));
Slope(1278,1)=Slope(1259,1);

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
%Slope - slope at current timestep
%slope - original slope
slope=Slope;
%Slope=repmat(slope',timesteps,1);

Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
STTime=Length./Velocity;%seconds, travel time through each link

lnkvol(1:timesteps,1:LinkNum)=NaN;
Lp=0.4; %porosity
%Lp=0.13+0.21/(((D*1000)+0.002)^0.21);%Wu and Wang 2006

%Ppulse(1:timesteps,1:LinkNum)=NaN;
%Ppvol(1:timesteps,1:LinkNum)=NaN;

%% Assign inputs
BMS_Inputs

%% Map input structure
% so inputs can always be added on top of what is in the system

P_idx_IN=P_idx;
P_loc_IN=P_loc;
P_storage_IN=P_storage;
P_vol_IN=P_vol;

P_idx=cell(timesteps,LinkNum);
P_loc=cell(timesteps,LinkNum);
P_storage=cell(timesteps,LinkNum);
P_vol=cell(timesteps,LinkNum);

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
    BMS_CapacitySlope

    Velocity=0.05./sqrt(9.81)./1.65./1.65./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta.*0.175;%m/s, realtime
    STTime=Length./Velocity;%seconds, travel time through each link

    for i = 1:LinkNum %step through each link
        
        % if current link is empty go to the next link
        if isempty(P_loc{t,i})
            continue
        end
        
        if Lake(i) % if a lake, do not move sediment
            continue
        end
        
        %capacity and slope here? no, needs to be above i loop because
        %slopes are updated from more than just its own link
        
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
            
            % Enter long term storage here?
            
            % Move parcels downstream
            sm=[];%time to move through current and downstream links
            csm=[];%cumulative time to move through current and ds links
            sm=STTime(i,1)*(1-P_loc{t,i}(1,p));%seconds, time to move out of current link
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
                if Lake(ii) %if move into a lake
                    sm=cat(1,sm,NaN);
                    csm=cat(1,csm,sum(sm));                   
                   break 
                end
                sm=cat(1,sm,STTime(ii,1));%add time through next link
                csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
            end
            % update parcel location
            if ~isnan(csm(end)) %check to make sure parcel is still in the system
                pi=Connect(i,length(csm)); %parcel link
                if length(csm)==1 %still in same link
                    pl=(STTime(i,1)*P_loc{t,i}(1,p)+dt)/STTime(i,1);%update location from current ploc
                else %moved to a ds link
                    pl=1-((csm(end)-dt)/STTime(pi,1));%update location from beginning of link
                    %note csm(end)-dt computes time remaining to move
                    %through the rest of the link
                end
                              
                P_loc{t+1,pi}=cat(2,P_loc{t+1,pi},pl);
                P_idx{t+1,pi}=cat(2,P_idx{t+1,pi},P_idx{t,i}(1,p));
                P_storage{t+1,pi}=cat(2,P_storage{t+1,pi},0);%activate
                P_vol{t+1,pi}=cat(2,P_vol{t+1,pi},P_vol{t,i}(1,p));
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

% lnkvol(t,:)=cellfun(@sum,P_vol(t,:));
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
seddepth=lnkvol./repmat(Length',timesteps,1)./repmat(B',timesteps,1);%m

%sum(sum(cellfun(@sum,P_storage),2))
%
% save('BE_NHD_MartinLakes4_BRU_X.mat','Ppulse','Ppvol','seddepth',...
%     'lnkvol','Elev','OutVol','OutArrival','L_arrival','-v7.3');

