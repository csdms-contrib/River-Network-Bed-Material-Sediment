%% Inputs
% This function assigns/generates inputs.

% Jon Czuba
% February 17, 2015

%% Initialize Variables
% Parcel properties for 
P_idx=cell(timesteps,LinkNum);% index of each parcel (Source)
P_loc=cell(timesteps,LinkNum);% location of each parcel within a link (Loc)
P_storage=cell(timesteps,LinkNum);% parcel in storage (inactive)

%P_d=0.0001;%m, parcel grain size
P_d=cell(timesteps,LinkNum);%m, parcel grain size
P_vol=cell(timesteps,LinkNum);%m3, parcel volume
%pvol=1;%m3
P_tt=cell(timesteps,LinkNum);%s, parcel travel time
%P_mass=cell(timesteps,LinkNum);%mg, parcel mass

%L_arrival=cell(LinkNum,1);%years, arrival time of parcels to each link

%%
% H2gage=0.1902.*170.^0.3897;
% Agage=3.4777e+08;%m2
% H2=(H2gage./(Agage.^0.4)).*(usarea).^0.4;
% tausgage=H2gage.*Slope(519)./1.65./0.064;
% 
% initD50=H2.*Slope./1.65./tausgage;
% %%
% DbinCDF=[2,4,8,16,32,64,128,256,512,1024];%mm
% phiDbinCDF=-log10(DbinCDF)./log10(2);
% mu=-log10(64)./log10(2);
% sig=-1;
% Fcdf=normcdf(-phiDbinCDF,-mu,-sig);
% Fpsd=diff(normcdf(-phiDbinCDF,-mu,-sig));
% Dpsd=2.^-(phiDbinCDF(1:end-1)-0.5);
% Dpsd=Dpsd./1000;%m
% 
% %%
% DbinCDF=[2,4,8,16,32,64,128,256,512,1024];%mm
% phiDbinCDF=-log10(DbinCDF)./log10(2);
% mu=-log10(initD50.*1000)./log10(2);
% sig=-1;
% Fcdf=normcdf(repmat(-phiDbinCDF,LinkNum,1),repmat(-mu,1,length(DbinCDF)),repmat(-sig,LinkNum,length(DbinCDF)));
% Fpsd=diff(Fcdf,1,2);
% Dpsd=2.^-(phiDbinCDF(1:end-1)-0.5);
% Dpsd=Dpsd./1000;%m
% 
% Fpsd=Fpsd.*0.6;%reduce by 40% to accomodate sand
% Fpsd=[repmat(0.4,LinkNum,1),Fpsd];
% Dpsd=[0.002,Dpsd];%m
% Fpsd(:,1)=Fpsd(:,1)+(1-sum(Fpsd,2));
% 
% Fcdf2=cumsum(Fpsd,2);
% 
% % %% compute Dxx
% DbinCDF=[2,4,8,32,64,128,256,512];%mm
% Fcdf=[];
% Dxx=0.5;
% DXX(1:LinkNum,1)=NaN;
% for i=1:size(Fcdf,1)
%     for j=1:length(DbinCDF)-1
%         if Fcdf(i,j)<Dxx && Fcdf(i,j+1)>=Dxx
%             DXX(i,1)=exp(interp1(Fcdf(i,j:j+1),log(DbinCDF(1,j:j+1)),Dxx));
%         end
%     end
% end
% DXX(isnan(DXX))=0.002;
% 
% %%
% theta=2.*DXX./1000;
% theta(theta<0.1)=0.1;
% capacity=Btmax.*Length.*theta;%m3
% 
% %%
% % figure
% % %plot(DbinCDF,Fcdf)
% % plot(Dpsd,Fpsd)
% % set(gca,'XScale','log')
% 
% %% Input subsurface layer t=0
% % Nisqually
% %Fpsd
% %Dpsd
% gsclass=length(Dpsd);
% 
% VolIn=Btmax.*Length.*theta;%m3
% 
% %subsurface
% VolInd=repmat(VolIn,1,gsclass).*Fpsd;
% pidx=1;
% for i = 1:length(VolInd)
%     for k = 1:gsclass
%         %j=7-k+1;
%         j=k;
%         if VolInd(i,j) == 0
%             continue
%         end
%         np=ceil(VolInd(i,j)./10);
%         pvol=VolInd(i,j)./np;
%         P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
%         P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
%         pidx=pidx+np;
%         P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
%         P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
%         P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
%         P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
%     end
% end
% 
% %surface
% VolInd=repmat(VolIn,1,gsclass).*Fpsd;
% for i = 1:length(VolInd)
%     for k = 1:gsclass
%         %j=7-k+1;
%         j=k;
%         if VolInd(i,j) == 0
%             continue
%         end
%         np=ceil(VolInd(i,j)./10);
%         pvol=VolInd(i,j)./np;
%         P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
%         P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
%         pidx=pidx+np;
%         P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
%         P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
%         P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
%         P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
%     end
% end
% %%
% pidxdf=pidx;
% %df input
% VolInd=repmat(VolIn,1,gsclass).*Fpsd;
% for i = 1:length(VolInd)
%     for k = 1:gsclass
%         %j=7-k+1;
%         j=k;
%         if VolInd(i,j) == 0
%             continue
%         end
%         np=ceil(VolInd(i,j)./10);
%         pvol=VolInd(i,j)./np;
%         P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
%         P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
%         pidx=pidx+np;
%         P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
%         P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
%         P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
%         P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
%     end
% end
% 



%% Input subsurface layer t=0
% Tushar
gsclass=8;
%Dpsd=[2,6,16,45,91,181,362]'./1000;
%Dpsd=[2,3,6,16,45,91,181]'./1000;
%Dpsd=[2,4,8,32,64,128,256]'./1000;%m
Fsspsd=[0.305,0.124,0.133,0.428,0.01,0,0,0]';
Fsfpsd=[0.046,0.024,0.048,0.274,0.278,0.166,0.1,0.064]';
Fdfpsd=[0.15,0.15,0.2,0.32,0.09,0.06,0.03,0]';
%Fsfpsd=[0.112,0.09,0.114,0.34,0.344,0,0]';
%Fdfpsd=[0.168,0.168,0.218,0.338,0.108,0,0]';
%Fdfpsd=[0.33,0.33,0.33,0,0,0,0]';

% gsclass=7;
% %Dpsd=[2,6,16,45,91,181,362]'./1000;
% %Dpsd=[2,3,6,16,45,91,181]'./1000;
Dpsd=[0.5,2.83,5.66,16,45.3,90.5,181,362]'./1000;%m
% Fsspsd=[0.1,0.03,0.05,0.15,0.12,0.18,0.37]';
% Fsfpsd=[0.1,0.03,0.05,0.15,0.12,0.18,0.37]';
% Fdfpsd=[0.1,0.03,0.05,0.15,0.12,0.18,0.37]';

VolIn=Btmax.*Length.*theta;%m3

%subsurface
VolInd=repmat(VolIn.*2,1,gsclass).*repmat(Fsspsd',LinkNum,1);
pidx=1;
for i = 1:length(VolInd)
    for k = 1:7
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

pidxss2=pidx;

%subsurface2
VolInd=repmat(VolIn,1,gsclass).*repmat(Fsfpsd',LinkNum,1);
for i = 1:length(VolInd)
    for k = 1:7
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

pidxsf=pidx; % surface parcels are from pidxsf:pidxdf-1

%surface
VolInd=repmat(VolIn,1,gsclass).*repmat(Fsfpsd',LinkNum,1);
for i = 1:length(VolInd)
    for k = 1:7
        %j=7-k+1;
        j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        %P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_loc{1,i}=cat(2,P_loc{1,i},rand(1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

% pidxdf=pidx;
% %df input
% VolInd=repmat(VolIn,1,gsclass).*repmat(Fsfpsd',LinkNum,1);
% for i = 1:length(VolInd)
%     for k = 1:7
%         %j=7-k+1;
%         j=k;
%         if VolInd(i,j) == 0
%             continue
%         end
%         np=ceil(VolInd(i,j)./10);
%         pvol=VolInd(i,j)./np;
%         P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
%         P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
%         pidx=pidx+np;
%         P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
%         P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
%         P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
%         P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
%     end
% end

%% Inputs from DF, Tushar
% DF_BLink=[];
% DF_VolIn=[];
% 
% VolInDF(1:LinkNum,1)=0;
% 
% for i=1:length(DF_VolIn)
%    VolInDF(DF_BLink(i)==GridID_Brendan,1)=DF_VolIn(i); 
% end

%load Tushar\TusharDFvols.mat
load Tushar\TusharDFvols2.mat

%%
pidxdf=pidx;

%df
VolInd=repmat(VolInDF,1,gsclass).*repmat(Fdfpsd',LinkNum,1);
for i = 1:length(VolInd)
    for k = 1:7
        j=7-k+1;
        %j=k;
        if VolInd(i,j) == 0
            continue
        end
        np=ceil(VolInd(i,j)./10);
        pvol=VolInd(i,j)./np;
        P_loc{1,i}=cat(2,P_loc{1,i},repmat(0,1,np));
        P_idx{1,i}=cat(2,P_idx{1,i},pidx+(0:1:np-1));
        pidx=pidx+np;
        P_storage{1,i}=cat(2,P_storage{1,i},zeros(1,np));%activate
        P_tt{1,i}=cat(2,P_tt{1,i},zeros(1,np));%all 0 sec
        P_vol{1,i}=cat(2,P_vol{1,i},repmat(pvol,1,np));
        P_d{1,i}=cat(2,P_d{1,i},repmat(Dpsd(j),1,np));
    end
end

%% random sort
for i=1:LinkNum
% sub
subidx=find(P_idx{1,i}<pidxsf);
n = length(subidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(subidx)=P_idx{1,i}(subidx(ridx));
    P_vol{1,i}(subidx)=P_vol{1,i}(subidx(ridx));
    P_d{1,i}(subidx)=P_d{1,i}(subidx(ridx)); 
end

% sub2
subidx2=find(and(P_idx{1,i}<pidxsf,P_idx{1,i}>=pidxss2));
n = length(subidx2);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(subidx2)=P_idx{1,i}(subidx2(ridx));
    P_vol{1,i}(subidx2)=P_vol{1,i}(subidx2(ridx));
    P_d{1,i}(subidx2)=P_d{1,i}(subidx2(ridx)); 
end

% sf
% sfidx=find(and(P_idx{1,i}<pidxdf,P_idx{1,i}>=pidxsf));
% n = length(sfidx);
% if n ~= 0
%     ridx=randperm(n); 
%     P_idx{1,i}(sfidx)=P_idx{1,i}(sfidx(ridx));
%     P_vol{1,i}(sfidx)=P_vol{1,i}(sfidx(ridx));
%     P_d{1,i}(sfidx)=P_d{1,i}(sfidx(ridx)); 
% end

% df
dfidx=find(P_idx{1,i}>=pidxdf);
n = length(dfidx);
if n ~= 0
    ridx=randperm(n); 
    P_idx{1,i}(dfidx)=P_idx{1,i}(dfidx(ridx));
    P_vol{1,i}(dfidx)=P_vol{1,i}(dfidx(ridx));
    P_d{1,i}(dfidx)=P_d{1,i}(dfidx(ridx)); 
end
end

%% Input to all links at only t=0
% input(1:LinkNum,:)=1;
% %input(Lake==0,:)=1;
% 
% inp=find(input==1);
% pidx=0;
% for in=1:length(inp)
%     P_idx{1,inp(in)}=inp(in);
%     P_loc{1,inp(in)}=0;
%     P_storage{1,inp(in)}=0;
%     P_vol{1,inp(in)}=0.0001;
%     pidx=pidx+1;
%     %P_mass{1,inp(in)}=25*Length(inp(in)).*B(inp(in)).*H(inp(in))*1000;%10^5;%mg
% end
% 
% clear input inp

%% Determine where and how much sediment from uplands, ravines, and bluffs enters the network
% %% Uplands
% UpNum = sum(Lake==0);
% Up_Link(1:UpNum,1)=GridID(Lake==0);
% Up_Loc(1:UpNum,1)=0;
% Up_In_Vol=[];
% %Up_In_Vol(1:UpNum,1)=NaN;
% 
% 
% for i=1:LinkNum
% % old code prior to having incremental areas computed with loaded data
% %     %determine incremental area, for now load BE inc area file first
% %     %compute Adj and Area_km then load BE_NHD_MartinLakes
% %     incarea(i,1) = Area_km(find(network(i).GridID==GridID_mrb));
% %     %scale factor for sublink on area
% %     sf = Length(i) / network(i).length_m;
% %     if sf > 1
% %         sf = 1;
% %     end
% %     Up_In_Vol(i,1) = sf * incarea * 10.^6 * 0.00002 / 1.8 * 0.1;%m3/yr sand
%     
% %Up_In_Vol = Area(Lake==0) .* 0.00002 ./ 1.8 .* 0.1;%m3/yr sand
% if Lake(i)==0
%     Up_In_Vol = cat(1,Up_In_Vol,...
%         Area(i,1) .* 0.00002 ./ 2.65 .* Up_sand_frac(i,1));%m3/yr sand
%         %        m2  *  Mg/m2/yr / Mg/m3 * %sand
% end
% 
% %    clear sf
% end
% % 2.14e+4 m3/yr input to the basin from uplands
% 
% %% Ravines
% RavNum = length(ravines);
% Rav_Link(1:RavNum,1)=NaN;
% Rav_Loc(1:RavNum,1)=NaN;
% Rav_In_Vol(1:RavNum,1)=NaN;
% 
% for i=1:RavNum
%     
%     Rav_Link(i,1) = ravines(i).FeatID;%set link ID where ravine inputs
%     x=network(Rav_Link(i,1)).X;
%     y=network(Rav_Link(i,1)).Y;
%     %determine vertext of link where ravine inputs
%     DT = delaunayTriangulation(x(1:end-1)',y(1:end-1)');
%     fidx = nearestNeighbor(DT,ravines(i).X,ravines(i).Y);
%     %fidx = intersect(find(ravines(i).X==x),find(ravines(i).Y==y));
%     %determine location along link where ravine inputs
%     Rav_Loc(i,1) = network(Rav_Link(i,1)).L(fidx);
%     
%     Rav_In_Vol(i,1) = ravines(i).Area * 0.0034 / 2.65 * 0.35;%m3/yr sand
%     %                      m2        * Mg/m2/yr / Mg/m3 * %sand
%     
%     clear x y DT fidx
%       
% end
% %9.15e+3 m3/yr input to the basin from ravines
% 
% %% Bluffs
% BlfNum = length(bluffs);
% Blf_Link(1:BlfNum,1)=NaN;
% Blf_Loc(1:BlfNum,1)=NaN;
% Blf_In_Vol(1:BlfNum,1)=NaN;
% 
% for i=1:BlfNum
%     
%     Blf_Link(i,1) = bluffs(i).FeatID;%set link ID where ravine inputs
%     x=network(Blf_Link(i,1)).X;
%     y=network(Blf_Link(i,1)).Y;
%     %determine vertext of link where ravine inputs
%     DT = delaunayTriangulation(x(1:end-1)',y(1:end-1)');
%     fidx = nearestNeighbor(DT,bluffs(i).X,bluffs(i).Y);
%     %fidx = intersect(find(bluffs(i).X==x),find(bluffs(i).Y==y));
%     %determine location along link where ravine inputs
%     Blf_Loc(i,1) = network(Blf_Link(i,1)).L(fidx);
%     
%     Blf_In_Vol(i,1) = bluffs(i).Bluff_avgE * bluffs(i).SA_AtanTht * 0.35 * 1.8 / 2.65;%m3/yr sand
%     %                          m/yr        *           m2       * %sand  * Mg/m3 / Mg/m3
%     
%     clear x y DT fidx
% end
% %1.015e+5 m3/yr input to the basin from bluffs
% 
% %%
% sedbudvolin=sum(Blf_In_Vol)+sum(Up_In_Vol)+sum(Rav_In_Vol)

%% Poisson inputs in continuous time for bluffs, ravines, and uplands
% % Poisson inputs in continuous time at upstream end of all links
% % variable parcel volume
% 
% % concatenate all inputs into one vector
% clear In_Vol In_Link In_Loc In_Type InNum
% In_Vol=cat(1,Up_In_Vol,Rav_In_Vol,Blf_In_Vol);
% In_Link=cat(1,Up_Link,Rav_Link,Blf_Link);
% In_Loc=cat(1,Up_Loc,Rav_Loc,Blf_Loc);
% In_Type=cat(1,ones(UpNum,1),2.*ones(RavNum,1),3.*ones(BlfNum,1));
% InNum=UpNum+RavNum+BlfNum;
% 
% % change the seed to give truer random numbers
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% mu=1; %years, mean of interarrival time, exp. dist., =1/lambda
% pidx=1; %set unique parcel index to 1 for first parcel
% for i=1:InNum %loop through input locations   
%     R = exprnd(mu); %exponential random variables, interarrival time, years
%     cR=R; %arrival time, years
%     while cR<tmax, %loop through inputs
%         % concatenate arrival times of parcels to each link
%         L_arrival{In_Link(i,1),1}=cat(2,L_arrival{In_Link(i,1),1},cR);
%         
%         % locate continuous time input at next timestep
%         add = ceil(cR./(daystp./365))+1; %arrival at discrete timestep
%         pdt=(time(add,1)-cR)*365*24*60*60; %parcel dt, convert years to seconds
%         
%         % Move parcels downstream
%         sm=[];%time to move through current and downstream links
%         csm=[];%cumulative time to move through current and ds links
%         sm=STTime(In_Link(i,1),1)*(1-In_Loc(i,1));%seconds, time to move out of current link
%         csm=sm;%seconds
%         ii=In_Link(i,1);
%         while csm(end)<=pdt %loop until cumulative time toward outlet exceeds partial timestep
%             ii=Connect(ii,2); %set ds link
%             if isnan(ii) %if downstream is the outlet
%                 %leave outlet before dt ends
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             if Lake(ii) %if move into a lake
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             sm=cat(1,sm,STTime(ii,1));%add time through next link
%             csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
%         end
%         % update parcel location
%         if ~isnan(csm(end)) %check to make sure parcel is still in the system
%             pi=Connect(In_Link(i,1),length(csm)); %parcel link
%             if length(csm)==1 %still in same link
%                 %pl=pdt/STTime(i,1);%update location from us end
%                 pl=(STTime(In_Link(i,1),1)*In_Loc(i,1)+pdt)/STTime(In_Link(i,1),1);%update location from input loc
%             else %moved to a ds link
%                 pl=1-((csm(end)-pdt)/STTime(pi,1));%update location from beginning of link
%                 %note csm(end)-dt computes time remaining to move
%                 %through the rest of the link
%             end
%             %determine number of parcels to break input volume into
%             %min vol of parcel is 1/2 min capacity
% %            expx=exprnd(1);%exponential multiplier
% %            np=ceil(In_Vol(i,1)*expx./(caplim/2));
%             np=ceil(In_Vol(i,1)./(caplim/2));
%             %determine volume of each of those parcels
%             pvol=In_Vol(i,1)./np;
%             P_loc{add,pi}=cat(2,P_loc{add,pi},repmat(pl,1,np));
%             P_idx{add,pi}=cat(2,P_idx{add,pi},pidx+(0:1:np-1));
%             pidx=pidx+np-1;
%             P_storage{add,pi}=cat(2,P_storage{add,pi},zeros(1,np));%activate
%             P_vol{add,pi}=cat(2,P_vol{add,pi},repmat(pvol,1,np));
%             
% %             P_loc{add,pi}=cat(2,P_loc{add,pi},pl);
% %             P_idx{add,pi}=cat(2,P_idx{add,pi},pidx);
% %             P_storage{add,pi}=cat(2,P_storage{add,pi},0);%activate
% %             P_vol{add,pi}=cat(2,P_vol{add,pi},In_Vol(i,1));
%         end
%         %assign arrival times
%         for ai=1:length(csm)-1 %as long the parcel moved to a new link
%             if isnan(Connect(In_Link(i,1),ai+1)) %left the system
%                 OutArrival=cat(1,OutArrival,...
%                     cR+csm(ai,1)/60/60/24/365);
%             else %arrival times at each link
%                 L_arrival{Connect(In_Link(i,1),ai+1),1}=cat(2,L_arrival{Connect(In_Link(i,1),ai+1),1},...
%                     cR+csm(ai,1)/60/60/24/365);
%             end
%         end
%      
%         pidx=pidx+1; %increment unique parcel index
%         
%         % next input
%         R = exprnd(mu);
%         cR = cR + R;
%     end
% end
% pidx
% clear R cR add sm csm ii pi pl ai np pvol

%% Poisson inputs in continuous time for bluffs, ravines, and uplands
% % Poisson inputs in continuous time at upstream end of all links
% % constant parcel volume version 1
% 
% Vp=10;%m3, parcel volume
% 
% % concatenate all inputs into one vector
% clear In_Vol In_Link In_Loc In_Type InNum
% clear In_Vol_all In_Link_all In_Loc_all In_Type_all InNum_all
% In_Vol_all=cat(1,Up_In_Vol,Rav_In_Vol,Blf_In_Vol);
% In_Link_all=cat(1,Up_Link,Rav_Link,Blf_Link);
% In_Loc_all=cat(1,Up_Loc,Rav_Loc,Blf_Loc);
% In_Type_all=cat(1,ones(UpNum,1),2.*ones(RavNum,1),3.*ones(BlfNum,1));
% InNum_all=UpNum+RavNum+BlfNum;
% 
% In_Vol=[];
% In_Link=[];
% In_Loc=[];
% In_Type=[];
% 
% In_Vol_all=round(In_Vol_all./10).*10;
% for i = 1:length(In_Vol_all)
%     np=In_Vol_all(i,1)./10;
%     In_Vol=cat(1,In_Vol,repmat(Vp,np,1));
%     In_Link=cat(1,In_Link,repmat(In_Link_all(i,1),np,1));
%     %In_Loc=cat(1,In_Loc,repmat(In_Loc_all(i,1),np,1));%all inputs at real location
%     In_Loc=cat(1,In_Loc,zeros(np,1));%all inputs at US end of link
%     In_Type=cat(1,In_Type,repmat(In_Type_all(i,1),np,1));
% end
% InNum=numel(In_Type);
% 
% clear In_Vol_all In_Link_all In_Loc_all In_Type_all InNum_all np
% %%
% % change the seed to give truer random numbers
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% mu=1; %years, mean of interarrival time, exp. dist., =1/lambda
% pidx=1; %set unique parcel index to 1 for first parcel
% for i=1:InNum %loop through input locations  
%     i
%     InNum
%     R = exprnd(mu); %exponential random variables, interarrival time, years
%     cR=R; %arrival time, years
%     while cR<tmax, %loop through inputs
%         % concatenate arrival times of parcels to each link
%         L_arrival{In_Link(i,1),1}=cat(2,L_arrival{In_Link(i,1),1},cR);
%         
%         % locate continuous time input at next timestep
%         add = ceil(cR./(daystp./365))+1; %arrival at discrete timestep
%         pdt=(time(add,1)-cR)*365*24*60*60; %parcel dt, convert years to seconds
%         
%         % Move parcels downstream
%         sm=[];%time to move through current and downstream links
%         csm=[];%cumulative time to move through current and ds links
%         sm=STTime(In_Link(i,1),1)*(1-In_Loc(i,1));%seconds, time to move out of current link
%         csm=sm;%seconds
%         ii=In_Link(i,1);
%         while csm(end)<=pdt %loop until cumulative time toward outlet exceeds partial timestep
%             ii=Connect(ii,2); %set ds link
%             if isnan(ii) %if downstream is the outlet
%                 %leave outlet before dt ends
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             if Lake(ii) %if move into a lake
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             sm=cat(1,sm,STTime(ii,1));%add time through next link
%             csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
%         end
%         % update parcel location
%         if ~isnan(csm(end)) %check to make sure parcel is still in the system
%             pi=Connect(In_Link(i,1),length(csm)); %parcel link
%             if length(csm)==1 %still in same link
%                 %pl=pdt/STTime(i,1);%update location from us end
%                 pl=(STTime(In_Link(i,1),1)*In_Loc(i,1)+pdt)/STTime(In_Link(i,1),1);%update location from input loc
%             else %moved to a ds link
%                 pl=1-((csm(end)-pdt)/STTime(pi,1));%update location from beginning of link
%                 %note csm(end)-dt computes time remaining to move
%                 %through the rest of the link
%             end
%             
% %             %determine number of parcels to break input volume into
% %             %min vol of parcel is 1/2 min capacity
% % %            expx=exprnd(1);%exponential multiplier
% % %            np=ceil(In_Vol(i,1)*expx./(caplim/2));
% %             np=ceil(In_Vol(i,1)./(caplim/2));
% %             %determine volume of each of those parcels
% %             pvol=In_Vol(i,1)./np;
% %             P_loc{add,pi}=cat(2,P_loc{add,pi},repmat(pl,1,np));
% %             P_idx{add,pi}=cat(2,P_idx{add,pi},pidx+(0:1:np-1));
% %             pidx=pidx+np-1;
% %             P_storage{add,pi}=cat(2,P_storage{add,pi},zeros(1,np));%activate
% %             P_vol{add,pi}=cat(2,P_vol{add,pi},repmat(pvol,1,np));
%             
%             P_loc{add,pi}=cat(2,P_loc{add,pi},pl);
%             P_idx{add,pi}=cat(2,P_idx{add,pi},pidx);
%             P_storage{add,pi}=cat(2,P_storage{add,pi},0);%activate
%             P_vol{add,pi}=cat(2,P_vol{add,pi},In_Vol(i,1));
%         end
%         %assign arrival times
%         for ai=1:length(csm)-1 %as long the parcel moved to a new link
%             if isnan(Connect(In_Link(i,1),ai+1)) %left the system
%                 OutArrival=cat(1,OutArrival,...
%                     cR+csm(ai,1)/60/60/24/365);
%             else %arrival times at each link
%                 L_arrival{Connect(In_Link(i,1),ai+1),1}=cat(2,L_arrival{Connect(In_Link(i,1),ai+1),1},...
%                     cR+csm(ai,1)/60/60/24/365);
%             end
%         end
%      
%         pidx=pidx+1; %increment unique parcel index
%         
%         % next input
%         R = exprnd(mu);
%         cR = cR + R;
%     end
% end
% pidx
% clear R cR add sm csm ii pi pl ai np pvol

%% add single pulse input
% tadd=2001;
% ladd=446;
% vadd=capacity(ladd).*4;
% 
% np=ceil(vadd./(caplim/2));
% %determine volume of each of those parcels
% pvol=vadd./np;
% 
% ppidxstart=pidx
% P_loc{tadd,ladd}=cat(2,P_loc{tadd,ladd},repmat(0,1,np));
% P_idx{tadd,ladd}=cat(2,P_idx{tadd,ladd},pidx+(0:1:np-1));
% pidx=pidx+np-1;
% P_storage{tadd,ladd}=cat(2,P_storage{tadd,ladd},zeros(1,np));%activate
% P_vol{tadd,ladd}=cat(2,P_vol{tadd,ladd},repmat(pvol,1,np));
% pidx
%% add pulse input everywhere
% tadd=2001;
% ppidxstart=pidx
% 
% for ladd=1:LinkNum;
%     
%     if Lake(ladd)>0
%         continue
%     end
%     
%     vadd=capacity(ladd).*2;
%     
%     np=ceil(vadd./(caplim/2));
%     %determine volume of each of those parcels
%     pvol=vadd./np;
%     
%     P_loc{tadd,ladd}=cat(2,P_loc{tadd,ladd},repmat(0,1,np));
%     P_idx{tadd,ladd}=cat(2,P_idx{tadd,ladd},pidx+(0:1:np-1));
%     pidx=pidx+np-1;
%     P_storage{tadd,ladd}=cat(2,P_storage{tadd,ladd},zeros(1,np));%activate
%     P_vol{tadd,ladd}=cat(2,P_vol{tadd,ladd},repmat(pvol,1,np));
% end
% pidx
%%
%% Input at bluffs and ravines
% % recurrent every year
% pidx=1; %set unique parcel index to 1 for first parcel
% %uplands
% for i=1:LinkNum
%     if ~Lake(i)
%         P_idx{1,i}=cat(2,P_idx{1,i},pidx); %each parcel a unique id
%         P_loc{1,i}=cat(2,P_loc{1,i},0); %each parcel input at US end
%         P_active{1,i}=cat(2,P_active{1,i},1); %each parcel active  
%         P_vol{1,i}=cat(2,P_vol{1,i},Up_In_Vol(i)); %input volume 
%     end
%     pidx=pidx+1; %increment unique parcel index
% end
% %ravines
% for i=1:Ravnum %loop through inputs
%     P_idx{1,Rav_Link(i)}=cat(2,P_idx{1,Rav_Link(i)},pidx); %each parcel a unique id
%     P_loc{1,Rav_Link(i)}=cat(2,P_loc{1,Rav_Link(i)},Rav_Loc(i)); %each parcel input at Loc
%     P_active{1,Rav_Link(i)}=cat(2,P_active{1,Rav_Link(i)},1); %each parcel active 
%     P_vol{1,Rav_Link(i)}=cat(2,P_vol{1,Rav_Link(i)},Rav_In_Vol(i)); %input volume 
%     pidx=pidx+1; %increment unique parcel index
% end
% %bluffs
% for i=1:Blfnum %loop through inputs
%     P_idx{1,Blf_Link(i)}=cat(2,P_idx{1,Blf_Link(i)},pidx); %each parcel a unique id
%     P_loc{1,Blf_Link(i)}=cat(2,P_loc{1,Blf_Link(i)},Blf_Loc(i)); %each parcel input at Loc
%     P_active{1,Blf_Link(i)}=cat(2,P_active{1,Blf_Link(i)},1); %each parcel active  
%     P_vol{1,Blf_Link(i)}=cat(2,P_vol{1,Blf_Link(i)},Blf_In_Vol(i)); %input volume 
%     pidx=pidx+1; %increment unique parcel index
% end
% pidx
% %% Repeat this set of inputs every year
% tinput = (1:1:floor(tmax))';
% tinidx(1:length(tinput),1)=NaN;
% for i=1:length(tinput)
%     tinidx(i)=find(time>tinput(i),1);
% end
% 
% for t=1:length(tinidx)
%     for i=1:LinkNum
%         P_idx{tinidx(t),i}=P_idx{1,i};%REPEATS SAME INDEX!!!
%         P_loc{tinidx(t),i}=P_loc{1,i};
%         P_active{tinidx(t),i}=P_active{1,i};
%         P_vol{tinidx(t),i}=P_vol{1,i};
%     end
% end
% clear tinput tinidx
%%
%% Poisson inputs in continuous time
% % Poisson inputs in continuous time at upstream end of all links
% % change the seed to give truer random numbers
% RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% 
% mu=1; %years, mean of interarrival time, exp. dist., =1/lambda
% pidx=1; %set unique parcel index to 1 for first parcel
% for i=1:LinkNum %loop through links
%     
%     if Lake(i)
%         continue
%     end
%     
%     R = exprnd(mu); %exponential random variables, interarrival time, years
%     cR=R; %arrival time, years
%     while cR<tmax, %loop through inputs 
%         % concatenate arrival times of parcels to each link
%         L_arrival{i,1}=cat(2,L_arrival{i,1},cR);
%         
%         % locate continuous time input at next timestep
%         add = ceil(cR./(daystp./365))+1; %arrival at discrete timestep
%         pdt=(time(add,1)-cR)*365*24*60*60; %parcel dt, convert years to seconds
%         
%         
%         % Move parcels downstream
%         sm=[];%time to move through current and downstream links
%         csm=[];%cumulative time to move through current and ds links
%         sm=STTime(i,1);%seconds, time to move out of current link
%         csm=sm;%seconds
%         ii=i;
%         while csm(end)<=pdt %loop until cumulative time toward outlet exceeds partial timestep
%             ii=Connect(ii,2); %set ds link
%             if isnan(ii) %if downstream is the outlet
%                 %leave outlet before dt ends
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             if Lake(ii) %if move into a lake
%                 sm=cat(1,sm,NaN);
%                 csm=cat(1,csm,sum(sm));
%                 break
%             end
%             sm=cat(1,sm,STTime(ii,1));%add time through next link
%             csm=cat(1,csm,sum(sm));%cumulative time to move through all subsequent links
%         end
%         % update parcel location
%         if ~isnan(csm(end)) %check to make sure parcel is still in the system
%             pi=Connect(i,length(csm)); %parcel link
%             if length(csm)==1 %still in same link
%                 pl=pdt/STTime(i,1);%update location from us end
%             else %moved to a ds link
%                 pl=1-((csm(end)-pdt)/STTime(pi,1));%update location from beginning of link
%                 %note csm(end)-dt computes time remaining to move
%                 %through the rest of the link
%             end
%             P_loc{add,pi}=cat(2,P_loc{add,pi},pl);
%             P_idx{add,pi}=cat(2,P_idx{add,pi},pidx);
%             P_storage{add,pi}=cat(2,P_storage{add,pi},0);%activate
% %            P_vol{add,pi}=cat(2,P_vol{add,pi},1);%parcel volume is 1 for all
%             P_vol{add,pi}=cat(2,P_vol{add,pi},Area_km(i,1).*sedbudvolin./max(usarea_km));
%         end
%         %assign arrival times
%         for ai=1:length(csm)-1 %as long the parcel moved to a new link
%             if isnan(Connect(i,ai+1)) %left the system
%                 OutArrival=cat(1,OutArrival,...
%                     cR+csm(ai,1)/60/60/24/365);
%             else %arrival times at each link
%                 L_arrival{Connect(i,ai+1),1}=cat(2,L_arrival{Connect(i,ai+1),1},...
%                     cR+csm(ai,1)/60/60/24/365);
%             end
%         end
%      
%         pidx=pidx+1; %increment unique parcel index
%         
%         % next input
%         R = exprnd(mu);
%         cR = cR + R;
%     end
% end
% pidx
% clear R cR add sm csm ii pi pl ai
