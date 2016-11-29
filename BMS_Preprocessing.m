%% BMS_Preprocessing
% This file contains useful code for preprocessing the network/lake
% shapefiles for use in the network-based, Bed-Material Sediment model code.

% Jon Czuba
% July 21, 2015

%% load shapefiles
%NHD BE
% network = shaperead('Map\MRB_NHD\BE_NHD_network_prj.shp');
% boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_diss_prj.shp');

%node = shaperead('Map\NodeXYElev.shp');
% Martin Lakes
network = shaperead('Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att.shp');
boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_diss_prj.shp');
%ravines = shaperead('Map\MRB_NHD\MartinLakes\GBER_ravines_merge_points_snapedit_BEnetMartin_noli.shp');
%bluffs = shaperead('Map\MRB_NHD\MartinLakes\GBER_bluffs_906_trim_points2_editsnap2_basinE_BEnetMartin_noli.shp');
lakes = shaperead('Map\MRB_NHD\MartinLakes\nhd_waterbody_netintsct_gt04km_poly2_prj.shp');
%% Martin Lakes LS
network = shaperead('Map\MRB_NHD\MartinLakes\LS_NHD_network_MartinLake_ds_prj_att.shp');
boundary = shaperead('Map\MRB_NHD\Catchment_NHD_MRB5_BE2_LS3_diss_prj.shp');
lakes = shaperead('Map\MRB_NHD\MartinLakes\nhd_waterbody_netintsct_gt04km_poly2_LS.shp');

%% Write Shapefiles
% for i=1:LinkNum
%    [network(i).seddepth] = seddepth(6001,i);
%    %[network(i).numfracocap] = numfracocap(i); 
%    %[network(i).exceedHs] = exceedHs(i);
%    %[network(i).wfpeak1] = lnkpl(i);
% end
% % Write Shapefile
% %shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_MartinLake_ds_prj_att_overcap.shp');
% shapewrite(network,'Map\MRB_NHD\MartinLakes\BE_NHD_network_stor6001seddepth.shp');

%%
% filename='Map\MRB_NHD\mnr_30m_w_prj_r2a.txt';
% %read grid
% [Z,R] = arcgridread(filename);

%% plot network
figure
box on
mapshow(network);
mapshow(boundary, 'FaceColor', 'none', 'EdgeColor', 'black');
mapshow(lakes, 'FaceColor', 'blue', 'EdgeColor', 'none');
%mapshow(bluffs);
%mapshow(ravines);
xlabel('Easting in meters')
ylabel('Northing in meters')
%%

DSlink(1:1360,1)=0;
for i=1339:1360
DSlink(i)=FeatID(find(dslink(i)==GridID,1),1);
end

%%
LinkNum=length(network);
GridID(1:LinkNum,1)=NaN;
%FromNode(1:LinkNum,1)=NaN;
ToNode(1:LinkNum,1)=NaN;
Length(1:LinkNum,1)=NaN;
Area_km(1:LinkNum,1)=NaN;
Slope(1:LinkNum,1)=NaN;
%Node(1:LinkNum+1,1)=NaN;
%NodeElev(1:LinkNum+1,1)=NaN;
usarea_km(1:LinkNum,1)=NaN;
%maxelev(1:LinkNum,1)=NaN;
%minelev(1:LinkNum,1)=NaN;
Lake(1:LinkNum,1)=NaN;

% GridID = [network.GridID]';
% FromNode = [network.from_node]';
% ToNode = [network.to_node]';
% Length = [network.Shape_Leng]';
% Area = [network.Shape_Area]';
% %Slope = [network.Slope]';
% Node = [node.Node]';
% NodeElev = [node.Elev]';
% 
% for i=1:LinkNum
%    Slope(i,1)= (NodeElev(find(FromNode(i,1)==Node))-...
%        NodeElev(find(ToNode(i,1)==Node)))./Length(i,1);
% end

% GridID = [network.Hydroseq]';
% FromNode = [network.UpHydroseq]';
% ToNode = [network.DnHydroseq]';
% Length = [network.LENGTHKM]';
% usarea_km = [network.TotDASqKM]';
% Slope = [network.SLOPE]';
% %Node = [node.Node]';
% %NodeElev = [node.Elev]';

% GridID = [network.GridID]';
% %FromNode = [network.UpHydroseq]';
% ToNode = [network.dslink]';
% Length = [network.length_m]';
% usarea_km = [network.usarea_km2]';
% Slope = [network.slope]';
% maxelev = [network.maxelevsm_]';%m
% minelev = [network.minelevsm_]';%m

GridID = [network.FeatID]';
ToNode = [network.FeatDSlnk]';
Length = [network.l_m]';
%usarea_km = [network.A_km2]';
Area_km = [network.a_km2]';
%Slope = [network.slope]';
mxelev = [network.mxelev_m]';%m
mnelev = [network.mnelev_m]';%m
Lake = [network.LakeID]';
%usarea=usarea_km.*10.^6;
Area=Area_km.*10.^6;

Slope=(mxelev-mnelev)./Length;
Slope(Slope<1e-4)=1e-4;

%% Determine outlet link ID
% for i=1:LinkNum
%     %find element of ToNode with one unique number
%     if numel(find(ToNode(i)==ToNode))==1
%         OutletLinkID=find(ToNode(i)==ToNode);
%         break
%     end
% end
%%
%OutletLinkID=find(ToNode==0);
%OutletLinkID_mrb=1802;%Cottonwood%1065;%5745;
OutletLinkID_mrbbe=1296;%MartinLakeLS
%OutletLinkID=1338;%MartinLakeBE
%% Reassign GridID, ToNode, OutletLinkID for subbasins
GridID_mrbbe=GridID;
ToNode_mrbbe=ToNode;
GridID=(1:LinkNum)';
OutletLinkID=find(GridID_mrbbe==OutletLinkID_mrbbe);
for i=1:LinkNum
    if i==OutletLinkID
        ToNode(i)=0;
    else
        ToNode(i)=find(ToNode_mrbbe(i)==GridID_mrbbe);
    end
end

%% Establish connectivity structure of each link to the outlet
% each row in Connect is the connectivity structure from that
% row to the outlet. the columns correspond to the subsequent
% links to connect to the outlet.
% Connect(1:LinkNum,1:LinkNum)=NaN;
% Connect(1:LinkNum,1)=1:LinkNum;
% for i=1:LinkNum
%     j=1;
%     while ToNode(Connect(i,j),1)~=ToNode(OutletLinkID)
%         Connect(i,j+1)=find(ToNode(Connect(i,j),1)==FromNode);
%         j=j+1;
%     end
% end
%
Connect(1:LinkNum,1:LinkNum)=NaN;
Connect(1:LinkNum,1)=1:LinkNum;
for i=1:LinkNum
    j=1;
    while ToNode(Connect(i,j),1)~=ToNode(OutletLinkID)
        Connect(i,j+1)=find(ToNode(Connect(i,j),1)==GridID);
        j=j+1;
    end
end
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

%%
% %% Remove Lakes from Sub
% % remove from Sub the contributions upstream of lakes for locations
% % downstream of lakes
% for i=1:LinkNum
%     if Lake(i)~=0 %for lake links
%         for j=1:LinkNum
%             if Sub(i,j) 
%                 if j~=i
%                 Sub(logical(Sub(:,i)),j)=0;
%                 end
%             end
%         end
%     end
% end
% 
% %% compute PSWF
% %same as PSWFi{OutletLinkID,1}
% PSWFout(1:LinkNum,1)=NaN;
% WFout(1:LinkNum,1)=NaN;
% %Remout(1:LinkNum,1)=NaN;
% for i=1:LinkNum
%     lastConn=find(Connect(i,:)==OutletLinkID);
%     PSWFout(i,1)=sum(STTime(Connect(i,1:lastConn),1))./60./60./24./365;%years to outlet 
%     WFout(i,1)=sum(Length(Connect(i,1:lastConn),1))./1000;%km to outlet
%     %Remout(i,1)=sum(Rem(Connect(i,1:lastConn),1));%kg removed from input to outlet
% end
% clear lastConn
% 
% %% compute max PSWF for each link
% PSWFi=cell(LinkNum,1);%years, PSWF for each link
% %WFi=cell(LinkNum,1);%km, WF for each link
% %WFidx=cell(LinkNum,1);%km, WF for each link
% for i=1:LinkNum
%     uslnks=find(Sub(:,i)==1);
%     for ii=1:length(uslnks)    
%         lastConn=find(Connect(uslnks(ii),:)==i);
%         trt=sum(STTime(Connect(uslnks(ii),1:lastConn),1))./60./60./24./365;%years to outlet
%         PSWFi{i,1}=cat(2,PSWFi{i,1},trt);%years
%         %dst=sum(Length(Connect(uslnks(ii),1:lastConn),1))./1000;%km to outlet
%         %WFi{i,1}=cat(2,WFi{i,1},dst);%km
%         %WFidx{i,1}=cat(2,WFidx{i,1},uslnks(ii));%index corresponding to distance
%     end 
% end
% 
% PSWFmx=cellfun(@max,PSWFi);
% clear uslnks ii lastConn trt dst
% 
% %hist(PSWFi{OutletLinkID,1},50)

%% Sum upstream area
usarea(1:LinkNum,1)=NaN;
for i=1:LinkNum
    usarea(i,1)=sum(Sub(:,i).*Area);
end
usarea_km=usarea./10.^6;
%give Sub NaN for now for plotting below
%Sub(Sub==0)=NaN;
%Sub(isnan(Sub))=0;

%% Construct Adjaceny matrix
Adj(1:LinkNum,1:LinkNum)=0;
%establish DS connection
for i=1:LinkNum
    if isnan(Connect(i,2))
        continue
    end
    Adj(Connect(i,1),Connect(i,2))=1;
end
%establish US connection
%Adj=Adj+Adj';
%Adj(Adj==0)=NaN;
%Adj(isnan(Adj))=0;

%% Determine incremental area
Area_km(1:LinkNum,1)=NaN;
for i=1:LinkNum
Area_km(i,1)=usarea_km(i)-sum(usarea_km(Adj(:,i)>0));
end

%% Stream Order
% %first order only
% %Ord1=find(nansum(Sub,1)==1)';
% Ord1(1:LinkNum,:)=0;
% for i=1:LinkNum
%    if isempty(find(i==ToNode,1)) 
%        Ord1(i)=1;
%    end
% end

%% Lake Properties
%%
% Lke_ID=[lakes.LakeID]';
% Lke_Larea_km=[lakes.Shape_Area]'./10^6;
% %%
% Lke_ToNode(1:length(Lke_ID),1)=NaN;
% Lke_Warea_km(1:length(Lke_ID),1)=NaN;
% for i=1:length(Lke_ID)
%     %find GridID of links with same LakeID
%     ntlk=find(Lake==Lke_ID(i));
%     %find links those links connect to
%     nlcon=ToNode(ntlk);
%     for j=1:length(ntlk)
%         %find link which connects to only 1 ds link
%         if numel(find(nlcon(j)==nlcon))==1
%             Lke_ToNode(i,1)=nlcon(j);
%             Lke_Warea_km(i,1)=usarea_km(ntlk(j));
%             break
%         end
%     end
% end

%%
%% Determine incremental length spatially though links 2-D
%length through link
for i=1:LinkNum
    network(i).L=network(i).X;
end
for i=1:LinkNum
    network(i).L(1)=0;
    for j=2:length(network(i).X)-1
        network(i).L(j)=sqrt((network(i).X(j)-network(i).X(j-1)).^2+...
            (network(i).Y(j)-network(i).Y(j-1)).^2);
    end
end
for i=1:LinkNum
    for j=3:length(network(i).X)-1
        network(i).L(j)=network(i).L(j)+network(i).L(j-1);
    end
    network(i).L=network(i).L./nanmax(network(i).L);
end

%% Determine distance between 2 points
% su=50;%42  %54;%46%site id of us point
% lu=1119;%686;%link of us point
% sd=51;%43  %53;%45%site id of ds point
% ld=1113;%1054;%link of ds point
% 
% iu=find(lu==GridID_mrb);
% id=find(ld==GridID_mrb);
% 
% d=sum(Length(Connect(iu,2:find(Connect(iu,:)==id)-1)))+...
% (1-network(iu).L(find(points(42).X==network(iu).X))).*Length(iu)+...%us point
% network(id).L(93).*Length(id);%ds point
% %network(id).L(find(points(43).X==network(id).X)).*Length(id);%ds point


%%
%% NHD_GIS_preprocessing
% In preparing the NHD network for use with other Matlab codes, this code
% helps preprocess the NHD network in conjunction with GIS.
% 
% %% Multipoint nodes to single point nodes
% x=[];
% y=[];
% %paste x, y values
% %%
% xe=x;
% ye=y;
% X=[];
% Y=[];
% Xg=[];
% Yg=[];
% %keep only points that occur 1 (ends) or 3 (junctions) times, points that
% %occur 2 (along flowpath) times are to be removed.
% %%
% for i=1:length(xe)
%     
%     if isnan(xe(i))
%         continue
%     end
%     clear ind1 ind2
%     ind1=find(xe(i)==xe);
%     ind2=find(ye(i)==ye);
%     ind=[];
%     for j=1:length(ind1)
%         if ~isempty(find(ind1(j)==ind2))
%             ind=cat(1,ind,ind1(j));
%         end
%     end
%     
%     if length(ind)==2
%         xe(ind)=NaN;
%         ye(ind)=NaN;
%     elseif length(ind)==1
%         X=cat(1,X,xe(ind(1)));
%         Y=cat(1,Y,ye(ind(1)));
%         xe(ind)=NaN;
%         ye(ind)=NaN;
%     elseif length(ind)==3
%         X=cat(1,X,xe(ind(1)));
%         Y=cat(1,Y,ye(ind(1)));
%         xe(ind)=NaN;
%         ye(ind)=NaN;
%     else
%         Xg=cat(1,Xg,xe(ind(1)));
%         Yg=cat(1,Yg,ye(ind(1)));
%         xe(ind)=NaN;
%         ye(ind)=NaN;
%     end
% end
% 
% %% Flowline info to condense
% % paste into here
% %LENGTHKM=[];%line length in kilometers
% DivDASqKM=[];%upstream drainage area in kilometers
% %MAXELEVRAW=[];%cm
% %MINELEVRAW=[];%cm
% MAXELEVSMO=[];%cm
% MINELEVSMO=[];%cm
% %SLOPE=[];
% SLOPELENKM=[];%length in kilometers used to compute slope
% Hydroseq=[];%id of link for ds connections
% DnHydroseq=[];%downstream Hydroseq
% FID=[];%id of target network
% %%
% FID=FID+1;
% links=max(FID);
% %assigned to here
% lengthkm(1:links,1)=NaN;%line length in kilometers
% divdasqkm(1:links,1)=NaN;%upstream drainage area in kilometers
% maxelevraw(1:links,1)=NaN;%cm
% minelevraw(1:links,1)=NaN;%cm
% maxelevsmo(1:links,1)=NaN;%cm
% minelevsmo(1:links,1)=NaN;%cm
% slope(1:links,1)=NaN;
% slopelenkm(1:links,1)=NaN;%length in kilometers used to compute slope
% hydroseq(1:links,1)=NaN;%id of link for ds connections
% dnhydroseq(1:links,1)=NaN;%downstream Hydroseq
% fid(1:links,1)=NaN;%id of target network
% 
% GridID(1:links,1)=NaN;
% dslink(1:links,1)=NaN;
% 
% %%
% for i=1:max(FID)
%    
%     clear ind
%     ind=find(FID==i);
%     
%     fid(i)=i;
%     
%     %max
%     divdasqkm(i)=max(DivDASqKM(ind));
%     %maxelevraw(i)=max(MAXELEVRAW(ind));
%     maxelevsmo(i)=max(MAXELEVSMO(ind));
%     
%     %min
% %     if numel(ind)==1 && MINELEVRAW(ind)==-9998
% %         minelevraw(i)=MINELEVRAW(ind);
% %     else
% %         minelevraw(i)=min(MINELEVRAW(ind(MINELEVRAW(ind)>0)));
% %     end
%     if numel(ind)==1 && MINELEVSMO(ind)==-9998
%         minelevsmo(i)=MINELEVSMO(ind);
%     else
%         minelevsmo(i)=min(MINELEVSMO(ind(MINELEVSMO(ind)>0)));
%     end
%     
%     %sum
%     %lengthkm(i)=sum(LENGTHKM(ind));
%     slopelenkm(i)=sum(SLOPELENKM(ind));
%     
%     for j=1:length(ind)
%         if sum(find(Hydroseq(ind(j))==DnHydroseq(ind)))==0
%             hydroseq(i)=Hydroseq(ind(j));
%         end
%         if sum(find(DnHydroseq(ind(j))==Hydroseq(ind)))==0
%             dnhydroseq(i)=DnHydroseq(ind(j));
%         end
%     end  
%     
% end
% %%
% slope=((maxelevsmo-minelevsmo)./100)./(slopelenkm.*1000);
% GridID=fid;
% fid=fid-1;
% %%
% for i=1:links
%     if dnhydroseq(i)==510003427 %outlet
%         dslink(i)=0;
%     else
%         dslink(i)=find(hydroseq==dnhydroseq(i));
%     end
% end
% %%
% slopelenkm.*1000;%m
% slope(slope<1e-5)=1e-5;%consistent with NHD slopes
% maxelevsmo./100;%m
% minelevsmo./100;%m

%% Extra processing of network
% for i=1:1360
%     %mlGridID(i,1)=network(i).GridID;
%     %mlFeatLeng(i,1)=network(i).FeatLeng;
%     %leng=[network.l_m]';
%     %slope=[network.slope]';
%     lkfeatid=[network2.FeatID]';
%     lklake=[network2.Lake]';
%     lklakeid=[network2.LakeID]';
%     %sf(i,1)=mlFeatLeng(i,1)./sum(mlFeatLeng(find(mlGridID(i)==mlGridID),1));
%     [val idx]=sort(usarea(find(mlGridID(i)==mlGridID)));
%     sz=length(idx);
%     gd=find(mlGridID(i)==mlGridID);
%     idx=gd(idx);
%     if sz==1
%         mxelev(i,1)=network(i).maxelevsm_;
%         mnelev(i,1)=network(i).minelevsm_;
%     else
%         mxelev(idx(1),1)=network(i).maxelevsm_;
%         mnelev(idx(end),1)=network(i).minelevsm_;
%         for j=1:sz-1
%             mnelev(idx(j),1)=mxelev(idx(j))-slope(idx(j)).*leng(idx(j));
%             mxelev(idx(j+1),1)=mnelev(idx(j),1);   
%         end
%         %i=max(idx);
%     end
% end
% %
% %a=incarea.*sf;
% for i=1:1360
%    if lklake(i)==1
%        network(lkfeatid(i)).LakeID=lklakeid(i);
%    end
% end
% %
% for i = 1:1360
%     %[network(i).incarea] = incarea(i);
%     %[network(i).sf] = sf(i,1);
%     %[network(i).a_km2] = incarea(i).*sf(i,1);
%     %[network(i).l_m] = network(i).length_m.*sf(i,1);
%     %[network(i).A_km2] = usarea(i);
%     %[network(i).mxelev_m] = mxelev(i);
%     %[network(i).mnelev_m] = mnelev(i);
%     [network(i).LakeID] = 0;
% end
