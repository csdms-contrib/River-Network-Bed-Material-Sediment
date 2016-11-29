%% Width Function Analysis
% This program scales the width function based on the process operating on 
% the network. For instance sediment transport. This version has link 
% contributions at the US end of links

% Jon Czuba
% Univ. of Minnesota
% January 2, 2013

%%
load BE_NHD.mat
%%
% %% Establish connectivity structure of each link to the outlet
% % each row in Connect is the connectivity structure from that
% % row to the outlet. the columns correspond to the subsequent
% % links to connect to the outlet.
% Connect(1:LinkNum,1:LinkNum)=NaN;
% Connect(1:LinkNum,1)=1:LinkNum;
% for i=1:LinkNum
%     j=1;
%     while ToNode(1,Connect(i,j))~=ToNode(OutletLinkID)
%         Connect(i,j+1)=find(ToNode(1,Connect(i,j))==FromNode);
%         j=j+1;
%     end
% end

%% Compute width function for every link
% for every row of Connect i (every link in the basin),
% find which links connect through that point and mark these
% links j in the column of Sub
% diagonal = 1

% Sub(i,j)
% *each row i indicates which indicies that nodes connects
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

%% Sum upstream area
% usarea(1:LinkNum,1)=NaN;
% for i=1:LinkNum
%     %[r c]=find(Connect(i,1)==Connect);
%     %usarea(i,1)=sum(Area(1,r));
%     %give Sub zeros
%     usarea(i,1)=sum(Sub(:,i).*Area');
% end
% usarea_km=usarea./10.^6;
% %give Sub NaN for now for plotting below
% Sub(Sub==0)=NaN;
% %Sub(isnan(Sub))=0;

%% Slope
% %replace links with zero slope with the S~A relation for the basin
% Slope(Slope==0)=NaN;
% A=usarea_km(~isnan(Slope)')';
% Y=Slope(~isnan(Slope));
% coeff = polyfit(log10(A),log10(Y),1);
% d=10.^coeff(2);
% h=coeff(1);
% Slope(isnan(Slope))=d.*usarea_km(isnan(Slope)).^h;
Slope(Slope<1e-4)=1e-4;
%%
%using NHD
usarea=usarea';

Q=(1.3956e-5).*(usarea').^(0.7015);
B=(0.0238).*(usarea').^(0.3397);
H=(0.0029).*(usarea').^(0.2943);
U=(0.1976).*(usarea').^(0.0679);
%Cf=(0.6995).*(usarea').^(-0.1753);
%Slope=(0.2953).*(usarea').^(-0.3078);
Cf=9.81.*H.*Slope./(U.^2);
tausS=H.*Slope./1.65./0.0001;
%tausSS=Cf.*U.^2./9.81./1.65./0.0001;
D=0.4;%mm - sediment grain size
D=D./1000;%m
R=1.65;%submerged specific gravity of sediment
g=9.81;%m/s2 - acceleration due to gravity
theta=0.1;%intermittency factor

VelocityS=0.05./sqrt(g)./R./R./D.*U.^2.*H.^(1/2).*Slope.^(3/2)./theta;

%some parameters are hardcoded here
tausG=H.*Slope./1.65./0.01;
VelocityG=sqrt(1.65.*9.81).*0.01.^(3/2)./0.5.*4.*(tausG-0.0495).^(3/2);
VelocityG(tausG<0.0495)=0;
VelocityM=U;

%% Grainsize for full transport at Q2
% H=(0.0029).*(usarea).^(0.2943);
% Dgm=H.*Slope./1.65./0.0495;
% 
% tausG=H.*Slope./1.65./0.01;
% VelocityG=sqrt(1.65.*9.81).*0.01.^(3/2)./0.5.*4.*(tausG-0.0495).^(3/2);
% VelocityG(tausG<0.0495)=0;

%% Sum distances from each link to the outlet
%using NHD
Length=Length';
VelocityS=VelocityS';
VelocityM=VelocityM';
VelocityG=VelocityG';

lastConn(1:LinkNum,1)=NaN;
Dist(1:LinkNum,1)=NaN;
%TimeU(1:LinkNum,1)=NaN;
TimeM(1:LinkNum,1)=NaN;
TimeS(1:LinkNum,1)=NaN;
TimeG(1:LinkNum,1)=NaN;

for i=1:LinkNum
    
    lastConn(i,1)=find(Connect(i,:)==OutletLinkID);
    
    Dist(i,1)=sum(Length(1,Connect(i,1:lastConn(i,1))));
    %Dist2(i,1)=sum(Sub(i,:).*Length);
%     TimeU(i,1)=sum(Length(1,Connect(i,1:lastConn(i,1)))./...
%         (VelocityU(1,Connect(i,1:lastConn(i,1)))));
    TimeM(i,1)=sum(Length(1,Connect(i,1:lastConn(i,1)))./...
        (VelocityM(1,Connect(i,1:lastConn(i,1)))));
    TimeS(i,1)=sum(Length(1,Connect(i,1:lastConn(i,1)))./...
        (VelocityS(1,Connect(i,1:lastConn(i,1)))));
    TimeG(i,1)=sum(Length(1,Connect(i,1:lastConn(i,1)))./...
        (VelocityG(1,Connect(i,1:lastConn(i,1)))));
end

%% Identify links for highlighting during plotting
% %% Identify order of streams for width function plotting
% %Ord=1 are the leaves
% Ord=1;
% 
% %% Identify low order streams entering high order
% %stream order connectivity
% OrdConnect=Connect;
% OrdConnect(:,:)=NaN;
% for i = 1:LinkNum
%     for j=1:LinkNum
%         if ~isnan(Connect(i,j))
%             OrdConnect(i,j)=Order(1,Connect(i,j));
%         end
%     end
% end
% % 1-3 entering 6-7
% part(1:LinkNum,1:LinkNum)=NaN;
% oc=diff(OrdConnect,1,2);
% for i = 1:LinkNum
%     for j=1:LinkNum-1
%         if oc(i,j)>=3 && OrdConnect(i,j+1)==6
%             part(i,j)=1;
%         elseif oc(i,j)>=4 && OrdConnect(i,j+1)==7
%             part(i,j)=1;
%         end
%     end
% end
% %row indentifies the partition
% [row col]=find(part==1);
% 
% partidx(1:LinkNum,1)=0;
% partidx(row,1)=1;

% %% Identify sub basin for width function plotting
% %subloc=5113;
% subloc=5072;%blue earth
% %subloc=2933;%headwaters
% j=find(GridID==subloc);
% %j=OutletLinkID;
% 
% %% Identify upper peak links
% 
% ans1=TimeS/60/60/24/365/0.175>58.6;
% ans2=TimeS/60/60/24/365/0.175<67.7;
% upi=(ans1+ans2)==2;
% 
% %% Identify mainstem sources
% %Mgrid=[];
% Mind(1:max(size(Mgrid)),1)=NaN;
% for i=1:max(size(Mgrid))
%     Mind(i,1)=find(GridID==Mgrid(i));
% end
% 
% nM=(1:max(size(GridID)))';
% nM(Mind)=0;
% nMind=nM(nM>0);

%% Plot Slope WF
% sptile=prctile(Slope,75);
% Y=Dist(Slope>=sptile)./1000;
% figure
% hold on
% [xL]=Plot_Width_Function(Y,100,'Distance to the outlet, km','b',1);

%% Plot Width Functions for the basin
figure
%axes('FontSize',14)
hold on
%plot width function for entire basin
[xL]=BMS_Plot_Width_Function2(Dist./1000,50,'Distance to the outlet, km','b',1,0);
%% Highlight contributions from a sub-area
% % [xL]=Plot_Width_Function(Dist./1000,[50 150 250 350 450 550 650],'Distance to the outlet, km','b',1);
% % histc(Dist./1000,[0 100 200 300 400 500 600 700]);
% % hold on
% % plot(50,ans(1,1),'.m','MarkerSize',30)
% % plot(150,ans(2,1),'.b','MarkerSize',30)
% % plot(250,ans(3,1),'.g','MarkerSize',30)
% % plot(350,ans(4,1),'.y','MarkerSize',30)
% % plot(450,ans(5,1),'.r','MarkerSize',30)
% % plot(550,ans(6,1),'.r','MarkerSize',30)
% % plot(650,ans(7,1),'.k','MarkerSize',30)
% %ylim([0 700])
% %[xL]=Plot_Width_Function([1 2 2 3 3 3 3 4 4],[1 2 3 4],'Distance to the outlet','b',1);
% %[xL]=Plot_Width_Function([2 3 3.5 3.5 3.5 4 4.5 5 5],[1.3 2.6 4 5.3],'Distance to the outlet','b',1);
% % figure
% % hold on
% %area partition
% % Plot_Area_Function(Dist./1000,xL,'Distance to the outlet, km','b',Area,1);
% 
% 
% %plot width function for sub basin
% %Y=Dist(~isnan(Sub(:,j)))./1000;%Basin outlet
% %Y=(Dist(~isnan(Sub(:,j)))-Dist(Connect(j,2)))./1000;%Sub basin outlet
% %Y=(Dist(HW==1)-Dist(734))./1000;%Sub basin outlet
% 
% %plot order X links
% %Y=Dist(Order==Ord)./1000;
% 
% %plot 1-3 entering 6-7
% %Y=Dist(row)./1000;
% 
% %plot mainstem sources
% % Y=Dist(Mind)./1000;
% % 
% Y=Dist(subbasin==3)./1000;
% 
% [xL]=Plot_Width_Function(Y,xL,'Distance to the outlet, km','r',2);
% %[xL]=Plot_Width_Function(Y,30,'Distance to the outlet, km','b',1);

%% Plot TimeU Width Functions for the basin
% figure
% hold on
% %plot width function for entire basin
% [xT]=Plot_Width_Function(TimeU/60/60/24,100,'Time to the outlet, days','b');
% 
% %plot width function for sub basin
% Y=TimeU(~isnan(Sub(:,j)))/60/60/24;%Basin outlet
% %Y=(TimeU(~isnan(Sub(:,j)))-TimeU(Connect(j,2)))/60/60/24;%Sub basin outlet
% 
% %plot order X links
% %Y=TimeU(Order==Ord)/60/60/24;
% 
% %plot 1-3 entering 6-7
% %Y=TimeU(row)/60/60/24;
% 
% [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, days','r');

%% Plot TimeM Width Functions for the basin
% figure
% %axes('FontSize',14)
% hold on
% %plot width function for entire basin
% [xT]=Plot_Width_Function(TimeM/60/60,100,'Time to the outlet, hrs','b',1);
% 
% % figure
% % hold on
% %area partition
% %Plot_Area_Function(TimeM/60/60,xT,'Time to the outlet, hrs','b',Area,1);
% 
% % scale=ones(max(size(GridID)),1);
% % scale(Mind)=75/100*2427/369;
% % scale(nMind)=25/100*2427/2058;
% % Plot_Area_Function(TimeM/60/60,xT,'Time to the outlet, hrs','b',scale,1);
% % hold on
% % Plot_Area_Function(TimeM(Mind)/60/60,xT,'Time to the outlet, hrs','r',scale(Mind),2);
% % Plot_Area_Function(TimeM(nMind)/60/60,xT,'Time to the outlet, hrs','g',scale(nMind),2);
% 
% %plot width function for sub basin
% %Y=TimeM(~isnan(Sub(:,j)))/60/60;%Basin outlet
% %Y=(TimeM(~isnan(Sub(:,j)))-TimeM(Connect(j,2)))/60/60;%Sub basin outlet
% 
% %plot order X links
% % Ord=7;
% % Y=TimeM(Order==Ord)/60/60;
% 
% %plot 1-3 entering 6-7
% % Y=TimeM(row)/60/60;
% 
% %plot mainstem sources
% % Y=TimeM(Mind)/60/60;
% % 
% % [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, hrs','r',2);
% % 
% % Y=TimeM(nMind)/60/60;
% % 
% % [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, hrs','g',2);


%% Plot TimeS Width Functions for the basin
figure
%axes('FontSize',14)
hold on
%box on
%plot width function for entire basin
[xT]=BMS_Plot_Width_Function2(TimeS/60/60/24/365/0.175,50,'Time to the outlet, yrs','k',1,0);

%% Highlight sub-area
% %title(cat(2,'betaA = ',num2str(betaA),'; betaS = ',num2str(betaS)))
% %
% %[a,b]=Plot_Width_Function(T,100,'Time to the outlet, yrs','b',1);
% %
% % figure
% % hold on
% %area partition
% %Plot_Area_Function(TimeS/60/60/24/365/0.175,xT,'Time to the outlet, yrs','b',Area,1);
% 
% % scale=ones(max(size(GridID)),1);
% % scale(Mind)=95/100*2427/369;
% % scale(nMind)=5/100*2427/2058;
% % Plot_Area_Function(TimeS/60/60/24/365/0.175,xT,'Time to the outlet, yrs','b',scale,1);
% % hold on
% % Plot_Area_Function(TimeS(Mind)/60/60/24/365/0.175,xT,'Time to the outlet, yrs','r',scale(Mind),2);
% % Plot_Area_Function(TimeS(nMind)/60/60/24/365/0.175,xT,'Time to the outlet, yrs','g',scale(nMind),2);
% 
% % N(1:400,1:5)=NaN;
% % xT(1:400,1:5)=NaN;
% % [N(:,1),xT(:,1)]=Plot_Width_Function(Time(:,1),400,'Time to the outlet, yrs','r',1);
% % [N(:,2),xT(:,2)]=Plot_Width_Function(Time(:,2),xT(:,1),'Time to the outlet, yrs','g',2);
% % [N(:,3),xT(:,3)]=Plot_Width_Function(Time(:,3),xT(:,1),'Time to the outlet, yrs','c',2);
% % [N(:,4),xT(:,4)]=Plot_Width_Function(Time(:,4),xT(:,1),'Time to the outlet, yrs','b',2);
% % [N(:,5),xT(:,5)]=Plot_Width_Function(Time(:,5),xT(:,1),'Time to the outlet, yrs','m',2);
% 
% %amplification
% % inc=(xT(73)-xT(25));
% % Y=cat(1,TimeS/60/60/24/365/0.175,TimeS/60/60/24/365/0.175+inc);
% % [xT]=Plot_Width_Function(Y,148,'Time to the outlet, yrs','c',1);
% % Y=cat(1,TimeS(upi)/60/60/24/365/0.175,TimeS(~isnan(Sub(:,j)))/60/60/24/365/0.175+inc);
% % [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, yrs','c',2);
% 
% %plot width function for sub basin
% %Y=TimeS(~isnan(Sub(:,j)))/60/60/24/365/0.175;%+inc;%Basin outlet
% %Y=(TimeS(~isnan(Sub(:,j)))-TimeS(Connect(j,2)))/60/60/24/365/0.175;%Sub basin outlet
% %Y=(TimeS(HW==1)-TimeS(734))/60/60/24/365/0.175;%Sub basin outlet
% 
% %
% %plot order X links
% %Ord=7;
% %Y=TimeS(Order==Ord)/60/60/24/365/0.175;
% 
% %plot 1-3 entering 6-7
% % Y=TimeS(row)/60/60/24/365/0.175;
% 
% %plot mainstem sources
% % Y=TimeS(Mind)/60/60/24/365/0.175;
% % 
% Y=TimeS(subbasin==3)/60/60/24/365/0.175;
% 
% [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, yrs','r',2);
% %[xT]=Plot_Width_Function(Y,30,'Time to the outlet, yrs','b',1);
% % 
% % Y=TimeS(nMind)/60/60/24/365/0.175;
% % [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, yrs','g',2);
% 
% %plot upper peak
% %Y=TimeS(upi)/60/60/24/365/0.175;
% %[xT]=Plot_Width_Function(Y,xT,'Time to the outlet, yrs','g',2);

% %% Plot TimeG Width Functions for the basin
% figure
% hold on
% %plot width function for entire basin
% [xT]=Plot_Width_Function(TimeG/60/60/24/365/0.02,100,'Time to the outlet, yrs','b');
% 
% %plot width function for sub basin
% Y=TimeG(~isnan(Sub(:,j)))/60/60/24/365/0.02;%Basin outlet
% %Y=(TimeG(~isnan(Sub(:,j)))-TimeG(Connect(j,2)))/60/60/24/365/0.02;%Sub basin outlet
% 
% %plot order X links
% %Y=TimeG(Order==Ord)/60/60/24/365/0.02;
% 
% %plot 1-3 entering 6-7
% %Y=TimeG(row)/60/60/24/365/0.02;
% 
% [xT]=Plot_Width_Function(Y,xT,'Time to the outlet, yrs','r');

