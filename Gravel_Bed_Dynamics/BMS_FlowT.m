%% BMS_FlowT

%%
% Q=[];%m3/s
% YMD=[];
%%
% DN=datenum(YMD);
%% plot timeseries
% figure; hold on; box on
% %axes('YScale','log')
% plot(DN,Q)
% plot(DN,Qgage)
% datetick('x')
% ylabel('Daily streamflow, m^3/s')
%%
% figure; hold on; box on
% plot((1:length(Qgage))./365,Qgage)
% ylabel('Daily streamflow at gage, m^3/s')
% xlabel('Time since debris flow input, years')
% xlim([0 5])

%%
% Tushar
%load Tushar\TusharQ.mat
load Tushar\TusharQ2.mat
%DN=datenum(YMD);
%July 8, 2011
%19700
Qgage=Q;
% Qgage=Q(19700:end,1);
% YMD=YMD(19700:end,:);
% DN=DN(19700:end,1);
%1970-1990
% Qgage=Q(4537:11841,1);
% YMD=YMD(4537:11841,:);
% DN=DN(4537:11841,1);
% 
% Qgage=cat(1,Qgage,Qgage);

%Methow
% load Methow\MethowQ.mat
% %Oct 1, 1969 -- 3654
% %Sept 30, 1974 -- 5479
% %UNCOMMENT FOR 5 years
% % Qgage=Q(3654:5479,1);
% % YMD=YMD(3654:5479,:);
% % DN=DN(3654:5479,1);
% %Flow for 30 days
% Qgage=Q(3654:3683,1);
% YMD=YMD(3654:3683,:);
% DN=DN(3654:3683,1);

%Nisqually
% load C:\Users\jczuba\Documents\Projects\Nisqually\MATLAB\NisquallyQ.mat
% % 1986
% Qgage=Q(15921:16285,1);
% YMD=YMD(15921:16285,:);
% DN=DN(15921:16285,1);

% % CY [1986-2015], 30 years
% Qgage=Q(15921:26877,1);
% YMD=YMD(15921:26877,:);
% DN=DN(15921:26877,1);
%%
% Tushar
Bgage=5.1865.*Qgage.^0.1413;
Hgage=0.2426.*Qgage.^0.3951;
%Ugage=Qgage./Bgage./Hgage;
%Frgage=Ugage./sqrt(g.*Hgage);
Agage=4.248e8;%m2

% Methow
% Bgage=30.906.*Qgage.^0.1215;
% Hgage=0.2703.*Qgage.^0.3447;
% %Ugage=Qgage./Bgage./Hgage;
% %Frgage=Ugage./sqrt(g.*Hgage);
% Agage=4.5895e+9;%m2

% Nisqually
% Bgage=17.355.*Qgage.^0.1596;
% Hgage=0.1902.*Qgage.^0.3897;
% %Ugage=Qgage./Bgage./Hgage;
% %Frgage=Ugage./sqrt(g.*Hgage);
% Agage=3.4777e+08;%m2

%%
B(1:timesteps,1:LinkNum)=NaN;
H(1:timesteps,1:LinkNum)=NaN;

% Tushar
B=(repmat(20,1,LinkNum)./(Agage.^0.5)).*repmat(usarea',timesteps,1).^0.5;
H=(repmat(Hgage,1,LinkNum)./(Agage.^0.4)).*repmat(usarea',timesteps,1).^0.4;
B(B<2)=2;
Btmax=max(B,[],1)';
%B=repmat(Btmax',timesteps,1);

% Methow
% B=(repmat(Bgage,1,LinkNum)./(Agage.^0.5)).*repmat(usarea',timesteps,1).^0.5;
% H=(repmat(Hgage,1,LinkNum)./(Agage.^0.4)).*repmat(usarea',timesteps,1).^0.4;
% % B(B<2)=2;
% Btmax=max(B,[],1)';
% % B=repmat(Btmax',timesteps,1);

% figure
% plot(usarea_km,Btmax,'.b')
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% %%
% Qgage=QQ(qq);
% %Qgage=6.4;%160;%54;%6.4;%1.7;%m3/s
% %Qgage=gageDMYQ(dmyq,4);
% %Q=gageDMYQ(dmyq,4)./usarea(472,1).*usarea;
% %Q=400./usarea(472,1).*usarea;
% Q=Qgage./usarea(472,1).*usarea;
% 
% %% Assign Q to ensure continuity from tributaries
% hw=(sum(Sub,1)==1);%identifies headwater, first-order links
% Q(~hw)=NaN;
% q=Q;
% 
% tocomp=sum(isnan(Q));
% while tocomp>0
%     for i=1:LinkNum
% 
%         if isnan(Q(i)) &&  ~isnan(sum(Q(Connect(:,2)==i)))
%             %add incremenal Q via Area and upstream Qs
%             q(i)=Qgage./usarea(472,1).*Area(i);
%             Q(i)=q(i)+sum(Q(Connect(:,2)==i));
%         end
%         
%     end
%     tocomp=sum(isnan(Q));
% end
% 
% %% Assign B
% %assign B at gage based on rating curve
% B(1:LinkNum,1)=NaN;
% if Q(472)<Qbf
%     B(472)=a1.*Q(472).^b1;
% else
%     B(472)=a2.*Q(472).^b2;
% end
% 
% %scale B throughout the basin
% B=(B(472)./(usarea(472).^(0.5))).*usarea.^(0.5);
% 
% %% Determine U, H, and wetland hydraulics
% n=0.035; %Manning's roughness
% %n=nval(nn);
% g=9.81; %m/s2 - acceleration due to gravity
% 
% U=(1/n.*(Q./B).^(2/3).*Slope.^(1/2)).^(3/5);
% H=Q./U./B;
% Fr=U./sqrt(g.*H);