%% BMS_CapacitySlope
% This function determines the capacity of an individual link and updates
% slope accordingly.

% Jon Czuba
% February 17, 2015

%%
%cycle through each link
for i=1:LinkNum
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
            %set their travel time to 0
            P_tt{t,i}(1,1:exc)=0;
            
            %UPDATE Elevations beginning here
            %compute volume inactive
            vstor=sum(P_vol{t,i}(1:exc))./(1-Lp);%m3 vol in storage, porosity Lp=0.4
            
            %determine upstream links
            usid=find(Connect(:,2)==i);
            %include current link and upstream links that are not lakes
            elevid=cat(1,i,usid(Lake(usid)==0));
            %update elevation at upstream end of link
            %volume is placed at upstream end of link and along current link
            %length and upstream link lengths to compute new elev
            Elev(t,i)=mxelevmod(i,1)+(2.*vstor)./sum(B(t,elevid)'.*Length(elevid,1));
            %UPDATE Elevations complete
            
            clear cvol exc vstor usid elevid
        end
    end
end

%compute/update slope
for i=1:LinkNum
    if i==OutletLinkID
        Slope(i,1)=(Elev(t,i)-mnelev(i,1))./Length(i,1);
    else
        Slope(i,1)=(Elev(t,i)-Elev(t,Connect(i,2)))./Length(i,1);
    end
end
Slope(Slope<1e-4)=1e-4;

%%

for i=1:LinkNum
    if isempty(P_loc{t,i})
        continue
    end
actidx=P_storage{t,i}==0;
actvol=sum(P_vol{t,i}(actidx));

% Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
% (-log10(P_d{t,i}(actidx).*1000)./log10(2)))));
% if Dg(t,i)==0.5
%     Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
%         (log10(P_d{t,i}(actidx).*1000)./log10(2)))))./1000;
% end

% Dg(t,i)=((sum(P_vol{t,i}(actidx)./actvol.*...
% ((P_d{t,i}(actidx))))));

Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
    (log10(P_d{t,i}(actidx))./log10(2)))));

% Dg=(2.^(sum(Fdfpsd.*...
%     (log10(Dpsd)./log10(2)))));

actsandidx=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(1));
actsandvol=sum(P_vol{t,i}(actsandidx));
Fs=actsandvol./actvol;

taursg=rho.*R.*g.*Dg(t,i).*(0.021+0.015.*exp(-20.*Fs));

% % cumulative PSD
% Dft=[2,4,8,32,64,128,256]'./1000;%m
% cpsd(1,1:gsclass,1)=NaN;
% for j=1:gsclass
%     actjidx=and(P_storage{t,i}==0,P_d{t,i}<=Dpsd(j));
%     actftjvol=sum(P_vol{t,i}(actjidx));
%     cpsd(1,j)=actftjvol./actvol;
% end
% 
% % compute D84
% Dxx=0.84;
% D84=NaN;
% for j=1:gsclass-1
%     if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
%         D84=exp(interp1(cpsd(1,j:j+1),log(Dft(j:j+1,1)),Dxx));
%     end
% end
% %DXX(isnan(DXX))=Dpsd(1);
% 
% % Schneider 2016 WRR
% a1=6.5;
% a2=3.97;
% e=1.5;
% 
% num=a2.*(H(t,i)./D84).^(5/6);
% den=sqrt(a1.^2+a2.^2.*(H(t,i)./D84).^(5/3));
% 
% Sred=Slope(i,1).*(num./den).^e;
% tau=rho.*g.*H(t,i).*Sred;

tau=rho.*g.*H(t,i).*Slope(i,1);

for k=1:gsclass
    actidxj=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(k));
    actvolj=sum(P_vol{t,i}(actidxj));
    Fj=actvolj./actvol;
    
    b=0.67./(1+exp(1.5-Dpsd(k)./Dg(t,i)));
    
    %v moved out of k loop
    %tau=rho.*g.*H(t,i).*Slope(i,1);
    taurj=taursg.*(Dpsd(k)./Dg(t,i)).^b;
    tautaurj=tau./taurj;
    
    if tautaurj<1.35
        Wj=0.002.*tautaurj.^7.5;
    else
        Wj=14.*(1-0.894./sqrt(tautaurj)).^4.5;
    end
    
%     if k==1
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta(i)./Wj./tau.^(3/2)./Fs;
%     else
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta(i)./Wj./tau.^(3/2)./(1-Fs)./Fj;
%     end
    
%     if k==1
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./Fs;
%     else
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./(1-Fs)./Fj;
%     end
    
    P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./Fj;
    
    clear actidxj actvolj Fj b taurj tautaurj Wj
    
end

clear actidx actvol actsandidx actsandvol Fs taursg tau a1 a2 e num den Sred
end

%%
% limit=ceil(capacity./pvol);
% 
% % determine which inputs are active during current timestep
% n = cellfun(@length,P_loc(t,:));
% 
% excess=find((n'-limit)>0);
% 
% if ~isempty(excess)
%     for k = 1:numel(excess)
%         
%         % First In First Out
%         [srt, sind]=sort(P_loc{t,excess(k)},'descend');
%         %include storage set to 0; no storage set to 1
%         %P_active{t,excess(k)}(1,sind(1:limit(excess(k))))=0;
%         P_active{t,excess(k)}(1,sind(limit(excess(k))+1:end))=0;
%         
%         % First In Last Out
% %        P_active{t,excess(k)}(1,limit(1:excess(k)))=0;
%         
%         
%         clear srt sind
%     end
% end

