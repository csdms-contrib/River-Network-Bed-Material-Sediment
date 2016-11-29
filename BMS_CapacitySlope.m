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
            Elev(t,i)=mxelevmod(i,1)+(2.*vstor)./sum(B(elevid,1).*Length(elevid,1));
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

