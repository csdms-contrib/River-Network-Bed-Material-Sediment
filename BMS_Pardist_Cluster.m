%% _BMS_Pardist_Cluster
%tic
% initialize
steps=timesteps;
Clst_spdist(1:steps,1:LinkNum)=NaN;
Clst_parconc(1:steps,1:LinkNum)=NaN;
numclst(1:steps,1)=NaN;
maxclst_spdist(1:steps,1)=NaN;
medclst_spdist(1:steps,1)=NaN;
ntclst_spdist(1:steps,1)=NaN;
meanclst_spdist(1:steps,1)=NaN;
 
LCS_spdist(1:steps,1:LinkNum)=NaN;%the sum of this through time gives CPI
LCS_parconc(1:steps,1:LinkNum)=NaN;

inputs=pidx;
PIDlink(1:steps,1:inputs)=NaN;
%cont594(1:steps,1:inputs)=NaN;%parcel sources affecting CPI at 594

Loc=P_loc;

% parcel ID
PID=Loc;
PAR_spdist=Loc;
PAR_parconc=Loc;
for t=1:steps
    n=1;
    for i=1:LinkNum
        if ~isempty(PID{t,i})
            PID{t,i}= n:n+numel(PID{t,i})-1;
            PIDlink(t,n:n+numel(PID{t,i})-1)=i;
            PAR_spdist{t,i}=linspace(0,0,numel(PID{t,i}));
            PAR_parconc{t,i}=linspace(0,0,numel(PID{t,i}));
            n=n+numel(PID{t,i});          
        end
    end
end

% compute distribution of parcel lengths
pardist(1:steps,1:inputs)=NaN;
pdid(1:inputs,1:2)=NaN;
for t=1:steps
    %t
    %go through each link and determine distance between the particle(s)
    %in that link with others in the system
    p=1;
    clear pdid
    pdid(1:inputs,1:2)=NaN;

    for i=1:LinkNum
        if isempty(Loc{t,i})
            continue
        end
        [slink, slid]=sort(Loc{t,i});       
        for j=1:numel(Loc{t,i})
            %find nearest DS parcel
            pdid(p,1)=PID{t,i}(slid(j));
            if j<numel(Loc{t,i}) %if in same link
                %pardist is length between parcels in link
                %pardist(1,p)=(Loc{t,i}(1,j+1)-Loc{t,i}(1,j)).*Length(i);
                pardist(t,p)=(slink(1,j+1)-slink(1,j)).*Length(i);
                pdid(p,2)=PID{t,i}(slid(j+1));
                p=p+1;
            else %search next DS link
                %pardist is the distance through the rest of the link +
                %pardist(1,p)=(1-Loc{t,i}(1,j)).*Length(i);
                pardist(t,p)=(1-slink(1,j)).*Length(i);
                %next downstream link
                nds=Connect(i,2);
                %if at outlet end search empty
                if i==OutletLinkID
                    pardist(t,p)=NaN;
                    continue
                end
                while isempty(Loc{t,nds}) %cycle through empty links
                    %add to pardist empty link lengths
                    pardist(t,p)=pardist(t,p)+Length(nds);
                    %next downstream link
                    nds=Connect(nds,2);
                    %if at outlet end search empty
                    if isnan(nds)
                        pardist(t,p)=NaN;
                        break
                    end
                end
                %if at outlet end search empty
                if isnan(nds)
                    pardist(t,p)=NaN;
                    continue
                end
                %add to pardist distance in link up to this parcel
                %pardist(1,p)=pardist(1,p)+Loc{t,nds}(1,1).*Length(nds);
                %pardist(t,)=pardist(t,p)+min(Loc{t,i}).*Length(nds);
                [mloc, mlocid] = min(Loc{t,nds});
                pardist(t,p)=pardist(t,p)+mloc.*Length(nds);
                pdid(p,2)=PID{t,nds}(mlocid);
                p=p+1;
            end
        end
    end
    
    % associate parcels with clusters
    
    clear pdidth pardistth
    dthresh=3000;%exp(8);%m
    % find 2 parcel sets within threshold distance
    pdidth=pdid(pardist(t,:)<=dthresh,:);
    pardistth=pardist(t,pardist(t,:)<=dthresh)';
    Cpid(1:inputs,1:inputs)=NaN;
    Cspdist(1:inputs,1:inputs)=NaN;
    Csplink(1:inputs,1:LinkNum)=NaN;
    i=1;
    while ~isempty(pdidth)
        clstpid=[];
        cv2s=[];
        clstspdist=[];
        clstsplink=[];
        % cycle through all sets
        % assign first sets to a cluster
        clstpid=cat(2,clstpid,pdidth(1,:));
        cv2s=cat(2,cv2s,pdidth(1,:));
        clstspdist=cat(2,clstspdist,pardistth(1,1));
        clstsplink=cat(2,clstsplink,...
            Connect(PIDlink(t,pdidth(1,1)),...
            1:find(PIDlink(t,pdidth(1,2))==Connect(PIDlink(t,pdidth(1,1)),:))));
        
        % remove assigned set from remaining sets
        pdidth(1,:)=[];
        pardistth(1,:)=[];
        
        while ~isempty(cv2s)
            % cycle through cluster
            [r,c]=find(pdidth(:,:)==cv2s(1));
            if ~isempty(r)
                for k=1:numel(r)
                    % assign next set to a cluster
                    if c(k)==2
                        clstpid=cat(2,clstpid,pdidth(r(k),1));
                        cv2s=cat(2,cv2s,pdidth(r(k),1));
                    elseif c(k)==1
                        clstpid=cat(2,clstpid,pdidth(r(k),2));
                        cv2s=cat(2,cv2s,pdidth(r(k),2));
                    end
                    clstspdist=cat(2,clstspdist,pardistth(r(k),1));
                    clstsplink=cat(2,clstsplink,...
                        Connect(PIDlink(t,pdidth(r(k),1)),...
                        1:find(PIDlink(t,pdidth(r(k),2))==Connect(PIDlink(t,pdidth(r(k),1)),:))));
                 end
                % remove assigned set from remaining sets
                pdidth(r,:)=[];
                pardistth(r,:)=[];
            end
            cv2s(1)=[];
        end
        Cpid(i,1:numel(clstpid))=clstpid;
        Cspdist(i,1:numel(clstspdist))=clstspdist;
        Csplink(i,1:numel(clstsplink))=clstsplink;
        i=i+1;
    end
    
%     if t==2
%         a=1;
%     end
    % compute cluster statistics
    numclst(t,1)=sum(~isnan(Cspdist(:,1)));
    Clst_spdist(t,1:numclst(t,1))=nansum(Cspdist(1:numclst(t,1),:),2);
    Clst_parconc(t,1:numclst(t,1))=nansum(~isnan(Cspdist(1:numclst(t,1),:)),2)+1;
    maxclst_spdist(t,1)=nanmax(Clst_spdist(t,:));
    medclst_spdist(t,1)=prctile(Clst_spdist(t,Clst_spdist(t,:)>0),50);
    meanclst_spdist(t,1)=mean(Clst_spdist(t,Clst_spdist(t,:)>0),2);
    ntclst_spdist(t,1)=prctile(Clst_spdist(t,Clst_spdist(t,:)>0),90);
    clear dist didx
    [dist, didx]=sort(Clst_spdist(t,:)');
    for j=1:numclst(t,1)
        LCS_spdist(t,Csplink(didx(j),~isnan(Csplink(didx(j),:))))=...
            Clst_spdist(t,didx(j));
        LCS_parconc(t,Csplink(didx(j),~isnan(Csplink(didx(j),:))))=...
            Clst_parconc(t,didx(j));
if t==7
    a=1;
end
        for k=1:sum(~isnan(Cpid(j,:)))         
            PAR_spdist{t,PIDlink(t,Cpid(j,k))}...
                (find(PID{t,PIDlink(t,Cpid(j,k))}==Cpid(j,k),1,'first'))...
                =Clst_spdist(t,j);       
            PAR_parconc{t,PIDlink(t,Cpid(j,k))}...
                (find(PID{t,PIDlink(t,Cpid(j,k))}==Cpid(j,k),1,'first'))...
                =Clst_parconc(t,j);
        end
              
    end
%     if ~isempty(Cpid(find(sum(Csplink==594,2)>0),:))
%         cont594(t,:)=Cpid(find(sum(Csplink==594,2)>0),:);
%         for pp=1:sum(~isnan(cont594(t,:)))
%             cont594(t,pp)=PIDlink(t,cont594(t,pp));
%         end
%     end
end
clear PID PIDlink Cpid Cspdist Csplink clstpid clstspdist clstsplink didx dist pdid
%toc