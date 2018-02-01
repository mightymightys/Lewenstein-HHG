function [classevents,numevents] = findevents(class,minEphot,maxEphot,Ip,ti,fullpulse)

classevents=cell(2,1);

    maxweight = max(class(:,4)); %heighest weight of all trajectories. We may want to keep only those events that include trajs. with 
                                 %weights of at least ...1e-4 or so of this

   
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% cut at jumps in ionization time (and keep only events with more than
    %  4 points (so there's no confusion between dimensions when taking "length" later...)
    %cutindx = find(diff(class(:,2))>3*(ti(2)-ti(1))); 
    cutindx = find(diff(sign(round(diff(diff(class(:,2))))))==-2)+1;   %indices where jumps occur within the series of
    cutindx = [0; cutindx; length(class)];                             %ioniz. instants for which resollisions have been found
                                                                       %- this is where the diff. "recoll. events" are separated
    k=1;
    for n=2:length(cutindx)
        tmp = class(cutindx(n-1)+1:cutindx(n),:);
        if length(tmp)>4
            ttmpevents{k} = tmp;
            k=k+1;
        end
    end
          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% if one event is partly at the beginning of the scanned ioniz. times
    % and partly at the end (so it's cut in two by the (arbitrary) division
    % of the n-period long time window, then wrap it around
    if (~fullpulse==1)
        if ttmpevents{1}(1,2)==min(ti) && ttmpevents{k-1}(end,2)==max(ti)
           ttmpevents{1}(:,2:3)=ttmpevents{1}(:,2:3) + (max(ti)-min(ti));
           ttmpevents{k-1} = [ttmpevents{k-1};ttmpevents{1}];
           ttmpevents = ttmpevents(2:end);
        end
    end   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now, within each event go along the recollision energies and at local minima,
    %divide once more into events.
    k=1;
    for n=1:length(ttmpevents)
        %if max(ttmpevents{1}(:,4))>1e-5*maxweight
            minindcs = find(diff(sign(diff(ttmpevents{n}(:,1))))==+2); %indices for local minima
            if isempty(minindcs) 
               if max(ttmpevents{n}(:,4))>1e-4*maxweight
                tmpevents{k} =  ttmpevents{n};
                k=k+1;
               end
            else
                eventcutindcs = [0; minindcs; length(ttmpevents{n})];
                for m=2:length(eventcutindcs)
                    tmp = ttmpevents{n}(eventcutindcs(m-1)+1:eventcutindcs(m),:);
                    if max(tmp(:,4))>1e-4*maxweight
                        tmpevents{k} = tmp;
                        k=k+1;
                    end
                end
            end
    end

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% now go along the recollision energies and at local maxima, divide into long/short traj
    k=1;
    for n=1:length(tmpevents)
        maxindx = find(diff(sign(diff(tmpevents{n}(:,1))))==-2)+1; %indices for local maxima
            %there should only be one local max... because you split events
            %at every local minimum. only if there's a real 'inflection point' there may be a problem...
            %also, if there's no actual maximum (i.e. no diff=0 crossing) or only very few points
            %before or after the max then this is an event for which trajectories got too long or smthg
            %and they are just discarded
            % also later check that there are enough points for each trajectory otherwise also discard
        if ~isempty(maxindx)
            [cutoff,mm] = max(tmpevents{n}(maxindx,1));  %should there be >1 candidates, take the index with the max. recoll. energy
            m=maxindx(mm);
            %if (cutoff > minEphot-1.3*Ip +(maxEphot-minEphot)*.02) %only keep events with cutoff slightly higher than what you care for
            if m > 3 && m < length(tmpevents{n})-3 && (cutoff > minEphot-Ip) %only keep events with cutoff slightly higher than what you care for
                classevents{2,k}=tmpevents{n}(1:m-1,:);      %long traj
                classevents{1,k}=tmpevents{n}(m+1:end,:);    %short traj
                k=k+1;
            end
        end
    end
         
    numevents = k-1;