function [EmAmplshort_out,EmAmpllong_out] = f_killdivergence(emampls,qevents,Ephot) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WORKING VERSION
% We have to discard one class of trajectories beyond the cutoff.
% we'll compare the emampls going from maxEphot downwards - one will be
% huge, the other tiny. From the energy where the two cross on upwards, we'll kill the first. 
% well, it turns out, they don't necessarily cross... so we'll take the
% highest local minimum of their difference!

Estp = Ephot(2)-Ephot(1);
shortweight=emampls{1}.*conj(emampls{1});   longweight=emampls{2}.*conj(emampls{2});

   if mean(abs(qevents{1}(:,1)-qevents{2}(:,1))) < 1e-6
                %if fsolve fell actually on the same traj for long and short...
         if (mean(diff(shortweight))>0)
            emampls{1} = zeros(size(emampls{1})); 
         elseif (mean(diff(longweight))>0)
             emampls{2} = zeros(size(emampls{2})); 
         end
      EmAmplshort_out=emampls{1}; EmAmpllong_out=emampls{2}; return
   elseif mean(abs(qevents{1}(end-100:end,1)-qevents{2}(end-100:end,1)))<1e-6 && mean(diff(shortweight(end-100:end)))>0
        %if fsolve fell actually on the same traj for long and short only beyond the cutoff and the emapmls are diverging...
        cutind=find(abs(abs(qevents{1}(:,1))-abs(qevents{2}(:,1)))>1e-6,1,'last');
        for traj=1:2
          emampls{traj}(cutind:end) = zeros(size(emampls{traj}(cutind:end))); 
        end
        EmAmplshort_out=emampls{1}; EmAmpllong_out=emampls{2}; return
   elseif min(abs(qevents{1}(:,1)-qevents{2}(:,1))) > 1
                %if there's simply no crossing of the ioniz. times most likely when
                %fsolve fell on two different traj. families...
         if shortweight(end)>1e10*longweight(end) || shortweight(1)>1e10*longweight(1)
            emampls{1} = zeros(size(emampls{1}));
            if mean(diff(longweight(end-10:end)))>0;  emampls{2} = zeros(size(emampls{2})); end
         elseif longweight(end)>1e10*shortweight(end) || longweight(1)>1e10*shortweight(1)
            emampls{2} = zeros(size(emampls{2}));
            if mean(diff(shortweight(end-10:end)))>0;  emampls{1} = zeros(size(emampls{1})); end
         else
            if shortweight(end)>1e-2 || shortweight(1)>1e-2
                emampls{1} = zeros(size(emampls{1})); 
            elseif longweight(end)>1e-2 || longweight(1)>1e-2
                emampls{2} = zeros(size(emampls{2})); 
            end
         end
      EmAmplshort_out=emampls{1}; EmAmpllong_out=emampls{2}; return
   else
       %cutind = find(diff(sign(shortweight-longweight)), 1, 'last' )+1; %highest index (energy) where the traj.weights of short and long are equal
       cutind = find(diff(sign(diff(abs(shortweight-longweight))))==+2, 1, 'last' )+1; %highest index (energy) where the traj.weight difference
                                                                                       %of short and long goes through a local minimum
       if isempty(cutind), cutind=1; end
       %just to be safe, we also do the "old way", looking for a ti-crossing and then a "kink" in the ti-curves below this crossing.
       %out of the two indices found like that, we then take the higher one 
       %    xcutind = find(diff(sign(qevents{1}(:,1)-qevents{2}(:,1))),1);
       %                        %find the indx when the ionization times cross
       %    if isempty(xcutind),  xcutind = length(Ephot)-1;  end
       %    [~,cutind2] = min( diff( abs(qevents{1}(1:xcutind+1,1)-qevents{2}(1:xcutind+1,1)) ) );
       %cutind=max(cutind1,xcutind);
       %cutind=xcutind;
       if (sign(shortweight(find(~isinf(shortweight),1,'last'))-longweight(find(~isinf(longweight),1,'last')))==1)
           traj=1; othertraj=2;
       else
           traj=2; othertraj=1;
       end
       %emampls{traj}(cutind:end) = exp(-.0005*Estp*(0:length(Ephot)-cutind)'.^2) .*emampls{traj}(cutind:end);
       %emampls{traj}(cutind:end) = exp(-1*Estp*(0:length(Ephot)-cutind)'.^2) *emampls{traj}(cutind);
       emampls{traj}(cutind:end) = emampls{othertraj}(cutind:end)/abs(emampls{othertraj}(cutind))*abs(emampls{traj}(cutind));
%               while (max(diff(emampls{traj}(cutind:end).*conj(emampls{traj}(cutind:end)))>0))
%                  %as long as anywhere beyond cutind, the em.ampl. increases again,
%                  %we add another exp.drop from that point on
%                  %cutagainind = find((diff(emampls{traj}(cutind:end).*conj(emampls{traj}(cutind:end)))>0),1);
%                  cutagainind =1;
%                  %emampls{traj}(cutind+cutagainind-1:end) = exp(-.005*Estp*(0:length(Ephot)-cutind-cutagainind+1)'.^2) .* emampls{traj}(cutind+cutagainind-1:end);
%                  emampls{traj}(cutind+cutagainind-1:end) = zeros(size(emampls{traj}(cutind+cutagainind-1:end)));
%               end
   end
   
   for traj=1:2
    if mean(diff(emampls{traj}(end-10:end).*conj(emampls{traj}(end-10:end))))>0
    %if emampls{traj}(end).*conj(emampls{traj}(end))>emampls{traj}(end-1).*conj(emampls{traj}(end-1))
        %if we missed something somehow and the emampls still diverged in the end,
        %we have to kill off the whole trajectory... this is rude, but I don't know better...
       emampls{traj} = zeros(size(emampls{traj})); 
    end
    if max(emampls{traj}(:).*conj(emampls{traj}(:)))>1e-2
        %if we missed something somehow and the emampls diverge to huge
        %values anywhere we have to kill off the whole trajectory... this is rude, but I don't know better...
       emampls{traj} = zeros(size(emampls{traj})); 
       disp(['I had to kill off traj; number',num2str(traj),'because it diverged. Sorry.'])
    end
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check the low-energy end - when you have a local minimum splitting a
% continuous "event-curve" into two then you get the same effects with the
% em.amplitudes correctly falling exponentially to zero on one
% side, but diverging on the other. You have to kill of this divergence, too.

    for traj=1:2
        %if the traj.weight is max at lowest energy, and this value is >1000x higher than the max.weight of the other trajectory, then I
        %suppose that we have a divergence towards low energies.
        if traj==1
            if emampls{traj}(1)*conj(emampls{traj}(1))==max(emampls{traj}.*conj(emampls{traj}(:))) && emampls{traj}(1)*conj(emampls{traj}(1)) >1e3*max(emampls{2}.*conj(emampls{2})) 
                yes=1;
            else
                yes=0;
            end
        else
            if emampls{traj}(1)*conj(emampls{traj}(1))==max(emampls{traj}.*conj(emampls{traj}(:))) && emampls{traj}(1)*conj(emampls{traj}(1)) >1e3*max(emampls{1}.*conj(emampls{1})) 
                yes=1;
            else
                yes=0;
            end
        end
        if yes
        %if (diff(emampls{traj}(1:2)*conj(emampls{traj}(1:2)))<-1e-3) %this means that in the very beginning (lowest energy step),
                                                         %the em. amplitude decreases much too fast - i.e.
                                                         %it has diverged from some higher energy on downwards
        %if (emampls{traj}(1)*conj(emampls{traj}(1))>1e-2) %the emission amplitude is absurdly big at the lowest energy ("normal" values are 1e-5 or so)
             tmp = (diff(abs(diff((qevents{traj}(:,1)))))); % where this has peaks, the ti(Ephot)-curve has a kink
                                                             %- this happens at places where the divergences occur
                 ttmp = sort(tmp); tttmp = ttmp(end-20); %take the 20th biggest value of tmp
             tmpindcs = find(tmp>tttmp);
             ttmpindcs=find(diff(sign(diff(tmp(tmpindcs,1))))==-2);  %of the points in tmp where values are bigger than tttmp
                                                               %we pick out the local maxima
             tmpindcs = tmpindcs(ttmpindcs);
             for n=tmpindcs'
                while n>0 && diff(emampls{traj}(n:n+1).*conj(emampls{traj}(n:n+1)))>0
                    n=n-1;    %from the found local maxima, go downwards until the slope of the emampls is negative (i.e. increases towards lower energies)
                end
                if ( max(diff(emampls{traj}(1:n).*conj(emampls{traj}(1:n))))<0 )
                    cutind=n; %if before the indx with the "kink" in the ti curve, the em-ampl. decreases *everywhere* with Ephot,
                           %write down the indx, and move to the next one to test. Once we've tested all indcs n,
                           %we have the highest one written down in "cutind"
                end
             end
             %now we have an indx, "cut", before which we must kill off the diverging traj
             %emampls{traj}(1:cutind) = exp(-.005*Estp*(-cutind+1:0)'.^2) .*emampls{traj}(1:cutind);
             emampls{traj}(1:cutind) = exp(-1*Estp*(-cutind+1:0)'.^2) *emampls{traj}(cutind);
             %while (max(diff(emampls{traj}(1:cutind).*conj(emampls{traj}(1:cutind)))<0))
                 %as long as anywhere before cutind, the em.ampl. still drops,
                 %we add another exp.increase from0to1 from that point on
             %    cutagainind = find((diff(emampls{traj}(1:cutind).*conj(emampls{traj}(1:cutind)))<0), 1, 'last' );
             %    emampls{traj}(1:cutagainind) = exp(-.01*Estp*(-cutagainind:-1)'.^2) .* emampls{traj}(1:cutagainind);
             %end
        end
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output
EmAmplshort_out = emampls{1};
EmAmpllong_out  = emampls{2};
