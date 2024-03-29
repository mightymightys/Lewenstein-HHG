%clear all;
%close all;
pause on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation parameters

if ~exist('nuclfactor','var')
    nuclfactor =0;
end

   if Ip>minEphot
       minEphot=Ip+.1;
    end


I0=intensity / 3.5094452E16; % total intensity in a.u.

lambda=basewavelength*1e-9; %base-wavelength corresponding to base-freq. omega

tau=pulseduration *1000/24.2;  %FWHM duration of the laser pulse, in a.u.

t_0_au = t_0 *1000/24.2	 %relative delay of pulses in a.u.

tlim=pi/2 *tau /2 /acos(2^(-0.25)) ;  % limits of the time-windows, where the cos^2-enevelopes are zero

omega=2*pi*299792458/lambda *24.2e-18; %laser freq. in a.u.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% queries

if ~exist('choice','var')
    choice = mymenu('How many CPU cores do you want to use?',{'One','Two','Three','Dunno, you decide.'}, 2);
end
if choice == 1
    maxNumCompThreads(1);
%     if matlabpool('size') == 0 || matlabpool('size') == 1
%     else
%         matlabpool close force local 
%     end
 elseif choice == 2
     maxNumCompThreads(2);
%     if matlabpool('size') == 2
%     elseif matlabpool('size') == 0
%         matlabpool open local 2
%     else
%         matlabpool close force local 
%         matlabpool open local 2
%     end
 elseif choice == 3
     maxNumCompThreads(3);
%     if matlabpool('size') == 3
%     elseif matlabpool('size') == 0
%         matlabpool open local 3
%     else
%         matlabpool close force local 
%         matlabpool open local 3
%     end
 elseif choice == 4
      maxNumCompThreads('automatic') 
%     if matlabpool('size') == 4
%     elseif matlabpool('size') == 0
%         matlabpool open local 4
%     else
%         matlabpool close force local 
%         matlabpool open local 4
%     end
else
     maxNumCompThreads('automatic') 
%     if matlabpool('size') == 2
%     elseif matlabpool('size') == 0
%         matlabpool open local 2
%     else
%         matlabpool close force local 
%         matlabpool open local 2
%     end
 end
% clear choice
% pause(.5)
disp(['Using max. ',num2str(maxNumCompThreads),' cpu cores.']);
%saveimgs = mymenu('Do you want to save all the figures?',{'sure.','nope.'}, 1);
%pause(.5)

%[savename,savepath] = uiputfile('results/*.mat','When done, save results as:');
disp('Results will be saved in:')
disp(savepath)
disp(savename)


excursion=cell(2,2,length(CEP));
for nCEP=1:length(CEP)
tic

freqs=ffreqs(CEP(nCEP),phi2, t_0_au);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix up the functions for the set parameters (so as to avoid stupid global variables)
classtraj = @(ti,tt,aalpha)fclasstraj(ti,tt,aalpha,freqs,omega,I0,tau,tlim); 
Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim); 
vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 
psalpha = @(ta,tb)fpsalpha(ta,tb,I0,tau,omega,freqs,tlim); 
emissionampl = @(Ephotindx,Ephot,qevents,ps,acts)f_emissionamplitude(Ephotindx,Ephot,qevents,ps,acts,I0,tau,omega,freqs,Ip,tlim,dip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classical trajectories for initial guesses
    disp(['CEP = ',num2str(CEP(nCEP)),'pi']);
    disp('finding classical trajectories...')
%tic    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up time vectors - decide which part of the pulse you want to consider

    tstp=.1; % timestep in a.u. with which to scan the birth times and calc. trajectories
    
    [n,m]=size(freqs);
    if n>1
        freqssorted = sort(freqs(:,1));
        waveperiod = 1/(freqssorted(2)-freqssorted(1)) * 2*pi/omega;
    else waveperiod = 2*pi/omega/freqs(1,1);
    end
    %lowestperiod = 1*2*pi/omega/min(freqs(:,1)); %x*one period of the lowest freq.component
    maxduration = longestexcursion*waveperiod;
       if fullpulse==1
           if (t_0_au > 0)
           t=(-2*tlim/3:tstp:2*tlim/3+maxduration+t_0_au);
           end
           if (t_0_au <= 0)
           t=(-2*tlim/3+t_0_au:tstp:2*tlim/3+maxduration);
           end
           %t=(-tlim:tstp:tlim+maxduration);  %(until x*periods of the lowest freq. component
                                                %after the end of the pulse. this also effectively limits the 
                                                %excursion times to be
                                                %considered to x*periods of the lowest freq. component
       else
           if (~exist('Nperiods','var'))
               Nperiods =1;
           end
           t=(-0.5*Nperiods*waveperiod:tstp:0.5*Nperiods*waveperiod+maxduration);
       end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make vector with integrated vector potential (alpha) (to be used by classtraj)
        alpha=zeros(1,length(t));
         for n=2:length(t)
            alpha(n) = alpha(n-1) + trapz(fvecpot(t(n-1:n),I0,tau,omega,freqs,tlim))*tstp;
                                                            %here, when you advance in very small time steps,
                                                            %trapz is totally good enough and muuuch faster. Even just summing up would do...
         end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find recolliding trajectories
    
    ti=(min(t):tstp:max(t)-maxduration); %birth times to be scanned 
    
    % can't really initialize the results-vector 'class'... dunno how long
    % it'll be. But, it's important to at least set it to a single line of
    % four zeros, so you don't drag along recollisions from the former
    % iterations!
    
    %cut the time-indcs to go through into three slices for parallel work
    loopindcs{1}=1:floor(length(ti)/3);
    loopindcs{2}=floor(length(ti)/3)+1:2*floor(length(ti)/3);
    loopindcs{3}=2*floor(length(ti)/3)+1:length(ti);
    
  parfor l=1:3
    classresult{l}=zeros(1,4); 
    m=1;
    %for n=1:length(ti) 
    for n=loopindcs{l}
        % here, we will for every birth time in ti calculate the class.
        % trajectory for [ti,ti+1period of the lowest freq. comp.], and
        % look for zero crossings
        aalpha=alpha(n+20:length(t)-length(ti)+n)-alpha(n);  % take the piece of alpha for [ti+20*tstp,ti+x*period of the lowest freq. comp.]
                                                          % and subtract alpha(ti) to get \int_ti^t A(t')dt'
                                                          % (we start a bit after ti to avoid finding "too early" crossings right around
                                                          % the ioniz. instant
        tt=t(n+20:length(t)-length(ti)+n);          % take [ti+20*tstp,ti+x*period of the lowest freq. comp.]
        [~,Ti] = crossings(classtraj(ti(n),tt,aalpha),tt,0,'int');  %find zero-crossings of classical trajectory for birth time ti(n)
        if (~isempty(Ti)) && (.5*(vecpot(min(Ti))-vecpot(ti(n)))^2)>0.15 
                                  % if the the output of the zero-crossing search is NOT empty,
                                  % and the recol. energy is >x a.u. (at low energies, all this is inaccurate anyway)
         classresult{l}(m,3)= min(Ti);     % write down the time value when the zero-crossing (=recollision) occurs
         classresult{l}(m,2)= ti(n);       % also write down the associated ionization time
         classresult{l}(m,1)= .5*(vecpot(min(Ti))-vecpot(ti(n)))^2; % also write down the recollision energy =0.5*(A(tr)-A(ti))^2
         classresult{l}(m,4)= 1E15*(exp(-2*(2*Ip)^(3/2)./abs(3*Efield(ti(n)))).*((classresult{l}(m,3)-ti(n)).^(-5))); %write down a "trajectory weight"
                                                                                                    %(already including a certain penalty for longer excursion times)
         m=m+1;                   % increase index of the "results"-vector by one for next zero crossing
        end 
    end
  end
class = [classresult{1}; classresult{2}; classresult{3}];
clear m n aalpha tt I Ti dump classresult

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% split classical traj. in "recollision events" (=one "contious" series of recollisions, incl all corresponding pairs
% of short and long traj.), and then split within each event in short and long.
% keep only those events with cutoff > x a.u.
    [classevents,numevents] = findevents(class,minEphot,maxEphot,Ip,ti,fullpulse);
%toc    

for i=1:numevents
    tmp(i) = max(classevents{1,i}(:,1));
    %disp(i)
end
checkcutoff = max(tmp)+1.3*Ip;
if checkcutoff>maxEphot
    oldmaxEphot=maxEphot;
    maxEphot=(ceil(checkcutoff*27.2*1.1)) / 27.2;
    disp('*********** !!! ***********');
    disp(['I found a classical cutoff at ',num2str(checkcutoff*27.2),' eV.']);
    disp(['You set the max. photon energy to consider, maxEphot=',num2str(round(oldmaxEphot*27.2)),' eV,']);
    disp('lower than this, which will probably lead to problems in the calculation.')
    disp(['I thus take the liberty to increase maxEphot to ',num2str(maxEphot*27.2),' eV.'])
    disp('Deal with it.');
    disp('*********** !!! ***********');
    
end
clear tmp oldmaxEphot checkcutoff

%    figure
%      hold on;
%      for n=1:length(class)
%          plot(class(n,2)*24.2/1000,class(n,1)*27.2,class(n,3)*24.2/1000,class(n,1)*27.2)
%      end
%     for n=1:numevents
%          plot(classevents{1,n}(:,3)*24.2/1000,classevents{1,n}(:,1)*27.2,'.',classevents{2,n}(:,3)*24.2/1000,classevents{2,n}(:,1)*27.2,'.' );
%          [y,ind] = max(classevents{1,n}(:,1));
%          text(classevents{1,n}(ind,3)*24.2/1000, y*27.2+5, num2str(n), 'HorizontalAlignment','center');
%     end
%     for n=1:numevents
%          plot(classevents{1,n}(:,2)*24.2/1000,classevents{1,n}(:,1)*27.2,'.',classevents{2,n}(:,2)*24.2/1000,classevents{2,n}(:,1)*27.2,'.' );
%          [y,ind] = max(classevents{1,n}(:,1));
%          text(classevents{1,n}(ind,2)*24.2/1000, y*27.2+5, num2str(n), 'HorizontalAlignment','center');
%     end
%     plot(t*24.2/1000,20*Efield(t)/max(Efield(t)),'red')
%     
%     hold off;
%     set(gcf,'DefaultAxesFontSize',9)
%     set(gca,'FontSize',9)
%     ylabel('Recollision energy (eV)','HorizontalAlignment','center')
%     xlabel('Ionization/Recollision time (fs)')
%     box on
%     %axis square
%     

%     figure(2)
%     hold on;
%     for n=1:numevents
%          plot(classevents{1,n}(:,2)*24.2/1000,classevents{1,n}(:,1)*27.2,'.',classevents{2,n}(:,2)*24.2/1000,classevents{2,n}(:,1)*27.2,'.' );
%          [y,ind] = max(classevents{1,n}(:,1));
%          text(classevents{1,n}(ind,2)*24.2/1000, y*27.2+5, num2str(n), 'HorizontalAlignment','center');
%     end
%     plot(t*24.2/1000,20*Efield(t)/max(Efield(t)),'red')
%     
%     hold off;
%     set(gcf,'DefaultAxesFontSize',9)
%     set(gca,'FontSize',9)
%     ylabel('Recollision energy (eV)','HorizontalAlignment','center')
%     xlabel('Ionization time (fs)')
%     box on
% clear y ind 

pause(2) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% solve saddle point equations for the separated events, but with equidistant photon energies
% 
%     firstevent = input('First recollision event you want the quantum trajectories for? [1]: ');
%         if isempty(firstevent) || firstevent < 1 ||  firstevent > numevents
%             firstevent = 1;
%             disp('Alright, let''s say you want to start with the first one.');
%         end
%     lastevent = input(['Last recollision event you want the quantum trajectories for? [',num2str(numevents),']: ']);
%         if isempty(lastevent) || lastevent < firstevent ||  lastevent > numevents
%             lastevent = numevents;
%             disp('Alright, let''s say you want to go until the last one.');
%         end
firstevent=1;
lastevent=numevents;


if (~exist('Estp','var'))
    Estp = .5*2*pi/waveperiod;
end

Ephot = (minEphot:Estp:maxEphot)';

shiftscale=250/max(Efield(t))/tstp;

numevents=lastevent-firstevent+1;
qevents=cell(2,numevents);

disp(['finding quantum trajectories for ',num2str(numevents),' events...'])
%tic
    options=optimset('Display','off');
    parfor_progress(numevents);
    parfor nn=1:numevents;
         n=nn-1+firstevent; 
         if (max(classevents{1,n}(:,1))-min(classevents{1,n}(:,1))>10*Estp) && (max(classevents{2,n}(:,1))-min(classevents{2,n}(:,1))>10*Estp)
            [dump,Estartindx] = min(abs(Ephot - (max(classevents{1,n}(:,1))+1.3*Ip -10*Estp)));  %at this element of Ephot we'll start searching for saddle points, using the closest 
                                                %"classical solution" as initial guess. This energy is about x*Estp below the cutoff of this recollision event.   
            tishift=waveperiod/200;
            teshift=waveperiod/100;
         elseif (max(classevents{1,n}(:,1))-min(classevents{1,n}(:,1))>5*Estp) && (max(classevents{2,n}(:,1))-min(classevents{2,n}(:,1))>5*Estp)
             [dump,Estartindx] = min(abs(Ephot - (max(classevents{1,n}(:,1))+1.3*Ip -5*Estp))); 
             tishift=waveperiod/150;
             teshift=waveperiod/75;
         else
             [dump,Estartindx] = min(abs(Ephot - (max(classevents{1,n}(:,1))+1.3*Ip -0*Estp)));
             tishift=waveperiod/100;
             teshift=waveperiod/50;
         end
   % do the short trajectories                                        
        tmpshort=zeros(length(Ephot),5);
        for m=Estartindx:length(Ephot)
            f = @(x)SPeq(x,Ip,Ephot(m),I0,tau,omega,freqs,tlim);  %create handle 'f' to anonymous function that executes the saddle-point eq. function with
                                          %the parameters Ip and Ephot(m
            if m==Estartindx
                [dump,k] = min(abs(classevents{1,n}(:,1)+1.3*Ip -Ephot(m))); %find closest class. traj for the photon energy tmpEphot
                guessti=classevents{1,n}(k,2) +tishift;
                %guessti = interp1(classevents{1,n}(:,1)+1.3*Ip, classevents{1,n}(:,2), Ephot(m)) +tishift;   
                    tmp=mean(diff(abs(Efield(guessti-10*tstp:tstp:guessti+10*tstp))));%if this is positive, the Efield-max is at larger t,
                                                                                 %and the ioniz.time-guess should be increased, or vice versa
                    guessti=guessti+ tmp*shiftscale; %this shifts the guess by a reasonable amount into the direction of the field peak
                guesste=classevents{1,n}(k,3)-teshift;   
                %guesste=interp1(classevents{1,n}(:,1)+1.3*Ip, classevents{1,n}(:,3), Ephot(m))-teshift;
                x0=[guessti sqrt(2*Ip)/abs(Efield(guessti))  guesste  0];
            else
                x0=tmpshort(m-1,1:4);
            end
                tmpshort(m,1:4) = fsolve(f,x0,options);   %solve the saddle point eqs. (we get them with the right parameters via handle 'f')
                tmpshort(m,5) = tmpshort(m,3)-tmpshort(m,1); %excursion time, while we're at it...
        end
        %
        for m=Estartindx-1:-1:1
            f = @(x)SPeq(x,Ip,Ephot(m),I0,tau,omega,freqs,tlim); 
            x0=tmpshort(m+1,1:4);
            tmpshort(m,1:4)=fsolve(f,x0,options);
            if tmpshort(m,1)<tmpshort(m+1,1) && tmpshort(m,3)>tmpshort(m+1,3) && min(classevents{1,n}(:,1)+1.3*Ip+.2) < Ephot(m)
                x0=[tmpshort(m+1,1)+2*tishift tmpshort(m+1,2) tmpshort(m+1,3)-5*teshift tmpshort(m+1,4)];
                tmpshort(m,1:4)=fsolve(f,x0,options);
            end
            tmpshort(m,5)=tmpshort(m,3)-tmpshort(m,1);
        end
        
   % do the long trajectories
        tmplong=zeros(length(Ephot),5);
        for m=Estartindx:length(Ephot)
            f = @(x)SPeq(x,Ip,Ephot(m),I0,tau,omega,freqs,tlim);  %create handle 'f' to anonymous function that executes the saddle-point eq. function with
                                          %the parameters Ip and Ephot(m
            if m==Estartindx
                [dump,k] = min(abs(classevents{2,n}(:,1)+1.3*Ip -Ephot(m))); %find closest class. traj for the photon energy tmpEphot
                guessti=classevents{2,n}(k,2)-tishift;
                %guessti = interp1(classevents{2,n}(:,1)+1.3*Ip, classevents{2,n}(:,2), Ephot(m)) -tishift;
                    tmp=mean(diff(abs(Efield(guessti-10*tstp:tstp:guessti+10*tstp))));
                    guessti=guessti+ tmp*shiftscale; 
                guesste=classevents{2,n}(k,3)+teshift;   
                %guesste=interp1(classevents{2,n}(:,1)+1.3*Ip, classevents{2,n}(:,3), Ephot(m)) +teshift;
                x0=[guessti sqrt(2*Ip)/abs(Efield(guessti))  guesste  0];
            else
                x0=tmplong(m-1,1:4);
            end
                tmplong(m,1:4)=fsolve(f,x0,options);   %solve the saddle point eqs. (we get them with the right parameters via handle 'f')
                tmplong(m,5)=tmplong(m,3)-tmplong(m,1); %excursion time, while we're at it...
        end
        %
        for m=Estartindx-1:-1:1
            f = @(x)SPeq(x,Ip,Ephot(m),I0,tau,omega,freqs,tlim);     
            x0=tmplong(m+1,1:4);
            tmplong(m,1:4)=fsolve(f,x0,options);
            if tmplong(m,1)>tmplong(m+1,1) && tmplong(m,3)<tmplong(m+1,3) && min(classevents{2,n}(:,1)+1.3*Ip+.2) < Ephot(m)
                x0=[tmplong(m+1,1)-2*tishift tmplong(m+1,2) tmplong(m+1,3)+5*teshift tmplong(m+1,4)];
                tmplong(m,1:4)=fsolve(f,x0,options);
            end
            tmplong(m,5)=tmplong(m,3)-tmplong(m,1);
        end
     qevents(:,nn)={tmpshort; tmplong};
   parfor_progress;
   end
 parfor_progress(0);
 %toc

 clear tmpEphot guesste guessti h k m options tmpp tstp x0 firstevent lastevent;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate stationary momentum and (not yet Legendre-transformed!) action
%
% I have checked, all my imaginary parts are fine, they are like other
% peoples' (in particular Olga's). And my phases (and thus the momenta) are
% fine, too. 

%tic
    ps=cell(size(qevents));
    acts=cell(size(qevents));
    for n=1:numevents
        %short
        ts=qevents{1,n};
        for m=1:length(ts)   %you have to do a loop over the energies because in pstat you call psalpha,
                             %which only takes scalar arguments (integration limits)
            [Rps,Ips] = pstat(ts(m,:),I0,tau,omega,freqs,tlim);
            ps{1,n}(m,:) = [Rps,Ips];
            [Ract,Iact] = action(ps{1,n}(m,:),ts(m,:),Ip,I0,tau,omega,freqs,tlim,Ephot(m));
            acts{1,n}(m,:) = [Ract,Iact];
        end
        %long
        ts=qevents{2,n};
        for m=1:length(ts)
            [Rps,Ips] = pstat(ts(m,:),I0,tau,omega,freqs,tlim);
            ps{2,n}(m,:) = [Rps,Ips];
            [Ract,Iact] = action(ps{2,n}(m,:),ts(m,:),Ip,I0,tau,omega,freqs,tlim,Ephot(m));
            acts{2,n}(m,:) = [Ract,Iact];
        end
    end
 
clear Rps Ips Ract Iact
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate everything to get higher energy resolution (makes sense if
%  you have many events (e.g. a full pulse) and narrow harmonics

% the very practical 'interp' funtion tends to produce a "riply" line after
% the cutoff though... dunno why. Normal interpolation is fine though.  

if fullpulse==1;
    interpfactor = ceil(Estp * max(tau)) ; 
else
    interpfactor = ceil(Estp * (t(end)-t(1)));
end


if interpfactor>1 
    Ephotorig =  Ephot;
    Ephot = (minEphot:Estp/interpfactor:maxEphot)';
    
    iqevents1 = zeros(length(Ephot),5);
    iqevents2 = zeros(length(Ephot),5);
    for n=1:numevents
        for m=1:5
            iqevents1(:,m) = interp1(Ephotorig, qevents{1,n}(:,m), Ephot,'spline');
            iqevents2(:,m) = interp1(Ephotorig, qevents{2,n}(:,m), Ephot,'spline');
        end
        qevents{1,n} = iqevents1;
        qevents{2,n} = iqevents2;
        clear iqevents1 iqevents2
        for m=1:2 %short and long
            acts{m,n} = [interp1(Ephotorig,acts{m,n}(:,1),Ephot,'spline'), interp1(Ephotorig, acts{m,n}(:,2), Ephot,'spline')];  %[Re ,Im]
            ps{m,n}   = [interp1(Ephotorig, ps{m,n}(:,1), Ephot,'spline'), interp1(Ephotorig, ps{m,n}(:,2),Ephot,'spline')];  %[Re ,Im]
        end
    end
    %Ephot = interp(Ephot,interpfactor);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate "emission amplitude", i.e. the whole term under the sum for the dipole amplitude
  
    emampl_raw=cell(size(qevents));
    emampl=cell(size(qevents));
    
    for n=1:numevents;
            emampl_raw{1,n} = emissionampl(1:length(Ephot),Ephot,qevents{1,n},ps{1,n},acts{1,n}).';
            emampl_raw{2,n} = emissionampl(1:length(Ephot),Ephot,qevents{2,n},ps{2,n},acts{2,n}).';
            %the dot with the transpose is super-important!
            %without it, matlab would do a complex-conjugate-transpose!!
            [emampl{1,n},emampl{2,n}] = f_killdivergence(emampl_raw(:,n),qevents(:,n),Ephot);
    end

%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% put results together
    for n=1:numevents
      	results(:,n,nCEP) = {[qevents{1,n}(:,:), emampl{1,n}]; [qevents{2,n}(:,:),emampl{2,n}]; classevents{1,n}(:,:); classevents{2,n}(:,:)};
    end

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% plot results to check if everyhting goes well...
    close all;
    
    plotindcs = 1:interpfactor:length(Ephot); %plot dots only for these energy-indices, otherwise the figure gets �ber-heavy...
    scrsz = get(0,'ScreenSize');
    figure('OuterPosition',[100 scrsz(4)/2 2*scrsz(3)/3 scrsz(4)/2],'Renderer','painters')
         subplot('Position',[0.08 0.45 0.87 0.52])
         hold on;
          for event=1:numevents
           if (~isempty(results{1,event,nCEP})) && (~max(isnan(results{1,event,nCEP}(:,6)))==1)
               scatter(results{1,event,nCEP}(plotindcs,3)*24.2/1000, Ephot(plotindcs)*27.2, 10, log10(Ephot(plotindcs).^4.*real(results{1,event,nCEP}(plotindcs,6).*conj(results{1,event,nCEP}(plotindcs,6)))),'filled');
               scatter(results{2,event,nCEP}(plotindcs,3)*24.2/1000, Ephot(plotindcs)*27.2, 10, log10(Ephot(plotindcs).^4.*real(results{2,event,nCEP}(plotindcs,6).*conj(results{2,event,nCEP}(plotindcs,6)))),'filled');
               scatter(results{1,event,nCEP}(plotindcs,1)*24.2/1000, Ephot(plotindcs)*27.2, 10, log10(Ephot(plotindcs).^4.*real(results{1,event,nCEP}(plotindcs,6).*conj(results{1,event,nCEP}(plotindcs,6)))),'filled');
               scatter(results{2,event,nCEP}(plotindcs,1)*24.2/1000, Ephot(plotindcs)*27.2, 10, log10(Ephot(plotindcs).^4.*real(results{2,event,nCEP}(plotindcs,6).*conj(results{2,event,nCEP}(plotindcs,6)))),'filled');
               p=plot(class(:,2)*24.2/1000,(class(:,1)+Ip)*27.2,'.',class(:,3)*24.2/1000,(class(:,1)+Ip)*27.2,'.');
               set(p,'Color','black','MarkerSize',4)
              [y,ind] = max(results{3,event,nCEP}(:,1)+1.3*Ip);
              text(results{3,event,nCEP}(ind,2)*24.2/1000-.1, y*27.2+5, num2str(event), 'HorizontalAlignment','center');
           end
          end
              box on
              axis tight
              load('MyGermanColormap','mycmap')
                  set(gcf,'Colormap',mycmap)
              caxis([logcmax-5 logcmax])
              xlim([min(t) max(t)]*24.2/1000)
              ylim([Ip maxEphot]*27.2)
            ylabel('Photon energy (eV)','HorizontalAlignment','center')
            set(gca,'XTickLabel','')
            text(0.95, 0.9, ['CEP = ',num2str(CEP(nCEP),'%02.2f'),'\pi'], 'Units','normalized','HorizontalAlignment','right');

         subplot('Position',[0.08 0.1 0.87 .34])
              plot(t*24.2/1000,Efield(t),'red');
              line([min(t*24.2/1000) max(t*24.2/1000)],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
              xlabel('Time (fs)')
              ylabel('E-field (a.u.)','HorizontalAlignment','center')
              box on
              xlim([min(t) max(t)]*24.2/1000)
              ylim([-0.09 0.09])
    pause(2)


end

if (~exist(savepath, 'dir'))
       mkdir(savepath);
end
save([savepath,savename], 'results', 'Ephot', 'interpfactor', 'Ip', 'I0', 'ffreqs', 'CEP', 'phi2', 't_0', 'lambda', 'omega', 't', 'tau', 'tlim', 'savepath', 'savename','-v7.3');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is just for plotting fields as funtion of time
%      
%     EE=Efield(t);
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % vector with vector potential (A)
% 
%         AA=vecpot(t);
% 
%         %if you want to integrate E numerically to get A
%          AAA=zeros(1,length(t));
%          for tt=2:length(t)
%           AAA(tt)=AAA(tt-1)+trapz(t(tt-1:tt),EE(tt-1:tt));
%          end
% 
%       %  tic
%       %  t=(-tlim:tstp:tlim); 
%       %  vp=zeros(1,length(t));
%       %  for n=2:1:length(t)
%       %  vp(n) = vecpottrapz(t(n));
%       %  end
%       %  toc
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % vector with integrated vector potential (alpha)
% 
%         tic
%         aalpha=zeros(1,length(t));
%         for k=2:1:length(t)
%             aalpha(k) = psalpha(-tlim,t(k));
%         end
%         toc
% 
%     figure
%     plot(t,EE/max(EE),t,AA/max(AA),t,aalpha/max(aalpha))
