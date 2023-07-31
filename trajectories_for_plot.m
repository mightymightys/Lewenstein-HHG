clear all;
%close all;
pause on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation parameters

%ffreqs=@(phi,phi2)[1/3,  .2, pi+phi2*pi; 2/3, .4, phi*pi; 1, .4, phi*pi];  % realistic case with 3microns
%ffreqs=@(phi,phi2)[1, 1.2e14, pi]; % monochromatic standard case
%ffreqs=@(phi,phi2)[0.5, .3, 3.0; 1, .6, 0; 2, .048, -1.9; 3, .031, 2.2; 4, .017, 0.1]; %luke's perfect wave
%ffreqs=@(phi,phi2)[2/3, .64, 0.85*pi; 1, .28, 0.85*pi; 2, .08, 0.12*pi]; % Luke's optimization, w/o considering free electrons
%ffreqs=@(phi,phi2)[2/3, .5, 0.6*pi; 1, .44, 0.6*pi; 2, .06, 1.7*pi]; %Luke's optimization for 72.5eV, <4% free electrons (at 0.75e14 intensity)
ffreqs=@(phi,phi2)[2/3, 6e13, phi*pi; 1, 5.7e13, phi*pi; 2, 3e12, phi*2*pi+phi2*pi];  % November2012 recent exp. setup - with intensities as estimated for the paper
%ffreqs=@(phi,phi2)[2/3, 4e13, phi*pi; 1, 7.7e13, phi*pi; 2, 3e12, phi*2*pi+phi2*pi];
%ffreqs=@(phi,phi2)[2/3, 1/2, phi*pi; 1, 1/2, phi*pi];  % 2-color 1+1,5
%ffreqs=@(phi,phi2)[2/3, .358, phi*pi; 1, .597, phi*pi; 2, .045,phi*2*pi+phi2*pi];  % April2012 exp. setup


choice = 2;

minEphot=40 /27.2; 

longestexcursion = 1; %longgest excursion time to consider (in the first classical part) in wave-periods
Nperiods = 1; %over how many waveperiods do you want to calculate (this is where the ionization times are initially scanned over)
    
basewavelength =1030; %fundamental wavelength in nm
pulseduration = 200; %intensity-envelope FWHM of all colors, in fs; negative value makes the wave CW

intensity = 1;
%intensity = 1.2e14;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
                                 % so this is just a multiplier that could take into account that our gas medium is shifted from the
                                 % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.

CEP=0.85;
phi2=0.2;
    
Ip = 15.76 /27.2;  

nuclfactor =0;


I0=intensity / 3.5094452E16; % total intensity in a.u.

lambda=basewavelength*1e-9; %base-wavelength corresponding to base-freq. omega

tau=pulseduration *1000/24.2;  %FWHM duration of the laser pulse, in a.u.

tlim=pi/2 * tau /2 /acos(2^(-0.25));  % limit of the time-window, where the cos^2-enevelope is zero

omega=2*pi*299792458/lambda *24.2e-18; %laser freq. in a.u.

%disp('Results will be saved in:')
%disp(savepath)
%disp(savename)

freqs=ffreqs(CEP,phi2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fix up the functions for the set parameters (so as to avoid stupid global variables)
classtraj = @(ti,tt,aalpha)fclasstraj(ti,tt,aalpha,freqs,omega,I0,tau,tlim); 
Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim); 
vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 
psalpha = @(ta,tb)fpsalpha(ta,tb,I0,tau,omega,freqs,tlim); 
emissionampl = @(Ephotindx,Ephot,qevents,ps,acts)f_emissionamplitude(Ephotindx,Ephot,qevents,ps,acts,I0,tau,omega,freqs,Ip,tlim,dip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% classical trajectories

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set up time vectors - decide which part of the pulse you want to consider

    tstp=.25; % timestep in a.u. with which to scan the birth times and calc. trajectories
    
    [n,m]=size(freqs);
    if n>1
        freqssorted = sort(freqs(:,1));
        waveperiod = 1/(freqssorted(2)-freqssorted(1)) * 2*pi/omega;
    else waveperiod = 2*pi/omega/freqs(1,1);
    end
    %lowestperiod = 1*2*pi/omega/min(freqs(:,1)); %x*one period of the lowest freq.component
    maxduration = longestexcursion*waveperiod;
    t=(-0.5*Nperiods*waveperiod:tstp:0.5*Nperiods*waveperiod+maxduration);
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % make vector with integrated vector potential (alpha) (to be used by classtraj)
        alpha=zeros(1,length(t));
         for n=2:length(t)
             %alpha(n) = psalpha(-tlim,t(n));
             %alpha(n) = alpha(n-1)+ psalpha(t(n-1),t(n));  %this really gives the exact same thing
             alpha(n) = alpha(n-1) + trapz(fvecpot(t(n-1:n),I0,tau,omega,freqs,tlim))*tstp;
                                                            %here, when you advance in very small time steps,
                                                            %trapz is totally good enough and muuuch faster. Even just summing up would do...
         end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find recolliding trajectories
 tic   
    ti=(min(t):tstp:max(t)-maxduration); %birth times to be scanned 
    
    % can't really initialize the results-vector 'class'... dunno how long
    % it'll be. But, it's important to at least set it to a single line of
    % four zeros, so you don't drag along recollisions from the former
    % iterations!
    trajs=zeros(length(t),length(ti)+1);
    trajs(:,1)=t;
    Ekin=zeros(length(t),length(ti)+1);
    Ekin(:,1)=t;
        
    for n=1:length(ti) 
        shiftstart=0;
        aalpha=alpha(n+shiftstart:end)-alpha(n); 
        tt=t(n+shiftstart:end);
        trajs(n+shiftstart:end,1+n) = classtraj(ti(n),tt,aalpha);
        Ekin(n+shiftstart:end,1+n) = .5*(vecpot(tt)-vecpot(ti(n))).^2; 
    end
toc    
    
%% find recollisions    
tic
m=1;
     for n=1:length(ti) 
         % here, we will for every birth time in ti calculate the class.
         % trajectory for [ti,ti+1period of the lowest freq. comp.], and
         % look for zero crossings
         [I,dmp] = crossings(trajs(n+20:end,n+1),[],0);  %find zero-crossings of classical trajectory for birth time ti(n)
         I=I+n+20;
         if (~isempty(I)) && (.5*(vecpot(t(min(I)))-vecpot(ti(n)))^2)>(minEphot-Ip) && (t(min(I))-ti(n))<longestexcursion*waveperiod
                                   % if the the output of the zero-crossing search is NOT empty,
                                   % and the recol. energy is >x eV. (at low energies, all this is inaccurate anyway)
%          class(m,3)= min(Ti);     % write down the time value when the zero-crossing (=recollision) occurs
%          class(m,2)= ti(n);       % also write down the associated ionization time
%          class(m,1)= .5*(vecpot(min(Ti))-vecpot(ti(n)))^2; % also write down the recollision energy =0.5*(A(tr)-A(ti))^2
%          class(m,4)= 1E15*(exp(-2*(2*Ip)^(3/2)./abs(3*Efield(ti(n)))).*((class(m,3)-ti(n)).^(-5))); %write down a "trajectory weight"
         recollindsI(m)=n;
         recollindsR(m)=min(I);
         m=m+1;                   % increase index of the "results"-vector by one for next zero crossing
         end 
     end
%%
for i=1:length(recollindsI)     
     recollenergies(i) = Ekin(recollindsR(i),recollindsI(i)+1);
end
[cutoff,cutoffind] = max(recollenergies);

toc


%% plot 
%trajweight=exp(-2*(2*Ip)^(3/2)./abs(3*Efield(t(recollindsI)))).*(t(recollindsR)-t(recollindsI)).^(-5); 
trajweight=exp(-2*(2*Ip)^(3/2)./abs(3*Efield(t(recollindsI)))); 

stretchtx=1.5;
stretchtt=12;

figure
    scatter3(stretchtt*((-2:.01:8.3)+t(recollindsI(cutoffind))*24.2/1000),zeros(1,length((-2:.01:8.3))),(Efield((-2:.01:8.3)*1000/24.2+t(recollindsI(cutoffind)))+.2) *300, 4,[0 .6 0],'filled')
    hold on
    %plot3(t*24.2/1000,(Efield(t)/max(Efield(t)))*50,zeros(1,length(t)),'Color','red','LineWidth',3)
for i=cutoffind
    %for i=1:5:length(recollindsI)
    %A=(Ekin(recollindsR(i),recollindsI(i)+1)-(minEphot-Ip))/(cutoff-(minEphot-Ip));
    A=(t(recollindsI(i):recollindsR(i))-t(recollindsI(i))).^2.2'/5.5e4;
    %colors= [A A (1-A)*A];
    %colors= [A (1-A).*A (1-A).*A];
    %colors = [ones(size(A)) A A];
    colors = [1 0 0];
    %colors =(1-trajweight(i)/max(trajweight)) * [1 1 1];
    %linewidths = 5* trajweight(i)/max(trajweight);
    %linewidths = 2+ .0004*(t(recollindsI(i):recollindsR(i))-t(recollindsI(i))).^3;
    linewidths = 2+ 2.5*(t(recollindsI(i):recollindsR(i))-t(recollindsI(i)));
    %linewidths = 2;
    %plot3(t(recollindsI(i):recollindsR(i))*24.2/1000,-abs(trajs(recollinds
    %I(i):recollindsR(i),recollindsI(i)))*0.53,Ekin(recollindsI(i):recollindsR(i),recollindsI(i))*27.2,'Color', colors,'LineWidth',linewidths)
    %scatter3(t(recollindsI(i):recollindsR(i))*24.2/1000,-abs(trajs(recollindsI(i):recollindsR(i),recollindsI(i)))*0.53,Ekin(recollindsI(i):recollindsR(i),recollindsI(i))*27.2,linewidths, colors, 'filled')
    %bubbleplot3(stretchtt*t(recollindsI(i):recollindsR(i))*24.2/1000,stretchtx*(trajs(recollindsI(i):recollindsR(i),recollindsI(i)))*0.53,Ekin(recollindsI(i):recollindsR(i),recollindsI(i))*27.2,.02*linewidths, colors, .05)
    for j=recollindsI(i):recollindsR(i)
        bubbleplot3(stretchtt*t(j)*24.2/1000,stretchtx*(trajs(j,recollindsI(i)))*0.53,Ekin(j,recollindsI(i))*27.2,.02*linewidths(j-recollindsI(i)+1), colors, 0.05*(1-A(j-recollindsI(i)+1)))
    end
    %scatter3(t(recollindsI(i):recollindsR(i))*24.2/1000,zeros(size(trajs(recollindsI(i):recollindsR(i),recollindsI(i)))),Ekin(recollindsI(i):recollindsR(i),recollindsI(i))*27.2, 4, .4*[1 1 1], 'filled')
    scatter3(stretchtt*t(recollindsI(i):recollindsR(i))*24.2/1000,-stretchtx*abs(trajs(recollindsI(i):recollindsR(i),recollindsI(i)))*0.53,zeros(size(Ekin(recollindsI(i):recollindsR(i),recollindsI(i)))), 6, .5*[1 1 1],'filled')
    scatter3(stretchtt*(8.3+t(recollindsI(i))*24.2/1000)*ones(size(t(recollindsI(i):recollindsR(i)))),-stretchtx*abs(trajs(recollindsI(i):recollindsR(i),recollindsI(i)))*0.53,Ekin(recollindsI(i):recollindsR(i),recollindsI(i))*27.2, 6, .5*[1 1 1],'filled')
    plot3(stretchtt*t(recollindsI(i))*24.2/1000*ones(size(0:cutoff*27.2)),zeros(size(0:cutoff*27.2)),(0:cutoff*27.2),  '--','Color', .5*[1 1 1], 'LineWidth', 1);
    %plot3(stretchtt*t(recollindsI(i))*24.2/1000*ones(size(0:80)),zeros(size(0:80)),(0:80),  '--','Color', .5*[1 1 1], 'LineWidth', 1);
    plot3(stretchtt*t(recollindsR(i))*24.2/1000*ones(size(0:cutoff*27.2)),zeros(size(0:cutoff*27.2)),(0:cutoff*27.2),  '--','Color', .5*[1 1 1], 'LineWidth', 1);
    plot3(stretchtt*(t(recollindsR(i))*24.2/1000:.1:8.3+t(recollindsI(cutoffind))*24.2/1000),zeros(size((t(recollindsR(i))*24.2/1000):.1:8.3+t(recollindsI(cutoffind))*24.2/1000)),cutoff*27.2*ones(size((t(recollindsR(i))*24.2/1000):.1:8.3+t(recollindsI(cutoffind))*24.2/1000)),  '--','Color', .5*[1 1 1], 'LineWidth', 1);

    set(gca,'FontSize',14);
    xlabel('Time (fs)')
    zlabel('Ekin (eV)')
    ylabel('x (nm)')
    axis tight
    zlim([0,120])
    ylim(stretchtx*[-50,0])
    xlim(stretchtt*([-2,8.3]+t(recollindsI(i))*24.2/1000))
    set(gca,'YTick',stretchtx*(-50:25:0))
    set(gca,'YTickLabel',{'-5','-2.5','0'})
    set(gca,'XTick',stretchtt*((-2:2:8)+t(recollindsI(i))*24.2/1000))
    set(gca,'XTickLabel',{'-2','0','2','4','6','8'})
    
    set(gcf,'Renderer','opengl',...
              'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[9.5 6.5],'PaperPosition',[0.25 0.25 9 6]);
          set(gcf,'DefaultAxesFontSize',14);

    view(-40,17)
end
