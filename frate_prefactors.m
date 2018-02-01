function pre = frate_prefactors(Ephotindx,qevents,ps,traj,event,I0,tau,omega,freqs,Ip,tlim)

Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim); 
vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 

pps = ps{traj,event}(Ephotindx,1) + j*ps{traj,event}(Ephotindx,2);
ti = qevents{traj,event}(Ephotindx,1) + j*qevents{traj,event}(Ephotindx,2);
tr = qevents{traj,event}(Ephotindx,3) + j*qevents{traj,event}(Ephotindx,4);

detS= ((pps +vecpot(tr))*(pps +A(ti))/(tr-ti))^2 ...
    -( -2*(omega-Ip)/(tr-ti) + Efield(tr)*(pps +vecpot(tr)))...
    *( -2*Ip/(tr-ti) - Efield(ti)*(pps +vecpot(ti)));