function [Rps,Ips] = pstat(ts,I0,tau,omega,freqs,tlim) 

psalpha = @(ta,tb)fpsalpha(ta,tb,I0,tau,omega,freqs,tlim); 

Rti=ts(1); Iti=ts(2); Rtr=ts(3); Itr=ts(4);
Rtau=Rtr-Rti; Itau=Itr-Iti;

alphari=psalpha(Rti+1i*Iti,Rtr+1i*Itr);

Rps = -1/(Rtau^2 +Itau^2) * (Rtau*real(alphari) + Itau*imag(alphari)); 
Ips = 1/(Rtau^2 +Itau^2) * (Itau*real(alphari) - Rtau*imag(alphari));
