% x is the vector of (Re(ti),Im(ti),Re(tr),Im(tr))
% fsolve finds the vector that solves y=F(x)=0

function y = SPeq(x,Ip,Ephot,I0,tau,omega,freqs,tlim) 

vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 
psalpha = @(ta,tb)fpsalpha(ta,tb,I0,tau,omega,freqs,tlim); 

Rti=x(1); Iti=x(2); Rtr=x(3); Itr=x(4);
Rtau=Rtr-Rti; Itau=Itr-Iti;

alphari=psalpha(Rti+1i*Iti,Rtr+1i*Itr);

Ai=vecpot(Rti+1i*Iti);
Ar=vecpot(Rtr+1i*Itr);

Rps = -1/(Rtau^2 +Itau^2) * (Rtau*real(alphari) + Itau*imag(alphari)); 
Ips = 1/(Rtau^2 +Itau^2) * (Itau*real(alphari) - Rtau*imag(alphari));



y(1)  = (Rps + real(Ai))^2 - (Ips +imag(Ai))^2 + 2*Ip;
y(2)  = (Rps + real(Ai))*(Ips +imag(Ai));
%y(1)  = - (Ips +imag(Ai))^2 + 2*Ip;
%y(2)  = Rps + real(Ai);
y(3)  = (Rps + real(Ar))^2 - (Ips +imag(Ar))^2 + 2*Ip -2*Ephot;
y(4)  = (Rps + real(Ar))*(Ips +imag(Ar));