%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% again, integrating only the carrier. actually, this only seems to work
% for a monochromatic carrier... interesting.

function salpha = simplalpha(t) 
global I0 tau omega freqs

[m,n]=size(freqs);
salpha = sum((sqrt(freqs(:,2)*I0)*ones(1,length(t))./((freqs(:,1)*omega^2)*ones(1,length(t)))) .*(ones(m,1)*cos(t/tau*2*acos(2^(-0.25))).^2) .*(cos(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))),1);





