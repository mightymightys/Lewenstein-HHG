function EmAmpl = f_emissionamplitude(Ephotindx,Ephot,qevents,ps,acts,I0,tau,omega,freqs,Ip,tlim,dip)

Ephot=Ephot';

Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim); 
vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 

aacts = (acts(Ephotindx,1) + 1i*acts(Ephotindx,2)).';  %the dot with the transpose is super-important!
                                                       %without it, matlab would do a complex-conjugate-transpose!!
pps =   (  ps(Ephotindx,1) + 1i*  ps(Ephotindx,2)).';
ti = (qevents(Ephotindx,1) + 1i*qevents(Ephotindx,2)).';
tr = (qevents(Ephotindx,3) + 1i*qevents(Ephotindx,4)).';

kti = pps + vecpot(ti);
%kti = real(kti);
%kti = abs(kti);
ktr = pps + vecpot(tr);

if dip==1
        % flat
        di = 1i*ones(1,length(Ephotindx));
        dr = 1i*ones(1,length(Ephotindx));
    elseif dip==2
        % Gaussian
        di=1i*(pi*0.8*Ip)^(-0.75) * kti/0.8/Ip .* exp(-kti.^2/(2*0.8*Ip)); 
        dr=1i*(pi*0.8*Ip)^(-0.75) * ktr/0.8/Ip .* exp(-ktr.^2/(2*0.8*Ip));
    elseif dip==3
        % scaled hydrogren 1s
        di=1i*(2^3.5*(2*Ip)^1.25/pi) * kti./(kti.*conj(kti) +2*Ip).^3;
        dr=1i*(2^3.5*(2*Ip)^1.25/pi) * ktr./(ktr.*conj(ktr) +2*Ip).^3;
    elseif dip==4
        % scaled hydrogenic 3p
        Z = sqrt(2*Ip);
        % it's problematic to use this for the ioniz. dipole, because 1) it's huge at
        % very low k and 2) the numerical noise in the imaginary part of the kti (kti is supposed to be the same for all saddle points...) 
        % somehow translates into ridiculously huge spikes in this dipole. Dunno how this happpens.. but anyway:
        % the stuff is gone if you use abs(kti) here...
        dmp=kti;
        kti=abs(kti);
        di = -sqrt(3/(2*pi)) * (pi*2^6/(2*pi*Z)^(3/2)) * Z^5 * ((20*kti.^4 - (5*kti.^2 - Z^2 ).^2) ./(kti.^2 + Z^2).^5);
        kti=dmp;
        %di = 1i*ones(1,length(Ephotindx));
        dr = -sqrt(3/(2*pi)) * (pi*2^6/(2*pi*Z)^(3/2)) * Z^5 * ((20*ktr.^4 - (5*ktr.^2 - Z^2 ).^2) ./(ktr.^2 + Z^2).^5);
    else
        % flat
        di = 1i*ones(1,length(Ephotindx));
        dr = 1i*ones(1,length(Ephotindx));
end
    
detS= (ktr.*kti./(tr-ti)).^2 -( -2*(Ephot(Ephotindx)-Ip*ones(1,length(Ephotindx)))./(tr-ti) + Efield(tr).*ktr) .*( 2*Ip./(tr-ti) - Efield(ti).*kti);

EmAmpl = 1i*2*pi./sqrt(detS)...
        .* (pi./(1e-9*ones(1,length(Ephotindx))+ 1i*(tr-ti)/2)).^(1.5)...
        .* Efield(ti) .* di...
        .* conj(dr)...
        .* exp(1i*aacts + 1i*Ephot(Ephotindx).*tr);