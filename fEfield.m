%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the Efield in a.u. for any (cplx valued) time t (in a.u.)

function E = fEfield(t,I0,tau,omega,freqs,tlim) 

%E = E0*cos(t/tau*2*acos(2^(-0.25))).^2.*cos(omega*t);

[m,~]=size(freqs);
if min(tau)<0
    E=sum((sqrt(freqs(:,2)*I0)*ones(1,length(t))) .*1.*(cos(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))),1);
else
    if length(tau)==m
        envelopes = cos(1./tau *t *2*acos(2^(-0.25))).^2;
        % make sure that the colors with the shorter envelopes have zero
        % envelope outside their central cos^2-peak
        [~,tauindx] = sort(tau);
        for n=1:length(tauindx)-1
            envelopes(tauindx(n),[find( t< -tlim/tau(tauindx(end))*tau(tauindx(n))), find(t> tlim/tau(tauindx(end))*tau(tauindx(n)))] )=0;
        end 
        % sum up the fields
        E=sum((sqrt(freqs(:,2)*I0)*ones(1,length(t))) .*envelopes .*(cos(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))),1);
    else
        envelopes = (ones(m,1)*cos(t/max(tau)*2*acos(2^(-0.25))).^2);
        E=sum((sqrt(freqs(:,2)*I0)*ones(1,length(t))) .*envelopes .*(cos(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))),1);
    end
end


%just to be sure; in case t somehow reaches beyond the longest
%cos^2-envelope duration
indxs=[find(t<-tlim),find(t>tlim)];
E(indxs)=0;

% this representation (envelope-whether Gaussian or cos^2 or whatever- times carrier)
% is actually inaccurate for rather short pulses (30fs or shorter, say)
% because the spectrum of such a pulse extends down to extremely low
% frequencies and eventually acquires dc-components. these are unphysical
% because they don't propagate. Funny thing: these dc-components depend
% strongly on the CEP: they are much more pronounced for cos pulses (makes
% sense since their field is the most asymmetric)
% The more accurate representation would be to construct the pulse in the
% frequency domain and to sum up discrete frequency components (much like
% the modes in an oscillator cavity). 