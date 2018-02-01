%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the vector potential -\int_{-infty}^t E(t')dt' in a.u. for any
% (cplx valued) time t (in a.u.).
% but this actually gives problemes with the dc-component that the short
% pulses get when representing them with envelope-times-carrier.
% for long enough pulses (>30fs or so), this problem gets small. But then, 
% the vector potential can also be well approximated by just "integrating the
% Efield carrier and neglecting the envelope", i.e. by phase-shifting the 
% by Efield-carrier by pi/2 and just dividing the same envelope by omega.
% So... since for short envelopes, stuff is inaccurate anyway, I just stick
% to the simple thing and calculate A analytically neglecting the envelope...

function A = fvecpot(t,I0,tau,omega,freqs,tlim)

%if abs(tau)< 20 *1000/24.2;
%    A = -quadgk(@Efield,-tlim,t);
%else
    %A = -E0/omega * cos(t/tau*2*acos(2^(-0.25))).^2.*sin(omega*t);
%end

[m,~]=size(freqs);
if tau<0
    A= -sum( (sqrt(freqs(:,2)*I0)./(freqs(:,1)*omega)*ones(1,length(t))) .*1 .*(sin(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))) ,1);
else
  if length(tau)==m
      envelopes = cos(1./tau *t *2*acos(2^(-0.25))).^2;
        % make sure that the colors with the shorter envelopes have zero
        % envelope outside their central cos^2-peak
        [~,tauindx] = sort(tau);
        for n=1:length(tauindx)-1
            envelopes(tauindx(n),[find(t< -tlim/tau(tauindx(end))*tau(tauindx(n))), find(t> tlim/tau(tauindx(end))*tau(tauindx(n)))] )=0;
        end 
      % sum up the vector-potentials
    A= -sum( (sqrt(freqs(:,2)*I0)./(freqs(:,1)*omega)*ones(1,length(t))) .* envelopes .*(sin(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))) ,1);
  else
    A= -sum( (sqrt(freqs(:,2)*I0)./(freqs(:,1)*omega)*ones(1,length(t))) .*(ones(m,1)*cos(t/max(tau)*2*acos(2^(-0.25))).^2) .*(sin(freqs(:,1)*omega*t+(freqs(:,3))*ones(1,length(t)))) ,1);
  end
end

%just to be sure; in case t somehow reachesbeyond the longest
%cos^2-envelope duration
indxs=[find(t<-tlim),find(t>tlim)];
A(indxs)=0;
