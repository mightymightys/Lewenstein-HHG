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
    envelopes = 1;
else
  if length(tau)==m
      envelopes = cos((1./tau *t - (1./tau).*freqs(:,4)*ones(1,length(t))) *2*acos(2^(-0.25))).^2;
      %envelopes = cos(1./tau *t - 1./tau.*freqs(:,4)*ones(1,length(t)) *2*acos(2^(-0.25))).^2;
  else
    envelopes = cos((ones(m,1)*t-freqs(:,4)*ones(1,length(t)))./max(tau)*2*acos(2^(-0.25))).^2;
  end
  % make sure that the envelopes are zero outside their central cos^2-peak
  for n=1:m
    envelopes(n,[ find(t-freqs(n,4)< -tlim(n)), find(t-freqs(n,4)> tlim(n)) ]) = 0;
  end
end

% sum up the vector-potentials
A = -sum((sqrt(freqs(:,2)*I0)./(freqs(:,1)*omega)*ones(1,length(t))) .* envelopes...
    .*sin( freqs(:,1)*omega.*(ones(m,1)*t - freqs(:,4)*ones(1,length(t)))...
          +(freqs(:,3))*ones(1,length(t))...
          + freqs(:,5).*(ones(m,1)*t-freqs(:,4)*ones(1,length(t))).*conj(ones(m,1)*t-freqs(:,4)*ones(1,length(t)))...
          ) ,1);

% 
% %just to be sure; in case t somehow reaches beyond the longest
% %cos^2-envelope duration
% indxs=[find(t<-tlim),find(t>tlim)];
% A(indxs)=0;
