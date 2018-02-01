%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the vector potential \int_{-infty}^t E(t')dt' in a.u. for any *vector of* (cplx valued) times t (in a.u.)

function AA = vecpottrapz(t,tlim tstp)

AA = -quadgk(@Efield,-tlim,t);

% AA=zeros(1,length(t));
% for n=1:length(t)
%  tt=(-tlim:tstp:t(n));
%  EE=Efield(tt);
%  AA(n)=trapz(tt,EE);
% end