%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the integrated vector potential \int_{-infty}^t A(t')dt' =
% \int_{-infty}^t \int_{-infty}^t' E(t'')dt'' dt'
% in a.u. for any (cplx valued) time t (in a.u.)

function alpha = fpsalpha(ta,tb,I0,tau,omega,freqs,tlim)

if isfinite(ta) && isfinite(tb)
    alpha = quadgk(@(t)fvecpot(t,I0,tau,omega,freqs,tlim),ta,tb);  
    %only works if vecpot calculates A analytically (or accepts vector arguments)
else alpha=0;
end
                                    
% in case you get A by integrating E with quadgk, you can't use quad again
% here because that makes a quad(quad(..)) and quad passes vector arguments
% but only accepts scalars...   so you have to do E->trapz->A->trapz->alpha
% stp=2;
% tt=[-tlim:stp:t];
% EE=Efield(tt);
% 
% AA=zeros(1,length(tt));
% for n=2:length(tt)
%   AA(n)=AA(n-1)+trapz(tt(n-1:n),EE(n-1:n));
% end
% 
% alpha=trapz(tt,AA);
