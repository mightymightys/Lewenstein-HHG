%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the classical trajectory for any birth time ti (in a.u.), for
% one period of the lowest frequency field

function traj = fclasstraj(ti,tt,aalpha,freqs,omega,I0,tau,tlim) 

%stp=.5;
%tt=(ti+stp:stp:ti+2*pi/omega/min(freqs(:,1)));

%alpha=zeros(1,length(tt));
%alpha(1)=psalpha(ti,ti+stp);
%for n=2:length(tt)
%    alpha(n)=alpha(n-1)+psalpha(tt(n-1),tt(n));
%end
       
vecpot = @(t)fvecpot(t,I0,tau,omega,freqs,tlim); 

traj = aalpha - (tt-ti)*vecpot(ti);