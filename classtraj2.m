%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the classical trajectory for any birth time ti (in a.u.), for
% one period of the lowest frequency field
% *WITHOUT* having calculated the alphas beforehand (which is more
% efficient but not practical if you want to jsut compute a few trajs)

function traj = classtraj(ti) 
global freqs omega

stp=.5;
tt=(ti+stp:stp:ti+2*pi/omega/min(freqs(:,1)));

alpha=zeros(1,length(tt));
alpha(1)=psalpha(ti,ti+stp);
for n=2:length(tt)
    alpha(n)=alpha(n-1)+psalpha(tt(n-1),tt(n));
end
       
traj = alpha - (tt-ti)*vecpot(ti);