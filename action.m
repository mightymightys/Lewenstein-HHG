function [Ract,Iact] = action(ps,ts,Ip,I0,tau,omega,freqs,tlim,~) 

p = ps(1) + 1i*ps(2);
ti = ts(1) + 1i*ts(2);
tr = ts(3) + 1i*ts(4);

if (isfinite(ti)) && (isfinite(tr))
    act = - quadgk(@(t)0.5*(p + fvecpot(t,I0,tau,omega,freqs,tlim)).^2+Ip,ti,tr);% + Ephot*tr;
else act=1i*1e6;
end

Ract=real(act);
Iact=imag(act);