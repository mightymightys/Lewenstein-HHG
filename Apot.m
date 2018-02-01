clear all;

global E0;
E0=1;
lambda=800e-9; %wavelength
global omega;
omega=2*pi*299792458/lambda *24.2e-18; %laser freq. in a.u.
global tau;
tau=10 *1000/24.2;  %FWHM duration of the laser pulse, in a.u.


global tstp;
tstp= 1; %timestep in a.u.
global tlim;
tlim=pi/2 * tau /2 /acos(2^(-0.25));

t=[-tlim:tstp:tlim];

E=Efield(t);

tic
A=zeros(1,length(E));
k=2;
for tt=-tlim+tstp:tstp:tlim
A(k)=A(k-1)+quad(@Efield,tt-tstp,tt);
k=k+1;
end
toc

% tic
% AA=zeros(1,length(E));
% for tt=2:length(E)
% AA(tt)=A(tt-1)+trapz(t(tt-1:tt),E(tt-1:tt));
% end
% toc

% tic
% AAA=zeros(1,length(E));
% for k=1:1:length(t)
% AAA(k) = vecpot(t(k));
% end
% toc

figure
plot(t,E,t,A)