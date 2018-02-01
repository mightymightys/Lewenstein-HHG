clear all
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% waveform composition configuration
% matrix with 3 columns: frequency (multiple of omega0), rel. power, phase in rad)
    %ffreqs=@(phi,phi2)[0.5, .3, 3.0; 1, .6, 0; 2, .048, -1.9; 3, .031, 2.2; 4, .017, 0.1]; %luke's perfect wave
    %ffreqs=@(phi,phi2)[2/3, .45, -2.5; 1, .45, 0; 2, .1, 2.8];   % envisaged experiment, optimized by hand
    %ffreqs=@(phi,phi2)[2/3, .67, 0.88+pi; 1, .26, 0.05+pi; 2, .07, -0.2+pi]; % envisaged exp with Luke's optimization
    %ffreqs=@(phi,phi2)[1/3, .2, 1.3; 2/3, .35, -2.5; 1, .35, 0; 2, .1, 2.8];  %by-hand-opti incl. 3micron
    %ffreqs=@(phi,phi2)[.625, 1/3, pi/2; 1, 1/3, pi/2; 1.625, 1/3, pi/2];   % example for single atto-pulse in H.-C. Banulet et al., PRA 2010
%    ffreqs=@(phi,phi2)[1, 1, phi*pi]; % monochromatic standard case
    ffreqs=@(phi,phi2)[2/3, .358, phi*pi; 1, .597, phi*pi; 2, .045,phi*2*pi+phi2*pi];  % April2012 exp. setup
    %ffreqs=@(phi,phi2)[2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi; 2, 4.4e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup
    %ffreqs=@(phi,phi2)[1/3,  3e13, pi+phi2*pi; 2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi];  % realistic case with 3microns
    %ffreqs=@(phi,phi2)[2/3, .5, phi*pi; 1, .45, phi*pi; 2, .05, phi*2*pi+phi2*pi]; % 3-color-modif1
    %ffreqs=@(phi,phi2)[2/3, .45, phi*pi; 1, .45, phi*pi; 2, .1, phi*2*pi+phi2*pi]; % envisaged experiment, optimized by hand, "modif2"
    %ffreqs=@(phi,phi2)[2/3, .5, phi*pi; 1, .5, phi*pi];  % 2-color-1+1,5


%% 1
%clear all
%close all;
%ffreqs=@(phi,phi2)[2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi; 2, 4.4e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup
clearvars -except ffreqs;

    minEphot = 30 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
    maxEphot = 310 / 27.2; %largest photon energy
    Estp=0.2;
    
    firstevent=1;
    lastevent=2;
    
    close all;
    %choice = 1;

    dip=1;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p

    longestexcursion = 0.5; %longgest excursion time to consider (in the first classical part)
                          %in wave-periods
    %longestexcursion = 1;
    fullpulse = 0;      %simulate whole pulse, or only the central period
    spectraonly = 1;    %plot only the spectra in the end (=1), or also the times...
    shortonly =1;
    savefigs=0;
    
    basewavelength = 1030; %fundamental wavelength in nm
    pulseduration = 200; %intensity-envelope FWHM of all colors, in fs
    intensity = (0.55:0.0125:4)*1e14;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
                                 % so this is just a multiplier that could take into account that our gas medium is shifted from the
                                 % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.


    CEP=1.8;
    phi2=1.2;

    %Ip = 15.76 /27.2;  
    Ip = 21.56 /27.2;  

    savename = ['alphatest.mat'];
    savepath = ['/Users/stefan/Dropbox/lewenstein-hhg/results/alpha/'];

    
         logcmax=-5.5;
         logcmin=logcmax-4;
         lincmax=3;
         lincmin=0;  

    saddlepointSFAforalpha;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%results{traj,event,nI0}(Ephot,6))
close all
clear A

Eind=24;

for i=1:length(intensity)
    A(i)=results{2,1,i}(Eind,6);
end
%for i=29:length(intensity)
%    A(i)=results{1,2,i}(Eind,6);
%end

Aphase=unwrap(angle(A));

Aampl=A.*conj(A);

first=11;

p=polyfit(intensity(first:end),Aphase(first:end),1);
R=corrcoef(intensity(first:end),Aphase(first:end));

disp(['alpha = ',num2str(-p(1)*1e14),'e-14 cm^2/W, with R = ',num2str( R(1,2))])


figure
subplot(211)
    plot(intensity, Aphase,'--rs',intensity, p(1)*intensity+p(2))
    ylabel('Dipole phase (rad)','HorizontalAlignment','center')
    xlabel('Driver intensity (W/cm^2)','HorizontalAlignment','center')
    text(0.94, 0.85, ['E_{phot} = ',num2str(Ephot(Eind)*27.2),' eV'], 'Units','normalized','HorizontalAlignment','right');
    text(0.95, 0.70, ['\alpha = ',num2str(-p(1)*1e14),'e-14 cm^2/W'], 'Units','normalized','HorizontalAlignment','right');
subplot(212)
    plot(intensity, log10(Aampl),'--rs')
    ylabel('Dipole intensity (arb.u.)','HorizontalAlignment','center')
    xlabel('Driver intensity (W/cm^2)','HorizontalAlignment','center')



%%%%%%%%%%%%%%%%
%% 
for Eind=1:41% length(Ephot) %19

    for i=1:length(intensity)
        Ashort(i)=results{1,1,i}(Eind,6);
        Along(i)=results{2,1,i}(Eind,6);
    end
%for i=29:length(intensity)
%    A(i)=results{1,2,i}(Eind,6);
%end

Ashortphase=unwrap(angle(Ashort));
Alongphase=unwrap(angle(Along));
Ashortampl=Ashort.*conj(Ashort);
Alongampl=Along.*conj(Along);

first=11;

p=polyfit(intensity(first:end),Ashortphase(first:end),1);
alphashort(Eind)=-p(1)*1e14;

p=polyfit(intensity(first:end),Alongphase(first:end),1);
alphalong(Eind)=-p(1)*1e14;


%disp(['alpha = ',num2str(-p(1)*1e14),'e-14 cm^2/W, with R = ',num2str( R(1,2))])

end

figure
    plot(Ephot(1:Eind)*27.2, alphashort,'--rs',Ephot(1:Eind)*27.2, alphalong,'--rs')
    ylabel('\alpha   (1e-14 cm^2/W)','HorizontalAlignment','center')
    xlabel('Photon energy (eV)','HorizontalAlignment','center')
%    text(0.94, 0.85, ['E_{phot} = ',num2str(Ephot(Eind)*27.2),' eV'], 'Units','normalized','HorizontalAlignment','right');
%    text(0.95, 0.70, ['\alpha = ',num2str(-p(1)*1e14),'e-14 cm^2/W'], 'Units','normalized','HorizontalAlignment','right');



