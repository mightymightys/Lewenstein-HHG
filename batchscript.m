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
%    ffreqs=@(phi,phi2)[2/3, .358, phi*pi; 1, .597, phi*pi; 2, .045,phi*2*pi+phi2*pi];  % April2012 exp. setup
    %ffreqs=@(phi,phi2)[2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi; 2, 4.4e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup
    %ffreqs=@(phi,phi2)[1/3,  3e13, pi+phi2*pi; 2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi];  % realistic case with 3microns
    %ffreqs=@(phi,phi2)[2/3, .5, phi*pi; 1, .45, phi*pi; 2, .05, phi*2*pi+phi2*pi]; % 3-color-modif1
    %ffreqs=@(phi,phi2)[2/3, .45, phi*pi; 1, .45, phi*pi; 2, .1, phi*2*pi+phi2*pi]; % envisaged experiment, optimized by hand, "modif2"
    %ffreqs=@(phi,phi2)[2/3, .5, phi*pi; 1, .5, phi*pi];  % 2-color-1+1,5

%% test
clear all
close all;
%ffreqs=@(phi,phi2)[1/3,  3e13, pi+phi2*pi; 2/3, 12e13, phi*pi];  % signal+idler
%ffreqs=@(phi,phi2)[2/3, 12e13, phi*pi];  % signal only

clearvars -except ffreqs;

    minEphot = 50 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
    maxEphot = 350 / 27.2; %largest photon energy

    clearvars -except phi2 ffreqs minEphot maxEphot ; 
    close all;
    choice = 3;
    saveimgs =0;

    dip=3;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p

    longestexcursion = 1; %longgest excursion time to consider (in the first classical part)
                          %in wave-periods
    fullpulse = 1;      %simulate whole pulse, or only the central period
    spectraonly = 0;    %plot only the spectra in the end (=1), or also the times...
    
    basewavelength = 1030; %fundamental wavelength in nm
    pulseduration = 80; %intensity-envelope FWHM of all colors, in fs
    intensity = 1;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
                                 % so this is just a multiplier that could take into account that our gas medium is shifted from the
                                 % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.


    CEPstp=.05;
    CEP=0;
%    CEP=(0:CEPstp:2-CEPstp);
    phi2=0;

    Ip = 15.76 /27.2;   %argon
    %Ip = 21.565 /27.2;   %neon
    %Ip = 24.587 /27.2;   %helium

    savename = ['signalonly_12e13.mat'];
    savepath = ['C:\Users\haessler\Documents\matlab\lewenstein-hhg\results\2colour-shortpulses\signal+idler\signal-only\'];

        
    logcmax=-8;
    logcmin=logcmax-4;
    lincmax=2;
    lincmin=0;
    
    saddlepointSFA;

      %E-bands to integrate spectra over, to compare to ions
        E1=50; %in eV
        E2=250;
        deltaE=2.5;
   
        pltminEphot=50;
        pltmaxEphot=400;
        %pltminEphot=38;
        %pltmaxEphot=56;

        CEPperiods=2;
        CEPjitter=0;
    
    plotstuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% signal+idler sims
% 
% clear all
% close all;
% ffreqs=@(phi,phi2)[1/3,  4e13, pi+phi2*pi; 2/3, 12e13, phi*pi];  % signal+idler
% %ffreqs=@(phi,phi2)[2/3, 30e13, phi*pi];  % signal+idler
% 
% clearvars -except ffreqs;
% 
%     minEphot = 50 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
%     maxEphot = 350 / 27.2; %largest photon energy
% 
%     clearvars -except phi2 ffreqs minEphot maxEphot ; 
%     close all;
%     choice = 3;
%     saveimgs =0;
% 
%     dip=3;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p
% 
%     longestexcursion = 1; %longgest excursion time to consider (in the first classical part)
%                           %in wave-periods
%     fullpulse = 1;      %simulate whole pulse, or only the central period
%     spectraonly = 0;    %plot only the spectra in the end (=1), or also the times...
%     
%     basewavelength = 1030; %fundamental wavelength in nm
%     pulseduration = 80; %intensity-envelope FWHM of all colors, in fs
%     intensity = 1;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
%                                  % so this is just a multiplier that could take into account that our gas medium is shifted from the
%                                  % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.
% 
% 
%     CEPstp=.05;
% %    CEP=0;
%     CEP=(0:CEPstp:2-CEPstp);
%     phi2=0;
% 
%     Ip = 15.76 /27.2;   %argon
%     %Ip = 21.565 /27.2;   %neon
%     %Ip = 24.587 /27.2;   %helium
% 
%    tmp=ffreqs(0,0);
%    R=tmp(2,2)/tmp(1,2);
%     
%     savename = ['2-color_signal+idler_H1s-dip_CEP-scan_Ar.mat'];
%     savepath = ['C:\Users\haessler\Documents\matlab\lewenstein-hhg\results\2colour-shortpulses\signal+idler\',num2str(pulseduration),'fs-envelopes_R=',num2str(R),'\'];
%     clear tmp R;
% %     savename = ['test_Ar.mat'];
% %     savepath = ['C:\Users\haessler\Documents\matlab\lewenstein-hhg\results\2colour-shortpulses\signal+idler\test\'];
% 
%         
%     logcmax=-8.5;
%     logcmin=logcmax-4;
%     lincmax=1;
%     lincmin=0;
%     
%     saddlepointSFA;
% 
%       %E-bands to integrate spectra over, to compare to ions
%         E1=50; %in eV
%         E2=70;
%         %E1=40; %in eV
%         %E2=55;
%         deltaE=2.5;
% 
%         pltminEphot=50;
%         pltmaxEphot=350;
%         %pltminEphot=38;
%         %pltmaxEphot=56;
% 
%         CEPperiods=2;
%         CEPjitter=0;
%     
%     plotstuff
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Vienna exp sims
% %clear all
% %close all;
% ffreqs=@(phi,phi2)[2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi; 2, 4.4e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup
% clearvars -except ffreqs;
% 
%     minEphot = 38 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
%     maxEphot = 130 / 27.2; %largest photon energy
% 
% for phi2=0.1:0.1:0.9
% %for phi2=[0.0, 0.1, 0.2]
%     clearvars -except phi2 ffreqs minEphot maxEphot ; 
%     close all;
%     choice = 2;
% 
%     dip=1;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p
% 
%     longestexcursion = 0.5; %longgest excursion time to consider (in the first classical part)
%                           %in wave-periods
%     fullpulse = 0;      %simulate whole pulse, or only the central period
%     spectraonly = 1;    %plot only the spectra in the end (=1), or also the times...
%     shortonly =1;
%     savefigs=0;
%     
%     basewavelength = 1030; %fundamental wavelength in nm
%     pulseduration = 200; %intensity-envelope FWHM of all colors, in fs
%     intensity = 0.55*0.92*1;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
%                                  % so this is just a multiplier that could take into account that our gas medium is shifted from the
%                                  % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.
% 
% 
%     CEPstp=.05;
%     %CEP=0.6;
%     CEP=(0:CEPstp:2-CEPstp);
%     phi2=phi2;
% 
%     Ip = 15.76 /27.2;  
% 
%     savename = ['3-color_sim-exp-november_flat-dip_intensity-scaling_0.55_phi2=',num2str(phi2,'%02.2f'),'_Ar.mat'];
%     savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\3-color_sim-exp-november_flat-dip_intensity-scaling_0.55\phi2=',num2str(phi2,'%02.2f'),'_withCEPjitter\'];
% 
%     
%          logcmax=-5.5;
%          logcmin=logcmax-4;
%          lincmax=3;
%          lincmin=0;  
% 
%     saddlepointSFA;
%      
%         %E-bands to integrate spectra over, to compare to ions
%         E1=50; %in eV
%         E2=70;
%         %E1=40; %in eV
%         %E2=55;
%         deltaE=2.5;
%      
%         pltminEphot=45;
%         pltmaxEphot=75;
%         %pltminEphot=38;
%         %pltmaxEphot=56;
% 
%         CEPperiods=3;
% 
%         CEPjitter=1;
%     
%     plotstuff
% end

% %%
% ffreqs=@(phi,phi2)[2/3, 7.9e13, phi*pi; 1, 7.5e13, phi*pi; 2, 4.4e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup
% clearvars -except ffreqs;
% 
%     minEphot = 38 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
%     maxEphot = 180 / 27.2; %largest photon energy

% for phi2=0:0.1:0.9
% %for phi2=[0.0, 0.1, 0.2]
%     clearvars -except phi2 ffreqs minEphot maxEphot ; 
%     close all;
%     choice = 2;
% 
%     dip=1;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p
% 
%     longestexcursion = 0.5; %longgest excursion time to consider (in the first classical part)
%                           %in wave-periods
%     fullpulse = 0;      %simulate whole pulse, or only the central period
%     spectraonly = 1;    %plot only the spectra in the end (=1), or also the times...
%     shortonly =1;
%     savefigs=0;
%     
%     basewavelength = 1030; %fundamental wavelength in nm
%     pulseduration = 200; %intensity-envelope FWHM of all colors, in fs
%     intensity = 0.8*0.92*1;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
%                                  % so this is just a multiplier that could take into account that our gas medium is shifted from the
%                                  % focus, 8% loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.
% 
% 
%     CEPstp=.05;
%     %CEP=0.6;
%     CEP=(0:CEPstp:2-CEPstp);
%     phi2=phi2;
% 
%     Ip = 15.76 /27.2;  
% 
%     savename = ['3-color_sim-exp-november_flat-dip_intensity-scaling_0.8_phi2=',num2str(phi2,'%02.2f'),'_Ar.mat'];
%     savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\3-color_sim-exp-november_flat-dip_intensity-scaling_0.8\phi2=',num2str(phi2,'%02.2f'),'_withCEPjitter\'];
% 
%     
%          logcmax=-6;
%          logcmin=logcmax-4;
%          lincmax=8;
%          lincmin=0;  
% 
%     saddlepointSFA;
%      
%         %E-bands to integrate spectra over, to compare to ions
%         E1=50; %in eV
%         E2=70;
%         %E1=40; %in eV
%         %E2=55;
%         deltaE=2.5;
%      
%         pltminEphot=45;
%         pltmaxEphot=75;
%         %pltminEphot=38;
%         %pltmaxEphot=56;
% 
%         CEPperiods=3;
% 
%         CEPjitter=1;
%     
%     plotstuff
% end



