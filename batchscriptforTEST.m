%secondwavelength=[1200, 1300, 1400, 1500, 1545, 1600, 1700, 1800]
for basewavelength=1030
    clearvars -except basewavelength;

    %basewavelength = 1030; %fundamental wavelength in nm
    %pulseduration = [60;180]; %intensity-envelope FWHM of all colors, in fs
    pulseduration = 300;
    %intensity = 0.8*0.92;
    %intensity = 1.6e15;
    %intensity = 1;
    intensity = 1.2e14;  %total intensity in W*cm^-2 (I now add the individual peak intensities of each color into the ffreqs-matrix.
                         %        so this is just a multiplier that could take into account that our gas medium is shifted from the
                         %        focus, loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.

    CEPstp=.025;
    %CEPstp=.05;
    %CEPstp=.1;
    %CEP=0.6;
    %CEP=0.85;
    %CEP=0.75;
    CEP=0;
    %CEP=(-1:CEPstp:1-CEPstp);
    %CEP=(0:CEPstp:1-CEPstp);
    %CEP=(0.4:CEPstp:2.4-CEPstp);
    %phi2=phi2;
    phi2=0.2;
    
    
%ffreqs=@(phi,phi2)[2/3, 1/(1+colorratio), phi*pi; 1, colorratio/(1+colorratio), phi*pi];  % 2-color-1+1,5 
%ffreqs=@(phi,phi2)[2/3+0.002*1.5, 1/(1+colorratio), phi*pi; 1.002, colorratio/(1+colorratio), phi*pi];  % 2-color-1+1,5 
%ffreqs=@(phi,phi2)[2, 1/(1+colorratio), phi*pi+1.8*pi; 1, colorratio/(1+colorratio), 0+0.9*pi];  % 2-color-1+0,5 
%ffreqs=@(phi,phi2)[4/3, 1/(1+colorratio), phi*pi+1.8*pi; 2/3,colorratio/(1+colorratio), 0+0.9*pi];  % 2-color-1,5+0,75 
%ffreqs=@(phi,phi2)[2, 1/(1+colorratio), phi*pi; 2/3 ,colorratio/(1+colorratio), 0 ];  % 2-color-1,5+0,5 
%ffreqs=@(phi,phi2)[2/3, 6e13, phi*pi; 1, 5.7e13, phi*pi; 2, 3e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup - with intensities as estimated for the paper
%ffreqs=@(phi,phi2)[2/3, .358, phi*pi; 1, .597, phi*pi; 2, .045,phi*2*pi+phi2*pi];  % April2012 exp. setup
ffreqs=@(phi,phi2)[1, 1, phi*pi];
%ffreqs=@(phi,phi2)[1, 1, phi*pi];
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1), phi*pi; 1, 1/(R+1), phi*pi];  % 2-color 1+1,5
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, 1/3, phi*pi; 1, 18/30, phi*pi;  2, 2/30, 2*phi*pi+phi2*pi];  % 3-color 1+1,5+0,5 again

    minEphot = 40 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
    maxEphot = 120 / 27.2; %largest photon energy

    choice = 3;

    dip=1;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p

    %longestexcursion = .5; %longgest excursion time to consider (in the first classical part) in wave-periods
    longestexcursion = .75;
    Nperiods = 3; %over how many waveperiods do you want to calculate (this is where the ionization times are initially scanned over)
	fullpulse = 0;      %simulate whole pulse, or only the central period
    shortonly = 0;
    spectraonly = 0;    %plot only the spectra in the end (=1), or also the times...
    
    
    
    Ip = 15.76 /27.2;
    %Ip = 21.56 /27.2;  
    %Ip = 24.6 /27.2;

    %savename = ['2-color_1+1,5_ratio-',num2str(colorratio,'%02.2f'),'_single-period_CEP-scan_te-lim-to-5fs_Ar.mat'];
    %savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\2-color-1+1,5_ratio-',num2str(colorratio,'%02.2f'),'\'];
    %savename = ['2-color_1+0,5_ratio-',num2str(colorratio,'%02.2f'),'_single-period_CEP-scan_te-lim-to-3,33fs_Ar.mat'];
    %savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\2-color-1+0,5_ratio-',num2str(colorratio,'%02.2f'),'\'];
    %savename = ['1,5+0,75_ratio-',num2str(colorratio,'%02.2f'),'_single-period_CEP-scan_te-lim-to-3,33fs_Ar.mat'];
    %avepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\2-color\1,5+0,75_ratio-',num2str(colorratio,'%02.2f'),'\'];
    %savename = ['0,5+1,5_ratio-',num2str(colorratio,'%02.2f'),'_single-period_CEP-scan_te-lim-to-5fs_Ar.mat'];
    %savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\2-color\0,5+1,5_ratio-',num2str(colorratio,'%02.2f'),'\'];
    %savename = ['2-color_1+',num2str(secondwavelength/1000,'%02.2f'),'_',num2str(pulseduration,'%02.2f'),'fs-pulses_CEP-scan_Ar.mat'];
    savename = ['1-color_',num2str(basewavelength/1000,'%02.2f'),'_',num2str(pulseduration,'%02.2f'),'fs-pulses_CEP-scan_Ar.mat'];
    savepath = ['C:\Users\haessler\Dropbox\lewenstein-hhg\results\2colour-shortpulses\1micron+signal\'];
    
    logcmax=-5;
    logcmin=logcmax-4;
    lincmax=1e2;
    lincmin=0;
    
    %E-bands to integrate spectra over, to compare to ions
    E1=50; %in eV
    E2=100;
    deltaE=2.5;
     
    saddlepointSFA;
   
     %Ephot plotting limis in eV
     pltminEphot=minEphot*27.2;
     pltmaxEphot=maxEphot*27.2;
%  
% 	CEPperiods=1;
% 
% 	CEPjitter=0; %in rad FWHM
%     
    savefigs = 1;
    plotspecmaps = 0;
    
     compute_spectra
%     
      Ecenter = 90;
      Ewidth= 15;
      temporalprofile
      xlabel('Time (fs)')
      ylabel('CEP (\pi rad)')
      xlim([-15,15])
%      
     figure;
    area(time*24.2/1000, It,'FaceColor',[0.5,0.9,1]),
    xlim([-15,15])
    hold on;
    plot([-15,15], max(It)*[.1,.1],'LineStyle','--','Color','black')
    xlabel('Time (fs)')
    ylabel('Intensity (arb.u.)')
    
    figure; 
    plot(Ephot*27.2,(spectrashort(:,1)))
    hold on;
    plot(Ephot*27.2, 1e-6*filter ,'red');
    xlabel('Photon energy (eV)')
    ylabel('Intensity (arb.u.)')
    xlim([pltminEphot,pltmaxEphot])
    %%save([savepath,savename(1:end-4),'_atto-pulses.mat'], 'spectrashort', 'filter', 'It', 'Et', 'time')
%     
    end

 