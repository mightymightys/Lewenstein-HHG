%for secondwavelength=[1200, 1300, 1400, 1500, 1545, 1600, 1700, 1800]
%for phi2=(0.1:.1:.9)
for R=[0.06, 0.1, 0.25, 0.5]
 for secondwavelength=[1155, 1200 ,1257, 1300] %,1440, 1715, 1325, 1370, 1333] % adapt the cutoff!
    R
    clearvars -except secondwavelength R phi2; 

    basewavelength = 800; %fundamental wavelength in nm
    %basewavelength = 1030;
    %pulseduration = [100;180;180]; %intensity-envelope FWHM of the colors, in fs, in the same order as in ffreqs
    %pulseduration = [50;180;180];
    pulseduration = [100;35];   %Attention, don't make this longer than freqs
                               %(i.e. give more pulse durations than color-components because
                               %then the program will simply make all pulse as long as the max value
    %pulseduration = 25;
    %intensity = 0.8*0.92;
    %intensity = 1.6e15;
    %intensity = 1;
    intensity = 1e14;  %total intensity in W*cm^-2 (I may also add the individual peak intensities of each color into the ffreqs-matrix.
                         %        so this is just a multiplier that could take into account that our gas medium is shifted from the
                         %        focus, loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.

	%CEPstp=.005;
    CEPstp=.05;
    %CEPstp=.1;
    %CEP=0.6;
    %CEP=0.85;
    %CEP=0.45;
    %CEP=0.675;
    %CEP=1.775;
    %CEP=(-1:CEPstp:1-CEPstp);
    %CEP=(0:CEPstp:2-CEPstp);
    CEP=(0:CEPstp:1-CEPstp);
    %CEP=(0.44:CEPstp:0.46);
    %CEP=(0.4:CEPstp:2.4-CEPstp);
    phi2=phi2;
    %phi2=0.1;
     
    
%ffreqs=@(phi,phi2)[2/3, 1/(1+colorratio), phi*pi; 1, colorratio/(1+colorratio), phi*pi];  % 2-color-1+1,5 
%ffreqs=@(phi,phi2)[2/3+0.002*1.5, 1/(1+colorratio), phi*pi; 1.002, colorratio/(1+colorratio), phi*pi];  % 2-color-1+1,5 
%ffreqs=@(phi,phi2)[2, 1/(1+colorratio), phi*pi+1.8*pi; 1, colorratio/(1+colorratio), 0+0.9*pi];  % 2-color-1+0,5 
%ffreqs=@(phi,phi2)[4/3, 1/(1+colorratio), phi*pi+1.8*pi; 2/3,colorratio/(1+colorratio), 0+0.9*pi];  % 2-color-1,5+0,75 
%ffreqs=@(phi,phi2)[2, 1/(1+colorratio), phi*pi; 2/3 ,colorratio/(1+colorratio), 0 ];  % 2-color-1,5+0,5 
%ffreqs=@(phi,phi2)[2/3, 6e13, phi*pi; 1, 5.7e13, phi*pi; 2, 3e12,phi*2*pi+phi2*pi];  % November2012 recent exp. setup - with intensities as estimated for the paper
%ffreqs=@(phi,phi2)[2/3, .358, phi*pi; 1, .597, phi*pi; 2, .045,phi*2*pi+phi2*pi];  % April2012 exp. setup
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, 1, phi*pi];
%ffreqs=@(phi,phi2)[1, 1, phi*pi];
ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1),     phi*pi; 1, 1/(R+1),      phi*pi];  % 2-color 1+1,5
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1+0.1), phi2*pi; 1, 1/(R+1+0.1), phi2*pi; 2, 0.1/(R+1+0.1), 2*phi2*pi + phi*pi]; %3-color 1+1,5+0,5,
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1+0.1), phi*pi; 1, 1/(R+1+0.1), phi*pi; 2, 0.1/(R+1+0.1), 2*phi*pi + phi2*pi]; %3-color 1+1,5+0,5,
%ffreqs=@(phi,phi2)[basewavelength/secondwavelength, 1/3, phi*pi; 1, 18/30, phi*pi;  2, 2/30, 2*phi*pi+phi2*pi];  % 3-color 1+1,5+0,5 again

    minEphot = 30 /27.2; %smallest harmonic photon energy you care about (we won't even look for quantum trajs. for smaller energies than that)
    maxEphot = 75 /27.2; %largest photon energy

    choice = 3 ;

    dip=3;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p

 %longestexcursion = .5; %longgest excursion time to consider (in the first classical part) in wave-periods
 longestexcursion = 1;
    Nperiods = 1; %over how many waveperiods do you want to calculate (this is where the ionization times are initially scanned over)
	fullpulse = 1;      %simulate whole pulse, or only the central period
    shortonly = 1;
    spectraonly = 1;    %plot only the spectra in the end (=1), or also the times...
    
    
    
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
    %savename = ['1-color_',num2str(secondwavelength/1000,'%02.2f'),'_',num2str(pulseduration,'%02.2f'),'fs-pulses_CEP-scan_Ar.mat'];
    savename = ['2-color_',num2str(basewavelength),'+',num2str(secondwavelength),'_R=',num2str(R,'%02.2f'),'_',num2str(pulseduration(2)),'+',num2str(pulseduration(1)),'fs-pulses_',num2str(1e14,'%02.1e'),'_CEP-scan_Ar.mat'];
    %savename = ['3-color_',num2str(basewavelength),'+',num2str(secondwavelength),'+',num2str(1/2*basewavelength),'_phi2=',num2str(phi2),'pi_R=',num2str(R,'%02.2f'),'_',num2str(pulseduration(1)),'+',num2str(pulseduration(2)),'+',num2str(pulseduration(3)),'fs-pulses_CEP-scan_Ar.mat'];
    %savepath = ['C:\Users\haessler\ownCloud\PALM-2-color\'];  
    %savepath = ['C:\Users\haessler\ownCloud\projects\PALM-2-color_FAB1\simulations\'];  
    savepath = ['C:\Users\haessler\Documents\matlab\lewenstein-hhg\results\ATTOLAB\'];
    
    logcmax=-6;
    logcmin=logcmax-4;
    lincmax=1e2;
    lincmin=0;
    
    %E-bands to integrate spectra over, to compare to ions
    E1=30; %in eV
    E2=60;
    deltaE=2.5;
     
    saddlepointSFA;
   
    %Ephot plotting limis in eV
    pltminEphot=minEphot*27.2;
    pltmaxEphot=maxEphot*27.2;
 
	CEPperiods=1;

	CEPjitter=0; %in rad FWHM
    
    %savefigs = 1;
   plotspecmaps = 0;
    
    compute_spectra
%%
    Ecenter = 65;
     Ewidth= 10;
    temporalprofile
     xlabel('Time (fs)')
     ylabel('CEP (\pi rad)')
     xlim([-20,20])
     
%     figure;
%     area(time*24.2/1000, It,'FaceColor',[0.5,0.9,1]),
%     xlim([-20,20])
%     hold on;
%     plot([-20,20], max(It)*[.1,.1],'LineStyle','--','Color','black')
%     xlabel('Time (fs)')
%     ylabel('Intensity (arb.u.)')
%     
%     figure; imagesc(Ephot*27.2,CEP,spectrashort')
%     xlabel('Photon energy (eV)')
%     ylabel('CEP (\pi rad)')
%     xlim([pltminEphot,pltmaxEphot])
%     %%save([savepath,savename(1:end-4),'_atto-pulses.mat'], 'spectrashort', 'filter', 'It', 'Et', 'time')
    
    end
end    

currdir=pwd;
cd(savepath)
plot_temporalprofiles
cd(currdir)