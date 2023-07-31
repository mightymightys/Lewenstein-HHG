for R=0.25; %[0.06, 0.1, 0.25, 0.5] % R is the intensity ratio of the two (main) colors, i.e. fundamental and OPA
 for secondwavelength=1257; % [1155, 1200 ,1257, 1300] %,1440, 1715, 1325, 1370, 1333] % the OPA wavelength; increase "maxEphot" for longer wavelengths! 
    R
    clearvars -except secondwavelength R phi2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE THE DRIVING PULSE
    basewavelength = 800; %fundamental wavelength in nm
    
    pulseduration = [100;35];  %intensity-envelope FWHM of the colors (cos^2), in fs, in the same order as in ffreqs
                               %Attention, don't make this longer than freqs (i.e. give more pulse durations than color-components) because
                               %then the program will simply make all pulses as long as the max value (why the heck did I do that?)
                               
    intensity = 1e14;  %total intensity in W*cm^-2 (I may also add the individual peak intensities of each color into the ffreqs-matrix.
                       %        so this is just a multiplier that could take into account that our gas medium is shifted from the
                       %        focus, loss at the entrance window etc. the value is then converted to a.u. in the saddlepointSFA-script.

    CEPstp=.05;
    CEP=(0:CEPstp:1-CEPstp);  %CEPs to scan over
    CEP = 0;
    phi2=0;                   %rel. phase of a potential third color
   
    % the very important ffreqs-matrix!  Different rows are for different
    % color components, and the 3 columns are 1) carrier-frequency/fundamental, 2) fraction of "intensity", 3) phase.
    ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1),     phi*pi; 1, 1/(R+1),      phi*pi];  % 2-color 1+1,5
    %ffreqs=@(phi,phi2)[basewavelength/secondwavelength, R/(R+1+0.1), phi*pi; 1, 1/(R+1+0.1),  phi*pi; 2, 0.1/(R+1+0.1), 2*phi*pi + phi2*pi]; %3-color 1+1,5+0,5,

% LEWENSTEIN-CALCULATION PARAMETERS

    Ip = 15.76 /27.2; %your target's ionization potential in atomic units
    %Ip = 21.56 /27.2;  
    %Ip = 24.6 /27.2;

    fullpulse = 1;   %simulate whole pulse, or only the central period

    minEphot = 30 /27.2; %smallest harmonic photon energy (atomic units) you care about (we won't even look for quantum trajs. for smaller energies than that)
    maxEphot = 75 /27.2; %largest photon energy (atomic units)

    dip=3;  %what kind of dipole matrix elements to use: 1..flat, 2..Gaussian, 3..scaled-hydrogen-1s, 4..scaled-hydrogen-3p

    %longestexcursion = .5; %longgest excursion time to consider (in the first classical part) in wave-periods: 0.5 is good for 3 colors, otherwise do 1 period
    longestexcursion = 1;
    
    Nperiods = 1;    %over how many waveperiods do you want to calculate (this is where the ionization times are initially scanned over)
    shortonly = 1;   %only consider short trajectories and skip all steps treating the long ones
    spectraonly = 1; %plot only the spectra in the end (=1), or also the times etc...
    
    savename = ['2-color_',num2str(basewavelength),'+',num2str(secondwavelength),'_R=',num2str(R,'%02.2f'),'_',num2str(pulseduration(2)),'+',num2str(pulseduration(1)),'fs-pulses_',num2str(1e14,'%02.1e'),'_CEP-scan_Ar.mat'];
    savepath = ['C:\Users\haessler\Documents\matlab\lewenstein-hhg\results\ATTOLAB\'];
    
    logcmax=-6;         %colorscale limits for plots in logscale
    logcmin=logcmax-4;
    lincmax=1e2;        %colorscale limits for plots in linscale
    lincmin=0;
    
    saddlepointSFA;  % to the Lewenstein thing, returns a ton of complex valued ionization and recollision times and saddle point momenta
                     % and saves all of that stuff in the 'savefile' in 'savepath' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% you could stop here and do the analysis of the results later, or continue and have some pretty plots right away:

% CALCULATE HHG SPECTRA
    %Ephot plotting limis in eV
    pltminEphot=minEphot*27.2;
    pltmaxEphot=maxEphot*27.2;
 
	CEPperiods=1; % how many CEP periods to plot 5it may look prettier to plot several...)

	CEPjitter=0; % CEP jitter in rad FWHM to simulate (convolution before plotting)
    
    plotspecmaps = 0; % compute and plot the fancy but heavy spectral maps or not
    
    compute_spectra  % compute the spectra from the Lewenstein results 
    
% CALCULATE TEMPORAL PROFILES
    Ecenter = 65;   %center of spectral filter 
    Ewidth= 10;     % width of spectral filter
    specfilter=2;   % type of filter: 1..Gaussian, 2..8th-order-super-Gaussian 
    
    temporalprofile % calculates temp profiles and plots them as function of the CEP 
     xlabel('Time (fs)')
     ylabel('CEP (\pi rad)')
     xlim([-20,20])

% the following are only useful if you don't do a scan but calculate for a single driver waveform only 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "plot_temporalprofiles" is a script for treatment of the saved results
% file, i.e. it contains its own def. of spectral filter and then again
% computes_spectra and temp. profiles and all - you can use this to play
% with the analysis of saved results from the Lewenstein-model-calculation
currdir=pwd;
cd(savepath)
plot_temporalprofiles 
cd(currdir)