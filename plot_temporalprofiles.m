clear all;
close all;

savefigs=1;

Ecenter = 50;  %central energy in eV
Ewidth = 10; %width in eV

specfilter=2; %super-Gaussian

currdir=pwd;
%files=dir([currdir,'/*20fs-pulses_CEP=-0.20scan_Ar.mat']);
files=dir([currdir,'\*.mat']);

for i=1:length(files)
    files(i).name
    load([currdir,'/',files(i).name]);

    %cd '/Users/stefan/Documents/work/matlab/lewenstein-hhg'
    cd 'C:\Users\haessler\Documents\matlab\lewenstein-hhg'
    tau=tau(1);
    compute_spectra
    
    figure;
    imagesc(CEP,Ephot*27.2,spectrashort)
    set(gca, 'CLim', [0, .8*max(max(spectrashort))]);
    %set(gca, 'CLim', [0, 20e-8]);
    ylabel('Photon energy (eV)')
    xlabel('CEP (\pi rad)')
    set(gcf,'Renderer','painters',...
     'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[12.5 8.5],'PaperPosition',[0.25 0.25 12 8]);
    set(gcf,'DefaultAxesFontSize',9);
    set(gca,'FontSize',9);
    if savefigs==1
        print('-dpdf','-painters',[currdir,'\',savename(1:end-4),'_spectra.pdf'])
    end
    
    temporalprofile
    set(gca, 'CLim', [0, max(max(It))]);
    %set(gca, 'CLim', [0, 10e-10]);
    xlim([-20,20])
    xlabel('Time (fs)')
    ylabel('CEP (\pi rad)')
    
    set(gcf,'Renderer','painters',...
     'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[12.5 8.5],'PaperPosition',[0.25 0.25 12 8]);
    set(gcf,'DefaultAxesFontSize',9);
    set(gca,'FontSize',9);
    if savefigs==1
        print('-dpdf','-painters',[currdir,'\',savename(1:end-4),'_atto-pulses_',num2str(Ecenter),'pm',num2str(Ewidth/2),'eV.pdf'])
    end
    
    
    cd(currdir)
end

%%
clear A B
smoothfactor=100;
tmp=mean(spectrashort(:,7:13),2);
for i=1:length(Ephot)-smoothfactor
        smoothwidth=round(smoothfactor*i/length(Ephot));
        A(i) = mean(tmp(i:i+smoothwidth));
end

tmp=mean(spectrashort(:,1:7),2);
for i=1:length(Ephot)-smoothfactor
        smoothwidth=round(smoothfactor*i/length(Ephot));
        B(i) = mean(tmp(i:i+smoothwidth));
end


figure; 
plt=plot(Ephot(1:end-smoothfactor)*27.2,A,Ephot(1:end-smoothfactor)*27.2,B);
xlabel('Photon energy (eV)')
ylabel('Spectral intensity (arb.u.)')
plt(1).LineWidth=1.5;
plt(2).LineWidth=1.5;
%plt.Color=[1 0 0];

legend(['CEP averaged between ',num2str(CEP(7)),'-',num2str(CEP(13)),'\pi'],...
            ['CEP averaged between ',num2str(CEP(1)),'-',num2str(CEP(7)),'\pi']);







       