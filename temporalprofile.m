%close all;
%% you have to run compute_spectra before

if (~exist('specfilter', 'var'))
specfilter=2;
end

padto = length(Ephot)*25; %zero-pad to this many samples

%% filter spectra
if (~exist('Ecenter', 'var'))
Ecenter = 110;  %central energy in eV
end
if (~exist('Ewidth', 'var'))
Ewidth = 10; %width in eV
end

[~,E0indx]=min(abs(Ephot*27.2-Ecenter));

if specfilter==1 %if Gaussian filter
    filter = exp(-((1:length(Ephot))'-E0indx).^2 / (4* (Ewidth/27.2/Estp/2/sqrt(2*log(2)))^2));
elseif specfilter==2 %if Super-Gaussian filter
    filter = exp(-((1:length(Ephot))'-E0indx).^8 / (4* (Ewidth/27.2/Estp/2/sqrt(2*log(2)))^8));
else
    filter = 1;
end

[n,m]=size(spectamplshort);
sspectamplshort = spectamplshort .* (sqrt(filter)*ones(1,m));

%% zero-pad spectra towards zero frequency  
Estp=Ephot(2)-Ephot(1);
EEphot = [(0:Estp:(Ephot(1)-Estp))'; Ephot];
sspectamplshort = [zeros(length((0:Estp:Ephot(1)-Estp)),m); sspectamplshort];

%% zero-pad on the high freq. side for higher time resolution
    EEphot = [EEphot; (Ephot(end)+Estp:Estp:padto*Estp)'];
    sspectamplshort = [sspectamplshort; zeros(padto-length(sspectamplshort),m)];

%%
Ew = (EEphot.^2*ones(1,m)) .*sspectamplshort;
Iw=Ew.*conj(Ew);
Et = fftshift(fft(Ew,[],1),1)*Estp;   % this scaling by Estp is crucial so you can compare attopulses from diff.
                               % calcualtions with diff. Estp! (so Parseval's theorem holds: FT is unitary and preserves norm) 
It = Et .*conj(Et);
%%
freq = EEphot/ 2/pi;
Nyquist = 1/2/(freq(2)-freq(1));
tstp=2*Nyquist/(length(freq)-1);
time=(-Nyquist:tstp:Nyquist);

%%
[~,cutindx1]=min(abs(time+max(tau)));
[~,cutindx2]=min(abs(time-max(tau)));

time=time(cutindx1:cutindx2);
Et=Et(cutindx1:cutindx2,:);
It=It(cutindx1:cutindx2,:);
%%
figure;
ax = axes;
imagesc(time*24.2/1000,CEP,It')
set(ax,'YDir','normal')       
%xlim([min(t),max(t)]*24.2/1000)

return


%%
figure;
plot(time*24.2/1000, It)
[AX,H1,H2] = plotyy(time*24.2/1000, It(:,1), time*24.2/1000, Efield(time) ,'plot');

%Timelim = [min(t),min(t)+Nperiods*waveperiod]*24.2/1000+3;
Timelim = [min(t),max(t)]*24.2/1000;
set(AX(1),'Xlim', Timelim)
set(AX(2),'Xlim', Timelim) 
Yticks=(0:.5:2)*1e-12;
set(AX(1),'Ylim', [min(Yticks) max(Yticks)])
set(AX(1),'Ytick',Yticks)
Yticks=(-.1:.05:.1);
set(AX(2),'Ylim', [-0.15 0.15])
set(AX(2),'Ytick',Yticks)


set(get(AX(1),'Ylabel'),'String','Intensity (arb.u.)') 
set(get(AX(2),'Ylabel'),'String','Driving electric field (a.u.)') 
xlabel('Time (fs)','HorizontalAlignment','center')
box on
          
set(gcf,'Renderer','painters',...
  'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[12.5 8.5],'PaperPosition',[0.25 0.25 12 8]);
set(gcf,'DefaultAxesFontSize',9);
set(gca,'FontSize',9);
if exist([savepath,'tempprofile\'],'dir')==0
   mkdir([savepath,'tempprofile\']);
end

saveas(gcf,[[savepath,'tempprofile\'],'tempprofile_intensity_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.fig'],'fig')
print('-dpdf','-painters',[[savepath,'tempprofile\'],'tempprofile_intensity_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.pdf'])
 
%%
figure; plot(time*24.2/1000, real(Et(:,1)))
xlim(Timelim)

ylabel('Electric field (arb.u.)')
xlabel('Time (fs)','HorizontalAlignment','center')
box on
          
set(gcf,'Renderer','painters',...
  'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[12.5 8.5],'PaperPosition',[0.25 0.25 12 8]);
set(gcf,'DefaultAxesFontSize',9);
set(gca,'FontSize',9);
if exist([savepath,'tempprofile\'],'dir')==0
   mkdir([savepath,'tempprofile\']);
end

saveas(gcf,[[savepath,'tempprofile\'],'tempprofile_field_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.fig'],'fig')
print('-dpdf','-painters',[[savepath,'tempprofile\'],'tempprofile_field_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.pdf'])

%%
figure
%plot(EEphot*27.2,log10(Iw),Ephot*27.2,log10(spectrashort));%,Ephot*27.2,log10(specref))
%plotyy(Ephot*27.2, log10(spectrashort), Ephot*27.2, filter );
[AX,H1,H2] = plotyy(Ephot*27.2, log10(spectrashort(:,1)), Ephot*27.2, filter ,'plot');

set(AX(1),'Xlim', [min(Ephot) max(Ephot)]*27.2)
set(AX(2),'Xlim', [min(Ephot) max(Ephot)]*27.2) 

Yticks=(ceil(max(log10(spectrashort)))-3:1:ceil(max(log10(spectrashort))));
set(AX(1),'Ylim', [Yticks(1) Yticks(end)])
set(AX(1),'Ytick',Yticks)
set(AX(2),'Ylim', [0 1.1]) 

set(get(AX(1),'Ylabel'),'String','log_{10}(Intensity) (arb.u.)') 
set(get(AX(2),'Ylabel'),'String','Spectral filter') 
xlabel('Photon energy (eV)','HorizontalAlignment','center')
box on
          
set(gcf,'Renderer','painters',...
  'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[12.5 8.5],'PaperPosition',[0.25 0.25 12 8]);
set(gcf,'DefaultAxesFontSize',9);
set(gca,'FontSize',9);
if exist([savepath,'tempprofile\'],'dir')==0
   mkdir([savepath,'tempprofile\']);
end

saveas(gcf,[[savepath,'tempprofile\'],'HHG-spectrum_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.fig'],'fig')
print('-dpdf','-painters',[[savepath,'tempprofile\'],'HHG-spectrum_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.pdf'])
       