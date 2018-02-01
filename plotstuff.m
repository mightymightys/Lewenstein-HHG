%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT results contained in the data-file saved by 'saddlepointSFA'. 
%
close all;

if (~exist('logcmax','var'))
    logcmax=-6;
end
if (~exist('logcmin','var'))
    logcmin=logcmax-4;
end
if (~exist('lincmax','var'))
    lincmax=2;
end
if (~exist('shortlincmax','var'))
    shortlincmax=.05*lincmax;
end
if (~exist('lincmin','var'))
    lincmin=0;
end
% logcmax=-4;
% logcmin=logcmax-5;
% 
% lincmax=15;
% lincmin=0;



if (~exist('savename','var')) || (~exist('savepath','var'))
    [savename,savepath] = uiputfile('results/*.mat','When done, save figures as:');
    disp('Figures will be saved as:')
    disp(savepath)
    disp(['...',savename(1:end-4),'...png'])
end
   
if (~exist('savefigs','var'))
    savefigs = 0;
end

if (~exist('fullpulse','var'))
    fullpulse = 0;
end

if (~exist('spectraonly','var'))
    spectraonly = 1;
end

if (~exist('waveperiod','var'))
    freqs=ffreqs(0,0);
    freqssorted = sort(freqs(:,1));
    waveperiod = 1/(freqssorted(2)-freqssorted(1)) * 2*pi/omega;
end
    
if (~exist('Ephot','var'))
    minEphot = results{1,1,1}(1,7);   %for old save-files...
    maxEphot = results{1,1,1}(end,7);
    Estp = results{1,1,1}(2,7)-results{1,1,1}(1,7);
    Ephot=(minEphot:Estp:maxEphot)';
else
    Estp = Ephot(2)-Ephot(1);
    minEphot=Ephot(1);
    maxEphot=Ephot(end);
end

if (~exist('pltminEphot','var')) ||  (~exist('pltmaxEphot','var'))
    pltminEphot=min(Ephot)*27.2;
    pltmaxEphot=max(Ephot)*27.2;
end

if (~exist('interpfactor','var'))
    interpfactor = ceil(omega/4/Estp);
end
%[~,Imin]=min(abs(Ephot-pltminEphot/27.2));
%[~,Imax]=min(abs(Ephot-pltmaxEphot/27.2));
%plotindcs = Imin:interpfactor:Imax;   %plot dots only for these energy-indices, otherwise figrues gets super-slow and files gigantic
plotindcs = 1:round(interpfactor/2):length(Ephot);

if (~exist('dE','var'))
    dE=1;
end
if (~exist('E1','var'))
    E1=min(Ephot*27.2)+dE/2;
end

if (~exist('E1','var'))
    E2=max(Ephot*27.2)-dE/2;
end


if length(CEP)>1
    CEPstp=abs(CEP(2)-CEP(1));
end
if length(size(results))==3
    [k numevents m] = size(results);
elseif length(size(results))==2
    [k numevents] = size(results);
    m=1;
end

if ~(spectraonly==1)
        for n=1:length(CEP)
         figure
          hold on;
          for event=1:numevents
           if (isempty(results{1,event,n})) || (max(isnan(results{1,event,n}(:,6)))==1)
           else
           scatter(Ephot(plotindcs)*27.2,results{1,event,n}(plotindcs,5)*24.2/1000,3, log10(Ephot(plotindcs).^4.*results{1,event,n}(plotindcs,6).*conj(results{1,event,n}(plotindcs,6))),'filled');
           scatter(Ephot(plotindcs)*27.2,results{2,event,n}(plotindcs,5)*24.2/1000,3, log10(Ephot(plotindcs).^4.*results{2,event,n}(plotindcs,6).*conj(results{2,event,n}(plotindcs,6))),'filled');
           %[tmp,ind] = min(abs(Ephot*27.2 - 60));
           %y = results{1,event,n}(ind,5)*24.2/1000;
           %text(results{1,event,n}(ind,7)*27.2, y, [num2str(CEP(n)),'pi'], 'HorizontalAlignment','center');
           end
          end
          text(0.9, 0.9, ['CEP = ',num2str(CEP(n)),'pi'], 'Units','normalized','HorizontalAlignment','right');
          ylabel('Re(tr)-Re(ti) = excursion time (fs)')
          xlabel('Photon energy (eV)','HorizontalAlignment','center')
          box on
          axis tight
            load('MyGermanColormap','mycmap')
            set(gcf,'Colormap',mycmap)
            caxis([logcmin logcmax])
          set(gcf,'Renderer','painters',...
              'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[9.5 6.5],'PaperPosition',[0.25 0.25 9 6]);
          set(gcf,'DefaultAxesFontSize',9);
          set(gca,'FontSize',9);
          
            if exist([savepath,'excursiontime\'],'dir')==0
                mkdir([savepath,'excursiontime\']);
            end
            if savefigs ==1
            saveas(gcf,[[savepath,'excursiontime\'],'excursiontime_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r600',[[savepath,'excursiontime\'],'excursiontime_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.png'])
            end
            print('-dpdf','-painters',[[savepath,'excursiontime\'],'excursiontime_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.pdf'])
         
        end
    
        close all;
    
        for n=1:length(CEP)
         figure
         subplot('Position',[0.18 0.45 0.77 0.52])
         %subplot(3,1,[1 2]) 
         hold on;
          for event=1:numevents
           if (~isempty(results{1,event,n})) && (~max(isnan(results{1,event,n}(:,6)))==1)
               scatter(results{1,event,n}(plotindcs,3)*24.2/1000,Ephot(plotindcs)*27.2,3, log10(Ephot(plotindcs).^4.*results{1,event,n}(plotindcs,6).*conj(results{1,event,n}(plotindcs,6))),'filled');
               scatter(results{2,event,n}(plotindcs,3)*24.2/1000,Ephot(plotindcs)*27.2,3, log10(Ephot(plotindcs).^4.*results{2,event,n}(plotindcs,6).*conj(results{2,event,n}(plotindcs,6))),'filled');
               scatter(results{1,event,n}(plotindcs,1)*24.2/1000,Ephot(plotindcs)*27.2,3, log10(Ephot(plotindcs).^4.*results{1,event,n}(plotindcs,6).*conj(results{1,event,n}(plotindcs,6))),'filled');
               scatter(results{2,event,n}(plotindcs,1)*24.2/1000,Ephot(plotindcs)*27.2,3, log10(Ephot(plotindcs).^4.*results{2,event,n}(plotindcs,6).*conj(results{2,event,n}(plotindcs,6))),'filled');
           end
               %[tmp,ind] = min(abs(Ephot*27.2 - 60));
               %y = results{1,event,n}(ind,3)*24.2/1000;
               %text(results{1,event,n}(ind,7)*27.2, y, [num2str(CEP(n)),'pi'], 'HorizontalAlignment','center');
          end
              box on
              axis tight
              load('MyGermanColormap','mycmap')
                  set(gcf,'Colormap',mycmap)
              caxis([logcmin logcmax])
              xlim([min(t) max(t)]*24.2/1000)
              ylim([0.98*minEphot maxEphot]*27.2)
            set(gcf,'DefaultAxesFontSize',9);
            set(gca,'FontSize',9);
            %xlabel('Re(tr) or Re(ti) (fs)')
            ylabel('Photon energy (eV)','HorizontalAlignment','center')
            set(gca,'XTickLabel','')
            text(0.95, 0.9, ['CEP = ',num2str(CEP(n),'%02.2f'),'\pi'], 'Units','normalized','HorizontalAlignment','right');

         subplot('Position',[0.18 0.15 0.77 .29])
         %subplot(3,1,3)
              freqs=ffreqs(CEP(n),phi2);
              Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim);
              %tt=(min(results{2,1,n}(:,1))-40:max(results{2,numevents,n}(:,3))+40);
              %tt=t(1:round(length(t)/2));
              tt=t;
              plot(tt*24.2/1000,Efield(tt),'red');
              line([min(tt*24.2/1000) max(tt*24.2/1000)],[0 0],'LineStyle','--','Color',[.5 .5 .5]);
              xlabel('Time (fs)')
              ylabel('E-field (a.u.)','HorizontalAlignment','center')
              box on
              xlim([min(t) max(t)]*24.2/1000)
              ylim([-0.09 0.09])

          set(gcf,'Renderer','painters',...
              'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[9.5 6.5],'PaperPosition',[0.25 0.25 9 6]);
          set(gcf,'DefaultAxesFontSize',9);
          set(gca,'FontSize',9);
          
              if exist([savepath,'ti-and-tr\'],'dir')==0
                mkdir([savepath,'ti-and-tr\']);
              end
          
              %print('-dpng','-r600',[[savepath,'ti-and-tr\'],'ti-and-tr_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.png'])
              print('-dpdf','-painters',[[savepath,'ti-and-tr\'],'ti-and-tr_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.pdf'])
         
         subplot('Position',[0.18 0.45 0.77 0.52])
         %subplot(3,1,[1 2])
         hold on;
          for event=1:numevents
           if (~isempty(results{1,event,n})) && (~max(isnan(results{1,event,n}(:,6)))==1)
              plot(results{3,event,n}(:,3)*24.2/1000,(results{3,event,n}(:,1)+1.3*Ip)*27.2,'-',results{4,event,n}(:,3)*24.2/1000,(results{4,event,n}(:,1)+1.3*Ip)*27.2,'-');
              plot(results{3,event,n}(:,2)*24.2/1000,(results{3,event,n}(:,1)+1.3*Ip)*27.2,'-',results{4,event,n}(:,2)*24.2/1000,(results{4,event,n}(:,1)+1.3*Ip)*27.2,'-');
              [y,ind] = max(results{3,event,n}(:,1)+1.3*Ip);
              text(results{3,event,n}(ind,2)*24.2/1000, y*27.2+5, num2str(event), 'HorizontalAlignment','center');
           end
         end
          if savefigs == 1
              saveas(gcf,[[savepath,'ti-and-tr\'],'ti-and-tr_CEP=',num2str(CEP(n),'%02.2f'),'pi_',savename(1:end-4),'.fig'],'fig')
          end
        end
    
        close all;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate and plot spectra

%Ephot=interp(Ephot,4);
spectamplshort = zeros(length(Ephot),length(CEP));
spectampllong = zeros(length(Ephot),length(CEP));
spectamplmacro = zeros(length(Ephot),length(CEP));
for nCEP=1:length(CEP)    
    for n=1:numevents
        if (isempty(results{1,n,nCEP})) || (isempty(results{2,n,nCEP}))
        elseif  (max(isnan(results{1,n,nCEP}(:,6)))==1) || (max(isnan(results{2,n,nCEP}(:,6)))==1)
            disp(['Dropping event no.',num2str(n),' for CEP=',num2str(CEP(nCEP)),'pi.'])
        else
            %spectamplshort(:,nCEP) = spectamplshort(:,nCEP) + interp(results{1,n,nCEP}(:,6),4); 
            %spectampllong(:,nCEP) = spectampllong(:,nCEP) + interp(results{2,n,nCEP}(:,6),4);
            % directly interpolating in the complex plane like this does
            % not work - where it counts, the phase varies rapidly and the
            % thing just rotates around like crazy and if you don't have
            % enough samples in the first place, interpolating will not
            % help.  maybe if you interpolate modulus and (unwrapped) phase
            % individually? 
            % nah, actually,it's best to interpolate before the spectral
            % amplitudes are computed. there, things are still smooth.
            spectamplshort(:,nCEP) = spectamplshort(:,nCEP) + results{1,n,nCEP}(:,6);
            spectampllong(:,nCEP)  = spectampllong(:,nCEP) + results{2,n,nCEP}(:,6);
            penaltyexponent=4;
            spectamplmacro(:,nCEP) = spectamplmacro(:,nCEP) + (waveperiod/4)^penaltyexponent*...
                                        (((results{1,n,nCEP}(:,5)).^(-penaltyexponent)).*results{1,n,nCEP}(:,6)+...
                                         ((results{2,n,nCEP}(:,5)).^(-penaltyexponent)).*results{2,n,nCEP}(:,6));
        end
    end
end

if E1>E2
    EE1=E1;
    E1=E2;
    E2=EE1;
end
[~,E1indx]=min(abs(Ephot*27.2-E1));
[~,E2indx]=min(abs(Ephot*27.2-E2));
dE=ceil(deltaE/27.2/Estp);
if E1indx<dE+1
    E1indx=1+dE;
end
if E2indx>length(Ephot)-dE
    E2indx=E2indx-dE;
end
spectrashort=abs(Ephot*ones(1,length(CEP)).^4.*(spectamplshort).*conj(spectamplshort));
spectralong =abs(Ephot*ones(1,length(CEP)).^4.*(spectampllong).*conj(spectampllong));
spectraall  =abs(Ephot*ones(1,length(CEP)).^4.*(spectamplshort+spectampllong).*conj(spectamplshort+spectampllong));
spectramacro=abs(Ephot*ones(1,length(CEP)).^4.*(spectamplmacro).*conj(spectamplmacro));

if exist('CEPjitter','var') && CEPjitter>0
    [x,y]=meshgrid(1:length(Ephot),1:length(CEP));
    sigmaCEP = CEPjitter/pi /CEPstp / (2*sqrt(2*log(2))); %standard deviation of the CEP in pixels
    %filter=  exp(-((y-length(CEP)/2)/(0.04*length(CEP))).^2);
    filter =  exp( -(y-length(CEP)/2).^2  /(2* (sigmaCEP/2/pi)^(-2) ) );
    
    %tmp=max(max(spectramacro));
    spectramacro = abs(ifft2(ifftshift(fftshift(fft2(spectramacro)).*filter')));
    %spectramacro = spectramacro/max(max(spectramacro))*tmp*0.8;
    %tmp=max(max(spectraall));
    spectraall = abs(ifft2(ifftshift(fftshift(fft2(spectraall)).*filter')));
    %spectraall = spectraall/max(max(spectraall))*tmp*0.8;
    %tmp=max(max(spectrashort));    
    spectrashort = abs(ifft2(ifftshift(fftshift(fft2(spectrashort)).*filter')));
    %spectrashort = spectrashort/max(max(spectrashort))*tmp*0.8;
    %tmp=max(max(spectralong));
    spectralong = abs(ifft2(ifftshift(fftshift(fft2(spectralong)).*filter')));
    %spectralong = spectralong/max(max(spectralong))*tmp*0.8;
    % these 2d-fft with filter only in one dimension gives the exact same
    % results as filtering at each energy-value individually
end

sumspectrashort1=sum(spectrashort(E1indx-dE:E1indx+dE,:));
sumspectrashort2=sum(spectrashort(E2indx-dE:E2indx+dE,:));
sumspectralong1=sum(spectralong(E1indx-dE:E1indx+dE,:));
sumspectralong2=sum(spectralong(E2indx-dE:E2indx+dE,:));
sumspectraall1=sum(spectraall(E1indx-dE:E1indx+dE,:));
sumspectraall2=sum(spectraall(E2indx-dE:E2indx+dE,:));
sumspectramacro1=sum(spectramacro(E1indx-dE:E1indx+dE,:));
sumspectramacro2=sum(spectramacro(E2indx-dE:E2indx+dE,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also calculate the ADK ionization rates and integrate -> see the ionzation
% level vary with CEP
ions=zeros(size(CEP));
for n=1:length(CEP)
     freqs=ffreqs(CEP(n),phi2);
              Efield = @(t)fEfield(t,I0,tau,omega,freqs,tlim);
              tt=(-0.5*waveperiod:t(2)-t(1):0.5*waveperiod);
              %tt=(-tau:t(2)-t(1):tau);
                
            nstar=1/sqrt(2*Ip);
            lstar=nstar-1;
            l=1; m=0;
            A=2^(2*nstar)/(nstar*gamma(nstar+lstar+1)*gamma(nstar-lstar));
            B=(2*l+1)*factorial(l+abs(m))/1;
            w=A*B*Ip*(2*(2*Ip)^1.5 ./abs(Efield(tt))).^(2*nstar-abs(m)-1) .*  exp(-2*(2*Ip)^1.5 / 3 ./ abs(Efield(tt)));
 
     %ions(n)=sum( exp(-2*(2*Ip)^1.5 / 3 ./ abs(Efield(tt))));
     ions(n)=1-exp(-sum(w)*(t(2)-t(1))); %actually this gives almost the same as just integrating the ADK rate... and there's not visible difference to just leaving all the funny pre-factors out.
end



%figure
%    imagesc(CEP,Ephot*27.2,log(spectraall))
%    set(gca,'YDir','normal')

%figure
%    imagesc(CEP,Ephot*27.2,log(spectrashort))
%    set(gca,'YDir','normal')
    
if exist([savepath,'spectra\'],'dir')==0
    mkdir([savepath,'spectra\']);
end

plotspectra
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%and the same thing for phi2+pi
%%%%%%%%%%%%
% if you also want the spectra for phi2+pi, all you have to do is take the spectra for phi2 and shift the
% CEP dependence by pi.
%(shifting phi+1 by pi flips the 0.5 up-down. Additionally shifting the CEP by pi does an up-down-flip of
% the 1+1.5-combination, such that the total field is flipped up-down, and
% the spectrum is the same as before the flip...
if length(CEP)>1 
    spectrashort_orig = spectrashort;
    spectrashort = [spectrashort(:,round(1/CEPstp)+1:end),spectrashort(:,1:round(1/CEPstp))];
    spectralong_orig = spectralong;
    spectralong = [spectralong(:,round(1/CEPstp)+1:end),spectralong(:,1:round(1/CEPstp))];
    spectraall_orig = spectraall;
    spectraall = [spectraall(:,round(1/CEPstp)+1:end),spectraall(:,1:round(1/CEPstp))];
    spectramacro_orig = spectramacro;
    spectramacro = [spectramacro(:,round(1/CEPstp)+1:end),spectramacro(:,1:round(1/CEPstp))];
    ions_orig = ions;
    ions = [ions(:,round(1/CEPstp)+1:end), ions(:,1:round(1/CEPstp))];
    sumspectrashort1_orig=sumspectrashort1;
    sumspectrashort2_orig=sumspectrashort2;
    sumspectrashort1 = [sumspectrashort1(:,round(1/CEPstp)+1:end), sumspectrashort1(:,1:round(1/CEPstp))];
    sumspectrashort2 = [sumspectrashort2(:,round(1/CEPstp)+1:end), sumspectrashort2(:,1:round(1/CEPstp))];
    sumspectralong1_orig = sumspectralong1;
    sumspectralong2_orig = sumspectralong2;
    sumspectralong1 = [sumspectralong1(:,round(1/CEPstp)+1:end), sumspectralong1(:,1:round(1/CEPstp))];
    sumspectralong2 = [sumspectralong2(:,round(1/CEPstp)+1:end), sumspectralong2(:,1:round(1/CEPstp))];
    sumspectraall1_orig = sumspectraall1;
    sumspectraall2_orig = sumspectraall2;
    sumspectraall1 = [sumspectraall1(:,round(1/CEPstp)+1:end), sumspectraall1(:,1:round(1/CEPstp))];
    sumspectraall2 = [sumspectraall2(:,round(1/CEPstp)+1:end), sumspectraall2(:,1:round(1/CEPstp))];
    sumspectramacro1_orig = sumspectramacro1;
    sumspectramacro2_orig = sumspectramacro2;
    sumspectramacro1 = [sumspectramacro1(:,round(1/CEPstp)+1:end), sumspectramacro1(:,1:round(1/CEPstp))];
    sumspectramacro2 = [sumspectramacro2(:,round(1/CEPstp)+1:end), sumspectramacro2(:,1:round(1/CEPstp))];
    
    savename_orig=savename;
    rplcindx = regexp(savename, 'phi2=');
    if (~isempty(rplcindx))
        savename(rplcindx:rplcindx+8)=['phi2=',num2str(phi2+1,'%02.2f')];
        plotspectra;
    end
    
    savename = savename_orig;
    spectrashort = spectrashort_orig;
    spectralong = spectralong_orig;
    spectraall = spectraall_orig;
    spectramacro = spectramacro_orig;
    sumspectrashort1=sumspectrashort1_orig;
    sumspectralong1 = sumspectralong1_orig;
    sumspectraall1 = sumspectraall1_orig;
    sumspectramacro1 = sumspectramacro1_orig;
    sumspectrashort2=sumspectrashort2_orig;
    sumspectralong2 = sumspectralong2_orig;
    sumspectraall2 = sumspectraall2_orig;
    sumspectramacro2 = sumspectramacro2_orig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%in the case of a calculation for the full pulse, also save spectra with
%smoothing to remove harmonics.

return 
if fullpulse == 1
    
    [n m ] = size(spectraall);
    
    % smooth in CEP-direction to simulate CEP-jitter
    %for k=1:n
    %   sspectraall(k,:) = smooth(spectraall(k,:),5);
    %end

    %smooth in spectral direction to remove harmonics
    for k=1:m
       spectraall(:,k) = smooth(spectraall(:,k),round(2*omega/Estp));
       spectrashort(:,k) = smooth(spectrashort(:,k),round(2*omega/Estp));
       spectralong(:,k) = smooth(spectralong(:,k),round(2*omega/Estp));
    end
    
savepath =[savepath,'smoothed']; 
if exist([savepath,'spectra\'],'dir')==0
    mkdir([savepath,'spectra\']);
end

plotspectra   

if length(CEP)>1 
    % and the same stuff for phi2+pi
    spectrashort = [spectrashort(:,round(1/CEPstp)+1:end),spectrashort(:,1:round(1/CEPstp))];
    spectralong = [spectralong(:,round(1/CEPstp)+1:end),spectralong(:,1:round(1/CEPstp))];
    spectraall = [spectraall(:,round(1/CEPstp)+1:end),spectraall(:,1:round(1/CEPstp))];

    if (~isempty(rplcindx))
    savename(rplcindx:rplcindx+8)=['phi2=',num2str(phi2+1,'%02.2f')];
    plotspectra;
    end
end
end

