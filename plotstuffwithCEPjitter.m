%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT results contained in the data-file saved by 'saddlepointSFA'. 
%
savepath=savepath(42:end);

logcmax=-6;
logcmin=logcmax-5;

lincmax=1.2;
lincmin=0;



if (~exist('savename','var')) || (exist('savepath','var')==0)
    [savename,savepath] = uiputfile('results/*.mat','When done, save figures as:');
    disp('Figures will be saved as:')
    disp(savepath)
    disp(['...',savename(1:end-4),'...png'])
end
   
if (~exist('saveimgs','var'))
    saveimgs = 1;
end

if (~exist('fullpulse','var'))
    fullpulse = 1;
end

if (~exist('spectraonly','var'))
    spectraonly = 1;
end

[k numevents m] = size(results);
minEphot = results{1,1,1}(1,7);
maxEphot = results{1,1,1}(end,7);
Estp = results{1,1,1}(2,7)-results{1,1,1}(1,7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate and plot spectra
for n=1:length(CEP)
    ttmp(n,:)= size(results{1,1,n}(:,1));
end
tmp = max(ttmp(:,1));
Ephot=minEphot:Estp:minEphot+(tmp-1)*Estp;

spectamplshort = 1e-15* ones(tmp,length(CEP));
spectampllong = 1e-15* ones(tmp,length(CEP));
for nCEP=1:length(CEP)    
    for n=1:numevents
        if (isempty(results{1,n,nCEP})) || (isempty(results{2,n,nCEP}))
        elseif  (max(isnan(results{1,n,nCEP}(:,6)))==1) || (max(isnan(results{2,n,nCEP}(:,6)))==1)
            disp(['Dropping event no.',num2str(n),' for CEP=',num2str(CEP(nCEP)),'pi.'])
        else
            spectamplshort(1:length(results{1,n,nCEP}(:,6)),nCEP) = spectamplshort(1:length(results{1,n,nCEP}(:,6)),nCEP) + results{1,n,nCEP}(:,6);
            spectampllong(1:length(results{1,n,nCEP}(:,6)),nCEP) = spectampllong(1:length(results{1,n,nCEP}(:,6)),nCEP) + results{2,n,nCEP}(:,6);
        end
    end
end
spectrashort=abs(Ephot'*ones(1,length(CEP)).^4.*(spectamplshort).*conj(spectamplshort));
spectralong=abs(Ephot'*ones(1,length(CEP)).^4.*(spectampllong).*conj(spectampllong));
spectraall=abs(Ephot'*ones(1,length(CEP)).^4.*(spectamplshort+spectampllong).*conj(spectamplshort+spectampllong));

    

    spectraall=[spectraall,spectraall,spectraall];

    [n m ] = size(spectraall);
    
    % smooth in CEP-direction to simulate CEP-jitter
    for k=1:n
       spectraall(k,:) = smooth(spectraall(k,:),7);
    end

    spectraall = spectraall(:,m/3+1:2*m/3);
    
    %smooth in spectral direction to remove harmonics
   % for k=1:m
   %    spectraall(:,k) = smooth(spectraall(:,k),11);
   %    spectrashort(:,k) = smooth(spectrashort(:,k),11);
   %    spectralong(:,k) = smooth(spectralong(:,k),11);
   % end
    
savepath =[savepath,'smoothed']; 
if (~exist([savepath,'spectra\'],'dir'))
    mkdir([savepath,'spectra\']);
end

%plotspectra   
  
    % and the same stuff for phi2+pi
    CEPstp=CEP(2)-CEP(1);
    
    %spectrashort = [spectrashort(:,round(1/CEPstp)+1:end),spectrashort(:,1:round(1/CEPstp))];
    %spectralong = [spectralong(:,round(1/CEPstp)+1:end),spectralong(:,1:round(1/CEPstp))];
    spectraall = [spectraall(:,round(1/CEPstp)+1:end),spectraall(:,1:round(1/CEPstp))];

    rplcindx = regexp(savename, 'phi2=');
    savename(rplcindx:rplcindx+8)=['phi2=',num2str(phi2+1,'%02.2f')];

%    plotspectra;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
    surface([CEP,CEP+2,CEP+4,CEP+6],Ephot*27.2,log10([spectraall,spectraall,spectraall,spectraall]));
    view(2)
        set(gca,'FontSize',10)
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        %load('MyGermanColormap','mycmap')
        %set(gcf,'Colormap',mycmap)
        %caxis([log10(max(max(spectraall)))-4-.5, log10(max(max(spectraall)))-.5])
        caxis([logcmin logcmax]);
        set(gcf,'Renderer','zbuffer',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[9.5 6.5],'PaperPosition',[0.25 0.25 9 6]);
        set(gcf,'DefaultAxesFontSize',8);
        set(gca,'FontSize',8);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','log_{10}[ \omega^4 |d(\omega)|^2 ]    (arb.u.)');
        if saveimgs == 1
       %     saveas(gcf,[[savepath,'spectra\'],'log-spectra_CEP-scan_short+long-traj_',savename(1:end-4),'.fig'],'fig')
            print('-dpng','-r1200',[[savepath,'spectra\'],'log-spectra_CEP-scan_short+long-traj_',savename(1:end-4),'.png'])
        end
        
figure
    surface([CEP,CEP+2,CEP+4,CEP+6],Ephot*27.2,(1e6*[spectraall,spectraall,spectraall,spectraall]));
    view(2)
        set(gca,'FontSize',10)
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        %load('MyGermanColormap','mycmap')
        %set(gcf,'Colormap',mycmap)
        %caxis([0, 0.75*1e6*(max(max(spectraall)))])
        caxis([lincmin lincmax]);
        set(gcf,'Renderer','zbuffer',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[9.5 6.5],'PaperPosition',[0.25 0.25 9 6]);
        set(gcf,'DefaultAxesFontSize',8);
        set(gca,'FontSize',8);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','\omega^4 |d(\omega)|^2   (arb.u.)'); 
        if saveimgs == 1
       %     saveas(gcf,[[savepath,'spectra\'],'spectra_CEP-scan_short+long-traj_',savename(1:end-4),'.fig'],'fig')
            print('-dpng','-r1200',[[savepath,'spectra\'],'spectra_CEP-scan_short+long-traj_',savename(1:end-4),'.png'])
        end
        