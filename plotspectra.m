plotindcs = 1:1:length(Ephot);

pltCEP=[];
pltspectramacro=[];
pltspectraall=[];
pltspectrashort=[];
pltspectralong=[];
pltions=[];
pltsumspectrashort1=[];
pltsumspectralong1=[];
pltsumspectraall1=[];
pltsumspectramacro1=[];
pltsumspectrashort2=[];
pltsumspectralong2=[];
pltsumspectraall2=[];
pltsumspectramacro2=[];

if length(CEP)>1 
    for n=1:CEPperiods
        pltCEP = [pltCEP,CEP+(n-1)*(max(CEP)+CEPstp)];
        pltspectramacro = [pltspectramacro,spectramacro];
        pltspectraall = [pltspectraall,spectraall];
        pltspectrashort = [pltspectrashort,spectrashort];
        pltspectralong = [pltspectralong,spectralong];
        pltions = [pltions,ions];
        pltsumspectrashort1 = [pltsumspectrashort1,sumspectrashort1];
        pltsumspectralong1 = [pltsumspectralong1,sumspectralong1];
        pltsumspectraall1 = [pltsumspectraall1, sumspectraall1];
        pltsumspectramacro1 = [pltsumspectramacro1, sumspectramacro1];
        pltsumspectrashort2 = [pltsumspectrashort2, sumspectrashort2];
        pltsumspectralong2 = [pltsumspectralong2, sumspectralong2];
        pltsumspectraall2 = [pltsumspectraall2, sumspectraall2];
        pltsumspectramacro2 = [pltsumspectramacro2, sumspectramacro2];
    end
    pltspectramacro = pltspectramacro(plotindcs,:);
    pltspectraall = pltspectraall(plotindcs,:);
    pltspectrashort = pltspectrashort(plotindcs,:);
    pltspectralong = pltspectralong(plotindcs,:);
    pltEphot = Ephot(plotindcs);
    

% if exist('CEPjitter','var') && CEPjitter==1
%     [x,y]=meshgrid(1:length(plotindcs),1:length(pltCEP));
%     filter=  exp(-( (y-length(pltCEP)/2)/(0.15*length(pltCEP)/CEPperiods)).^2);
%     pltspectramacro = abs(ifft2(ifftshift(fftshift(fft2(pltspectramacro)).*filter')));
%     pltspectraall = abs(ifft2(ifftshift(fftshift(fft2(pltspectraall)).*filter')));
%     pltspectrashort = abs(ifft2(ifftshift(fftshift(fft2(pltspectrashort)).*filter')));
%     pltspectralong = abs(ifft2(ifftshift(fftshift(fft2(pltspectralong)).*filter')));
% end
% I add the jitter directly after the calculation in plotstuff so that the
% sums also include the averaging. It gives the same result.

else
    pltCEP=CEP;
    pltspectramacro = spectramacro(plotindcs,:);
    pltspectraall=spectraall(plotindcs,:);
    pltspectrashort=spectrashort(plotindcs,:);
    pltspectralong=spectralong(plotindcs,:);
    pltEphot = Ephot(plotindcs);
    pltions = ions;
    pltsumspectrashort1 = sumspectrashort1 ;
    pltsumspectralong1 = sumspectralong1;
    pltsumspectraall1 = sumspectraall1 ;
    pltsumspectramacro1 = sumspectramacro1;
    pltsumspectrashort2 = sumspectrashort2 ;
    pltsumspectralong2 = sumspectralong2;
    pltsumspectraall2 = sumspectraall2 ;
    pltsumspectramacro2 = sumspectramacro2;
end

%%
figure
 sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectrashort1/max(pltsumspectrashort1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectrashort2/max(pltsumspectrashort2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,log10(pltspectrashort));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2,log10(pltspectrashort));
    set(gca,'ydir','normal');
         set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
       %caxis([log10(max(max(spectrashort)))-4-0, log10(max(max(spectrashort)))-0])
       caxis([logcmin logcmax]);
       ylim([pltminEphot pltmaxEphot]);
       set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','log_{10}[ \omega^4 |d(\omega)|^2 ]    (arb.u.)');
        set(cb,'YTick', logcmin:(logcmax-logcmin)/4:logcmax)
        
        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 

        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'log-spectra_short-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'log-spectra_short-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'log-spectra_short-traj_',savename(1:end-4),'.pdf'])
        pause(2); close

%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectrashort1/max(pltsumspectrashort1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectrashort2/max(pltsumspectrashort2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,1e6*(pltspectrashort));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2,1e6*(pltspectrashort));
    set(gca,'ydir','normal');
         set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
       %caxis([0, 0.99*1e6*(max(max(spectrashort)))])
        caxis([lincmin shortlincmax]);
        ylim([pltminEphot pltmaxEphot]);
       set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','\omega^4 |d(\omega)|^2   (arb.u.)'); 
        set(cb,'YTick', lincmin:(shortlincmax-lincmin)/2:shortlincmax ,'YTickLabel',{'0','0.5','1'}) 

        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 

        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'spectra_short-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'spectra_short-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'spectra_short-traj_',savename(1:end-4),'.pdf'])
        pause(2); close

if exist('shortonly','var') && shortonly==1
    return
end
%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectramacro1/max(pltsumspectramacro1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectramacro2/max(pltsumspectramacro2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
        set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,log10(pltspectramacro));
    imagesc(pltCEP,pltEphot*27.2,log10(pltspectramacro));
    set(gca,'ydir','normal');
    %view(2)
        set(gca,'FontSize',11)
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        %caxis([log10(max(max(spectraall)))-4-.5,
        %log10(max(max(spectraall)))-.5])
        caxis([logcmin logcmax]);
        ylim([pltminEphot pltmaxEphot]);
       
        set(gcf,'DefaultAxesFontSize',10);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','log_{10}[ \omega^4 |d(\omega)|^2 ]    (arb.u.)');
        set(cb,'YTick', logcmin:(logcmax-logcmin)/4:logcmax)
        
        set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[.25 .25 11 8]);
      
        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
        set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)])
        set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02])

        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'log-spectra_all-traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r600',[[savepath,'spectra\'],'log-spectra_all_traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'log-spectra_all_traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.pdf'])
        
     %pause(2); close

%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectramacro1/max(pltsumspectramacro1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectramacro2/max(pltsumspectramacro2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
        set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,(1e5*pltspectramacro));
    imagesc(pltCEP,pltEphot*27.2,(1e5*pltspectramacro));
    set(gca,'ydir','normal');
        set(gca,'FontSize',11)
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        
        caxis([lincmin shortlincmax]);
        ylim([pltminEphot pltmaxEphot]);
        
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','\omega^4 |d(\omega)|^2  (arb.u.)'); 
        set(cb,'YTick', lincmin:(shortlincmax-lincmin)/2:shortlincmax ,'YTickLabel',{'0','0.5','1'}) 
        
        set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
      
     sp1pos = get(sp(1), 'Position');
     sp2pos = get(sp(2), 'Position');
        set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
        set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 

     if savefigs == 1
         saveas(gcf,[[savepath,'spectra\'],'spectra_all-traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.fig'],'fig')
         %print('-dpng','-r1200',[[savepath,'spectra\'],'spectra_all-traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.png'])
     end
     print('-dpdf','-painters',[[savepath,'spectra\'],'spectra_all_traj_penalty-exp=',num2str(penaltyexponent),'_',savename(1:end-4),'.pdf'])
     pause(2); close  
      
      
%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectraall1/max(pltsumspectraall1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectraall2/max(pltsumspectraall2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,log10(pltspectraall));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2,log10(pltspectraall));
    set(gca,'ydir','normal');
        set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        %caxis([log10(max(max(spectraall)))-4-.5, log10(max(max(spectraall)))-.5])
        caxis([logcmin logcmax]);
        ylim([pltminEphot pltmaxEphot]);
        set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','log_{10}[ \omega^4 |d(\omega)|^2 ]    (arb.u.)');
        set(cb,'YTick', logcmin:(logcmax-logcmin)/4:logcmax)
        
        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 

        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'log-spectra_short+long-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'log-spectra_short+long-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'log-spectra_short+long-traj_',savename(1:end-4),'.pdf'])
        pause(2); close
     
     
%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectraall1/max(pltsumspectraall1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectraall2/max(pltsumspectraall2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,(1e6*pltspectraall));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2,(1e6*pltspectraall));
    set(gca,'ydir','normal');
        set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        %caxis([0, 0.75*1e6*(max(max(spectraall)))])
        caxis([lincmin lincmax]);
        ylim([pltminEphot pltmaxEphot]);
        set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','\omega^4 |d(\omega)|^2   (arb.u.)'); 
        set(cb,'YTick', lincmin:(lincmax-lincmin)/4:lincmax ,'YTickLabel',{'0',num2str(0.25*lincmax/shortlincmax),num2str(0.5*lincmax/shortlincmax),num2str(0.75*lincmax/shortlincmax),num2str(lincmax/shortlincmax)}) 
        
        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 

        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'spectra_short+long-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'spectra_short+long-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'spectra_short+long-traj_',savename(1:end-4),'.pdf'])
        pause(2); close  


%%
figure
sp=[];
   sp(1) = subplot(4,1,1);
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectralong1/max(pltsumspectralong1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectralong2/max(pltsumspectralong2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
        
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2,log10(pltspectralong));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2,log10(pltspectralong));
    set(gca,'ydir','normal');
         set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        %caxis([log10(max(max(spectralong)))-4-.5, log10(max(max(spectralong)))-.5]) 
        caxis([logcmin logcmax]);
        ylim([pltminEphot pltmaxEphot]);
        set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','log_{10}[ \omega^4 |d(\omega)|^2 ]    (arb.u.)'); 
        set(cb,'YTick', logcmin:(logcmax-logcmin)/4:logcmax)
        
        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 
        
        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'log-spectra_long-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'log-spectra_long-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'log-spectra_long-traj_',savename(1:end-4),'.pdf'])
        pause(2); close

%%     
figure
sp=[];
   sp(1) = subplot(4,1,1); 
        plot(pltCEP,pltions/max(pltions),'Color','blue','Linewidth',1)
        hold on
        plot(pltCEP,pltsumspectralong1/max(pltsumspectralong1),'Color',[.9,0.1,0],'Linewidth',1)
        plot(pltCEP,pltsumspectralong2/max(pltsumspectralong2),'Color',[1,0.65,0],'Linewidth',1)
        ylim([0 1.1])
         set(gca,'YTick',[0 1])
        ylabel('Ions / HHG','FontSize',11);
        set(gca,'FontSize',11);
     
    sp(2)=subplot(4,1,[2 4]);
    %surface(pltCEP,pltEphot*27.2, 1e6*(pltspectralong));
    %view(2)
    imagesc(pltCEP,pltEphot*27.2, 1e6*(pltspectralong));
    set(gca,'ydir','normal');
         set(gca,'FontSize',11);
        xlabel('CEP (\pi rad)')
        ylabel('Photon energy (eV)')
        shading interp
        axis tight
        load('MyGermanColormap','mycmap')
        set(gcf,'Colormap',mycmap)
        %caxis([0, 0.75*1e6*(max(max(spectralong)))]) 
        caxis([lincmin lincmax]);
        ylim([pltminEphot pltmaxEphot]);
     set(gcf,'Renderer','painters',...
          'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperSize',[11.5 8.5],'PaperPosition',[0.25 0.25 11 8]);
        set(gcf,'DefaultAxesFontSize',11);
        set(gca,'FontSize',11);
        box on
        cb = colorbar('vert');  zlab = get(cb,'ylabel'); 
        set(zlab,'String','\omega^4 |d(\omega)|^2   (arb.u.)'); 
        set(cb,'YTick', lincmin:(lincmax-lincmin)/4:lincmax ,'YTickLabel',{'0',num2str(0.25*lincmax/shortlincmax),num2str(0.5*lincmax/shortlincmax),num2str(0.75*lincmax/shortlincmax),num2str(lincmax/shortlincmax)}) 

        sp1pos = get(sp(1), 'Position');
        sp2pos = get(sp(2), 'Position');
           set(sp(1), 'Position', [sp1pos(1)-0.02 sp1pos(2)+0.02 sp2pos(3)-0.05 sp1pos(4)]) 
           set(sp(2), 'Position', [sp2pos(1)-0.02 sp2pos(2)+0.02 sp2pos(3)-0.05 sp2pos(4)-0.02]) 
        
        if savefigs == 1
            saveas(gcf,[[savepath,'spectra\'],'spectra_long-traj_',savename(1:end-4),'.fig'],'fig')
            %print('-dpng','-r1200',[[savepath,'spectra\'],'spectra_long-traj_',savename(1:end-4),'.png'])
        end
        print('-dpdf','-painters',[[savepath,'spectra\'],'spectra_long-traj_',savename(1:end-4),'.pdf'])
        pause(2); close  
