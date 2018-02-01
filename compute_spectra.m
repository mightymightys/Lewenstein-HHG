if (~exist('waveperiod','var'))
    freqs=ffreqs(0,0);
    if length(freqs(:,1))>1
        freqssorted = sort(freqs(:,1));
        waveperiod = 1/(freqssorted(2)-freqssorted(1)) * 2*pi/omega;
    else
        waveperiod = 1/freqs(:,1) * 2*pi/omega;
    end
end

if length(size(results))==3
    [k numevents m] = size(results);
elseif length(size(results))==2
    [k numevents] = size(results);
    m=1;
end

if (~exist('dE','var'))
    dE=1;
end
if (~exist('E1','var'))
    E1=min(Ephot*27.2)+dE/2;
end

if (~exist('E2','var'))
    E2=max(Ephot*27.2)-dE/2;
end

if (~exist('deltaE','var'))
    deltaE=1;
end

Estp=Ephot(2)-Ephot(1);

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
            % directly interpolating in the complex plane likes this does
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
if E2indx<dE+1
    E2indx=1+dE;
end
if E1indx>length(Ephot)-dE
    E1indx=E1indx-dE;
end
if E2indx>length(Ephot)-dE
    E2indx=E2indx-dE;
end
%%
spectrashort=abs(Ephot.^4*ones(1,length(CEP)).*(spectamplshort).*conj(spectamplshort));
spectralong =abs(Ephot.^4*ones(1,length(CEP)).*(spectampllong).*conj(spectampllong));
spectraall  =abs(Ephot.^4*ones(1,length(CEP)).*(spectamplshort+spectampllong).*conj(spectamplshort+spectampllong));
spectramacro=abs(Ephot.^4*ones(1,length(CEP)).*(spectamplmacro).*conj(spectamplmacro));
%%
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
     ions(n)=sum( exp(-2*(2*Ip)^1.5 / 3 ./ abs(Efield(tt))));
end
