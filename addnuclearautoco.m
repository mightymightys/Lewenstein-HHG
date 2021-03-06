
if length(size(results))==3
    [k numevents m] = size(results);
elseif length(size(results))==2
    [k numevents] = size(results);
    m=1;
end

for nCEP=1:length(CEP)    
    for n=1:numevents
        if      (isempty(results{1,n,nCEP})) || (isempty(results{2,n,nCEP}))
        elseif  (max(isnan(results{1,n,nCEP}(:,6)))==1) || (max(isnan(results{2,n,nCEP}(:,6)))==1)
                disp(['Dropping event no.',num2str(n),' for CEP=',num2str(CEP(nCEP)),'pi.'])
        else
        
        %short
            emampl = results{1,n,nCEP}(:,6);
            exctau = results{1,n,nCEP}(:,5);

        %     % model the nuclear autocorrelation funtion as a cos^2 with the nucl.
        %     % oscillation period and with a phase that grows lin. with excusrion time
        %     nuclperiod = 8; %period of nuclear oscillation in a.u. (e.g. for NO2 it's about 8fs)
        %     phaseprop= 0.1; %proportionality factor for phase in rad/fs

            emampl = emampl .* cos(pi/(nuclperiod*1000 /24.2) *exctau).^2 .*exp(1i*phaseprop/1000*24.2 *exctau);

            results{1,n,nCEP}(:,6) = emampl;
        %long
            emampl = results{2,n,nCEP}(:,6);
            exctau = results{2,n,nCEP}(:,5);

            emampl = emampl .* cos(pi/(nuclperiod*1000 /24.2) *exctau).^2 .*exp(1i*phaseprop/1000*24.2 *exctau);

            results{2,n,nCEP}(:,6) = emampl;
        end
    
    end
end