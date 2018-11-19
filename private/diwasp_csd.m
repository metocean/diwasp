function [S, f] = diwasp_csd(x,y,nfft,fs)
%Diwasp cross spectral density.
%If CPSD is available from Matlab Signal Processing Toolbox, then use that.
%Otherwise this function will calc cross-spectra density using inbuild FFT function
%
%[Pxy, f] = diwasp_csd(x,y,nfft,fs)
%
%

%Make a windowed estimate of CSD
hann=0.5*(1-cos(2*pi*(1:nfft/2)/(nfft+1)));
win = [hann hann(end:-1:1)];
    
if exist('cpsd')==2
    [S, f] = cpsd(y,x,win,0,nfft,fs);
else
    nw = length(win);
    nseg=fix(length(x)/nw);
    S = zeros(nfft,1);
    for iseg=0:nseg-1
        ind=nw*iseg+[1:nw];
        xw = win'.*x(ind);
        yw = win'.*y(ind);
        Px = fft(xw,nfft);
        Py = fft(yw,nfft);
        Pxy = Py.*conj(Px);
        S = S + Pxy;
    end
    nfac=(fs*nseg*norm(win)^2);
    S=[S(1); 2*S(2:nfft/2); S(nfft/2+1)]/nfac;
    f=(fs/nfft)*[0:nfft/2]';
end
