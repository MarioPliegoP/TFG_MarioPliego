%% SIMULATION

out=sim('ADC_model_jitter_metaestab.slx');
assignin('base','out',out)

Xtim= FS/(2^B-1) * out.do.signals.values; 


%% USE OF SPECT_METRICS
flag_plot=1;
flag_meas=1;

% Xtim:         Temporal sequence
% fs:           Sampling frequency (delta_t=1/fs)
% flag_plot:    Flag to control spectrum plot -> PSD [dB] vs. f [Hz]
%               0 - Don't plot / 1 - Do plot
% flag_meas:    Flag to control computation of metrics from spectrum
%               Psignal, SNR, SNDR, HDs, THD, ENOB
%               0 - Don't compute / 1 - Do compute
% fin_aprox:    Input tone frequency (for input signal identification)
% span_bins:    Number of bins around central bin for power computation

N=2^floor(log2(length(Xtim))); % Vector length. Number of points rounded to power of 2.

Xtim=Xtim(end-N+1:end);
% time=time(end-N+1:end);
% freq=freq(end-N+1:end);

% 
W= 	hanning(N);

xw=W.*Xtim;
pxx=abs(fft(xw/sqrt(N))).^2;

pxx=pxx/norm(W)^2; % Factor 2 for the single-side representation
psd=2*pxx(1:ceil(N/2));	% Factor 2 for the single-side representation

fbin=fs/N;
f=(0:fbin:fs); f=f(1:ceil(N/2));

if flag_meas==1
    %% Computation of Psig, Pnoise, PHD, ...
    bin_sig=find((fin<f),1);
    sig_span=bin_sig-span_bins:1:bin_sig+span_bins;
    Psig=sum(psd(sig_span));
    dc_span=1:span_bins;
    Pdc=sum(psd(dc_span));
    nhar=floor(fs/2/fin); % Number of harmonics laying in the band
    Phar=zeros(1,nhar-1);
    for i=2:nhar
        bin_hari=find((i*fin<f),1);
        hari_span=bin_hari-span_bins:1:bin_hari+span_bins;
        Phar(i-1)=sum(psd(hari_span));        
    end
    Phartot=sum(Phar);

    Pnoise=sum(psd)-Psig-Pdc-Phartot;
    SNR=10*log10(Psig/Pnoise);
    SNDR=10*log10(Psig/(Pnoise+Phartot));
    ENOB=(SNDR-1.76)/6.02;
end
 