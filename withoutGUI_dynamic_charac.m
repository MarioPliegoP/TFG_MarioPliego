%% DYNAMIC CHARACTERIZATION
% Mario Pliego Padilla 08/02/2024.
%% DESCRIPTION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes the frequency spectrum of a sinusoidal input.  
% It extracts THD, SNR, SNDR, ENOB, Psig, Pnoise, PHD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIMULATION

out=sim('ADC_model_jitter_metaestab.slx');

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

if flag_plot==1
    plot(f,10*log10(psd)); 
    grid("on"); 
    hold("on")
    xlabel('Frequency [Hz]'); 
    ylabel('PSD [dB]');
    axis('tight')
end

if flag_meas==1
    %% Computation of Psig, Pnoise, PHD, ...
    bin_sig=find((fin<f),1);
    sig_span=bin_sig-span_bins:1:bin_sig+span_bins;
    Psig=sum(psd(sig_span));
    if flag_plot==1
        plot(f(sig_span),10*log10(psd(sig_span)),'k');
        text(f(bin_sig),10*log10(psd(bin_sig)),' SIG=H1');
    end
    dc_span=1:span_bins;
    Pdc=sum(psd(dc_span));
    if flag_plot==1
        plot(f(dc_span),10*log10(psd(dc_span)),'r');
        text(f(1),10*log10(psd(1)),' DC');
    end
    nhar=min(floor(fs/2/fin),5); % Number of harmonics laying in the band
    Phar=zeros(1,nhar-1);
    for i=2:nhar
        bin_hari=find((i*fin<f),1);
        hari_span=bin_hari-span_bins:1:bin_hari+span_bins;
        Phar(i-1)=sum(psd(hari_span));
        if flag_plot==1
            plot(f(hari_span),10*log10(psd(hari_span)),'r');
            text(f(bin_hari),10*log10(psd(bin_hari)),sprintf(' H%i',i));
        end
    end
    Phartot=sum(Phar);

    Pnoise=sum(psd)-Psig-Pdc-Phartot;
    SNR=10*log10(Psig/Pnoise);
    SNDR=10*log10(Psig/(Pnoise+Phartot));
    ENOB=(SNDR-1.76)/6.02;

 
    % for i=2:nhar
    %     fprintf('\n            HD%i = %.2f dBc / A%i = %.2f V = %.2e Vrms', i, 10*log10(Phar(i-1)/Psig),i,sqrt(2*Phar(i-1)),sqrt(Phar(i-1)));
    % end

    if flag_plot==1
        title(sprintf('SPECTRUM (2^{%i}-point FFT)',log2(N)));                                                                                                                                                      % Cambiado el +1
    end
end


%% Representation without spect_metrics

% % t=1/fs*(0:(N-1))';

% t=1/fs*(0:(length(Xtim)-1))';
% 
% figure; plot(t,Xtim); grid on; hold on;
% title('TEMPORAL SEQUENCE');
% xlabel('Time, s'); ylabel('X(t)');
% axis([0 t(end) 1.1*min(Xtim) 1.1*max(Xtim)]);
% shg

% % W= hanning(N);

% W = hanning(length(Xtim));
% xw=W.*Xtim;
% 
% figure; plot(t,xw); grid on; hold on;
% title('TEMPORAL SEQUENCE');
% xlabel('Time, s'); ylabel('X(t)');
% axis([0 t(end) 1.1*min(Xtim) 1.1*max(Xtim)]);
% shg





