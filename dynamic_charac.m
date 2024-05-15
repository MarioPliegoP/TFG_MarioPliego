%% DYNAMIC CHARACTERIZATION
% Mario Pliego Padilla 08/02/2024.
%% DESCRIPTION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script computes the frequency spectrum of a sinusoidal input.  
% It extracts THD, SNR, SNDR, ENOB, Psig, Pnoise, PHD.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SIMULATION
tic
out=sim('ADC_model_jitter_metaestab_v2.slx');
toc
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
if app.Windowfunction.Value== "Hanning"
    W= 	hanning(N);
elseif app.Windowfunction.Value == "Blackman-Harris"
    W = blackmanharris(N);
elseif app.Windowfunction.Value == "Bartlett"
    W = bartlett(N);
elseif app.Windowfunction.Value == "Chebyshev"
    W = chebwin(N);
elseif app.Windowfunction.Value == "Hamming"
    W = hamming(N);
elseif app.Windowfunction.Value == "Rectangular"
    W = rectwin(N);
end

xw=W.*Xtim;
pxx=abs(fft(xw/sqrt(N))).^2;

pxx=pxx/norm(W)^2; % Factor 2 for the single-side representation
psd=2*pxx(1:ceil(N/2));	% Factor 2 for the single-side representation

fbin=fs/N;
f=(0:fbin:fs); f=f(1:ceil(N/2));

if flag_plot==1
    plot(app.spectfig,f,10*log10(psd)); 
    grid(app.spectfig,"on"); 
    hold(app.spectfig,"on")
    xlabel(app.spectfig,'Frequency [Hz]'); 
    ylabel(app.spectfig,'PSD [dB]');
    axis(app.spectfig,'tight')
end

if flag_meas==1
    %% Computation of Psig, Pnoise, PHD, ...
    bin_sig=find((fin<f),1);
    sig_span=bin_sig-span_bins:1:bin_sig+span_bins;
    Psig=sum(psd(sig_span));
    if flag_plot==1
        plot(app.spectfig,f(sig_span),10*log10(psd(sig_span)),'k');
        text(app.spectfig,f(bin_sig),10*log10(psd(bin_sig)),' SIG=H1');
    end
    dc_span=1:span_bins;
    Pdc=sum(psd(dc_span));
    if flag_plot==1
        plot(app.spectfig,f(dc_span),10*log10(psd(dc_span)),'r');
        text(app.spectfig,f(1),10*log10(psd(1)),' DC');
    end
    nhar=min(floor(fs/2/fin),5); % Number of harmonics laying in the band
    Phar=zeros(1,nhar-1);
    for i=2:nhar
        bin_hari=find((i*fin<f),1);
        hari_span=bin_hari-span_bins:1:bin_hari+span_bins;
        Phar(i-1)=sum(psd(hari_span));
        if flag_plot==1
            plot(app.spectfig,f(hari_span),10*log10(psd(hari_span)),'r');
            text(app.spectfig,f(bin_hari),10*log10(psd(bin_hari)),sprintf(' H%i',i));
        end
    end
    Phartot=sum(Phar);

    Pnoise=sum(psd)-Psig-Pdc-Phartot;
    SNR=10*log10(Psig/Pnoise);
    SNDR=10*log10(Psig/(Pnoise+Phartot));
    ENOB=(SNDR-1.76)/6.02;

    app.text1.Text=sprintf('\n SIGNAL:    Psig = %.2f dB / A1 = %.2f V = %.2e Vrms', 10*log10(Psig), sqrt(2*Psig), sqrt(Psig));
    app.text2.Text=sprintf('\n HARMONICS: Phar = %.2f dB / THD = %.2f dB', 10*log10(Phartot), 10*log10(Phartot/Psig));
    % for i=2:nhar
    %     fprintf('\n            HD%i = %.2f dBc / A%i = %.2f V = %.2e Vrms', i, 10*log10(Phar(i-1)/Psig),i,sqrt(2*Phar(i-1)),sqrt(Phar(i-1)));
    % end
    if length(Phar) >= 1
        app.text3.Text=sprintf('\n            HD2 = %.2f dBc / A2 = %.2f V = %.2e Vrms', 10*log10(Phar(2-1)/Psig),sqrt(2*Phar(2-1)),sqrt(Phar(2-1)));
    else
        app.text3.Text= '';
    end
    if length(Phar) >= 2
        app.text4.Text=sprintf('\n            HD3 = %.2f dBc / A3 = %.2f V = %.2e Vrms', 10*log10(Phar(3-1)/Psig),sqrt(2*Phar(3-1)),sqrt(Phar(3-1)));
    else
        app.text4.Text = '';
    end
    if length(Phar) >=3
        app.text5.Text=sprintf('\n            HD4 = %.2f dBc / A4 = %.2f V = %.2e Vrms', 10*log10(Phar(4-1)/Psig),sqrt(2*Phar(4-1)),sqrt(Phar(4-1)));
    else
        app.text5.Text='';
    end
    if length(Phar) >=4
    app.text6.Text=sprintf('\n            HD5 = %.2f dBc / A5 = %.2f V = %.2e Vrms', 10*log10(Phar(5-1)/Psig),sqrt(2*Phar(5-1)),sqrt(Phar(5-1)));
    else
    app.text6.Text='';
    end
    
    app.text7.Text=sprintf('\n NOISE:     Pnoise  = %.2f dB = %.2e Vrms \n', 10*log10(Pnoise),sqrt(Pnoise));

    app.text8.Text=sprintf('\n SNR  = %.2f dB', 10*log10(Psig/Pnoise));
    app.text9.Text=sprintf('\n SNDR = %.2f dB', 10*log10(Psig/(Pnoise+Phartot)));
    app.text10.Text=sprintf('\n ENOB = %.2f bit \n', ENOB);

    if flag_plot==1
        title(app.spectfig,sprintf('SPECTRUM (2^{%i}-point FFT)',log2(N)));                                                                                                                                                      % Cambiado el +1
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











%% -------------------------------------------------- END  -------------------------------------------------- %%