%==========================================================================
%
% Ex3c_MDCT_colored.m
%
% Exercício 3 - Variante II: Quantização com RUÍDO COLORIDO no domínio da MDCT
%
% Based on the provided qzmdct.m skeleton.
%
% The quantisation step is INDIVIDUALISED for each MDCT coefficient
% so that the SMR is UNIFORM across the spectrum (colored noise shaped
% to follow the signal spectrum).
%
%   Per coefficient k:  SMR = mdct[k]² / noise_k = SMR_lin (constant)
%   => noise_k = mdct[k]² / SMR_lin
%   => Δ_k = sqrt(12 * mdct[k]² / SMR_lin) = |mdct[k]| * sqrt(12/SMR_lin)
%
% M.EEC045 - Codificação de Informação Multimédia - FEUP 2025/2026
%==========================================================================

close all;
clear all;

% -----------------------------------------------------------------------
% Configuration
% -----------------------------------------------------------------------
SMR_dB = 10.0;    % <-- desired SMR in dB (try 5, 10, 15, 20)
plots  = 0;       % 1 = frame-by-frame plots, 0 = batch

inpfile = 'sound.wav';
outfile = 'out_MDCT_colored.wav';

% -----------------------------------------------------------------------
% Read WAV and write temp PCM
% -----------------------------------------------------------------------
[datain, FS] = audioread(inpfile);
datain = datain(:,1);
nread  = length(datain);
fprintf('Sampling frequency: %d Hz\n', FS);

bits = 16;
fid = fopen('tmpinpaudio.pcm','w');
fwrite(fid, datain*(2^(bits-1)-1), 'short');
fclose(fid);

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
N  = 1024;
N2 = N / 2;
win = sin(pi/N * ([0:N-1] + 0.5));

idata   = zeros(1, N);
odata   = zeros(1, N);
osdata  = zeros(1, N2);
tmpdata = zeros(1, N2);
fdata   = zeros(1, N);
ofdata  = zeros(1, N);

direxp = exp(-1i*pi*(0:N-1)/N);
invexp = exp( 1i*pi*(0:N-1)/N);
cosargmdct = cos(pi*(0.5+(0:N2-1))*(1.0+N2)/N);
sinargmdct = sin(pi*(0.5+(0:N2-1))*(1.0+N2)/N);

SMR_lin    = 10^(SMR_dB / 10);
step_scale = sqrt(12 / SMR_lin);   % Δ_k = |mdct[k]| * step_scale

% -----------------------------------------------------------------------
% Processing loop
% -----------------------------------------------------------------------
freqregion = 1:N2;

fidr = fopen('tmpinpaudio.pcm', 'r');
fidw = fopen('tmpoutaudio.pcm', 'w');

[data, nread_blk] = fread(fidr, N2, 'short');
idata(1, N2+1:N) = data(1:N2, 1).';

saved_signal_mdct = [];
saved_noise_mdct  = [];
saved_frame_done  = false;
SAVE_FRAME = 5;

k = 0;

while nread_blk == N2

    k = k + 1;

    if plots
        figure(1);
        plot(0:N-1, idata(1:N));
        xlabel('Samples'); ylabel('Amplitude');
        title(sprintf('Frame %d – time signal', k));
    end

    idata = idata .* win;
    fdata = idata .* direxp;
    odft  = fft(fdata);

    % --- ODFT -> MDCT ---
    mdct = real(odft(1:N2)) .* cosargmdct + imag(odft(1:N2)) .* sinargmdct;

    % =============================================================
    % Variant II: COLORED NOISE quantization of MDCT coefficients
    % =============================================================
    % Per-coefficient step proportional to |mdct[k]|
    Delta_k  = abs(mdct) * step_scale;   % vector of length N2
    mdct_q   = zeros(1, N2);

    for ki = 1:N2
        if Delta_k(ki) > 0
            mdct_q(ki) = Delta_k(ki) * floor(mdct(ki)/Delta_k(ki) + 0.5);
        else
            mdct_q(ki) = 0;   % zero coefficient stays zero
        end
    end

    % --- Save one frame for illustration ---
    if k == SAVE_FRAME && ~saved_frame_done
        saved_signal_mdct = mdct;
        saved_noise_mdct  = mdct_q - mdct;
        saved_frame_done  = true;
    end

    % =============================================================
    % IMDCT -> IODFT
    % =============================================================
    fdata(1:N2) = 2 * (mdct_q .* cosargmdct + 1i * mdct_q .* sinargmdct);

    fdata(N:-1:N2+1) = conj(fdata(1:N2));
    ofdata = ifft(fdata);
    ofdata = ofdata .* invexp;
    odata  = real(ofdata);
    odata  = odata .* win;

    tmpdata(1:N2) = floor(0.5 + osdata + odata(1:N2));
    osdata = odata(N2+1:N);
    fwrite(fidw, tmpdata(1:N2), 'short');

    idata(1, 1:N2) = data(1:N2, 1).';
    [data, nread_blk] = fread(fidr, N2, 'short');
    if nread_blk < N2
        data(nread_blk+1:N2, 1) = zeros(N2-nread_blk, 1);
    end
    idata(1, N2+1:N) = data(1:N2, 1).';

    if plots; pause; end
end

fclose(fidr);
fclose(fidw);

% -----------------------------------------------------------------------
% Convert output PCM to WAV
% -----------------------------------------------------------------------
fidr = fopen('tmpinpaudio.pcm','r');
fidw = fopen('tmpoutaudio.pcm','r');
[data1, ~] = fread(fidr, 'short');
[data2, ~] = fread(fidw, 'short');
fclose(fidr);
fclose(fidw);

audiowrite(outfile, data2/(2^(bits-1)-1), FS);

% -----------------------------------------------------------------------
% SNR calculation (with system delay alignment)
% -----------------------------------------------------------------------
shift = N2;
neval = length(data2) - shift;
if neval > length(data1); neval = length(data1); end

difsignal = data1(1:neval) - data2(1+shift:shift+neval);
Psignal   = sum(data1(1:neval).^2);
Pnoise    = sum(difsignal.^2);

fprintf('\n=== MDCT Colored Noise Quantization ===\n');
if Pnoise == 0
    fprintf('SNR = infinity (perfect reconstruction)\n');
    SNR_measured = Inf;
else
    SNR_measured = 10*log10(Psignal / Pnoise);
    fprintf('Target SMR : %.1f dB\n', SMR_dB);
    fprintf('Measured SNR: %.2f dB\n', SNR_measured);
end

% -----------------------------------------------------------------------
% Coding noise time-domain plot
% -----------------------------------------------------------------------
figure(10);
plot((0:neval-1)*1000/FS, difsignal);
xlabel('Time (ms)'); ylabel('Amplitude');
title('MDCT Colored Noise – Coding noise (time domain)');
grid on;

% -----------------------------------------------------------------------
% Spectrum illustration
% -----------------------------------------------------------------------
if ~isempty(saved_signal_mdct)
    figure(11);
    freqs = FS/N * (0:N2-1);
    plot(freqs, 20*log10(eps + abs(saved_signal_mdct)), 'b', 'LineWidth', 1.2);
    hold on;
    plot(freqs, 20*log10(eps + abs(saved_noise_mdct)),  'r', 'LineWidth', 1.2);
    hold off;
    xlabel('Frequency (Hz)'); ylabel('Spectral Density (dB)');
    title(sprintf('MDCT Colored Noise – Frame %d Spectra (SMR=%.0f dB, SNR≈%.1f dB)', ...
                  SAVE_FRAME, SMR_dB, SNR_measured));
    legend('Signal MDCT','Noise MDCT');
    grid on;
end

fprintf('\nOutput written to: %s\n', outfile);
fprintf('Done.\n');
