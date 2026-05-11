%==========================================================================
%
% Ex2b_DFT_colored.m
%
% Exercício 2 - Variante II: Quantização com RUÍDO COLORIDO no domínio da DFT
%
% The quantisation step is INDIVIDUALISED for each spectral coefficient
% so that the Signal-to-Mask Ratio is UNIFORM across the spectrum.
%
% Per bin:  SMR = |X[k]|² / noise_power_k  = constant = SMR_lin
%   => noise_power_k = |X[k]|² / SMR_lin
%   => Δ_k = sqrt(12 * noise_power_k) = |X[k]| * sqrt(12 / SMR_lin)
%
% The quantisation noise spectrum thus follows (is shaped like) the
% signal spectrum — hence "colored" noise.
%
% M.EEC045 - Codificação de Informação Multimédia - FEUP 2025/2026
%==========================================================================

close all;
clear all;

% -----------------------------------------------------------------------
% Configuration
% -----------------------------------------------------------------------
SMR_dB = 10.0;    % <-- desired SMR in dB (try 5, 10, 15, 20)
plots  = 1;       % 1 = show plots, 0 = skip

inpfile = 'sound.wav';
outfile = 'out_DFT_colored.wav';

% -----------------------------------------------------------------------
% Read WAV
% -----------------------------------------------------------------------
[datain, FS] = audioread(inpfile);
datain = datain(:,1);
nread  = length(datain);
fprintf('Sampling frequency: %d Hz\n', FS);

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
N  = 1024;
N2 = N / 2;
win = sin(pi/N * (0:N-1).');

dataout = zeros(size(datain));

SMR_lin = 10^(SMR_dB / 10);

% Pre-compute the per-bin step scale factor
%   Δ_k = |X[k]| * sqrt(12 / SMR_lin)
step_scale = sqrt(12 / SMR_lin);

% -----------------------------------------------------------------------
% Overlap-add loop
% -----------------------------------------------------------------------
regiaofreq = 1:N2+1;

saved_signal_spectrum = [];
saved_noise_spectrum  = [];
saved_frame_done      = false;
SAVE_FRAME = 5;

frame = 1;
while (frame+1)*N2 < nread

    varre   = 1 + (frame-1)*N2 : (frame+1)*N2;
    tmpdata = datain(varre);

    if plots
        figure(1);
        plot(0:N-1, tmpdata);
        xlabel('Samples'); ylabel('Amplitude');
        title(sprintf('Frame %d – time signal (before window)', frame));
        drawnow;
    end

    tmpdata = tmpdata .* win;
    fdata   = fft(tmpdata);
    magnitude = abs(fdata);

    % -----------------------------------------------------------
    % Variant II: COLORED NOISE quantization
    %   Per-bin step  Δ_k = |X[k]| * step_scale
    %   (gives uniform SMR across all bins)
    % -----------------------------------------------------------

    % Quantise each bin's real and imaginary part independently
    % with its own step size
    Delta_k  = magnitude * step_scale;   % length-N vector
    Delta_k(magnitude == 0) = 0;         % avoid NaN for zero bins

    fdata_q = zeros(N, 1);
    for k = 1:N
        if Delta_k(k) > 0
            fdata_q(k) = Delta_k(k) * (floor(real(fdata(k))/Delta_k(k) + 0.5) + ...
                              1i * floor(imag(fdata(k))/Delta_k(k) + 0.5));
        else
            fdata_q(k) = 0;   % zero-magnitude bin
        end
    end

    % -----------------------------------------------------------
    % Noise spectrum
    % -----------------------------------------------------------
    noise_spec = fdata_q - fdata;

    if frame == SAVE_FRAME && ~saved_frame_done
        saved_signal_spectrum = fdata;
        saved_noise_spectrum  = noise_spec;
        saved_frame_done      = true;
    end

    if plots
        figure(2);
        plot(FS/N*(regiaofreq-1), 20*log10(eps + abs(fdata(regiaofreq))), 'b');
        hold on;
        plot(FS/N*(regiaofreq-1), 20*log10(eps + abs(fdata_q(regiaofreq))), 'r--');
        hold off;
        xlabel('Frequency (Hz)'); ylabel('Spectral density (dB)');
        title(sprintf('Frame %d – Blue: original | Red: quantised (colored)', frame));
        legend('Original','Quantised');
        drawnow;
    end

    % Enforce Hermitian symmetry
    fdata_q(N:-1:N2+2) = conj(fdata_q(2:N2));

    tmpdata = real(ifft(fdata_q));
    tmpdata = tmpdata .* win;
    dataout(varre) = dataout(varre) + tmpdata;

    if plots
        figure(1); hold on;
        plot(0:N-1, tmpdata, 'r');
        hold off;
        drawnow;
    end

    frame = frame + 1;
end

% -----------------------------------------------------------------------
% Write output WAV
% -----------------------------------------------------------------------
dataout = max(-1, min(1, dataout));
audiowrite(outfile, dataout, FS);

% -----------------------------------------------------------------------
% SNR calculation
% -----------------------------------------------------------------------
neval     = min(length(datain), length(dataout));
difsignal = datain(1:neval) - dataout(1:neval);
Psignal   = sum(datain(1:neval).^2);
Pnoise    = sum(difsignal.^2);

if Pnoise == 0
    fprintf('\nSNR = infinity (perfect reconstruction)\n');
    SNR_measured = Inf;
else
    SNR_measured = 10*log10(Psignal / Pnoise);
    fprintf('\n=== DFT Colored Noise Quantization ===\n');
    fprintf('Target SMR : %.1f dB\n', SMR_dB);
    fprintf('Measured SNR: %.2f dB\n', SNR_measured);
end

% -----------------------------------------------------------------------
% Spectrum illustration: signal vs. noise (frame SAVE_FRAME)
% -----------------------------------------------------------------------
if ~isempty(saved_signal_spectrum)
    figure(3);
    freqs = FS/N * (regiaofreq - 1);
    plot(freqs, 20*log10(eps + abs(saved_signal_spectrum(regiaofreq))), 'b', 'LineWidth', 1.2);
    hold on;
    plot(freqs, 20*log10(eps + abs(saved_noise_spectrum(regiaofreq))),  'r', 'LineWidth', 1.2);
    hold off;
    xlabel('Frequency (Hz)'); ylabel('Spectral Density (dB)');
    title(sprintf('DFT Colored Noise – Frame %d Spectra (SMR=%.0f dB, SNR≈%.1f dB)', ...
                  SAVE_FRAME, SMR_dB, SNR_measured));
    legend('Signal spectrum','Noise spectrum');
    grid on;
end

fprintf('\nOutput written to: %s\n', outfile);
fprintf('Done.\n');
