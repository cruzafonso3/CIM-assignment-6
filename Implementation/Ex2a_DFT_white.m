%==========================================================================
%
% Ex2a_DFT_white.m
%
% Exercício 2 - Variante I: Quantização com RUÍDO BRANCO no domínio da DFT
%
% Uniform quantization step is the SAME for all spectral coefficients
% (real and imaginary parts) within a frame.
% The step is derived from the desired SMR (Signal-to-Mask Ratio)
% applied to the average power of the frame (sum of squared magnitudes
% of all spectral coefficients).
%
% M.EEC045 - Codificação de Informação Multimédia - FEUP 2025/2026
%==========================================================================

close all;
clear all;

% -----------------------------------------------------------------------
% Configuration
% -----------------------------------------------------------------------
SMR_dB = 10.0;    % <-- desired SMR in dB (try 5, 10, 15, 20)
                  %     SMR ≈ SNR after processing (iterate if needed)
plots  = 1;       % 1 = show plots, 0 = skip (faster batch runs)

inpfile = 'sound.wav';
outfile = 'out_DFT_white.wav';

% -----------------------------------------------------------------------
% Read WAV
% -----------------------------------------------------------------------
[datain, FS] = audioread(inpfile);   % normalized to [-1, 1]
datain = datain(:,1);                % mono
nread  = length(datain);
fprintf('Sampling frequency: %d Hz\n', FS);

% -----------------------------------------------------------------------
% Parameters
% -----------------------------------------------------------------------
N  = 1024;
N2 = N / 2;
% Analysis/synthesis sine window (as in original code)
win = sin(pi/N * (0:N-1).');

% Output accumulator
dataout = zeros(size(datain));

% -----------------------------------------------------------------------
% Convert SMR_dB to linear ratio
%   Quantisation step Δ for a uniform mid-rise quantiser:
%   noise power ≈ Δ²/12
%   signal power ≈ P_frame (sum |X[k]|²)
%   SMR = P_signal / P_noise  =>  Δ = sqrt(12 * P_frame / SMR_linear)
%
%   BUT here P_frame is the sum of |X[k]|² over all N coefficients.
%   The real & imaginary parts of X[k] each have variance P_frame / N.
%   We quantise real and imaginary parts with the same step Δ.
% -----------------------------------------------------------------------
SMR_lin = 10^(SMR_dB / 10);

% -----------------------------------------------------------------------
% Overlap-add loop
% -----------------------------------------------------------------------
regiaofreq = 1:N2+1;   % positive-frequency indices for display

% Storage for one illustrative frame (saved for spectrum plot)
saved_signal_spectrum = [];
saved_noise_spectrum  = [];
saved_frame_done      = false;
SAVE_FRAME = 5;        % save this frame index for illustration

frame = 1;
while (frame+1)*N2 < nread

    % --- extract segment ---
    varre   = 1 + (frame-1)*N2 : (frame+1)*N2;
    tmpdata = datain(varre);

    if plots
        figure(1);
        plot(0:N-1, tmpdata);
        xlabel('Samples'); ylabel('Amplitude');
        title(sprintf('Frame %d – time signal (before window)', frame));
        drawnow;
    end

    % --- apply window ---
    tmpdata = tmpdata .* win;

    % --- FFT ---
    fdata     = fft(tmpdata);
    magnitude = abs(fdata);
    phase_sig = angle(fdata);

    % -----------------------------------------------------------
    % Variant I: WHITE NOISE quantization
    %   Single quantisation step for all spectral coefficients
    % -----------------------------------------------------------

    % Frame power = sum of squared magnitudes of all spectral bins
    P_frame = sum(magnitude.^2);

    % Noise power budget for desired SMR
    P_noise_budget = P_frame / SMR_lin;

    % Each spectral sample (real + imag) contributes Δ²/12 to noise.
    % Total noise power = N * 2 * (Δ²/12)  [N bins, 2 components each]
    % => P_noise_budget = N * Δ²/6
    % => Δ = sqrt(6 * P_noise_budget / N)
    if P_frame > 0
        Delta = sqrt(6 * P_noise_budget / N);
    else
        Delta = 0;
    end

    % Quantise real and imaginary parts with step Delta
    if Delta > 0
        fdata_q = Delta * (floor(real(fdata)/Delta + 0.5) + ...
                       1i * floor(imag(fdata)/Delta + 0.5));
    else
        fdata_q = fdata;   % zero-power frame: pass through
    end

    % -----------------------------------------------------------
    % Compute noise spectrum for this frame (for illustration)
    % -----------------------------------------------------------
    noise_spec = fdata_q - fdata;

    % Save one representative frame
    if frame == SAVE_FRAME && ~saved_frame_done
        saved_signal_spectrum = fdata;
        saved_noise_spectrum  = noise_spec;
        saved_frame_done      = true;
    end

    % -----------------------------------------------------------
    % Show spectra before / after quantisation
    % -----------------------------------------------------------
    if plots
        figure(2);
        plot(FS/N*(regiaofreq-1), 20*log10(eps + abs(fdata(regiaofreq))), 'b');
        hold on;
        plot(FS/N*(regiaofreq-1), 20*log10(eps + abs(fdata_q(regiaofreq))), 'r--');
        hold off;
        xlabel('Frequency (Hz)'); ylabel('Spectral density (dB)');
        title(sprintf('Frame %d – Blue: original | Red: quantised', frame));
        legend('Original','Quantised');
        drawnow;
    end

    % -----------------------------------------------------------
    % Enforce Hermitian symmetry for real IFFT output
    % -----------------------------------------------------------
    fdata_q(N:-1:N2+2) = conj(fdata_q(2:N2));

    % --- IFFT ---
    tmpdata = real(ifft(fdata_q));

    % --- apply window ---
    tmpdata = tmpdata .* win;

    % --- overlap-add ---
    dataout(varre) = dataout(varre) + tmpdata;

    if plots
        figure(1); hold on;
        plot(0:N-1, tmpdata, 'r');
        hold off;
        legend('Original','Reconstructed');
        drawnow;
    end

    frame = frame + 1;
end

% -----------------------------------------------------------------------
% Write output WAV
% -----------------------------------------------------------------------
% Clip to [-1, 1] to avoid audiowrite warnings
dataout = max(-1, min(1, dataout));
audiowrite(outfile, dataout, FS);

% -----------------------------------------------------------------------
% SNR calculation (time domain)
% -----------------------------------------------------------------------
neval      = min(length(datain), length(dataout));
difsignal  = datain(1:neval) - dataout(1:neval);
Psignal    = sum(datain(1:neval).^2);
Pnoise     = sum(difsignal.^2);

if Pnoise == 0
    fprintf('\nSNR = infinity (perfect reconstruction)\n');
else
    SNR_measured = 10*log10(Psignal / Pnoise);
    fprintf('\n=== DFT White Noise Quantization ===\n');
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
    title(sprintf('DFT White Noise – Frame %d Spectra (SMR=%.0f dB, SNR≈%.1f dB)', ...
                  SAVE_FRAME, SMR_dB, SNR_measured));
    legend('Signal spectrum','Noise spectrum');
    grid on;
end

fprintf('\nOutput written to: %s\n', outfile);
fprintf('Done.\n');
