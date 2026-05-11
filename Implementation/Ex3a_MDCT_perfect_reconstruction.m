%==========================================================================
%
% Ex3a_MDCT_perfect_reconstruction.m
%
% Exercício 3 - Verificação da RECONSTRUÇÃO PERFEITA da cadeia MDCT
%
% Runs qzmdct.m processing chain WITHOUT any quantisation (SMR = Inf)
% and confirms the output equals the input (SNR = infinity).
%
% M.EEC045 - Codificação de Informação Multimédia - FEUP 2025/2026
%==========================================================================

close all;
clear all;

inpfile = 'sound.wav';
outfile = 'out_MDCT_perf_recon.wav';

[datain, FS] = audioread(inpfile);
datain = datain(:,1);
fprintf('Sampling frequency: %d Hz\n', FS);

bits = 16;
fid = fopen('tmpinpaudio.pcm','w');
fwrite(fid, datain*(2^(bits-1)-1), 'short');
fclose(fid);

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

fidr = fopen('tmpinpaudio.pcm','r');
fidw = fopen('tmpoutaudio.pcm','w');

[data, nread_blk] = fread(fidr, N2, 'short');
idata(1, N2+1:N) = data(1:N2, 1).';

k = 0;

while nread_blk == N2
    k = k + 1;
    idata = idata .* win;
    fdata = idata .* direxp;
    odft  = fft(fdata);

    % MDCT
    mdct = real(odft(1:N2)) .* cosargmdct + imag(odft(1:N2)) .* sinargmdct;

    % NO quantisation  <-- this is the perfect reconstruction test
    mdct_q = mdct;

    % IMDCT
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
end

fclose(fidr);
fclose(fidw);

fidr = fopen('tmpinpaudio.pcm','r');
fidw = fopen('tmpoutaudio.pcm','r');
[data1, ~] = fread(fidr, 'short');
[data2, ~] = fread(fidw, 'short');
fclose(fidr);
fclose(fidw);

audiowrite(outfile, data2/(2^(bits-1)-1), FS);

shift = N2;
neval = length(data2) - shift;
if neval > length(data1); neval = length(data1); end

difsignal = data1(1:neval) - data2(1+shift:shift+neval);
Psignal   = sum(data1(1:neval).^2);
Pnoise    = sum(difsignal.^2);

figure(1);
plot((0:neval-1)*1000/FS, difsignal);
xlabel('Time (ms)'); ylabel('Amplitude');
title('MDCT Perfect Reconstruction Test – Coding noise');
grid on;

fprintf('\n=== MDCT Perfect Reconstruction Test ===\n');
if Pnoise == 0
    fprintf('Result: SNR = INFINITY  -->  PERFECT RECONSTRUCTION CONFIRMED!\n');
else
    SNR = 10*log10(Psignal/Pnoise);
    fprintf('Result: SNR = %.2f dB\n', SNR);
    fprintf('(Non-zero noise is due to integer rounding in PCM write/read,\n');
    fprintf(' which is normal and expected at high SNR.)\n');
end
fprintf('Done.\n');
