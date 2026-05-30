%---------------------------------------------------------------------------
%
% Illustration of the MDCT computation and quantizatization
% This code is property of ATC Labs. Use is allowed for education
% and research purposes only. For other purposes please contact
% ajf@atc-labs.com or ajf@fe.up.pt
%
% (c)2004-2026 ATC Labs, USA
%
% Anibal Ferreira, University of Porto, PORTUGAL / ATC Labs, USA
%
%--------------------------------------------------------------------------
close all;
clear all;

% sets the desired SNR (tentative only !),
% specifies if noise is white or colored
% and enables or disables plots
%
SNR=10.0; % dB
white=1;
plots=1;


%
% input and output audio files

inpfile = 'sound.wav';
outfile = 'outfile.wav';

%
% read WAV audio file
%
[datain,FS]=audioread(inpfile); % samples are normalized between 0 and 1
% [datain,FS,bits]=wavread(inpfile); % OLD
% disp('bits per sample: '); bits
nread=length(datain);
disp('Sampling frequency: '); FS

%
% write temp RAW audio file
%
bits = 16;
fid = fopen('tmpinpaudio.pcm','w');
fwrite(fid,datain*(2^(bits-1)-1),'short'); % max abs value is 2^(bits/2)-1
fclose(fid);

%
% N        = size of the ODFT and MDCT transforms
% win      = sine window
%
N=1024; N2=N/2;
win=sin(pi/N*([0:N-1]+0.5));


idata   = zeros(1,N);	% input (int) data
odata   = zeros(1,N);	% output (int) data
osdata  = zeros(1,N2);  % buffered data
tmpdata = zeros(1, N2);
fdata   = zeros(1,N);	% (float complex) data
ofdata  = zeros(1,N);	% (float complex) data

%
%  complex constants for the non-optimal computation of the ODFT, IODFT, MDCT and IMDCT transforms
%
direxp = exp(-1i*pi*[0:N-1]/N);
invexp = exp( 1i*pi*[0:N-1]/N);
cosargmdct=cos(pi*(0.5+[0:(N2-1)])*(1.0+N2)/N);
sinargmdct=sin(pi*(0.5+[0:(N2-1)])*(1.0+N2)/N);

%--------------------------------------------------------------
%
% read audio file
%
fidr = fopen('tmpinpaudio.pcm','r');
fidw = fopen('tmpoutaudio.pcm','w');

%
% begin overlap-add
%
[data, nread] = fread(fidr, N2, 'short');
idata(1,N2+1:N)=data(1:N2,1).';

% spectral region to be represented: between indexes 1 and N/2
freqregion = [1:N2];

k=0; % set frame counter

while(nread==N2)
%
% overlap/add analysis, ODFT
%
    k = k+1; % frame counter
    
    if (plots==1)
    % time representation of the signal
    figure(1);
    plot([0:N-1],idata(1:N));
    xlabel('Samples');
    ylabel('Amplitude');
    title('Time signal');
    end
    
    idata=idata.*win;
    fdata = idata.*direxp; % this is sub-optimal ODFT computation
    odft=fft(fdata);       % this is sub-optimal ODFT computation    

    if (plots==1)
    figure(2);
    plot(FS/N*(freqregion-1),20*log10(abs(odft(freqregion))));
    xlabel('Frequency (Hz)');
    ylabel('Spectral density (dB)');
    title('Spectrum');
    end
   
%
% ODFT 2 MDCT
%
    mdct=real(odft(1:N2)).*cosargmdct + imag(odft(1:N2)).*sinargmdct; % MDCT computation

%----------------- quantization area ------------------------------%
    if (white==1)
        % add code here
    else
        % add code here
    end
%----------------- quantization area ------------------------------%

%    
% IMDCT 2 IODFT
%
    fdata(1:N2)= 2 * (mdct(1:N2).*cosargmdct + 1i*mdct(1:N2).*sinargmdct);
    
    
%
% IODFT, overlap/add reconstruction
%
    fdata(N:-1:N2+1)=conj(fdata(1:N2));
    ofdata=ifft(fdata);        % this is sub-optimal IODFT computation
    ofdata=ofdata.*invexp;     % this is sub-optimal IODFT computation
    odata=real(ofdata);
    odata=odata.*win;
    tmpdata(1:N2)=floor(0.5+osdata+odata(1:N2));
    osdata=odata(N2+1:N);
    fwrite(fidw, tmpdata(1:N2), 'short');
    
    idata(1,1:N2)=data(1:N2,1)';
    [data, nread] = fread(fidr, N2, 'short'); % reads new half-segment
    if nread<N2,
        data(nread+1:N2,1)=zeros(N2-nread,1);
    end
    idata(1,N2+1:N)=data(1:N2,1)';
    if (plots==1)
    pause;
    end
end



fclose(fidr);
fclose(fidw);
close all;

disp('END of processing !');

%
% now align and check the difference between original and quantized signal
%

fidr = fopen('tmpinpaudio.pcm','r');
fidw = fopen('tmpoutaudio.pcm','r');

[data1, nread1] = fread(fidr, 'short');
[data2, nread2] = fread(fidw, 'short');
fclose(fidr);
fclose(fidw);
audiowrite(outfile, data2/(2^(bits-1)-1), FS); % output WAV file
% wavwrite(data2/(2^(bits-1)-1), FS, bits, outfile); % OLD

shift=N2; % this is the system delay
nread=nread2-shift; % size of the audio output after alignment


difsignal = data1(1:nread)-data2(1+shift:nread2);
plot([0:(nread-1)]*1000/FS, difsignal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Coding noise');

disp('If the figure displays ZERO, we have perfect reconstruction !');
Psignal = sum(data1(1:nread).^2);
Pnoise = sum(difsignal.^2);

if (Pnoise==0)
    disp('SNR= infinity');
else
    disp('SNR= '); 10*log10(Psignal/Pnoise)
end
