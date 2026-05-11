%
% Este conjunto de comandos Matlab faz parte das aulas
% da UC CODIFICAÇÃO DE INFORMAÇÃO MULTIMÉDIA
% (CIM), integrada no Mestrado  em
% Engenharia Electrotécnica e Computadores (M.EEC), da
% Faculdade de Engenharia da Universidade do Porto (FEUP).
%
% É permitida a utilização por outras pessoas que não
% frequentem a disciplina, unicamente para objectivos académicos
% ou auto-estudo e desde que seja feita menção, em qualquer trabalho
% publicado ou exposto publicamente, da fonte: "Elementos de Estudo
% da Unidade Curricular CIM do curso M.EEC (FEUP)".
%
% Qualquer outra utilização deverá ser objecto de um pedido
% de autorização dirigido por escrito ao autor e que só
% será válido se houver uma resposta EXPLÍCITA e favorável.
%
% (c) Aníbal Ferreira; ajf(at)fe.up.pt
%

% objectivo: análise, modificação no domínio das frequências e síntese de
% um sinal áudio



% input audio file (raw PCM)
%
inpfile = 'sound.wav';
outfile = 'newsound.wav';

%
% N        = comprimento da transformada DFT (calculada através da FFT)
% N/2      = sobreposição entre tramas, i.e, overlap é de 50%
% win      = janela de análise (qual a sua utilidade ?)
%
N=1024; N2=N/2;
win=sin(pi/N*[0:(N-1)].');

%
% lê ficheiro áudio
%
[datain,FS]=audioread(inpfile); % NOTA: dados surgem normalizados entre 0 e 1
% [datain,FS,bits]=wavread(inpfile); % OLD
% disp('bits por amostra: '); bits

nread=length(datain);
disp('Frequência de amostragem: '); FS

% vector de saída
dataout = zeros(size(datain));

%
% NOTA: alternativa de leitura de ficheiro áudio em PCM "raw"
% (ou seja sem cabeçalho identificador), a leitura é binária
% e através de "shorts" já que se sabe que cada amostra está
% representada numa palavra representando um inteiro de 16 bits,
% como esta palavra inclui sinal, a gama possível é [-32768, 32767]
% a gama normalizada correpondente é [-1.0, 1.0[
%
%fid = fopen(inpfile,'r');
%[data, nread] = fread(fid,'short');
%fclose(fid);



tmpdata = zeros(N,1);	%% input (int) data

% NOTA: usamos para visualizar só os pontos
% entre 1 e N/2+1 da análise de Fourier (calculada
% através da FFT) já que os restantes são uma
% repeticão destes, porquê ?
regiaofreq = [1:N2+1];

frame=1;

while((frame+1)*N2 < nread),

%    
% NOTA: a instrução seguinte é usada para
% delimitar o segmento que se pretende seleccionar

   varre = [1+(frame-1)*N2:(frame+1)*N2];

   tmpdata=datain(varre);

% representa sinal temporal
   figure(1);
   plot([0:N-1],tmpdata(1:N));
   xlabel('Amostras temporais');
   ylabel('Amplitude normalizada');
   title('Sinal nos Tempos');
   
%
% FFT
%
   tmpdata=tmpdata.*win;
   fdata=fft(tmpdata);
   magnitude=eps+abs(fdata);
   fase=phase(fdata);

   figure(2);
   plot(FS/N*(regiaofreq-1),20*log10(magnitude(regiaofreq)));
   xlabel('Frequência');
   ylabel('Densidade Espectral (dB)');
   title('Sinal nas Frequências');
   
   
%   aqui tem lugar a modificação espectral

%   ruído de quantização fica com espectro plano
%   ...


%  ruído de quantização fica com espectro semelhante ao do sinal (SMR constante)
%  ...



%
% IFFT
%
%    fdata=magnitude.*exp(i*fase);
   
   figure(2);
   hold on;
   plot(FS/N*(regiaofreq-1),20*log10(abs(fdata(regiaofreq))),'r');
   hold off;

%    fdata(1)=magnitude(1)*sign(exp(i*fase(1))); % porquê ?
%    fdata(N2+1)=magnitude(N2+1)*sign(exp(i*fase(N2+1))); % porquê ?
   fdata(N:-1:N2+2)=conj(fdata(2:N2)); % porquê ?
   
   tmpdata=real(ifft(fdata));
   tmpdata=tmpdata.*win;
   dataout(varre)=dataout(varre)+tmpdata;
   
   figure(1);
   hold on;
   plot([0:N-1],tmpdata(1:N),'r');
   hold off;

%    pause;
   frame=frame+1

end

% grava sinal de saída
audiowrite(outfile, dataout, FS);
% wavwrite(dataout, FS, bits, outfile); % OLD

% calculo do SNR: de modo identico ao utilizado em qzmdct.m

disp('Fim de processamento.');