%
% Este conjunto de comandos Matlab faz parte das aulas
% da UC CODIFICAÇÃO DE INFORMAÇÃO MULTIMÉDIA (CIM)
%

% =============================================================
% CONFIGURAÇÕES DA QUANTIZAÇÃO (PONTO 2)
% =============================================================
% Escolhe qual a variante a utilizar:
% 1 -> Variante I (Ruído Branco / Passo de quantização igual)
% 2 -> Variante II (Ruído Colorido / Passo de quantização individual)
TIPO_VARIANTE = 2; 

% Ajusta o SMR para atingir os diferentes SNR Globais (5, 10, 15 ou 20 dB)
% SMR_dB = 16.18; % SNR 20dB Variante I
% SMR_dB = 9.5; % SNR 15dB Variante I
% SMR_dB = 1.84; % SNR 10dB Variante I
% SMR_dB = -6.8; % SNR 5dB Variante I

% SMR_dB = 19.9; % SNR 20dB Variante II
% SMR_dB = 15.79; % SNR 15dB Variante II
% SMR_dB = 10.3; % SNR 10dB Variante II
 SMR_dB = 7.43; % SNR 5dB Variante II
% =============================================================

% input audio file (raw PCM)
inpfile = 'sound.wav';
outfile = 'newsound.wav';

% N        = comprimento da transformada DFT (calculada através da FFT)
% N/2      = sobreposição entre tramas, i.e, overlap é de 50%
N=1024; N2=N/2;
win=sin(pi/N*[0:(N-1)].');

% lê ficheiro áudio
[datain,FS]=audioread(inpfile); 
nread=length(datain);
disp('Frequência de amostragem: '); disp(FS);

% vector de saída
dataout = zeros(size(datain));
tmpdata = zeros(N,1);	
regiaofreq = [1:N2+1];

frame=1;

while((frame+1)*N2 < nread)
    
   % delimitar o segmento que se pretende seleccionar
   varre = [1+(frame-1)*N2:(frame+1)*N2];
   tmpdata=datain(varre);
   
   % ==========================================
   % GRÁFICO 1 - TEMPO (SINAL ORIGINAL)
   % ==========================================
   figure(1);
   p_orig_t = plot([0:N-1], tmpdata(1:N), 'b', 'LineWidth', 0.8);
   xlabel('Amostras temporais', 'FontSize', 11, 'FontWeight', 'bold');
   ylabel('Amplitude normalizada', 'FontSize', 11, 'FontWeight', 'bold');
   title('Sinal Temporal: Original vs Reconstruído', 'FontSize', 13, 'FontWeight', 'bold');
   axis tight;
   
   % FFT
   tmpdata=tmpdata.*win;
   fdata=fft(tmpdata);
   magnitude=eps+abs(fdata);
   fase=angle(fdata); 
   
   % ==========================================
   % GRÁFICO 2 - FREQUÊNCIAS (SINAL ORIGINAL)
   % ==========================================
   figure(2);
   p_orig_f = plot(FS/N*(regiaofreq-1), 20*log10(magnitude(regiaofreq)), 'b', 'LineWidth', 0.8);
   xlabel('Frequência (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
   ylabel('Densidade Espectral (dB)', 'FontSize', 11, 'FontWeight', 'bold');
   title(sprintf('Comparação de Espetros (Variante %d)', TIPO_VARIANTE), 'FontSize', 13, 'FontWeight', 'bold');
   axis tight;
   
   % -------------------------------------------------------------
   % AQUI TEM LUGAR A MODIFICAÇÃO ESPECTRAL
   % -------------------------------------------------------------
   SMR_lin = 10^(SMR_dB / 10);
   
   % GUARDAR O SINAL ORIGINAL ANTES DE QUANTIZAR PARA CALCULAR O ERRO!
   fdata_original = fdata; 
   
   if TIPO_VARIANTE == 1
       % --- VARIANTE I (Ruído Branco) ---
       P_media = mean(abs(fdata(regiaofreq)).^2);
       delta = sqrt((12 * P_media) / SMR_lin);
       
       if delta > 0
           real_q = round(real(fdata(regiaofreq)) / delta) * delta;
           imag_q = round(imag(fdata(regiaofreq)) / delta) * delta;
           fdata(regiaofreq) = real_q + 1i * imag_q;
       end
       
   elseif TIPO_VARIANTE == 2
       % --- VARIANTE II (Ruído Colorido) ---
       P_k = abs(fdata(regiaofreq)).^2;
       delta_k = sqrt((12 .* P_k) ./ SMR_lin);
       delta_k(delta_k == 0) = eps; % Proteção divisão por zero
       
       real_q = round(real(fdata(regiaofreq)) ./ delta_k) .* delta_k;
       imag_q = round(imag(fdata(regiaofreq)) ./ delta_k) .* delta_k;
       fdata(regiaofreq) = real_q + 1i * imag_q;
   end
   % -------------------------------------------------------------
   
   % ==========================================
   % GRÁFICO 2 - FREQUÊNCIAS (SINAL QUANTIZADO E ERRO)
   % ==========================================
   figure(2); hold on;
   
   % TRUQUE PARA O GRÁFICO: Criar um "chão" no -Infinito (-100 dB)
   mag_quantizada = abs(fdata(regiaofreq));
   mag_quantizada(mag_quantizada < 1e-5) = 1e-5; % Valores nulos são forçados a -100 dB
   
   % CÁLCULO DO ERRO ESPETRAL (Original - Quantizado)
   erro_espetral = abs(fdata_original(regiaofreq) - fdata(regiaofreq));
   erro_espetral(erro_espetral < 1e-5) = 1e-5; % Aplicar o chão também ao erro
   
   % Plot do Sinal Quantizado (Vermelho)
   p_quant_f = plot(FS/N*(regiaofreq-1), 20*log10(mag_quantizada), 'r', 'LineWidth', 0.8);
   
   % Plot do Erro de Quantização (Amarelo Dourado para visibilidade)
   p_erro_f = plot(FS/N*(regiaofreq-1), 20*log10(erro_espetral), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 0.8);
   
   legend([p_orig_f, p_quant_f, p_erro_f], 'Sinal Original', 'Sinal Quantizado', 'Erro de Quantização', 'Location', 'northeast', 'FontSize', 10);
   
   % Limitar o eixo Y para o gráfico não descer além do nosso "chão" artificial
   ylim([-100, max(20*log10(magnitude(regiaofreq)))+5]); 
   hold off;
   
   % IFFT
   fdata(N:-1:N2+2)=conj(fdata(2:N2)); 
   
   tmpdata=real(ifft(fdata));
   tmpdata=tmpdata.*win;
   dataout(varre)=dataout(varre)+tmpdata;
   
   % ==========================================
   % GRÁFICO 1 - TEMPO (SINAL RECONSTRUÍDO)
   % ==========================================
   figure(1); hold on;
   p_recon_t = plot([0:N-1], tmpdata(1:N), 'r', 'LineWidth', 0.8);
   legend([p_orig_t, p_recon_t], 'Original (Janelado)', 'Reconstruído', 'Location', 'northeast', 'FontSize', 10);
   hold off;
   
   frame=frame+1;
end

% grava sinal de saída
audiowrite(outfile, dataout, FS);

% =============================================================
% CÁLCULO DO SNR E GRÁFICO DO RUÍDO FINAL
% =============================================================
tamanho_valido = length(dataout); 
difsignal = datain(1:tamanho_valido) - dataout(1:tamanho_valido);

Psignal = sum(datain(1:tamanho_valido).^2);
Pnoise = sum(difsignal.^2);

% Calcular primeiro o SNR para o colocar no título do gráfico
if (Pnoise == 0)
    snr_str = 'Infinito (Perfeito)';
    disp('SNR = infinity (Reconstrução Perfeita!)');
else
    snr_val = 10 * log10(Psignal / Pnoise);
    snr_str = sprintf('%.2f dB', snr_val);
    disp(['SNR Global atingido: ', snr_str]);
end

% ==========================================
% GRÁFICO 3 - RUÍDO GLOBAL NOS TEMPOS
% ==========================================
figure(3);
plot([0:(tamanho_valido-1)]*1000/FS, difsignal, 'Color', [0.850 0.325 0.098], 'LineWidth', 0.8);
xlabel('Tempo (s)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Amplitude do Erro', 'FontSize', 11, 'FontWeight', 'bold');
title(sprintf('Ruído de Codificação (Variante %d) | SNR: %s', TIPO_VARIANTE, snr_str), 'FontSize', 13, 'FontWeight', 'bold');
axis tight;

disp('Fim de processamento.');