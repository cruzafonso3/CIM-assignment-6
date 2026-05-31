%---------------------------------------------------------------------------
% GERAÇÃO AUTOMÁTICA EM LOTE - PONTO 3 (MDCT)
%---------------------------------------------------------------------------
close all; clear all;

% =============================================================
% MATRIZ DE TESTES (AJUSTA OS SMR AQUI!)
% =============================================================
% Formato: [Variante, SMR_dB_usado, SNR_Alvo_Desejado]
% ATENÇÃO: Os valores de SMR abaixo são EXEMPLOS! Vais ter de
% os afinar até o output da consola dar perto dos 20, 15, 10 e 5 dB.
casos_teste = [
    1,  16,33, 20; % Variante 1 (Branco),   Alvo: 20dB
    1,  9.6, 15; % Variante 1 (Branco),   Alvo: 15dB
    1,   2, 10; % Variante 1 (Branco),   Alvo: 10dB
    1,   -6.6,  5; % Variante 1 (Branco),   Alvo: 5dB
    2,  9.965, 20; % Variante 2 (Colorido), Alvo: 20dB
    2,  9.375, 15; % Variante 2 (Colorido), Alvo: 15dB
    2,  8.41, 10; % Variante 2 (Colorido), Alvo: 10dB
    2,   6.95,  5  % Variante 2 (Colorido), Alvo: 5dB
];
% =============================================================

% input audio file
inpfile = 'sound.wav';
[datain,FS] = audioread(inpfile); 
bits = 16;

% N = size of the ODFT and MDCT transforms
N=1024; N2=N/2;
win=sin(pi/N*([0:N-1]+0.5));
direxp = exp(-1i*pi*[0:N-1]/N);
invexp = exp( 1i*pi*[0:N-1]/N);
cosargmdct=cos(pi*(0.5+[0:(N2-1)])*(1.0+N2)/N);
sinargmdct=sin(pi*(0.5+[0:(N2-1)])*(1.0+N2)/N);
freqregion = [1:N2];

disp('========================================');
disp(' A INICIAR GERAÇÃO AUTOMÁTICA (MDCT)...');
disp('========================================');

for teste = 1:size(casos_teste, 1)
    
    TIPO_VARIANTE = casos_teste(teste, 1);
    SMR_dB        = casos_teste(teste, 2);
    SNR_ALVO      = casos_teste(teste, 3);
    white         = (TIPO_VARIANTE == 1); % white=1 (Branco), white=0 (Colorido)
    
    disp(['A processar Variante ', num2str(TIPO_VARIANTE), ' - Alvo: ', num2str(SNR_ALVO), ' dB (SMR = ', num2str(SMR_dB), ') ...']);
    
    % Inicializar/Resetar variáveis do Overlap-Add
    idata   = zeros(1,N);	
    odata   = zeros(1,N);	
    osdata  = zeros(1,N2);  
    tmpdata = zeros(1, N2);
    fdata   = zeros(1,N);	
    ofdata  = zeros(1,N);
    
    % Regenerar o ficheiro PCM temporário de input para cada iteração
    fid = fopen('tmpinpaudio.pcm','w');
    fwrite(fid, datain*(2^(bits-1)-1), 'short'); 
    fclose(fid);
    
    fidr = fopen('tmpinpaudio.pcm','r');
    fidw = fopen('tmpoutaudio.pcm','w');
    
    [data, nread] = fread(fidr, N2, 'short');
    idata(1,N2+1:N)=data(1:N2,1).';
    
    while(nread==N2)
        idata = idata.*win;
        fdata = idata.*direxp; 
        odft  = fft(fdata);       
       
        % ODFT 2 MDCT
        mdct = real(odft(1:N2)).*cosargmdct + imag(odft(1:N2)).*sinargmdct; 
        mdct_original = mdct; % Guardar para o gráfico
    
        % --- QUANTIZAÇÃO ---
        SMR_lin = 10^(SMR_dB / 10);
        if (white)
            P_media = mean(mdct.^2);
            delta = sqrt((12 * P_media) / SMR_lin);
            if delta > 0, mdct = round(mdct / delta) * delta; end
        else
            P_k = mdct.^2;
            delta_k = sqrt((12 .* P_k) ./ SMR_lin);
            delta_k(delta_k == 0) = eps; 
            mdct = round(mdct ./ delta_k) .* delta_k;
        end
    
        % IMDCT 2 IODFT e Reconstrução
        fdata(1:N2) = 2 * (mdct(1:N2).*cosargmdct + 1i*mdct(1:N2).*sinargmdct);
        fdata(N:-1:N2+1) = conj(fdata(1:N2));
        ofdata = ifft(fdata);        
        ofdata = ofdata.*invexp;     
        odata = real(ofdata).*win;
        
        tmpdata(1:N2) = floor(0.5 + osdata + odata(1:N2));
        osdata = odata(N2+1:N);
        fwrite(fidw, tmpdata(1:N2), 'short');
        
        idata(1,1:N2) = data(1:N2,1)';
        [data, nread] = fread(fidr, N2, 'short');
        if nread < N2
            data(nread+1:N2,1) = zeros(N2-nread,1);
        end
        idata(1,N2+1:N) = data(1:N2,1)';
    end
    
    fclose(fidr);
    fclose(fidw);
    
    % ==========================================
    % GERAR GRÁFICO DO ESPETRO (ÚLTIMA FRAME)
    % ==========================================
    figure(1);
    p_orig = plot(FS/N*(freqregion-1), 20*log10(abs(mdct_original)+eps), 'b', 'LineWidth', 0.8); hold on;
    mag_quantizada = abs(mdct); mag_quantizada(mag_quantizada < 1e-5) = 1e-5; 
    erro_espetral = abs(mdct_original - mdct); erro_espetral(erro_espetral < 1e-5) = 1e-5; 
    
    p_quant = plot(FS/N*(freqregion-1), 20*log10(mag_quantizada), 'r', 'LineWidth', 0.8);
    p_erro = plot(FS/N*(freqregion-1), 20*log10(erro_espetral), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 0.8);
    
    xlabel('Frequência (Hz)', 'FontWeight', 'bold'); ylabel('Densidade Espectral (dB)', 'FontWeight', 'bold');
    title(sprintf('Espetro MDCT (Var. %d) | Alvo: %ddB', TIPO_VARIANTE, SNR_ALVO), 'FontWeight', 'bold');
    legend([p_orig, p_quant, p_erro], 'Original', 'Quantizado', 'Erro (Ruído)', 'Location', 'northeast');
    ylim([-100, max(20*log10(abs(mdct_original)+eps))+5]); hold off;
    
    % Guardar gráfico com nome preparado para o LaTeX
    nome_figura = sprintf('espetro_mdct_v%d_%ddb.png', TIPO_VARIANTE, SNR_ALVO);
    saveas(figure(1), nome_figura);
    
    % ==========================================
    % VERIFICAR O ÁUDIO E CALCULAR SNR FINAL
    % ==========================================
    fidr = fopen('tmpinpaudio.pcm','r');
    fidw = fopen('tmpoutaudio.pcm','r');
    [data1, nread1] = fread(fidr, 'short');
    [data2, nread2] = fread(fidw, 'short');
    fclose(fidr);
    fclose(fidw);
    
    % Gravar o Áudio Final
    nome_audio = sprintf('audio_mdct_v%d_%ddb.wav', TIPO_VARIANTE, SNR_ALVO);
    audiowrite(nome_audio, data2/(2^(bits-1)-1), FS); 
    
    shift = N2; 
    nread_final = nread2 - shift; 
    difsignal = data1(1:nread_final) - data2(1+shift:nread2);
    
    Psignal = sum(data1(1:nread_final).^2);
    Pnoise = sum(difsignal.^2);
    
    if (Pnoise == 0)
        disp(' -> Feito! SNR = INFINITO (Reconstrução Perfeita)');
    else
        snr_final = 10*log10(Psignal/Pnoise);
        disp([' -> Feito! SNR Final = ', num2str(snr_final), ' dB']);
    end
    
    % Fechar gráficos para poupar memória na próxima ronda
    close all; 
end

disp('========================================');
disp(' TUDO CONCLUÍDO COM SUCESSO!');
disp('========================================');