% =========================================================================
% CIM - GERADOR DE GRÁFICOS DE APRECIAÇÃO SUBJETIVA ITU-R (PREMIUM)
% IDENTIDADE VISUAL: Marcadores preenchidos, grelha translúcida, box on
% =========================================================================
close all; clear all;

% Eixo X comum (Valores de SNR Global Alvo)
snr_vals = [5, 10, 15, 20];

% -------------------------------------------------------------------------
% INTRODUÇÃO DE DADOS (Insiram aqui as notas médias da dupla)
% -------------------------------------------------------------------------
% Ponto 2: Dados da DFT 
dft_var1 = [1.0, 2.0, 3.0, 4.0]; 
dft_var2 = [2.0, 3.0, 4.0, 4.5]; 

% Ponto 3: Dados da MDCT 
mdct_var1 = [1.5, 2.5, 3.5, 4.2]; 
mdct_var2 = [2.5, 3.8, 4.5, 4.8]; 

% Configuração padrão da escala do Eixo Y da norma ITU-R BS.1116
labels_itu = {
    '1. Muito incomodativo', ...
    '2. Incomodativo', ...
    '3. Ligeiramente incomodativo', ...
    '4. Percetível, não incomod.', ... 
    '5. Impercetível'
};

% Cores Oficiais Modernas do MATLAB
c_blue   = [0.0000 0.4470 0.7410];
c_orange = [0.8500 0.3250 0.0980];
c_green  = [0.4660 0.6740 0.1880];
c_purple = [0.4940 0.1840 0.5560];

% =========================================================================
% IMAGEM 1: APRECIAÇÃO SUBJETIVA - APENAS DFT
% =========================================================================
figure('Position', [100, 100, 700, 400]); 
p1 = plot(snr_vals, dft_var1, '-o', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_blue, 'MarkerFaceColor', c_blue); hold on;
p2 = plot(snr_vals, dft_var2, '-s', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_orange, 'MarkerFaceColor', c_orange);

ylim([0.8, 5.8]); yticks(1:5); yticklabels(labels_itu);
xlim([4, 21]); xticks(snr_vals);

xlabel('SNR Global (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('Apreciação Subjetiva: DFT', 'FontSize', 12, 'FontWeight', 'bold');
legend([p1, p2], 'Var. I (Ruído Branco)', 'Var. II (Ruído Colorido)', 'Location', 'northwest');

% O toque de mestre visual: Grelha translúcida e caixa fechada
set(gca, 'FontSize', 10, 'FontName', 'Helvetica', 'Box', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.25);
grid on; hold off;

saveas(gcf, 'Fig1_Apreciacao_Subjetiva_DFT.png');

% =========================================================================
% IMAGEM 2: APRECIAÇÃO SUBJETIVA - APENAS MDCT
% =========================================================================
figure('Position', [150, 150, 700, 400]);
p3 = plot(snr_vals, mdct_var1, '-^', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_green, 'MarkerFaceColor', c_green); hold on;
p4 = plot(snr_vals, mdct_var2, '-d', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_purple, 'MarkerFaceColor', c_purple);

ylim([0.8, 5.8]); yticks(1:5); yticklabels(labels_itu);
xlim([4, 21]); xticks(snr_vals);

xlabel('SNR Global (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('Apreciação Subjetiva: MDCT', 'FontSize', 12, 'FontWeight', 'bold');
legend([p3, p4], 'Var. I (Ruído Branco)', 'Var. II (Ruído Colorido)', 'Location', 'northwest');

set(gca, 'FontSize', 10, 'FontName', 'Helvetica', 'Box', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.25);
grid on; hold off;

saveas(gcf, 'Fig2_Apreciacao_Subjetiva_MDCT.png');

% =========================================================================
% IMAGEM 3: COMPARAÇÃO GLOBAL
% =========================================================================
figure('Position', [200, 200, 750, 450]); 
% Linhas tracejadas nas Variantes I e sólidas nas Variantes II para melhor distinção a preto e branco
g1 = plot(snr_vals, dft_var1, '--o', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_blue, 'MarkerFaceColor', c_blue); hold on;
g2 = plot(snr_vals, dft_var2, '-s', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_orange, 'MarkerFaceColor', c_orange);
g3 = plot(snr_vals, mdct_var1, '--^', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_green, 'MarkerFaceColor', c_green);
g4 = plot(snr_vals, mdct_var2, '-d', 'LineWidth', 1.5, 'MarkerSize', 7, 'Color', c_purple, 'MarkerFaceColor', c_purple);

ylim([0.8, 6.2]); yticks(1:5); yticklabels(labels_itu);
xlim([4, 21]); xticks(snr_vals);

xlabel('SNR Global (dB)', 'FontSize', 11, 'FontWeight', 'bold');
title('Comparação Global: DFT vs MDCT', 'FontSize', 12, 'FontWeight', 'bold');
legend([g1, g2, g3, g4], ...
    'DFT - Var I (Branco)', 'DFT - Var II (Colorido)', ...
    'MDCT - Var I (Branco)', 'MDCT - Var II (Colorido)', ...
    'Location', 'northwest', 'NumColumns', 2); 

set(gca, 'FontSize', 10, 'FontName', 'Helvetica', 'Box', 'on', 'GridLineStyle', '--', 'GridAlpha', 0.25);
grid on; hold off;

saveas(gcf, 'Fig3_Comparacao_Global.png');
disp('-> Gráficos Premium gerados com sucesso!');