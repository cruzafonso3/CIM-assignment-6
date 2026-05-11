%==========================================================================
%
% Ex2_3_batch_runner.m
%
% Batch runner – gera resultados para todas as combinações:
%   Transform: DFT, MDCT
%   Variant:   White noise (I), Colored noise (II)
%   SNR targets: ~5, 10, 15, 20 dB
%
% For each combination, the function sweeps SMR_dB values until the
% measured SNR is within ±1 dB of the target, then records the result.
%
% At the end it prints a summary table and plots the ITU-R subjective
% score versus SNR.
%
% USAGE: Place sound.wav in the MATLAB working directory, then run this
%        script.  Each variant writes its own output WAV file.
%
% M.EEC045 - Codificação de Informação Multimédia - FEUP 2025/2026
%==========================================================================

close all;
clear all;

inpfile = 'sound.wav';

% -----------------------------------------------------------------------
% ITU-R BS.1116 subjective impairment scale (fill these in after listening)
% -----------------------------------------------------------------------
% 5.0 = Imperceptible
% 4.0 = Perceptible, but not annoying
% 3.0 = Slightly annoying
% 2.0 = Annoying
% 1.0 = Very annoying

% -----------------------------------------------------------------------
% Targets
% -----------------------------------------------------------------------
SNR_targets = [5, 10, 15, 20];   % dB

% -----------------------------------------------------------------------
% Helper: run DFT white quantization for given SMR, return measured SNR
% -----------------------------------------------------------------------

    function SNR_meas = run_DFT_white(inpfile, SMR_dB, outfile)
        [datain, FS] = audioread(inpfile);
        datain = datain(:,1);
        N  = 1024; N2 = N/2;
        win = sin(pi/N * (0:N-1).');
        SMR_lin = 10^(SMR_dB/10);
        dataout = zeros(size(datain));
        frame = 1;
        nread = length(datain);
        while (frame+1)*N2 < nread
            varre   = 1+(frame-1)*N2 : (frame+1)*N2;
            tmpdata = datain(varre) .* win;
            fdata   = fft(tmpdata);
            P_frame = sum(abs(fdata).^2);
            if P_frame > 0
                Delta  = sqrt(6 * P_frame / (SMR_lin * N));
                fdata_q = Delta*(floor(real(fdata)/Delta+0.5) + 1i*floor(imag(fdata)/Delta+0.5));
            else
                fdata_q = fdata;
            end
            fdata_q(N:-1:N2+2) = conj(fdata_q(2:N2));
            tmpdata = real(ifft(fdata_q)) .* win;
            dataout(varre) = dataout(varre) + tmpdata;
            frame = frame+1;
        end
        dataout = max(-1, min(1, dataout));
        audiowrite(outfile, dataout, FS);
        neval = min(length(datain), length(dataout));
        d = datain(1:neval) - dataout(1:neval);
        Ps = sum(datain(1:neval).^2);
        Pn = sum(d.^2);
        if Pn==0; SNR_meas=Inf; else; SNR_meas=10*log10(Ps/Pn); end
    end

    function SNR_meas = run_DFT_colored(inpfile, SMR_dB, outfile)
        [datain, FS] = audioread(inpfile);
        datain = datain(:,1);
        N  = 1024; N2 = N/2;
        win = sin(pi/N * (0:N-1).');
        SMR_lin = 10^(SMR_dB/10);
        step_scale = sqrt(12/SMR_lin);
        dataout = zeros(size(datain));
        frame = 1;
        nread = length(datain);
        while (frame+1)*N2 < nread
            varre   = 1+(frame-1)*N2 : (frame+1)*N2;
            tmpdata = datain(varre) .* win;
            fdata   = fft(tmpdata);
            mag     = abs(fdata);
            Delta_k = mag * step_scale;
            fdata_q = zeros(N,1);
            for ki=1:N
                if Delta_k(ki)>0
                    fdata_q(ki) = Delta_k(ki)*(floor(real(fdata(ki))/Delta_k(ki)+0.5)+...
                                               1i*floor(imag(fdata(ki))/Delta_k(ki)+0.5));
                end
            end
            fdata_q(N:-1:N2+2) = conj(fdata_q(2:N2));
            tmpdata = real(ifft(fdata_q)) .* win;
            dataout(varre) = dataout(varre) + tmpdata;
            frame = frame+1;
        end
        dataout = max(-1, min(1, dataout));
        audiowrite(outfile, dataout, FS);
        neval = min(length(datain), length(dataout));
        d = datain(1:neval) - dataout(1:neval);
        Ps = sum(datain(1:neval).^2);
        Pn = sum(d.^2);
        if Pn==0; SNR_meas=Inf; else; SNR_meas=10*log10(Ps/Pn); end
    end

    function SNR_meas = run_MDCT_white(inpfile, SMR_dB, outfile)
        [datain, FS] = audioread(inpfile);
        datain = datain(:,1);
        bits = 16;
        fid = fopen('tmp_in.pcm','w');
        fwrite(fid, datain*(2^(bits-1)-1), 'short');
        fclose(fid);
        N=1024; N2=N/2;
        win=sin(pi/N*([0:N-1]+0.5));
        direxp=exp(-1i*pi*(0:N-1)/N);
        invexp=exp( 1i*pi*(0:N-1)/N);
        cosargmdct=cos(pi*(0.5+(0:N2-1))*(1.0+N2)/N);
        sinargmdct=sin(pi*(0.5+(0:N2-1))*(1.0+N2)/N);
        SMR_lin=10^(SMR_dB/10);
        idata=zeros(1,N); osdata=zeros(1,N2); tmpdata=zeros(1,N2); fdata=zeros(1,N);
        fidr=fopen('tmp_in.pcm','r'); fidw=fopen('tmp_out.pcm','w');
        [data, nb]=fread(fidr,N2,'short');
        idata(1,N2+1:N)=data(1:N2,1).';
        while nb==N2
            id=idata.*win;
            fd=id.*direxp; od=fft(fd);
            mdct=real(od(1:N2)).*cosargmdct+imag(od(1:N2)).*sinargmdct;
            P=sum(mdct.^2);
            if P>0 && ~isinf(SMR_dB)
                D=sqrt(12*P/(SMR_lin*N2));
                mdct_q=D*floor(mdct/D+0.5);
            else
                mdct_q=mdct;
            end
            fdata(1:N2)=2*(mdct_q.*cosargmdct+1i*mdct_q.*sinargmdct);
            fdata(N:-1:N2+1)=conj(fdata(1:N2));
            of=ifft(fdata).*invexp;
            od2=real(of).*win;
            tmpdata(1:N2)=floor(0.5+osdata+od2(1:N2));
            osdata=od2(N2+1:N);
            fwrite(fidw,tmpdata(1:N2),'short');
            idata(1,1:N2)=data(1:N2,1).';
            [data,nb]=fread(fidr,N2,'short');
            if nb<N2; data(nb+1:N2,1)=zeros(N2-nb,1); end
            idata(1,N2+1:N)=data(1:N2,1).';
        end
        fclose(fidr); fclose(fidw);
        fr=fopen('tmp_in.pcm','r'); fw=fopen('tmp_out.pcm','r');
        [d1,~]=fread(fr,'short'); [d2,~]=fread(fw,'short');
        fclose(fr); fclose(fw);
        audiowrite(outfile,d2/(2^(bits-1)-1),FS);
        sh=N2; nev=min(length(d1),length(d2)-sh);
        dif=d1(1:nev)-d2(1+sh:sh+nev);
        Ps=sum(d1(1:nev).^2); Pn=sum(dif.^2);
        if Pn==0; SNR_meas=Inf; else; SNR_meas=10*log10(Ps/Pn); end
    end

    function SNR_meas = run_MDCT_colored(inpfile, SMR_dB, outfile)
        [datain, FS] = audioread(inpfile);
        datain = datain(:,1);
        bits = 16;
        fid = fopen('tmp_in.pcm','w');
        fwrite(fid, datain*(2^(bits-1)-1), 'short');
        fclose(fid);
        N=1024; N2=N/2;
        win=sin(pi/N*([0:N-1]+0.5));
        direxp=exp(-1i*pi*(0:N-1)/N);
        invexp=exp( 1i*pi*(0:N-1)/N);
        cosargmdct=cos(pi*(0.5+(0:N2-1))*(1.0+N2)/N);
        sinargmdct=sin(pi*(0.5+(0:N2-1))*(1.0+N2)/N);
        SMR_lin=10^(SMR_dB/10);
        step_scale=sqrt(12/SMR_lin);
        idata=zeros(1,N); osdata=zeros(1,N2); tmpdata=zeros(1,N2); fdata=zeros(1,N);
        fidr=fopen('tmp_in.pcm','r'); fidw=fopen('tmp_out.pcm','w');
        [data, nb]=fread(fidr,N2,'short');
        idata(1,N2+1:N)=data(1:N2,1).';
        while nb==N2
            id=idata.*win;
            fd=id.*direxp; od=fft(fd);
            mdct=real(od(1:N2)).*cosargmdct+imag(od(1:N2)).*sinargmdct;
            Dk=abs(mdct)*step_scale;
            mdct_q=zeros(1,N2);
            for ki=1:N2
                if Dk(ki)>0; mdct_q(ki)=Dk(ki)*floor(mdct(ki)/Dk(ki)+0.5); end
            end
            fdata(1:N2)=2*(mdct_q.*cosargmdct+1i*mdct_q.*sinargmdct);
            fdata(N:-1:N2+1)=conj(fdata(1:N2));
            of=ifft(fdata).*invexp;
            od2=real(of).*win;
            tmpdata(1:N2)=floor(0.5+osdata+od2(1:N2));
            osdata=od2(N2+1:N);
            fwrite(fidw,tmpdata(1:N2),'short');
            idata(1,1:N2)=data(1:N2,1).';
            [data,nb]=fread(fidr,N2,'short');
            if nb<N2; data(nb+1:N2,1)=zeros(N2-nb,1); end
            idata(1,N2+1:N)=data(1:N2,1).';
        end
        fclose(fidr); fclose(fidw);
        fr=fopen('tmp_in.pcm','r'); fw=fopen('tmp_out.pcm','r');
        [d1,~]=fread(fr,'short'); [d2,~]=fread(fw,'short');
        fclose(fr); fclose(fw);
        audiowrite(outfile,d2/(2^(bits-1)-1),FS);
        sh=N2; nev=min(length(d1),length(d2)-sh);
        dif=d1(1:nev)-d2(1+sh:sh+nev);
        Ps=sum(d1(1:nev).^2); Pn=sum(dif.^2);
        if Pn==0; SNR_meas=Inf; else; SNR_meas=10*log10(Ps/Pn); end
    end

% -----------------------------------------------------------------------
% Sweep SMR to hit each SNR target (±1 dB tolerance)
% -----------------------------------------------------------------------
fprintf('\n========================================================\n');
fprintf(' Batch processing – finding SMR values for target SNRs\n');
fprintf('========================================================\n\n');

% Results tables  [target_SNR, SMR_used, measured_SNR]
res_DFT_w   = zeros(length(SNR_targets), 3);
res_DFT_c   = zeros(length(SNR_targets), 3);
res_MDCT_w  = zeros(length(SNR_targets), 3);
res_MDCT_c  = zeros(length(SNR_targets), 3);

for ti = 1:length(SNR_targets)
    tgt = SNR_targets(ti);
    fprintf('--- Target SNR = %d dB ---\n', tgt);

    % Initial SMR guess: start with SMR = target SNR (rough estimate)
    % then iterate if needed
    for variant = 1:4
        SMR_try = tgt;     % first guess
        for iter = 1:8
            switch variant
                case 1
                    out = sprintf('out_DFT_w_%ddB.wav', tgt);
                    snr_m = run_DFT_white(inpfile, SMR_try, out);
                case 2
                    out = sprintf('out_DFT_c_%ddB.wav', tgt);
                    snr_m = run_DFT_colored(inpfile, SMR_try, out);
                case 3
                    out = sprintf('out_MDCT_w_%ddB.wav', tgt);
                    snr_m = run_MDCT_white(inpfile, SMR_try, out);
                case 4
                    out = sprintf('out_MDCT_c_%ddB.wav', tgt);
                    snr_m = run_MDCT_colored(inpfile, SMR_try, out);
            end
            fprintf('  [V%d] SMR=%.1f dB  ->  SNR=%.2f dB\n', variant, SMR_try, snr_m);
            if abs(snr_m - tgt) < 1.0; break; end
            % Adjust SMR: if SNR too high, increase SMR (less noise)
            SMR_try = SMR_try + (tgt - snr_m) * 0.8;
        end
        switch variant
            case 1; res_DFT_w(ti,:)  = [tgt, SMR_try, snr_m];
            case 2; res_DFT_c(ti,:)  = [tgt, SMR_try, snr_m];
            case 3; res_MDCT_w(ti,:) = [tgt, SMR_try, snr_m];
            case 4; res_MDCT_c(ti,:) = [tgt, SMR_try, snr_m];
        end
    end
end

% -----------------------------------------------------------------------
% Summary table
% -----------------------------------------------------------------------
fprintf('\n\n');
fprintf('=================================================================\n');
fprintf(' RESULTS SUMMARY TABLE\n');
fprintf('=================================================================\n');
fprintf('%-10s | %-12s | %-12s | %-12s | %-12s\n', ...
        'Target SNR', 'DFT White', 'DFT Colored', 'MDCT White', 'MDCT Colored');
fprintf('%s\n', repmat('-',1,67));
for ti = 1:length(SNR_targets)
    fprintf('%-10d | %-12.2f | %-12.2f | %-12.2f | %-12.2f\n', ...
            SNR_targets(ti), ...
            res_DFT_w(ti,3), res_DFT_c(ti,3), ...
            res_MDCT_w(ti,3), res_MDCT_c(ti,3));
end
fprintf('=================================================================\n');
fprintf('(All values are MEASURED SNR in dB)\n\n');

% -----------------------------------------------------------------------
% ITU-R BS.1116 subjective score placeholder
% Fill in after listening tests. Scale: 1=Very annoying ... 5=Imperceptible
% -----------------------------------------------------------------------
% Example placeholder values (replace with your own after listening):
itu_DFT_w   = [1.5, 2.5, 3.5, 4.5];
itu_DFT_c   = [2.0, 3.0, 4.0, 4.8];
itu_MDCT_w  = [1.5, 2.5, 3.5, 4.5];
itu_MDCT_c  = [2.0, 3.0, 4.0, 4.8];

% -----------------------------------------------------------------------
% ITU-R subjective score vs SNR plot
% -----------------------------------------------------------------------
snr_vals = SNR_targets;

figure(20);
plot(snr_vals, itu_DFT_w,  'b-o', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(snr_vals, itu_DFT_c,  'b--s','LineWidth', 1.5, 'MarkerSize', 8);
plot(snr_vals, itu_MDCT_w, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(snr_vals, itu_MDCT_c, 'r--s','LineWidth', 1.5, 'MarkerSize', 8);
hold off;
xlabel('SNR (dB)');
ylabel('ITU-R BS.1116 Subjective Score');
yticks(1:5);
yticklabels({'1 – Very annoying','2 – Annoying','3 – Slightly annoying',...
             '4 – Perceptible, not annoying','5 – Imperceptible'});
ylim([0.5 5.5]);
xlim([3 22]);
title('Subjective Quality vs SNR – ITU-R BS.1116 Scale');
legend('DFT White noise','DFT Colored noise','MDCT White noise','MDCT Colored noise',...
       'Location','southeast');
grid on;

fprintf('Batch complete.  Edit itu_* arrays above with your listening scores.\n');
