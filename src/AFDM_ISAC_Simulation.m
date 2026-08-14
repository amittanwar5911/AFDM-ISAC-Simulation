clc; clear; close all; rng(2025);

%% =====================================================================
%  SECTION 1 — SYSTEM PARAMETERS
% ======================================================================

N         = 512;
Nsym      = 32;
M         = 16;
bps       = log2(M);

SNRdB     = 15;
SNR_range = -5:5:30;
Nmc       = 20;          % Monte Carlo trials per SNR point for BER

fc        = 24e9;
DeltaF    = 22.729e3;
Fs        = N * DeltaF;
Ts        = 1 / Fs;

alphamax  = 2;
kv        = 4;
c1        = (2*(alphamax + kv) + 1) / (2*N);
c2        = 1 / (128*N);

Ncp       = 64;
TAFDM     = (N + Ncp) / Fs;
c_light   = 3e8;

% Pilot structure: comb, 1 per Np subcarriers
Np        = 8;
pilot_idx = 1:Np:N;
data_idx  = setdiff(1:N, pilot_idx);
Npilots   = numel(pilot_idx);
Ndata     = numel(data_idx);

% Spectral efficiency
SE_gross  = bps;                             % bits/subcarrier without pilots
SE_net    = bps * (Ndata/N);                 % bits/subcarrier with pilot overhead
SE_bpshz  = SE_net / (1 + Ncp/N);           % bits/s/Hz (CPP overhead included)

% Processing gain
Gp_dB     = 10*log10(N * Nsym);

fprintf('=================================================================\n');
fprintf('  AFDM-ISAC Simulation v3  —  arXiv:2208.13430\n');
fprintf('=================================================================\n');
fprintf('N=%d | Nsym=%d | %d-QAM | Ncp=%d\n', N, Nsym, M, Ncp);
fprintf('c1=%.6f | c2=%.3e\n', c1, c2);
fprintf('Fs=%.2f MHz | fc=%.0f GHz\n', Fs/1e6, fc/1e9);
fprintf('Pilots: %d/%d (%.1f%% overhead)\n', Npilots, N, 100*Npilots/N);
fprintf('Spectral efficiency: %.3f bits/s/Hz\n', SE_bpshz);
fprintf('Processing gain: %.1f dB\n\n', Gp_dB);

%% =====================================================================
%  SECTION 2 — QAM SYMBOL GENERATION + PILOT INSERTION
% ======================================================================

total_bits = Ndata * Nsym * bps;
tx_bits    = randi([0 1], total_bits, 1);
data_syms  = qammod(tx_bits, M, 'InputType','bit', 'UnitAveragePower',true);
data_mat   = reshape(data_syms, Ndata, Nsym);

X               = zeros(N, Nsym);
X(pilot_idx, :) = 1 + 0j;        % Known pilots
X(data_idx,  :) = data_mat;

fprintf('[TX] Bits: %d | Pilots: %d | Data subcarriers: %d\n\n', ...
        total_bits, Npilots, Ndata);

%% =====================================================================
%  SECTION 3 — AFDM MODULATION (IDAFT, Eq. 1)
% ======================================================================

s_blocks = zeros(N, Nsym);
for k = 1:Nsym
    s_blocks(:,k) = idaft_col_fn(X(:,k), c1, c2, N);
end

%% =====================================================================
%  SECTION 4 — CPP INSERTION
% ======================================================================

s_tx_mat = zeros(N + Ncp, Nsym);
for k = 1:Nsym
    blk = s_blocks(:,k);
    s_tx_mat(:,k) = [blk(end-Ncp+1:end); blk];
end
s_tx = s_tx_mat(:);
s_mat = reshape(s_tx, N+Ncp, Nsym);

fprintf('[TX] Signal: %d samples (%.3f ms)\n\n', numel(s_tx), numel(s_tx)*Ts*1e3);

%% =====================================================================
%  SECTION 5 — TARGET PARAMETERS
% ======================================================================

R_true = [300, 200, 120];
v_true = [ 20,   8,  -8];
P      = numel(R_true);

tau = 2*R_true / c_light;
li  = round(tau / Ts);
fd  = 2*v_true*fc / c_light;
fi  = fd * Ts;
h_i = (randn(1,P) + 1j*randn(1,P)) / sqrt(2*P);

fprintf('[RADAR] %d Targets:\n', P);
for p = 1:P
    fprintf('  T%d: R=%3.0f m | v=%+3.0f m/s | l=%d samp | fd=%.1f Hz | nu=%.3f*Df\n', ...
            p, R_true(p), v_true(p), li(p), fd(p), fd(p)/DeltaF);
end
fprintf('\n');

if any(li >= Ncp)
    error('CPP too short! Max delay=%d >= Ncp=%d.', max(li), Ncp);
end

%% =====================================================================
%  SECTION 6 — RADAR ECHO CHANNEL (Eq. 3)
% ======================================================================

L     = (N + Ncp) * Nsym;
n_vec = (0:L-1).';

r_echo = zeros(L, 1);
for p = 1:P
    s_shift = zeros(size(s_mat));
    for k = 1:Nsym
        s_shift(:,k) = circshift(s_mat(:,k), li(p));
    end
    r_echo = r_echo + h_i(p) * (s_shift(:) .* exp(1j*2*pi*fi(p)*n_vec));
end
r_radar = awgn(r_echo, SNRdB, 'measured');
fprintf('[RADAR RX] SNR = %d dB\n\n', SNRdB);

%% =====================================================================
%  SECTION 7 — COMMUNICATIONS CHANNEL (doubly-selective, separate)
% ======================================================================

Pc   = 3;
li_c = [1, 3, 6];
fi_c = [0.05, -0.12, 0.08] * DeltaF * Ts;
h_c  = (randn(1,Pc) + 1j*randn(1,Pc)) / sqrt(2*Pc);

r_comms_echo = zeros(L, 1);
for p = 1:Pc
    s_shift_c = zeros(size(s_mat));
    for k = 1:Nsym
        s_shift_c(:,k) = circshift(s_mat(:,k), li_c(p));
    end
    r_comms_echo = r_comms_echo + h_c(p) * ...
        (s_shift_c(:) .* exp(1j*2*pi*fi_c(p)*n_vec));
end
r_comms = awgn(r_comms_echo, SNRdB, 'measured');
fprintf('[COMMS CH] %d paths | max delay=%d samp | SNR=%d dB\n\n', Pc, max(li_c), SNRdB);

%% =====================================================================
%  SECTION 8 — CPP REMOVAL
% ======================================================================

% Radar path
rx_mat_r   = reshape(r_radar, N+Ncp, Nsym);
Y_blocks_r = rx_mat_r(Ncp+1:end, :);

% Comms path
rx_mat_c   = reshape(r_comms, N+Ncp, Nsym);
Y_blocks_c = rx_mat_c(Ncp+1:end, :);
Y_daft_c   = daft_colwise_fn(Y_blocks_c, c1, c2, N);

%% =====================================================================
%  SECTION 9 — PILOT-AIDED LS + ZF EQUALISATION
% ======================================================================

% LS at pilot positions
H_pilots = Y_daft_c(pilot_idx, :) ./ X(pilot_idx, :);

% Linear interpolation to all subcarriers
H_full = zeros(N, Nsym);
for k = 1:Nsym
    H_full(:,k) = interp1(pilot_idx.', real(H_pilots(:,k)), (1:N).', ...
                           'linear','extrap') + ...
               1j*interp1(pilot_idx.', imag(H_pilots(:,k)), (1:N).', ...
                           'linear','extrap');
end

% ZF equalisation
Xhat_data  = Y_daft_c(data_idx, :) ./ (H_full(data_idx, :) + 1e-10);

% EVM
X_data_ref = X(data_idx, :);
EVM_rms    = 100 * sqrt(mean(abs(Xhat_data(:) - X_data_ref(:)).^2)) / ...
                   sqrt(mean(abs(X_data_ref(:)).^2));

fprintf('[COMMS] Pilot-aided LS+ZF:\n  EVM = %.2f%%\n', EVM_rms);

%% =====================================================================
%  SECTION 10 — QAM DEMODULATION + BER / SER
% ======================================================================

rx_bits_hat = qamdemod(Xhat_data(:), M, 'OutputType','bit', 'UnitAveragePower',true);
[~, BER_op] = biterr(tx_bits, rx_bits_hat);

rx_sym_hat  = qamdemod(Xhat_data(:), M, 'UnitAveragePower',true);
X_data_vec  = X_data_ref(:);
tx_sym_idx  = qamdemod(X_data_vec, M, 'UnitAveragePower',true);
SER_op      = mean(rx_sym_hat ~= tx_sym_idx);

fprintf('  BER = %.4e | SER = %.4e\n\n', BER_op, SER_op);

%% =====================================================================
%  SECTION 11 — MONTE CARLO BER vs SNR SWEEP
%  AFDM vs OFDM vs OTFS over doubly-selective fading channel.
%  Nmc trials per SNR point for statistical reliability.
%  FIX: pilot subcarriers used for channel estimation in BER loop —
%       no longer using circular estimate that returns x_s directly.
% ======================================================================

fprintf('[BER] Monte Carlo sweep: %d SNR pts x %d trials...\n', ...
        numel(SNR_range), Nmc);

BER_afdm = zeros(size(SNR_range));
BER_ofdm = zeros(size(SNR_range));
BER_otfs = zeros(size(SNR_range));

% Fixed pilot pattern for BER loop (Np=8, same as main sim)
p_idx_ber = 1:Np:N;
d_idx_ber = setdiff(1:N, p_idx_ber);
Nd_ber    = numel(d_idx_ber);

for si = 1:numel(SNR_range)
    snr = SNR_range(si);
    err_a = 0; err_o = 0; err_t = 0;
    nb    = 0;

    for mc = 1:Nmc
        % Fresh random data bits
        b_ber   = randi([0 1], Nd_ber*bps, 1);
        x_data  = qammod(b_ber, M, 'InputType','bit', 'UnitAveragePower',true);

        % Build full symbol vector with pilots
        x_full              = zeros(N, 1);
        x_full(p_idx_ber)   = 1 + 0j;
        x_full(d_idx_ber)   = x_data;

        n_s = (0:N+Ncp-1).';

        %--- AFDM ---
        s_a   = idaft_col_fn(x_full, c1, c2, N);
        cpp_a = [s_a(end-Ncp+1:end); s_a];
        r_a   = zeros(N+Ncp, 1);
        for pp = 1:Pc
            r_a = r_a + h_c(pp)*circshift(cpp_a,li_c(pp)) .* ...
                  exp(1j*2*pi*fi_c(pp)*n_s);
        end
        r_a   = awgn(r_a, snr, 'measured');
        y_a   = r_a(Ncp+1:end);
        Y_a   = daft_colwise_fn(y_a, c1, c2, N);
        % Pilot-aided LS estimate
        H_a   = Y_a(p_idx_ber) ./ x_full(p_idx_ber);
        H_a_f = interp1(p_idx_ber.', real(H_a), (1:N).','linear','extrap') + ...
             1j*interp1(p_idx_ber.', imag(H_a), (1:N).','linear','extrap');
        xh_a  = Y_a(d_idx_ber) ./ (H_a_f(d_idx_ber) + 1e-10);
        b_rx_a = qamdemod(xh_a, M, 'OutputType','bit', 'UnitAveragePower',true);
        err_a  = err_a + sum(b_ber ~= b_rx_a);

        %--- OFDM ---
        s_o   = ifft(x_full) * sqrt(N);
        cp_o  = [s_o(end-Ncp+1:end); s_o];
        r_o   = zeros(N+Ncp, 1);
        for pp = 1:Pc
            r_o = r_o + h_c(pp)*circshift(cp_o,li_c(pp)) .* ...
                  exp(1j*2*pi*fi_c(pp)*n_s);
        end
        r_o   = awgn(r_o, snr, 'measured');
        y_o   = r_o(Ncp+1:end);
        Y_o   = fft(y_o) / sqrt(N);
        H_o   = Y_o(p_idx_ber) ./ x_full(p_idx_ber);
        H_o_f = interp1(p_idx_ber.', real(H_o), (1:N).','linear','extrap') + ...
             1j*interp1(p_idx_ber.', imag(H_o), (1:N).','linear','extrap');
        xh_o  = Y_o(d_idx_ber) ./ (H_o_f(d_idx_ber) + 1e-10);
        b_rx_o = qamdemod(xh_o, M, 'OutputType','bit', 'UnitAveragePower',true);
        err_o  = err_o + sum(b_ber ~= b_rx_o);

        %--- OTFS (correct: delay-Doppler grid, ISFFT + Heisenberg/Wigner pair) ---
        % OTFS frame: N_o delay bins x M_o Doppler bins, N_o*M_o = N
        % Using N_o=N/4, M_o=4 for a minimal 2D grid (fair single-symbol equiv.)
        N_o = N/4; M_o = 4;   % delay bins x Doppler bins
        % Reshape data into delay-Doppler grid
        x_dd_2d = reshape(x_full, N_o, M_o);   % N_o x M_o
        % ISFFT: IDFT along Doppler axis (cols), DFT along delay axis (rows)
        x_tf_2d = fft(ifft(x_dd_2d, [], 2), [], 1);  % N_o x M_o time-freq
        % Heisenberg: serialize to time domain
        s_otfs  = x_tf_2d(:);            % N x 1 time-domain OTFS signal
        cp_t    = [s_otfs(end-Ncp+1:end); s_otfs];
        r_t     = zeros(N+Ncp, 1);
        for pp = 1:Pc
            r_t = r_t + h_c(pp)*circshift(cp_t,li_c(pp)) .* ...
                  exp(1j*2*pi*fi_c(pp)*n_s);
        end
        r_t    = awgn(r_t, snr, 'measured');
        y_t    = r_t(Ncp+1:end);
        % Wigner: reshape and SFFT (DFT along Doppler, IDFT along delay)
        y_tf_2d = reshape(y_t, N_o, M_o);
        Y_dd_2d = ifft(fft(y_tf_2d, [], 1), [], 2);  % back to delay-Doppler
        Y_dd    = Y_dd_2d(:);
        % Pilot-aided LS in delay-Doppler domain
        H_t    = Y_dd(p_idx_ber) ./ x_full(p_idx_ber);
        H_t_f  = interp1(p_idx_ber.', real(H_t), (1:N).','linear','extrap') + ...
              1j*interp1(p_idx_ber.', imag(H_t), (1:N).','linear','extrap');
        xh_t   = Y_dd(d_idx_ber) ./ (H_t_f(d_idx_ber) + 1e-10);
        b_rx_t = qamdemod(xh_t, M, 'OutputType','bit', 'UnitAveragePower',true);
        err_t  = err_t + sum(b_ber ~= b_rx_t);

        nb = nb + numel(b_ber);
    end

    BER_afdm(si) = err_a / nb;
    BER_ofdm(si) = err_o / nb;
    BER_otfs(si) = err_t / nb;
end

% Theoretical AWGN reference
EbN0_lin   = 10.^(SNR_range/10) / bps;
BER_theory = (4/bps)*(1-1/sqrt(M))*qfunc(sqrt(3*bps*EbN0_lin/(M-1)));

fprintf('[BER] Complete.\n\n');

figure('Name','BER vs SNR','Color','w','Position',[50 50 750 520]);
semilogy(SNR_range, BER_afdm,   'b-o','LineWidth',2,'MarkerSize',7, ...
         'DisplayName','AFDM (proposed)');
hold on;
semilogy(SNR_range, BER_ofdm,   'r-s','LineWidth',2,'MarkerSize',7, ...
         'DisplayName','OFDM (baseline)');
semilogy(SNR_range, BER_otfs,   'm-^','LineWidth',2,'MarkerSize',7, ...
         'DisplayName','OTFS (baseline)');
semilogy(SNR_range, BER_theory, '--','Color',[0.2 1.0 0.4],'LineWidth',2.0, ...
         'DisplayName',sprintf('Theoretical %d-QAM AWGN',M));
grid on; xlabel('SNR (dB)'); ylabel('BER');
title(sprintf('BER vs SNR — %d-QAM, N=%d, Doubly-Selective Fading (%d-path)',M,N,Pc));
legend('Location','southwest'); ylim([1e-5 1]); xlim([SNR_range(1) SNR_range(end)]); box on;
text(SNR_range(end)-8, 0.3, sprintf('Nmc=%d trials/pt',Nmc), ...
     'FontSize',9,'Color',[0.4 0.4 0.4]);

%% =====================================================================
%  SECTION 12 — BER vs NORMALISED DOPPLER  (key diversity comparison)
%  Shows AFDM maintaining BER across all Doppler values while OFDM
%  collapses near integer Doppler shifts (nu = 0.5, 1.0, 1.5...).
% ======================================================================

fprintf('[BER-Doppler] Sweeping normalised Doppler nu = 0..2...\n');

nu_range  = 0:0.1:2.0;        % Normalised Doppler nu = fd/DeltaF
BER_afdm_nu = zeros(size(nu_range));
BER_ofdm_nu = zeros(size(nu_range));
snr_fixed   = 10;              % Fixed SNR for this sweep

for ni = 1:numel(nu_range)
    nu_i  = nu_range(ni);
    fi_nu = nu_i * DeltaF * Ts;    % Normalised Doppler for single path

    err_a_nu = 0; err_o_nu = 0; nb_nu = 0;

    for mc = 1:Nmc
        b_nu  = randi([0 1], Nd_ber*bps, 1);
        x_nu  = zeros(N,1);
        x_nu(p_idx_ber) = 1;
        x_nu(d_idx_ber) = qammod(b_nu, M, 'InputType','bit', 'UnitAveragePower',true);
        n_nu  = (0:N+Ncp-1).';

        % Single-path channel with varying Doppler (delay=2 samples)
        li_nu = 2;

        % AFDM
        s_a_nu  = idaft_col_fn(x_nu, c1, c2, N);
        cpp_a_nu = [s_a_nu(end-Ncp+1:end); s_a_nu];
        r_a_nu  = circshift(cpp_a_nu, li_nu) .* exp(1j*2*pi*fi_nu*n_nu);
        r_a_nu  = awgn(r_a_nu, snr_fixed, 'measured');
        y_a_nu  = r_a_nu(Ncp+1:end);
        Y_a_nu  = daft_colwise_fn(y_a_nu, c1, c2, N);
        H_an    = Y_a_nu(p_idx_ber) ./ x_nu(p_idx_ber);
        H_an_f  = interp1(p_idx_ber.', real(H_an),(1:N).','linear','extrap') + ...
               1j*interp1(p_idx_ber.', imag(H_an),(1:N).','linear','extrap');
        xh_an   = Y_a_nu(d_idx_ber) ./ (H_an_f(d_idx_ber) + 1e-10);
        b_an    = qamdemod(xh_an, M, 'OutputType','bit', 'UnitAveragePower',true);
        err_a_nu = err_a_nu + sum(b_nu ~= b_an);

        % OFDM
        s_o_nu  = ifft(x_nu)*sqrt(N);
        cp_o_nu = [s_o_nu(end-Ncp+1:end); s_o_nu];
        r_o_nu  = circshift(cp_o_nu, li_nu) .* exp(1j*2*pi*fi_nu*n_nu);
        r_o_nu  = awgn(r_o_nu, snr_fixed, 'measured');
        y_o_nu  = r_o_nu(Ncp+1:end);
        Y_o_nu  = fft(y_o_nu)/sqrt(N);
        H_on    = Y_o_nu(p_idx_ber) ./ x_nu(p_idx_ber);
        H_on_f  = interp1(p_idx_ber.', real(H_on),(1:N).','linear','extrap') + ...
               1j*interp1(p_idx_ber.', imag(H_on),(1:N).','linear','extrap');
        xh_on   = Y_o_nu(d_idx_ber) ./ (H_on_f(d_idx_ber) + 1e-10);
        b_on    = qamdemod(xh_on, M, 'OutputType','bit', 'UnitAveragePower',true);
        err_o_nu = err_o_nu + sum(b_nu ~= b_on);

        nb_nu = nb_nu + numel(b_nu);
    end

    BER_afdm_nu(ni) = err_a_nu / nb_nu;
    BER_ofdm_nu(ni) = err_o_nu / nb_nu;
end

fprintf('[BER-Doppler] Complete.\n\n');

figure('Name','BER vs Doppler','Color','w','Position',[780 50 720 480]);
semilogy(nu_range, BER_afdm_nu, 'b-o','LineWidth',2,'MarkerSize',6, ...
         'DisplayName','AFDM (proposed)');
hold on;
semilogy(nu_range, BER_ofdm_nu, 'r-s','LineWidth',2,'MarkerSize',6, ...
         'DisplayName','OFDM (baseline)');
% Mark integer Doppler shifts where OFDM is worst
for nu_int = 0.5:0.5:2
    xline(nu_int,'--','Color',[1.0 0.85 0.0],'LineWidth',1.2);
end
grid on; xlabel('Normalised Doppler shift \nu = f_d/\Deltaf'); ylabel('BER');
title(sprintf('BER vs Normalised Doppler — %d-QAM, SNR=%d dB',M,snr_fixed));
legend('Location','northwest'); ylim([1e-4 1]); box on;
text(0.55, 5e-4, 'ICI peaks', 'FontSize',8, 'Color',[1.0 0.85 0.0]);

%% =====================================================================
%  SECTION 13 — CONSTELLATION DIAGRAM
% ======================================================================

Nplot       = min(3000, numel(Xhat_data));
idx_p       = randperm(numel(Xhat_data), Nplot);
X_data_plot = X_data_ref(:);

figure('Name','Constellation','Color','k','Position',[780 560 820 400]);
subplot(1,2,1);
ax1 = gca; ax1.Color = [0.08 0.08 0.08]; ax1.XColor = 'w'; ax1.YColor = 'w';
ax1.GridColor = [0.4 0.4 0.4]; ax1.GridAlpha = 0.5;
plot(real(X_data_plot(idx_p)), imag(X_data_plot(idx_p)), ...
     '.','Color',[0.0 1.0 1.0],'MarkerSize',6); axis([-1.5 1.5 -1.5 1.5]); grid on;
title(sprintf('TX: %d-QAM (DAFT domain)',M),'Color','w');
xlabel('In-Phase (I)','Color','w'); ylabel('Quadrature (Q)','Color','w');

subplot(1,2,2);
ax2 = gca; ax2.Color = [0.08 0.08 0.08]; ax2.XColor = 'w'; ax2.YColor = 'w';
ax2.GridColor = [0.4 0.4 0.4]; ax2.GridAlpha = 0.5;
plot(real(Xhat_data(idx_p)), imag(Xhat_data(idx_p)), ...
     '.','Color',[1.0 0.4 0.1],'MarkerSize',6); axis([-1.5 1.5 -1.5 1.5]); grid on;
title(sprintf('RX: Pilot LS+ZF (EVM=%.1f%%)',EVM_rms),'Color','w');
xlabel('In-Phase (I)','Color','w'); ylabel('Quadrature (Q)','Color','w');
sgtitle('Constellation Diagram','FontWeight','bold','Color','w');

%% =====================================================================
%  SECTION 14 — PHYSICAL SENSING AXES
% ======================================================================

range_bin_m = c_light / (2*Fs);
ranges_plot = (0:Ncp-1) * range_bin_m;
fd_axis     = (-floor(Nsym/2):ceil(Nsym/2)-1) / (Nsym*TAFDM);
vel_axis    = fd_axis * c_light / (2*fc);

fprintf('[SENSING] Range resolution    : %.3f m\n',  range_bin_m);
fprintf('[SENSING] Velocity resolution : %.3f m/s\n', vel_axis(2)-vel_axis(1));
fprintf('[SENSING] Max range (Ncp)     : %.1f m\n',  (Ncp-1)*range_bin_m);
fprintf('[SENSING] Max velocity        : +/-%.2f m/s\n', max(abs(vel_axis)));
fprintf('[SENSING] Processing gain     : %.1f dB\n\n', Gp_dB);

%% =====================================================================
%  SECTION 15 — TIME-DOMAIN FCCR  (Sec. III-A)
% ======================================================================

fprintf('[FCCR] Computing Range-Doppler map...\n');

Rf           = fft(Y_blocks_r, [], 1);
Sf           = fft(s_blocks,   [], 1);
Z_fccr       = ifft(Rf .* conj(Sf), [], 1);
RD_fccr_full = fft(Z_fccr, [], 2);
RD_fccr      = abs(RD_fccr_full(1:Ncp, :));
RD_fccr_show = fftshift(RD_fccr, 2);

%% =====================================================================
%  SECTION 16 — DAFT-DOMAIN SENSING  (Sec. III-B, Eq. 10–11)
% ======================================================================

fprintf('[DAFT] Computing Range-Doppler map...\n');

Y_daft_r  = daft_colwise_fn(Y_blocks_r, c1, c2, N);
p_vec     = (0:N-1).';
RD_daft   = zeros(Ncp, Nsym);
peak_vals = zeros(1, Ncp);

FYr = fft(Y_daft_r, [], 1);
FX  = fft(X,        [], 1);

for l = 0:Ncp-1
    Ll      = exp(1j * 2*pi/N * l * p_vec);
    corr_f  = (Ll .* FYr) .* conj(FX);
    Wl_time = ifft(corr_f, [], 1);
    Wl      = fftshift(fft(Wl_time, [], 2), 2);
    RD_daft(l+1, :) = max(abs(Wl), [], 1);
    peak_vals(l+1)  = max(abs(Wl(:)));
end

%% =====================================================================
%  SECTION 17 — CFAR DETECTION ON DAFT RDM
%  Cell-Averaging CFAR: compares each cell to average of surrounding
%  guard+training cells. Gives clean automated target detection.
% ======================================================================

fprintf('[CFAR] Applying CA-CFAR detector...\n');

Pfa        = 1e-3;           % Target probability of false alarm
Ntrain     = 8;              % Training cells each side (range & vel)
Nguard     = 2;              % Guard cells each side
% Correct 2D CA-CFAR training cell count:
% Total window = (2*Ntrain+2*Nguard+1)^2, minus guard+CUT = (2*Nguard+1)^2
Ncfar      = (2*Ntrain+2*Nguard+1)^2 - (2*Nguard+1)^2;
alpha_cfar = Ncfar * (Pfa^(-1/Ncfar) - 1);  % CFAR threshold multiplier

RD_daft_cfar = fftshift(RD_daft, 2);
RDM_power    = RD_daft_cfar.^2;
[nR, nV]     = size(RDM_power);
CFAR_mask    = false(nR, nV);

for ri = Ntrain+Nguard+1 : nR-Ntrain-Nguard
    for vi = Ntrain+Nguard+1 : nV-Ntrain-Nguard
        % Training region (exclude guard)
        r_lo_t = ri-Ntrain-Nguard; r_hi_t = ri+Ntrain+Nguard;
        v_lo_t = vi-Ntrain-Nguard; v_hi_t = vi+Ntrain+Nguard;
        region = RDM_power(r_lo_t:r_hi_t, v_lo_t:v_hi_t);
        % Zero out guard region
        region(Ntrain+1:Ntrain+2*Nguard+1, Ntrain+1:Ntrain+2*Nguard+1) = 0;
        noise_est = sum(region(:)) / Ncfar;
        threshold = alpha_cfar * noise_est;
        if RDM_power(ri, vi) > threshold
            CFAR_mask(ri, vi) = true;
        end
    end
end

% Connected components: extract detected target locations
CC = bwconncomp(CFAR_mask);
detected_R = zeros(1, CC.NumObjects);
detected_V = zeros(1, CC.NumObjects);
for t = 1:CC.NumObjects
    [ri_c, vi_c] = ind2sub([nR nV], CC.PixelIdxList{t});
    detected_R(t) = ranges_plot(round(mean(ri_c)));
    detected_V(t) = vel_axis(round(mean(vi_c)));
end

fprintf('  CFAR detected %d object(s):\n', CC.NumObjects);
for t = 1:CC.NumObjects
    fprintf('    R=%.1f m | v=%.2f m/s\n', detected_R(t), detected_V(t));
end
fprintf('\n');

%% =====================================================================
%  SECTION 18 — PSLR COMPUTATION
% ======================================================================

fprintf('[PSLR] Computing...\n');
guard_r   = 2; guard_v = 2;
PSLR_fccr = compute_pslr_fn(RD_fccr_show, R_true, v_true, ...
                              ranges_plot, vel_axis, guard_r, guard_v);
PSLR_daft = compute_pslr_fn(RD_daft_cfar, R_true, v_true, ...
                              ranges_plot, vel_axis, guard_r, guard_v);
fprintf('  PSLR-FCCR = %.2f dB\n  PSLR-DAFT = %.2f dB\n\n', PSLR_fccr, PSLR_daft);

%% =====================================================================
%  SECTION 19 — IMAGE SNR vs INPUT SNR
% ======================================================================

fprintf('[ImgSNR] Computing...\n');
SNR_img_range = -20:5:20;
ImgSNR_fccr   = zeros(size(SNR_img_range));
ImgSNR_daft   = zeros(size(SNR_img_range));

for si = 1:numel(SNR_img_range)
    snr_i = SNR_img_range(si);
    re_i  = zeros(L,1);
    for p = 1:P
        s_sh = zeros(size(s_mat));
        for k = 1:Nsym, s_sh(:,k) = circshift(s_mat(:,k), li(p)); end
        re_i = re_i + h_i(p)*(s_sh(:).*exp(1j*2*pi*fi(p)*n_vec));
    end
    r_i  = awgn(re_i, snr_i, 'measured');
    Yb_i = reshape(r_i, N+Ncp, Nsym); Yb_i = Yb_i(Ncp+1:end,:);

    % FCCR image SNR
    Rf_i = fft(Yb_i,[],1);
    RD_i = abs(fftshift(fft(ifft(Rf_i.*conj(Sf),[],1),[],2),2));
    RD_i = RD_i(1:Ncp,:);
    [~,ri1] = min(abs(ranges_plot-R_true(1)));
    [~,vi1] = min(abs(vel_axis-v_true(1)));
    pk_f = RD_i(ri1,vi1)^2;
    msk  = true(size(RD_i));
    msk(max(1,ri1-2):min(end,ri1+2),max(1,vi1-2):min(end,vi1+2)) = false;
    ImgSNR_fccr(si) = 10*log10(pk_f/(mean(RD_i(msk).^2)+eps));

    % DAFT image SNR (at known delay l=li(1))
    Yd_i = daft_colwise_fn(Yb_i, c1, c2, N);
    FY_i = fft(Yd_i,[],1);
    Ll1  = exp(1j*2*pi/N*li(1)*p_vec);
    Wl1  = fftshift(fft(ifft((Ll1.*FY_i).*conj(FX),[],1),[],2),2);
    pk_d = max(abs(Wl1(:)))^2;
    ImgSNR_daft(si) = 10*log10(pk_d/(mean(abs(Wl1(:)).^2)+eps));
end
fprintf('[ImgSNR] Done.\n\n');

figure('Name','Image SNR','Color','w','Position',[50 50 680 470]);
plot(SNR_img_range, ImgSNR_fccr, 'b-o','LineWidth',2,'MarkerSize',7, ...
     'DisplayName','AFDM-FCCR');
hold on;
plot(SNR_img_range, ImgSNR_daft, 'r-^','LineWidth',2,'MarkerSize',7, ...
     'DisplayName','AFDM-DAFT-domain');
grid on; xlabel('Input SNR (dB)'); ylabel('Image SNR (dB)');
title('Image SNR vs Input SNR  (reproduces paper Fig. 4)');
legend('Location','northwest'); box on;

%% =====================================================================
%  SECTION 20 — PARAMETER EXTRACTION (Eq. 12–16)
% ======================================================================

fprintf('[DAFT] Extracting target parameters...\n');
threshold         = 0.25 * max(peak_vals);
estimated_targets = [];

for l = 0:Ncp-1
    if peak_vals(l+1) > threshold
        Ll   = exp(1j*2*pi/N*l*p_vec);
        Wl   = fftshift(fft(ifft((Ll.*FYr).*conj(FX),[],1),[],2),2);
        [~,max_idx]    = max(abs(Wl(:)));
        [p_bar, k_bar] = ind2sub([N,Nsym], max_idx);

        l_hat     = l;
        R_hat     = l_hat * range_bin_m;
        alpha_hat = mod(2*N*c1*l_hat - (1-p_bar), N);
        if alpha_hat > N/2, alpha_hat = alpha_hat - N; end
        b_hat    = (k_bar - 1 - Nsym/2) / Nsym;
        ratio    = N/(N+Ncp);
        beta_min = ceil( ratio*(alpha_hat-b_hat) - 0.5*ratio );
        beta_max = floor(ratio*(alpha_hat-b_hat) + 0.5*ratio );

        if beta_min == beta_max
            beta_hat = beta_max;
        else
            p_prev = mod(p_bar-2,N)+1; p_next = mod(p_bar,N)+1;
            if abs(Wl(p_prev,k_bar)) > abs(Wl(p_next,k_bar))
                beta_hat = beta_min;
            else
                beta_hat = beta_max;
            end
        end

        a_hat  = ratio*beta_hat + ratio*b_hat - alpha_hat;
        nu_hat = alpha_hat + a_hat;
        fd_hat = nu_hat * DeltaF;
        v_hat  = fd_hat * c_light / (2*fc);
        estimated_targets(end+1,:) = [R_hat, v_hat, l_hat, fd_hat]; %#ok
    end
end

%% =====================================================================
%  SECTION 21 — ESTIMATION ERROR + DOPPLER AMBIGUITY TABLE
% ======================================================================

fprintf('\n--- Ground Truth vs Estimates ---\n');
fprintf('  %-6s %-9s %-10s %-9s %-10s %-9s %-9s\n', ...
        'Target','R_true','v_true','R_est','v_est','dR(m)','dv(m/s)');

R_err = zeros(1,P); V_err = zeros(1,P);
R_est_all = zeros(1,P); v_est_all = zeros(1,P);

for p = 1:P
    if ~isempty(estimated_targets)
        [~,best]     = min(abs(estimated_targets(:,1) - R_true(p)));
        R_err(p)     = abs(estimated_targets(best,1) - R_true(p));
        V_err(p)     = abs(estimated_targets(best,2) - v_true(p));
        R_est_all(p) = estimated_targets(best,1);
        v_est_all(p) = estimated_targets(best,2);
        fprintf('  T%-5d %-9.1f %-10.1f %-9.1f %-10.2f %-9.2f %-9.2f\n', ...
                p, R_true(p), v_true(p), R_est_all(p), v_est_all(p), ...
                R_err(p), V_err(p));
    end
end
fprintf('\n  Mean dR=%.3f m  |  Mean dv=%.3f m/s\n', mean(R_err), mean(V_err));

v_max_td   = 1/(2*TAFDM) * c_light/(2*fc);
v_max_daft = (alphamax+kv+0.5) * DeltaF * c_light/(2*fc);
fprintf('\n--- Doppler Ambiguity (paper Table II) ---\n');
fprintf('  FCCR max velocity : +/-%.2f m/s\n', v_max_td);
fprintf('  DAFT max velocity : +/-%.2f m/s  (%.1fx extension)\n\n', ...
        v_max_daft, v_max_daft/v_max_td);

% Error bar chart
dv_res = vel_axis(2)-vel_axis(1);
figure('Name','Estimation Errors','Color','w','Position',[50 560 820 420]);
subplot(1,2,1);
bar(1:P, R_err, 'FaceColor',[0.2 0.55 0.9]); grid on; hold on;
yline(range_bin_m,'r--','LineWidth',1.5,'DisplayName','\DeltaR resolution');
xlabel('Target'); ylabel('Range Error (m)'); title('Range Estimation Error');
xticks(1:P); xticklabels(arrayfun(@(i)sprintf('T%d',i),1:P,'UniformOutput',false));
legend('show','Location','northwest');

subplot(1,2,2);
bar(1:P, V_err, 'FaceColor',[0.9 0.4 0.2]); grid on; hold on;
yline(abs(dv_res),'r--','LineWidth',1.5,'DisplayName','\Deltav resolution');
xlabel('Target'); ylabel('Velocity Error (m/s)'); title('Velocity Estimation Error');
xticks(1:P); xticklabels(arrayfun(@(i)sprintf('T%d',i),1:P,'UniformOutput',false));
legend('show','Location','northwest');
sgtitle(sprintf('Estimation Errors  |  PSLR-FCCR=%.1f dB  |  PSLR-DAFT=%.1f dB', ...
        PSLR_fccr, PSLR_daft),'FontWeight','bold');

%% =====================================================================
%  SECTION 22 — RANGE-DOPPLER MAP VISUALISATION
% ======================================================================

[VV, RR]   = meshgrid(vel_axis, ranges_plot);
[~,kidxT]  = arrayfun(@(v) min(abs(vel_axis-v)), v_true);
ridxT      = min(max(round(R_true/range_bin_m)+1,1), Ncp);

% FCCR 2D
figure('Name','RDM FCCR 2D','Color','w','Position',[780 50 720 500]);
imagesc(vel_axis, ranges_plot, mag2db_norm_fn(RD_fccr_show));
axis xy; colormap turbo; colorbar;
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title(sprintf('FCCR Range-Doppler Map  (Sec. III-A)  —  PSLR=%.1f dB',PSLR_fccr));
hold on;
scatter(v_true, R_true, 140,'r','filled','MarkerEdgeColor','k','LineWidth',1.5,...
        'DisplayName','True targets');
legend('Location','northoutside'); grid on; box on;

% FCCR 3D
figure('Name','RDM FCCR 3D','Color','w','Position',[780 560 720 480]);
Zf = mag2db_norm_fn(RD_fccr_show);
surf(VV, RR, Zf, 'EdgeColor','none'); view(40,30);
xlabel('Velocity (m/s)'); ylabel('Range (m)'); zlabel('Magnitude (dB)');
title('3D Range-Doppler: FCCR'); colormap turbo; colorbar; grid on; box on; hold on;
zT = arrayfun(@(ii) Zf(ridxT(ii),kidxT(ii)),1:P);
scatter3(v_true, R_true, zT, 120,'r','filled','MarkerEdgeColor','k','LineWidth',1.5);

% DAFT-domain 2D with CFAR overlay
figure('Name','RDM DAFT 2D','Color','w','Position',[50 50 720 500]);
imagesc(vel_axis, ranges_plot, mag2db_norm_fn(RD_daft_cfar));
axis xy; colormap turbo; colorbar;
xlabel('Velocity (m/s)'); ylabel('Range (m)');
title(sprintf('DAFT-Domain RDM  (Sec. III-B)  —  PSLR=%.1f dB',PSLR_daft));
hold on;
scatter(v_true, R_true, 140,'r','filled','MarkerEdgeColor','k','LineWidth',1.5,...
        'DisplayName','True targets');
% CFAR detections
if ~isempty(detected_R)
    scatter(detected_V, detected_R, 120,'g','filled','MarkerEdgeColor','k', ...
            'LineWidth',1.2,'DisplayName','CFAR detections');
end
legend('Location','northoutside'); grid on; box on;

% DAFT-domain 3D
figure('Name','RDM DAFT 3D','Color','w','Position',[50 560 720 480]);
Zd = mag2db_norm_fn(RD_daft_cfar);
surf(VV, RR, Zd, 'EdgeColor','none'); view(40,30);
xlabel('Velocity (m/s)'); ylabel('Range (m)'); zlabel('Magnitude (dB)');
title('3D Range-Doppler: DAFT-Domain'); colormap turbo; colorbar; grid on; box on; hold on;
zDT = arrayfun(@(ii) Zd(ridxT(ii),kidxT(ii)),1:P);
scatter3(v_true, R_true, zDT, 120,'r','filled','MarkerEdgeColor','k','LineWidth',1.5);

%% =====================================================================
%  SECTION 23 — PILOT STRUCTURE VISUALISATION
% ======================================================================

figure('Name','Pilot Structure','Color','w','Position',[50 50 700 420]);
imagesc(1:Nsym, 1:N, abs(X));
colormap hot; colorbar;
xlabel('Symbol index (slow time)'); ylabel('Subcarrier index (DAFT domain)');
title(sprintf('DAFT-Domain Frame: Pilots (white) + Data  |  Overhead=%.1f%%', ...
      100*Npilots/N));
axis tight; box on;

%% =====================================================================
%  SECTION 24 — SUMMARY
% ======================================================================

fprintf('=================================================================\n');
fprintf('  AFDM-ISAC v3 — Complete\n');
fprintf('  Waveform: AFDM | N=%d | Nsym=%d | %d-QAM | CPP=%d\n', N,Nsym,M,Ncp);
fprintf('  ----- Communications -----\n');
fprintf('    BER  = %.4e | SER = %.4e\n', BER_op, SER_op);
fprintf('    EVM  = %.2f%%\n', EVM_rms);
fprintf('    Spectral Efficiency = %.3f bits/s/Hz\n', SE_bpshz);
fprintf('  ----- Sensing -----\n');
fprintf('    PSLR-FCCR = %.2f dB | PSLR-DAFT = %.2f dB\n', PSLR_fccr, PSLR_daft);
fprintf('    Mean dR   = %.3f m  | Mean dv   = %.3f m/s\n', mean(R_err), mean(V_err));
fprintf('    CFAR detections: %d\n', CC.NumObjects);
fprintf('    Max vel FCCR: +/-%.2f m/s | DAFT: +/-%.2f m/s\n', v_max_td, v_max_daft);
fprintf('  ----- Figures -----\n');
fprintf('    1. BER vs SNR  (AFDM / OFDM / OTFS, %d-trial MC)\n', Nmc);
fprintf('    2. BER vs Normalised Doppler\n');
fprintf('    3. Constellation diagram\n');
fprintf('    4. Image SNR vs Input SNR\n');
fprintf('    5. FCCR RDM 2D + 3D\n');
fprintf('    6. DAFT RDM 2D + CFAR overlay + 3D\n');
fprintf('    7. Estimation error bar charts\n');
fprintf('    8. Pilot structure\n');
fprintf('=================================================================\n');

%% =========================================================================
%  LOCAL FUNCTIONS  (must be at end of file in MATLAB)
% ==========================================================================

function s = idaft_col_fn(x, c1, c2, N)
% IDAFT  Eq.(1): s[n] = exp(j2pi*c1*n^2) * IFFT{exp(j2pi*c2*m^2)*x[m]} * sqrt(N)
n = (0:N-1).'; m = (0:N-1).';
s = exp(1j*2*pi*c1*(n.^2)) .* (ifft(exp(1j*2*pi*c2*(m.^2)) .* x) * sqrt(N));
end

function Y = daft_colwise_fn(S, c1, c2, N)
% DAFT  Eq.(5): Y[m] = exp(-j2pi*c2*m^2) * FFT{exp(-j2pi*c1*n^2)*s[n]} / sqrt(N)
n = (0:N-1).'; m = (0:N-1).';
Y = (fft(exp(-1j*2*pi*c1*(n.^2)) .* S, [], 1) / sqrt(N)) .* ...
     exp(-1j*2*pi*c2*(m.^2));
end

function ydb = mag2db_norm_fn(Y)
Yn = Y ./ (max(Y(:)) + eps);
ydb = 20*log10(Yn + eps);
end

function pslr = compute_pslr_fn(RDM, true_ranges_m, true_vels_ms, ...
                                  range_ax, vel_ax, guard_r, guard_v)
% PSLR: average over all targets of peak/sidelobe ratio in dB
pslr_vals = zeros(1, numel(true_ranges_m));
for t = 1:numel(true_ranges_m)
    [~,ri] = min(abs(range_ax - true_ranges_m(t)));
    [~,vi] = min(abs(vel_ax   - true_vels_ms(t)));
    peak_power = abs(RDM(ri,vi))^2;
    mask = true(size(RDM));
    mask(max(1,ri-guard_r):min(end,ri+guard_r), ...
         max(1,vi-guard_v):min(end,vi+guard_v)) = false;
    sidelobe_power = mean(abs(RDM(mask)).^2);
    pslr_vals(t) = 10*log10(peak_power / (sidelobe_power + eps));
end
pslr = mean(pslr_vals);
end