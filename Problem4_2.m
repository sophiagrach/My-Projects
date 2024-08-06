close all
clearvars
%% Read the audio signal x and the sampling rate fs
[x, fs] = audioread('/Users/sophiagrach/Documents/MATLAB/bison_and_marmot 3.wav');
[n, ~] = size(x);
t = (1:n)'*(1/fs);
%% Plot the left (x(:,1)) and right channel (x(:,2)) of the audio signal
figure(1)
subplot(2,1,1);
plot(t, x(:, 1))
grid on;
title("The bison and the marmot: left channel");
xlabel("Time [s]");
axis tight;
subplot(2,1,2);
plot(t, x(:, 2))
grid on;
title("The bison and the marmot: right channel");
xlabel("Time [s]");
axis tight;
saveas(gcf,"bison_and_marmot","pdf");
%% Plot the DFT of the left and right channel of the audio signal
dt = 1 / fs;
fHz = [0:n/2-1, n/2, -n/2+1:-1]' * (fs / n); % freq vector
fkHz = fHz / 1000;
fft_l = fft(x(:, 1));
fft_r = fft(x(:, 2));
abs_fft_l = abs(fft_l) * dt;
abs_fft_r = abs(fft_r) * dt;
figure(2)
subplot(2,1,1);
plot(fkHz, abs_fft_l);
grid on;
title("The bison and the marmot: left channel frequency");
xlabel("Frequency [kHz]");  
axis tight;
subplot(2,1,2);
plot(fkHz, abs_fft_r);
grid on;
title("The bison and the marmot: right channel frequency");
xlabel("Frequency [kHz]");
axis tight;
saveas(gcf,"bison_and_marmot_frequency","pdf");
%% Set-up the low-pass filter, separate bison from marmot
% d) - h)
f_cut = 1.3 * 1000;
% bison
idx_b = find(abs(fHz) <= f_cut);
fft_b_l = zeros(n, 1);
fft_b_r = zeros(n, 1);
fft_b_l(idx_b) = fft_l(idx_b);
fft_b_r(idx_b) = fft_r(idx_b);
% marmot
idx_m = find(abs(fHz) > f_cut);
fft_m_l = zeros(n, 1);
fft_m_r = zeros(n, 1);
fft_m_l(idx_m) = fft_l(idx_m);
fft_m_r(idx_m) = fft_r(idx_m);
% plots
figure(3)
subplot(2,1,1);
plot(fkHz, abs(fft_b_l)*dt);
grid on;
title("The bison: left channel frequency");
xlabel("Frequency [kHz]");  
axis tight;
subplot(2,1,2);
plot(fkHz, abs(fft_m_l)*dt);
grid on;
title("The marmot: left channel frequency");
xlabel("Frequency [kHz]");  
axis tight;
% inverse ft, revover sounds
x_b_l = ifft(fft_b_l);
x_b_r = ifft(fft_b_r);
x_m_l = ifft(fft_m_l);
x_m_r = ifft(fft_m_r);
% write the the recovered signals to .wav files
audiowrite("bison_recovered_CT.wav", [x_b_l, x_b_r], fs);
audiowrite("marmot_recovered_CT.wav",[x_m_l, x_m_r], fs);
% plots
figure(4)
subplot(2,1,1);
plot(fkHz, x_b_l)
title("American bison - snorts: left channel (recovered)");
xlabel("Time [s]");
axis tight;
subplot(2,1,2);
plot(fkHz, x_b_r)
title("American bison - snorts: right channel (recovered)");
xlabel("Time [s]");
axis tight;
saveas(gcf,"bison_recovered","pdf");
% Plot the left and right channel of the audio signal containing the marmot
figure(5)
subplot(2,1,1);
plot(fkHz, x_m_l)
title("Alpine marmot - close-up shrill call: left channel (recovered)");
xlabel("Time [s]");
axis tight;
subplot(2,1,2);
plot(fkHz, x_m_r)
title("Alpine marmot - close-up shrill call: right channel (recovered)");
xlabel("Time [s]");
axis tight;
saveas(gcf,"marmot_recovered","pdf");