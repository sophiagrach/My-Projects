%Problem 3
close all
clear vars
% fourier transform defs
% h(w) = int[h(t)*exp(-iwt)*dt]
% h(t) = 1/2pi * int[ H(w) * exp(iwt) dw ] - inverse
close all;
clearvars;
N = 256;
dx = 2 * pi / N;
x = (0:N-1)' * dx;
w_n = (0:N-1)'; % freq
sn = 0.1;
dz = sn*randn(N,1);
y = sin(5*x);
z = y + dz;
figure(1);
title('y and z');
plot(x, y, x, z);
grid on;
% fourier transform / spectrum of y
y_w = fft(y);
y_abs_dx_plot = abs(y_w) * dx / (2*pi);  %  plot amplitudes (normalized)
z_w = fft(z);
z_abs_dx_plot = abs(z_w) * dx / (2*pi);  %  plot amplitudes (normalized)
figure(2);
stem(w_n, z_abs_dx_plot); %
title('spectrum')
hold on
stem(w_n, y_abs_dx_plot); % freq at 5 and 256-5, 256-5 represents freq -5
text(50, 0.45,  'freq 5 and 256-5. Interpret freq 256-5 as -5')
grid on;
z_w_2 = thresh(z_w, 2);
z_w_3 = thresh(z_w, 3);
z_w_4 = thresh(z_w, 4);
s_2 = ifft(z_w_2);
s_3 = ifft(z_w_3);
s_4 = ifft(z_w_4);
figure(3)
plot(x, z, x, s_2, x, s_3, x, s_4); % cannot see well
grid on
legend('z', 's_2', 's_3', 's_4')
figure(4)
subplot(4, 1, 1);
plot(x, z); grid on;
title('z')
subplot(4, 1, 2);
plot(x, s_2); grid on;
title('s_2')
subplot(4, 1, 3);
plot(x, s_3); grid on;
title('s_3')
subplot(4, 1, 4);
plot(x, s_4); grid on;
title('s_4')
% zoom into filtered s_4
figure(5)
plot(x, s_4, x, y); grid on;
legend('s_4', 'y')
% conclusions
% original y gets progreessively better extracted with the increased threshold
% (in other words noise better filtered as we increase the threshold)
% (s2, s3, s4 are progressively closer to the original y)
% however original y gets distorted - as one can see from fig 5 -
% s_4 > 1 at peak 2, 4 and less s4 < 1, i.e. we don't fully recover
% original y



