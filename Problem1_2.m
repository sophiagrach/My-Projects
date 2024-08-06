% Problem1
%a 
t = [1900; 1910; 1920; 1930; 1940; 1950; 1960; 1970; 1980; 1990];
T = length(t);
A = [ones(T,1), t, t.^2, t.^3];


b = [75.995; 91.972; 105.711; 123.203; 131.669; 150.697; 179.323; 203.212; 226.505; 249.633];

%b
% A'Ax = A'b
x = (A' * A)\ A' * b;

tt = (1900:1990)';
pt = x(1) + x(2) * tt +  x(3) * tt.^2 + x(4) * tt.^3;

figure(1)
plot(tt, pt, t, b, '-o');

%c
% Prove that x'A'Ax>= for any x
% Say Ax = y
% y'y = norm(y) >= for any y
% This implies that x'A'Ax>= for any x
rank(A)
% rank(A) = 4
L = chol(A' * A, 'lower');
z = A' * b ; 
y = L \ z; 
x2 = L' \ y; 
 

%d
[Q,R] = qr(A);
%Rx = Q'b
bQ = Q'*b;
x3 = R(1:4, :) \ bQ(1:4);

figure(4)
tt = (1900:1990)';
pz = x3(1) + x3(2) * tt +  x3(3) * tt.^2 + x3(4) * tt.^3;
plot(tt, pz, t, b, '-o')

%Problem2
%a
t = [1900; 1910; 1920; 1930; 1940; 1950; 1960; 1970; 1980; 1990];
T = length(t);
x = (t - 1945) / 45; 

p0 = x.^0;
p1 = x;
p2 = 0.5 * (3*x.^2 - 1);
p3 = 0.5 * (5*x.^3 - 3*x);

P = [ones(T, 1), p1, p2, p3];

%The resulting polynomials are orthogonal because the dot products of the
%vectors are all 0

%b
% P'Py= P'b
y = (P' * P) \(P' *b);
% y = [152.4013 86.6472 12.5160 0.2035]'
% y is much larger than x
tt = ((1900-1945)/45:(1990-1945)/45)';
py = y(1) + y(2) * tt +  y(3) * tt.^2 + y(4) * tt.^3;

figure(2)
plot(tt, py, t, b, '-o')

%c
[Q,R] = qr(P);
bQ = Q'*b;
y2 = R(1:4, :) \ bQ(1:4);

% y2 = [152.4013 86.6472 12.516 0.2035]'
tt = ((1900-1945)/45:(1990-1945)/45)';
px = y2(1) + y2(2) * tt +  y2(3) * tt.^2 + y2(4) * tt.^3;

figure(3)
plot(tt, px, t, b, '-o')

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
