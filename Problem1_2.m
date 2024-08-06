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

