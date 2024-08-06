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
