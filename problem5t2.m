% Problem5
clear vars;
close all;

% a : A = Q0 * A1 * Q0'
% If A1 = R0 * Q0 then R0 =A1 * Q0'
% If A = Q0 * R0, then if we substitute R0,
% the resulting equation is :
% A = Q0 * A1 * Q0'

% b : A = (Q0*Q1) * A2 * (Q0*Q1)'
% If A2 = R1 * Q1, then R1 = A2 * Q1' 
% Therefore if we substitute R1, 
% A1 = Q1 * A2 * Q1' 
% If A1 = R0 * Q0, then Q1 * A2 * Q1' = R0 * Q0 
% and R0 = Q1 * A2 * Q1' * Q0' 
% If we substitute R0 in the equation for A, 
% A = Q0 * R0 = Q0 * Q1 *A2 * Q1' * Q0' 
% Therefore A = (Q0 * Q1) * A2 * (Q0 * Q1)' 

% c : Q0*Q1 is orthogonal
% Q0 is orthogonal by definition of QR factorization
% Q1 is orthogonal by definition of QR factorization
% Q0'*Q0 = I Q1'*Q1 = I
% (Q0 * Q1)' * Q0 * Q1 = Q1' * Q0' * Q0 * Q1
%  = Q1' * I * Q1 = Q1' * Q1 = I 
% Therefore Q0 * Q1 is orthogonal

% d : A, A1, A2 all have the same eigenvalues
% From the resulting transforms in a, b, and c, we know that matrices 
% A, A1, and A2 are similar. If the 3 matrices are similar, then they must
% have the same eigenvalues 

% e & f

fprintf('M1:\n')
M1 = [1 -2 8; 
      7 -7 6; 
      5 7 -8];

tol = 1e-5;
N = 100;
A = M1;
for ii = 1:N
    [Q, R] = qr(A);
    A = R * Q;
    if abs(A(2, 1)) <= tol && ...
            abs(A(3, 1)) <= tol && ...
            abs(A(3, 2)) <= tol
        disp('converged')
        break
    end
end

%converged
%number of iterations 97
%eigen values:
%-13.051159
%-6.767949
%5.819108

fprintf('M2:\n')
M2 = [ 4 -2 3 -7; 1 2 6 8; 8 5 1 -5; -5 8 -5 3];
tol = 1e-5;
N = 300;
A = M2;
for ii = 1:N
    [Q, R] = qr(A);
    A = R * Q;
    if abs(A(2, 1)) <= tol && ...
            abs(A(3, 1)) <= tol && ...
            abs(A(3, 2)) <= tol && ...
            abs(A(4, 1)) <= tol && ...
            abs(A(4, 2)) <= tol && ...
            abs(A(4, 3)) <= tol
        disp('converged')
        break
    end
end

% converged
%number of iterations 107
%eigen values:
%13.829866
%-10.245535
%8.916750
%-2.501082
num_iterations = ii;
if ii == N
    disp('did not converge')
end
eigen_values = diag(A);
fprintf('number of iterations %d\n', num_iterations)
fprintf('eigen values:\n')
fprintf('%f\n', eigen_values)

