% test for getting an idea of the lambda_min of CC^T, where C is the
% concatenation matrix resulting from the newton sketch method

clear all

n = 4;
s=2;
C1 = zeros(s,n);
C1(1) = 1;
C2 = zeros(s,n);
C2(2) = 1;
C3 = zeros(s,n);
C3(3) = 1;
C4 = zeros(s,n);
C4(4) = 1;

C = C1'*C1+C2'*C2+C3'*C3+C4'*C4;
C1'*C1