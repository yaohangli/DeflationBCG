clear all
close all
clc

%% setup
% rng(114773)
n = 2000;
A = rand(n);
A = A'*A;

[rows, cols] = size(A);

s = 20;

[realU,realD,realV] = svd(full(A));
realD = diag(realD);

M = eye(rows);

Omega = randn(rows,s-1);
b = rand(rows,1);
B = [b,Omega];


iters = 600;
cols_W = s-1;
errR = zeros(10000,4);
%% BCG
nstep =  0;
Xj = B;      % initial guess
Rj = B - A*Xj;  % residual space
Zj = M*Rj;     % preconditioned residual space
Pj = orth(Zj);   % search space
%
for i = 1:iters
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_bcg(A,B,Pj,Zj,Rj,Xj,M);   % a single step in Block CG
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,1) = norm(RjTemp)/norm(b);
end

%%  deflated BCG
nstep = 0;

Xj = B;
Rj = B - A*Xj;  % residual space

W = realU(:, end-cols_W+1:end);
Xj = Xj + W*((W'*A*W)\(W'*Rj));
Rj = B - A*Xj;
Zj = M*Rj;
Pj = orth(Zj - W*((W'*A*W)\(W'*A*Zj)));

for i = 1:iters
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W);
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,2) = norm(RjTemp)/norm(b);
end


%%  rand BCG
nstep = 0;

Xj = B;
Rj = B - A*Xj;  % residual space

[tempU,tempD,tempV] = svd(B);
W = tempU(:, end-cols_W+1:end);
Xj = Xj + W*((W'*A*W)\(W'*Rj));
Rj = B - A*Xj;
Zj = M*Rj;
Pj = orth(Zj - W*((W'*A*W)\(W'*A*Zj)));

for i = 1:iters
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W);
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,3) = norm(RjTemp)/norm(b);
end


%% adaptive BCG

nstep =  0;
Xj = B;      % initial guess
Rj = B - A*Xj;  % residual space
Zj = M*Rj;     % preconditioned residual space
Pj = orth(Zj);   % search space
%
for i = 1:iters/4
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_bcg(A,B,Pj,Zj,Rj,Xj,M);   % a single step in Block CG
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,4) = norm(RjTemp)/norm(b);
end

Q = orth(Xj(:,end-s+2:end));
B2 = [b,Q];
% Xj = B2;      % initial guess
Xj =[Xj(:,1),Q];
Rj = B2 - A*Xj;  % residual space

W = Q;
Xj = Xj + W*((W'*A*W)\(W'*Rj));
Rj = B2 - A*Xj;
Zj = M*Rj;
Pj = orth(Zj - W*((W'*A*W)\(W'*A*Zj)));

for i = 1:iters/4
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W);
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,4) = norm(RjTemp)/norm(b);
end

Q = orth(Xj(:,end-s+2:end));
B2 = [b,Q];
% Xj = B2;      % initial guess
Xj =[Xj(:,1),Q];
Rj = B2 - A*Xj;  % residual space

W = Q;
Xj = Xj + W*((W'*A*W)\(W'*Rj));
Rj = B2 - A*Xj;
Zj = M*Rj;
Pj = orth(Zj - W*((W'*A*W)\(W'*A*Zj)));

for i = 1:iters/4
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W);
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,4) = norm(RjTemp)/norm(b);
end

Q = orth(Xj(:,end-s+2:end));
B2 = [b,Q];
% Xj = B2;      % initial guess
Xj =[Xj(:,1),Q];
Rj = B2 - A*Xj;  % residual space

W = Q;
Xj = Xj + W*((W'*A*W)\(W'*Rj));
Rj = B2 - A*Xj;
Zj = M*Rj;
Pj = orth(Zj - W*((W'*A*W)\(W'*A*Zj)));

for i = 1:iters/4
    i
    nstep = nstep +1 ;
    [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W);
    RjTemp = b - A*Xj(:,1);  % residual space
    errR(nstep,4) = norm(RjTemp)/norm(b);
end

%%
figure
semilogy(errR);