function [Pj,Zj,Rj,Xj] = singleStep_deflation_bcg(A,Pj,Zj,Rj,Xj,M,W)

[n,r] = size(Pj);
% disp(['The size of Pj: ','[',num2str(n),',',num2str(r),']']);
Qj = A*Pj;

alphaj = (Pj'*Qj)\(Pj'*Rj);

Xj1 = Xj + Pj*alphaj;
Rj1 = Rj - Qj*alphaj;
Rj1 = Rj1 - W*((W'*W)\(W'*Rj1));  % deflation

Zj1 = M*Rj1;

betaj = -(Pj'*Qj)\(Qj'*Zj1);
% betaj = (Rj'*Zj)\(Rj1'*Zj1);


Pj1 = orth(Zj1 + Pj*betaj - W*((W'*A*W)\(W'*A*Zj1)));

Xj = Xj1;
Rj = Rj1;
Zj = Zj1;
Pj = Pj1;