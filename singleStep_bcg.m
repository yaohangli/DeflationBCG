function [Pj,Zj,Rj,Xj] = singleStep_bcg(A,B,Pj,Zj,Rj,Xj,M)

[n,r] = size(Pj);
% disp(['The size of Pj: ','[',num2str(n),',',num2str(r),']']);
Qj = A*Pj;

alphaj = (Pj'*Qj)\(Pj'*Rj);

Xj1 = Xj + Pj*alphaj;
Rj1 = Rj - Qj*alphaj;
% Rj1 = B - A*Xj1;
Zj1 = M*Rj1;

betaj = -(Pj'*Qj)\(Qj'*Zj1);
% betaj = (Rj'*Zj)\(Rj1'*Zj1);

Pj1 = orth(Zj1 + Pj*betaj);

Xj = Xj1;
Rj = Rj1;
Zj = Zj1;
Pj = Pj1;