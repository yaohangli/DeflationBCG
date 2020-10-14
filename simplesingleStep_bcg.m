function [Pj,Zj,Rj] = simplesingleStep_bcg(A,Pj,Zj,Rj)

[n,r] = size(Pj);
disp(['The size of Pj: ','[',num2str(n),',',num2str(r),']']);
Qj = A*Pj;

alphaj = (Pj'*Qj)\(Pj'*Rj);

Rj1 = Rj - Qj*alphaj;
disp(Rj);
disp(Qj*alphaj);
disp(Rj-Qj*alphaj);
disp(Rj1'*Rj)
% Rj1 = B - A*Xj1;
Zj1 = Rj1;

betaj = -(Pj'*Qj)\(Qj'*Zj1);
% betaj = (Rj'*Zj)\(Rj1'*Zj1);
% disp(-Pj'*Qj);
% disp(Qj'*Zj1);
% disp(Zj1 + Pj*betaj);

Pj1 = orth(Zj1 + Pj*betaj);

Rj = Rj1;
Zj = Zj1;
Pj = Pj1;