function [softNullSumRate,softNulluLRate,softNulldLRate] = cal_SoftNull_sumRate(GmaDFdd,HdBtFdd,HbrBtFdd,GmaUFdd,HbrUFdd,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Md,K,rho)
%   Author  : Wei Zhang (wzhang@fudan.edu.cn)
%   Date    : 02/02/23

alpha=1-rho;
Mr = size(HbrBtFdd,1);
%% DL
[~,S,V]=svd(HbrBtFdd);
sVec = diag(S);
[~,idx] = sort(sVec, 'ascend');
Psic = V(:,idx(1:Md));
Hsic = GmaDFdd*HdBtFdd*Psic;

tmp = (Hsic*Hsic')^(-1);
tmpDiag = diag(tmp);
tmpTrace = trace(tmp);
softNulldLRate = K*log2(dLTxPow/usrFdNoisePow+tmpTrace) - sum(log2(K*tmpDiag));
Rsd = diag((dLTxPow+tmpTrace*usrFdNoisePow)./(K*tmpDiag))-usrFdNoisePow;
%% UL
Pd = Hsic'/(Hsic*Hsic');
P = Psic*Pd;
Hu = HbrUFdd*GmaUFdd;
Hsi = HbrBtFdd*P;
Ryb = txUsrPow*(Hu*Hu')+(Hsi*Rsd*Hsi')+ bsFdNoisePow*eye(Mr);
RybDiag = diag(Ryb);
Rnq = rho*alpha*diag(RybDiag);
Qb = Rnq+alpha^2*bsFdNoisePow*eye(Mr);
softNulluLRate = real(log2(det(eye(K)+alpha^2*txUsrPow*(Hu'/Qb*Hu))));
softNullSumRate = softNulldLRate + softNulluLRate;

end