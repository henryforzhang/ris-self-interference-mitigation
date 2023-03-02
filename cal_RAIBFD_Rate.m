function [raibfdSumRate,uLRate,dLRate] = cal_RAIBFD_Rate(HrBtFd,HdRFd,HdBtFd,HbrBtFd,GmaDFd,HbrRFd,HrUFd,HbrUFd,GmaUFd,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Psic,D,K,rho)
%   Author  : Wei Zhang (wzhang@fudan.edu.cn)
%   Date    : 02/02/23

alpha = 1-rho;
Mr = size(HbrBtFd,1);
%% DL
Hd = GmaDFd*(HdRFd*D*HrBtFd+HdBtFd);
Hsic = Hd*Psic;
tmp = (Hsic*Hsic')^(-1);
tmpDiag = diag(tmp);
tmpTrace = trace(tmp);
dLRate = K*log2(dLTxPow/usrFdNoisePow+tmpTrace)-sum(log2(K*tmpDiag));
Rsd = diag((dLTxPow+tmpTrace*usrFdNoisePow)./(K*tmpDiag))-usrFdNoisePow; % DL power allocation

%% UL
Pd = Hsic'*tmp;
P = Psic*Pd; % DL precoding
Hu = (HbrRFd*D*HrUFd+HbrUFd)*GmaUFd;
Gsi = HbrRFd*D*HrBtFd+HbrBtFd;
Hsi = Gsi*P;
Ryb = txUsrPow*(Hu*Hu')+(Hsi*Rsd*Hsi')+ bsFdNoisePow*eye(Mr);
RybDiag = diag(Ryb);
Rnq = rho*alpha*diag(RybDiag);
Qb = Rnq+alpha^2*bsFdNoisePow*eye(Mr);
uLRate = real(log2(det(eye(K)+alpha^2*txUsrPow*(Hu'/Qb*Hu))));

raibfdSumRate = dLRate + uLRate;
end