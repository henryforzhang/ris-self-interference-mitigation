function P = Fd_SIC_precoding_opt_Algo(Harr,Hrat,Had,Md,D)
G = Harr*D*Hrat + Had;
R = G'*G;
[U,S,~] = eig(R);
SDiag = diag(S);
[~,idx] = sort(real(SDiag),'ascend');
P = U(:,idx(1:Md));
end