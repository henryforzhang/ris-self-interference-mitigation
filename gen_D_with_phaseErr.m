function DErr = gen_D_with_phaseErr(D,sgmaP)
[Mris] = size(D,1);
phaseErr = (2*rand(1,Mris)-1)*sgmaP;
DErr = D.*diag(exp(1j*phaseErr));
end