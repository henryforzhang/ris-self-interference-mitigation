function [D, P] = moBased_SIC_Algo(HbrR,HrBt,HbrBt,Md,phi,b)
%% algorithm
D = diag(exp(1j*phi));
if b < inf
    D = quantize_RIS_element(D,b);
end
Mr = size(HbrR,1);
iterNum = 1e5;
test = zeros(iterNum,1);

for ii = 1 : iterNum
    %% precoding matrix
    P = Fd_SIC_precoding_opt_Algo(HbrR,HrBt,HbrBt,Md,D);
    testBf = norm((HbrR*D*HrBt+HbrBt)*P,'fro')^2/(Md*Mr);
    test(ii) = testBf;
    %% d (manifold)
    A = HrBt*P;
    B = HbrBt*P;
    C = krb(A.',HbrR);
    bVec = B(:);
    [DInf,cost] = my_manopt_method(bVec,C,D);
    if b < inf
        Dquant = quantize_RIS_element(DInf,b);
        testAf = norm((HbrR*Dquant*HrBt+HbrBt)*P,'fro')^2/(Md*Mr);
        if testAf <= testBf
            D = Dquant;
        end
    else
        D = DInf;
    end
    if (ii > 5 && abs(test(ii)-test(ii-5))/test(ii-5) < 1e-2) || abs(test(ii)) < 1e-15
        break;
    end
end
% figure; semilogy(test(1:ii))
end