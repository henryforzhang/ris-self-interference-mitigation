function [iraibfdSumRate,ulRate,dlRate] = cal_IRAIBFD_sumRate(GmaDFd,HrBtFd,HdRFd,HdBtFd,dLTxPow,usrFdNoisePow,GmaUFd,HbrRFd,HrUFd,HbrUFd,txUsrPow,bsFdNoisePow,phi)
%   Author  : Wei Zhang (wzhang@fudan.edu.cn)
%   Date    : 02/02/23

d = exp(1j*phi);
dLen = length(d);

Mr = size(HbrRFd,1);
Kd = size(HdBtFd,1);

alphaU = txUsrPow/bsFdNoisePow;
alphaD = dLTxPow/usrFdNoisePow;

%% Manifold Optimization
iterNum = 1e3;
negSumRateTmp = zeros(iterNum,1);
ulRateSet = zeros(iterNum,1);
dlRateSet = zeros(iterNum,1);
[negSumRateTmp(1),ulRateSet(1),dlRateSet(1)] = neg_sum_rate_func(d);
initT = 0;
for ii = 2 : iterNum
    if ii == 2
        oldGrad = neg_sum_rate_der(d);
        dir = -oldGrad; % update direction
    end
    df0 = real(oldGrad'*dir);
    [stepSize,newd,initT] = back_tracking_line_search(d,df0,dir,initT);
    [negSumRateTmp(ii),ulRateSet(ii),dlRateSet(ii)] = neg_sum_rate_func(newd);
    dirTransp = proj(newd,dir);
    %%
    newGrad = neg_sum_rate_der(newd);
    betaRcg = norm(newGrad)^2/norm(oldGrad)^2;
    dir = -newGrad+betaRcg*dirTransp;
    %%
    d = newd;
    oldGrad = newGrad;
    if ii > 5 && abs(negSumRateTmp(ii)-negSumRateTmp(ii-5))/abs(negSumRateTmp(ii)) < 1e-2
        break;
    end
end
negSumRateSet = negSumRateTmp(1:ii);
% figure; plot(negSumRateSet);
iraibfdSumRate = -negSumRateSet(ii);
ulRate = ulRateSet(ii);
dlRate = dlRateSet(ii);



    function [negSumRate,dlRate,ulRate] = neg_sum_rate_func(d)
        D = diag(d);
        Hu = (HbrRFd*D*HrUFd+HbrUFd)*GmaUFd;
        ulRate = real(log2(det(eye(Mr)+alphaU*(Hu*Hu'))));
        
        Hd = GmaDFd*HdRFd*D*HrBtFd+GmaDFd*HdBtFd;
        tmp = (Hd*Hd')^(-1);
        tmpDiag = diag(tmp);
        tmpTrace = trace(tmp);
        dlRate = (Kd*log2(alphaD+tmpTrace)-sum(log2(Kd*tmpDiag)));
        negSumRate = -(dlRate + ulRate);
    end
    function negSumRateDerd = neg_sum_rate_der(d)
        D = diag(d);
        Hu = HbrRFd*D*HrUFd*GmaUFd+HbrUFd*GmaUFd;
        Hd = GmaDFd*HdRFd*D*HrBtFd+GmaDFd*HdBtFd;
        tmp = (Hd*Hd')^(-1);
        tmpTrace = trace(tmp);
        tmpDiag = diag(tmp);
        
        Fd = sum(((HrBtFd*(Hd)'*tmp).*(tmp*GmaDFd*HdRFd).'),2);
        ulRateDer = 1/(log(2))*conj(alphaU*(diag(HrUFd*GmaUFd/(eye(Kd)+alphaU*(Hu'*Hu))*Hu'*HbrRFd)));
        dlRateDer = 1/(log(2))*conj(- Kd*Fd/(alphaD+tmpTrace)+...
            sum(((HrBtFd*(Hd)'*tmp).*(tmp*GmaDFd*HdRFd).')./(ones(dLen,1)*tmpDiag.'),2));
        enegSumRateDerd = -(ulRateDer+dlRateDer);
        negSumRateDerd = proj(d,enegSumRateDerd);
    end
    function y = retr(z, v, t)
        if nargin <= 2
            t = 1.0;
        end
        y = z+t*v;
        y = y ./ abs(y);
    end
    function y = proj(z,u)
        y = u - real( conj(u) .* z ) .* z;
    end
    function [stepSize,nextd,initT] = back_tracking_line_search(d,df0,direction,initT)
        alpha = 0.5;
        beta = 0.5;
        dirNorm = norm(direction);
        if initT == 0
            t = 1/dirNorm;
        else
            t = initT;
        end
        [tmp2,~,~] = neg_sum_rate_func(d);
        costEva = 1;
        while 1
            nextd = retr(d, direction, t);
            [tmp1,~,~] = neg_sum_rate_func(nextd);
            if tmp1 > tmp2 + alpha*t*df0
                if t < 1e-5
                    nextd = d;
                    break;
                end
                t = beta*t;
            else
                break;
            end
            costEva = costEva+1;
        end
        stepSize = t*dirNorm;
        switch costEva
            case 1
                % If things go well, push your luck.
                initT = 2*t;
            case 2
                % If things go smoothly, try to keep pace.
                initT = t;
            otherwise
                initT = 0;
        end
                   
    end
end