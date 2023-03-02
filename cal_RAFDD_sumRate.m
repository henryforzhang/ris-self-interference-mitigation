function [rafddSumRate,rafddulRate,rafdddlRate] = cal_RAFDD_sumRate(GmaDFdd,HrBtFdd,HdRFdd,HdBtFdd,dLTxPow,usrFddNoisePow,GmaUFdd,HbrRFdd,HrUFdd,HbrUFdd,txUsrPow,bsFddNoisePow,phi)
%   Author  : Wei Zhang (wzhang@fudan.edu.cn)
%   Date    : 02/02/23

d = exp(1j*phi);
dLen = length(d);

Mr = size(HbrRFdd,1);
K = size(HdBtFdd,1);

alphaU = txUsrPow/bsFddNoisePow;
alphaD = dLTxPow/usrFddNoisePow;

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
    if ii > 2 && abs(negSumRateTmp(ii)-negSumRateTmp(ii-2))/abs(negSumRateTmp(ii)) < 1e-2
        break;
    end
end
negSumRateSet = negSumRateTmp(1:ii);
rafddSumRate = -1/2*negSumRateSet(ii);
rafddulRate = 1/2*ulRateSet(ii);
rafdddlRate = 1/2*dlRateSet(ii);



    function [negSumRate,dlRate,ulRate] = neg_sum_rate_func(d)
        D = diag(d);
        Hu = (HbrRFdd*D*HrUFdd+HbrUFdd)*GmaUFdd;
        ulRate = real(log2(det(eye(Mr)+alphaU*(Hu*Hu'))));
        
        Hd = GmaDFdd*HdRFdd*D*HrBtFdd+GmaDFdd*HdBtFdd;
        tmp = (Hd*Hd')^(-1);
        tmpDiag = diag(tmp);
        tmpTrace = trace(tmp);
        dlRate = (K*log2(alphaD+tmpTrace)-sum(log2(K*tmpDiag)));
        negSumRate = -(dlRate + ulRate);
    end
    function negSumRateDerd = neg_sum_rate_der(d)
        D = diag(d);
        Hu = HbrRFdd*D*HrUFdd*GmaUFdd+HbrUFdd*GmaUFdd;
        Hd = GmaDFdd*HdRFdd*D*HrBtFdd+GmaDFdd*HdBtFdd;
        tmp = (Hd*Hd')^(-1);
        tmpTrace = trace(tmp);
        tmpDiag = diag(tmp);
        
        Fd = sum(((HrBtFdd*(Hd)'*tmp).*(tmp*GmaDFdd*HdRFdd).'),2);
        ulRateDer = 1/(log(2))*conj(alphaU*(diag(HrUFdd*GmaUFdd/(eye(K)+alphaU*(Hu'*Hu))*Hu'*HbrRFdd)));
        dlRateDer = 1/(log(2))*conj(- K*Fd/(alphaD+tmpTrace)+...
            sum(((HrBtFdd*(Hd)'*tmp).*(tmp*GmaDFdd*HdRFdd).')./(ones(dLen,1)*tmpDiag.'),2));
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