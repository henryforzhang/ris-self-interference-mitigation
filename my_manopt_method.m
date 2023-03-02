function [D,cost] = my_manopt_method(bVec,C,D)
%% manifold optimization
d = diag(D);
iterNum = 1e3;
costFunc = zeros(iterNum,1);
[costFunc(1)] = cost_func(d);
initT = 0;
for ii = 2 : iterNum
    %%
    if ii == 2
        oldGrad = cost_func_der(d);
        dir = -oldGrad; % update direction
    end
    df0 = real(oldGrad'*dir);
    [stepSize,newf,initT] = back_tracking_line_search(d,df0,dir,initT);
    costFunc(ii) = cost_func(newf);
    dirTransp = proj(newf,dir);
    %%
    newGrad = cost_func_der(newf);
    betaRcg = norm(newGrad)^2/norm(oldGrad)^2;
    dir = -newGrad+betaRcg*dirTransp;
    %%
    d = newf;
    oldGrad = newGrad;
    if ii > 5 && abs(costFunc(ii)-costFunc(ii-5))/costFunc(ii-5) < 1e-2
        break;
    end
end
D = diag(d);
cost = costFunc(1:ii);
% figure; semilogy(cost)

    function y = cost_func(fbreve)
        yTmp = bVec+C*fbreve;
        y = yTmp'*yTmp;
    end
    function y = cost_func_der(fbreve)
        ey = C'*(bVec+C*fbreve);
        y = proj(fbreve,ey);
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
    function [stepSize,nextf,initT] = back_tracking_line_search(f,df0,direction,initT)
        alpha = 0.5;
        beta = 0.5;
        dirNorm = norm(direction);
        if initT == 0
            t = 1/dirNorm;
        else
            t = initT;
        end
        tmp2 = cost_func(f);
        costEva = 1;
        while 1
            nextf = retr(f, direction, t);
            tmp1 = cost_func(nextf);
            if tmp1 > tmp2 + alpha*t*df0
                if t < 1e-5
                    nextf = f;
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
