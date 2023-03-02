function [duSet,ddSet] = gen_usr_posi(K,lowBd,upBd)
% lowBd = 100;
% upBd = 140;
y = lowBd + (upBd - lowBd)*rand(K,2);
duSet = y(:,1);
ddSet = y(:,2);
end