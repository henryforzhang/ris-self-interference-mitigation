function [HrU,HdR,HdBt,HbrU] = gen_Saleh_Valenzuela_model(K,MrisSet,MtSet,Nc,Nray)
%% user to RIS
numRis = MrisSet(1)*MrisSet(2);
HrU = zeros(numRis,K);
for kk = 1 : K
    HrU(:,kk) = channel_generation(MrisSet,Nc,Nray);
end
%% RIS to user
HdR = zeros(K,numRis);
for kk = 1 : K
    tmp = channel_generation(MrisSet,Nc,Nray);
    HdR(kk,:) = tmp';
end
%% Tx to user
numAnte = MtSet(1)*MtSet(2);
HdBt = zeros(K,numAnte);
for kk = 1 : K
    tmp = channel_generation(MtSet,Nc,Nray);
    HdBt(kk,:) = tmp';
end
%% user to Rx
HbrU = zeros(numAnte,K);
for kk = 1 : K
    HbrU(:,kk) = channel_generation(MtSet,Nc,Nray);
end
end