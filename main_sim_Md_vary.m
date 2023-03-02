%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [4,4;8,8;16,16];
MrisLen = size(MrisSet,1);
%% Bs setting
MrSet = [8,1;4,2];
MtSet = [8,1;4,2];
mrMtLen = size(MtSet,1);
K = 3;
Ku = K;
Kd = K;
bAdc = 12;
rho = pi*sqrt(3)/2*2^(-2*bAdc);
%% ibfd RIS to RX and Tx antenna channel
fFd = 2.4e9;
c = 3e8;
lmaFd = c/fFd;
lowBd = 100;
upBd = 140;
[duSet,ddSet] = gen_usr_posi(K,lowBd,upBd);
betaULDbFd = gen_pathloss(duSet,fFd);
betaDLDbFd = gen_pathloss(ddSet,fFd);
GmaUFd = gen_gma(betaULDbFd);
GmaDFd = gen_gma(betaDLDbFd);
%% FD noise power
bsFdNoisePowdBm = -96; % dBm
bsFdNoisePow = dBm2pow(bsFdNoisePowdBm);
usrFdNoisePowdBm = -96; % dBm
usrFdNoisePow = dBm2pow(usrFdNoisePowdBm);
%% channel
dRxTxIbfd = 6*lmaFd/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRisIbfd = lmaFd/2;
dRisShiftibfd = dRxTxIbfd/2;%lma/2*MrisSet(2)/2;
bPs = Inf;
Nc = 3;Nray = 5;
%% tx power
dLTxPowdBm = 30; % dBm
dLTxPow = dBm2pow(dLTxPowdBm);
txUsrPowDbm = 10; %dBm
txUsrPow = dBm2pow(txUsrPowDbm);
%% Full Duplex
MdSet = [3,4,5,6,7,8];
MdSetLen = length(MdSet);
raibfdsumRateSet = zeros(MrisLen,mrMtLen,MdSetLen);
raibfduLRateSet = zeros(MrisLen,mrMtLen,MdSetLen);
raibfddLRateSet = zeros(MrisLen,mrMtLen,MdSetLen);
innernMonte = 2e1;
outernMonte = 5e1;

for kk = 1 : MrisLen
    MrisSize = MrisSet(kk,:);
    Mris = MrisSize(1)*MrisSize(2);
    for ii = 1 : mrMtLen
        MtSize = MtSet(ii,:);
        MrSize = MrSet(ii,:);
        Mr = MrSize(1)*MrSize(2);
        [HrBtFd,HbrRFd,HbrBtFd] = gen_nearFieldChan0516(MtSize,MrSize,MrisSize,dRxTxIbfd,dRxRisIbfd,dRisShiftibfd,fFd);
        for oo = 1 : outernMonte
            [HrU,HdR,HdBt,HbrU] = gen_Saleh_Valenzuela_model(K,MrisSize,MtSize,Nc,Nray);
            for jj = 1 : MdSetLen
                Md = MdSet(jj);
                raibfdsumRateTmp = 0; raibfduLRateTmp = 0; raibfddLRateTmp = 0;
                %                 softNullsumRateTmp = 0; softNulluLRateTmp = 0; softNulldLRateTmp = 0;
                for nn = 1 : innernMonte
                    if (nn == 1 || mod(nn, 5) == 0) && (oo == 1 || mod(oo,5) == 0)
                        fprintf(['\n',num2str(MrisSize(1)),'x',num2str(MrisSize(2)),' RIS, ', num2str(MrSize(1)),'x',num2str(MrSize(2)),...
                            ' Array, ', 'Md = ',num2str(Md),', nn = ',num2str(nn),', oo = ',num2str(oo), ', ',datestr(now,"HH:MM"), '\n']);
                    end
                    phi = rand(Mris,1)*2*pi;
                    [D, Psic] = moBased_SIC_Algo(HbrRFd,HrBtFd,HbrBtFd,Md,phi,bPs);
                    [raibfdsumRateTmp1,raibfduLRateTmp1,raibfddLRateTmp1] = cal_RAIBFD_Rate(HrBtFd,HdR,HdBt,HbrBtFd,GmaDFd,HbrRFd,HrU,...
                        HbrU,GmaUFd,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Psic,D,K,rho);
                    raibfdsumRateTmp = raibfdsumRateTmp + raibfdsumRateTmp1; raibfduLRateTmp = raibfduLRateTmp + raibfduLRateTmp1; raibfddLRateTmp = raibfddLRateTmp + raibfddLRateTmp1;
                    %                      [softNullSumRateTmp,softNulluLRateTmp,softNulldLRateTmp] = cal_SoftNull_sumRate(GmaDFd,HdBt,HbrBtFd,GmaUFd,HbrU,txUsrPow,...
                    %                          bsFdNoisePow,dLTxPow,usrFdNoisePow,Md,K,rho);
                    %                      softNullsumRateTmp = softNullsumRateTmp + softNullSumRateTmp; softNulluLRateTmp = softNulluLRateTmp + softNulluLRateTmp; softNulldLRateTmp = softNulldLRateTmp + softNulldLRateTmp;
                end
                raibfdsumRateSet(kk,ii,jj) = raibfdsumRateSet(kk,ii,jj)+raibfdsumRateTmp/innernMonte;
                raibfduLRateSet(kk,ii,jj) = raibfduLRateSet(kk,ii,jj)+raibfduLRateTmp/innernMonte;
                raibfddLRateSet(kk,ii,jj) = raibfddLRateSet(kk,ii,jj)+raibfddLRateTmp/innernMonte;
            end
        end
    end
end
raibfdsumRateSetAve = raibfdsumRateSet/oo;
raibfduLRateSetAve = raibfduLRateSet/oo;
raibfddLRateSetAve = raibfddLRateSet/oo;

%% UL
plt = {'k-o','k--o';'r-x','r--x';'b-*','b--*'};
figure;
for kk = 1 : MrisLen
    for ii = 1 : mrMtLen
        plot(MdSet,squeeze(raibfduLRateSetAve(kk,ii,:)),plt{kk,ii},'linewidth',2,'markersize',8)
        hold on
    end
end
grid on
xlabel('Number of DL effectve antenna')
ylabel('UL Rate (bps/Hz)');
legend('ULA, 4x4 RIS','URA, 4x4 RIS',...
    'ULA, 8x8 RIS','URA, 8x8 RIS',...
    'ULA, 16x16 RIS','URA, 16x16 RIS');
%% DL
plt = {'k-o','k--o';'r-x','r--x';'b-*','b--*'};
figure;
for kk = 1 : MrisLen
    for ii = 1 : mrMtLen
        plot(MdSet,squeeze(raibfddLRateSetAve(kk,ii,:)),plt{kk,ii},'linewidth',2,'markersize',8)
        hold on
    end
end
grid on
xlabel('Number of DL effectve antenna')
ylabel('DL Rate (bps/Hz)');
legend('ULA, 4x4 RIS','URA, 4x4 RIS',...
    'ULA, 8x8 RIS','URA, 8x8 RIS',...
    'ULA, 16x16 RIS','URA, 16x16 RIS');
%% Sumrate
plt = {'k-o','k--o';'r-x','r--x';'b-*','b--*'};
figure;
for kk = 1 : MrisLen
    for ii = 1 : mrMtLen
        plot(MdSet,squeeze(raibfdsumRateSetAve(kk,ii,:)),plt{kk,ii},'linewidth',2,'markersize',8)
        hold on
    end
end
grid on
xlabel('Number of DL effectve antenna')
ylabel('Sum Rate (bps/Hz)');
legend('ULA, 4x4 RIS','URA, 4x4 RIS',...
    'ULA, 8x8 RIS','URA, 8x8 RIS',...
    'ULA, 16x16 RIS','URA, 16x16 RIS');















