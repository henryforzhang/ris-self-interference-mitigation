clear;

%% RIS setting
MrisSet = [16,16];
Mris = MrisSet(1)*MrisSet(2);
%% Bs setting
MrSet = [8,1];
Mr = MrSet(1)*MrSet(2);
% MrSetMax = max(MrSet);
MtSet = [8,1];
Mt = MtSet(1)*MtSet(2);
bAdc = 12;
rho = pi*sqrt(3)/2*2^(-2*bAdc);
bPsSet = [1 2 3 4 5 6 Inf];
bPsSetLen = length(bPsSet);
phsErrSet = [0,5,10,15,20,25,30];
phsErrSetLen = length(phsErrSet);
K = 3;
Ku = K;
Kd = K;
Md = 8;
% user position
lowBd = 100;
upBd = 140;
[duSet,ddSet] = gen_usr_posi(K,lowBd,upBd);
% duSet = 100*ones(K,1);
% ddSet = 500*ones(K,1);

%% tx power
dLTxPowdBm = 30; % dBm
dLTxPow = dBm2pow(dLTxPowdBm);
txUsrPowDbm = 10; %dBm
txUsrPow = dBm2pow(txUsrPowDbm);

%% FD noise power
bsFdNoisePowdBm = -96; % dBm
bsFdNoisePow = dBm2pow(bsFdNoisePowdBm);
usrFdNoisePowdBm = -96; % dBm
usrFdNoisePow = dBm2pow(usrFdNoisePowdBm);
%% ibfd RIS to RX and Tx antenna channel
fFd = 2.4e9;
c = 3e8;
lmaFd = c/fFd;
betaULDbFd = gen_pathloss(duSet,fFd);
betaDLDbFd = gen_pathloss(ddSet,fFd);
GmaUFd = gen_gma(betaULDbFd);
GmaDFd = gen_gma(betaDLDbFd);
dRxTxFd = 6*lmaFd/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRisFd = lmaFd/2;
dRisShiftFd = dRxTxFd/2;%lma/2*MrisSet(2)/2;
[HrBtFd,HbrRFd,HbrBtFd] = gen_nearFieldChan0516(MtSet,MrSet,MrisSet,dRxTxFd,dRxRisFd,dRisShiftFd,fFd);
Nc = 3;Nray = 5;
%% Algorithm Test
raibfdSumRate = zeros(bPsSetLen,phsErrSetLen);raibfduLRate = zeros(bPsSetLen,phsErrSetLen);raibfddLRate = zeros(bPsSetLen,phsErrSetLen);
outnMonte = 5e1;
innerLp = 2e1;
for nn = 1 : outnMonte
    [HrU,HdR,HdBt,HbrU] = gen_Saleh_Valenzuela_model(K,MrisSet,MtSet,Nc,Nray);
    tmp7 = zeros(bPsSetLen,phsErrSetLen);tmp8 = zeros(bPsSetLen,phsErrSetLen);tmp9 = zeros(bPsSetLen,phsErrSetLen);
    for ii = 1 : innerLp
        phi = rand(Mris,1)*2*pi;
        raibfdSumRateTmp = 0;raibfduLRateTmp = 0; raibfddLRateTmp = 0;
        for pp = 1 : bPsSetLen
            bPs = bPsSet(pp);
            %% RAIBFD
            [D, Psic] = moBased_SIC_Algo(HbrRFd,HrBtFd,HbrBtFd,Md,phi,bPs);
            for ee = 1 : phsErrSetLen
                if (nn == 1 || mod(nn, 5) == 0) &&  (ii == 1 || mod(ii, 5) == 0)
                    fprintf(['\n',num2str(MrisSet(1)),'x',num2str(MrisSet(2)),' RIS, ', num2str(MrSet(1)),'x',num2str(MrSet(2)),...
                        ' Array, ', 'OuterLoop: ',num2str(nn),', InnerLoop: ',num2str(ii),', ',datestr(now,"HH:MM"), '\n']);
                end
                sgmaP = phsErrSet(ee)/180*pi;
                DErr = gen_D_with_phaseErr(D,sgmaP);
                [raibfdSumRateTmp,raibfduLRateTmp,raibfddLRateTmp] = cal_RAIBFD_Rate(HrBtFd,HdR,HdBt,HbrBtFd,GmaDFd,HbrRFd,HrU,...
                    HbrU,GmaUFd,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Psic,DErr,K,rho);
                tmp7(pp,ee) = tmp7(pp,ee) + raibfdSumRateTmp; tmp8(pp,ee) = tmp8(pp,ee) + raibfduLRateTmp; tmp9(pp,ee) = tmp9(pp,ee) + raibfddLRateTmp;
            end
        end
    end
    raibfdSumRate = raibfdSumRate+tmp7/innerLp; raibfduLRate = raibfduLRate+tmp8/innerLp; raibfddLRate = raibfddLRate+tmp9/innerLp;
end
raibfdSumRateAve=raibfdSumRate/nn;raibfduLRateAve=raibfduLRate/nn; raibfddLRateAve = raibfddLRate/nn;
%% plot
plt1 = {'k-o','k-x','k-*','k-^','k-v','k->','k-<'};
plt2 = 'r-';
plt3 = 'b-';
plt4 = 'm-';
%% Sum-rate
figure;
for ii = 1:bPsSetLen
    plot(phsErrSet,raibfdSumRateAve(ii,:),plt1{ii},'lineWidth',2,'markersize',8);hold on
end
set(gca,'FontSize',12);
xlabel('Phase Error: \sigma_p (^\circ)','FontSize', 14);
ylabel('Sum-Rate (bps/Hz)','FontSize',14);
legend('b=1','b=2','b=3','b=4','b=5','b=6','b=Inf','FontSize',12);
grid on




