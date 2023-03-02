clear;
% basic setting
%% RIS setting
MrisSet = [16,16];
Mris = MrisSet(1)*MrisSet(2);
%% Bs setting
MrSet = [8,1];
Mr = MrSet(1)*MrSet(2);
% MrSetMax = max(MrSet);
MtSet = [8,1];
Mt = MtSet(1)*MtSet(2);
bAdcSet = [8,9,10,11,12];
bAdcSetLen = length(bAdcSet);
bPsSet = [1 2 3 4 5 6 Inf];
bPsSetLen = length(bPsSet);

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
%% HD noise power (3dB smaller than Fd)
bsHdNoisePowdBm = -99; % dBm
bsHdNoisePow = dBm2pow(bsHdNoisePowdBm);
usrHdNoisePowdBm = -99; % dBm
usrHdNoisePow = dBm2pow(usrHdNoisePowdBm);

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
%% fdd RIS to RX and Tx antenna channel
fuHd = 1760e6; % MHz
fdHd = 1855e6; % MHz
lmadHd = c/fdHd;
lmauHd = c/fuHd;
betaULDbHd = gen_pathloss(duSet,fuHd);
betaDLDbHd = gen_pathloss(ddSet,fdHd);
GmaUHd = gen_gma(betaULDbHd);
GmaDHd = gen_gma(betaDLDbHd);
dRxTxHd = 6*lmadHd/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRisHd = lmadHd/2;
dRisShiftHd = dRxTxHd/2;
[HrBtHd,HbrRHd,HbrBtHd] = gen_nearFieldChan0516(MtSet,MrSet,MrisSet,dRxTxHd,dRxRisHd,dRisShiftHd,fdHd);

Nc = 3;Nray = 5;
%% Algorithm Test
raibfdSumRate = zeros(bAdcSetLen,bPsSetLen);raibfduLRate = zeros(bAdcSetLen,bPsSetLen);raibfddLRate = zeros(bAdcSetLen,bPsSetLen);
iraibfdSumRate = 0;iraibfduLRate = 0; iraibfddLRate = 0;
rafddSumRate = 0;rafdduLRate = 0; rafdddLRate = 0;
softNullSumRate = zeros(bAdcSetLen,1); softNulluLRate = zeros(bAdcSetLen,1); softNulldLRate = zeros(bAdcSetLen,1);
outnMonte = 5e1;
innerLp = 2e1;
for nn = 1 : outnMonte
    if nn == 1 || mod(nn,5) == 0
        fprintf(['\nThe Monte Carol Run: %d, ',datestr(now,"HH:MM"),'\n'], nn); 
    end
    [HrU,HdR,HdBt,HbrU] = gen_Saleh_Valenzuela_model(K,MrisSet,MtSet,Nc,Nray);
    tmp1 = 0;tmp2 = 0;tmp3 = 0; tmp4 = 0;tmp5 = 0;tmp6 = 0;
    tmp7 = zeros(bAdcSetLen,bPsSetLen);tmp8 = zeros(bAdcSetLen,bPsSetLen);tmp9 = zeros(bAdcSetLen,bPsSetLen);
    tmp10= zeros(bAdcSetLen,1); tmp11= zeros(bAdcSetLen,1); tmp12= zeros(bAdcSetLen,1);
    for ii = 1 : innerLp
        phi = rand(Mris,1)*2*pi;
        raibfdSumRateTmp = 0;raibfduLRateTmp = 0; raibfddLRateTmp = 0;
        %% I-RAIBFD
        [iraibfdSumRateTmp,iraibfduLRateTmp,iraibfddLRateTmp] = cal_IRAIBFD_sumRate(GmaDFd,HrBtFd,HdR,HdBt,dLTxPow,usrFdNoisePow,GmaUFd,HbrRFd,HrU,HbrU,txUsrPow,bsFdNoisePow,phi);
        tmp1 = tmp1 + iraibfdSumRateTmp; tmp2 = tmp2 + iraibfduLRateTmp; tmp3 = tmp3+ iraibfddLRateTmp;
        %% RAFDD
        [rafddSumRateTmp,rafdduLRateTmp,rafdddLRateTmp] = cal_RAFDD_sumRate(GmaDHd,HrBtHd,HdR,HdBt,dLTxPow,usrHdNoisePow,GmaUHd,HbrRHd,HrU,HbrU,txUsrPow,bsHdNoisePow,phi);
        tmp4 = tmp4 + rafddSumRateTmp; tmp5 = tmp5 + rafdduLRateTmp; tmp6 = tmp6 + rafdddLRateTmp;
        for bb = 1 : bAdcSetLen
            bAdc = bAdcSet(bb);
            rho = pi*sqrt(3)/2*2^(-2*bAdc);
            for pp = 1 : bPsSetLen
                bPs = bPsSet(pp);
                %% RAIBFD
                [D, Psic] = moBased_SIC_Algo(HbrRFd,HrBtFd,HbrBtFd,Md,phi,bPs);
                [raibfdSumRateTmp,raibfduLRateTmp,raibfddLRateTmp] = cal_RAIBFD_Rate(HrBtFd,HdR,HdBt,HbrBtFd,GmaDFd,HbrRFd,HrU,...
                    HbrU,GmaUFd,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Psic,D,K,rho);
                tmp7(bb,pp) = tmp7(bb,pp) + raibfdSumRateTmp; tmp8(bb,pp) = tmp8(bb,pp) + raibfduLRateTmp; tmp9(bb,pp) = tmp9(bb,pp) + raibfddLRateTmp;
            end
            %% SoftNull
            [softNullSumRateTmp,softNulluLRateTmp,softNulldLRateTmp] = cal_SoftNull_sumRate(GmaDFd,HdBt,HbrBtFd,GmaUFd,HbrU,txUsrPow,bsFdNoisePow,dLTxPow,usrFdNoisePow,Md,K,rho);
            tmp10(bb) = tmp10(bb)+softNullSumRateTmp; tmp11(bb) = tmp11(bb)+softNulluLRateTmp; tmp12(bb) = tmp12(bb)+softNulldLRateTmp;
%             softNullSumRate(bb) = softNullSumRate(bb)+softNullSumRateTmp; softNulluLRate(bb) = softNulluLRate(bb)+softNulluLRateTmp; softNulldLRate(bb) = softNulldLRate(bb)+softNulldLRateTmp;
        end
    end
    raibfdSumRate = raibfdSumRate+tmp7/innerLp; raibfduLRate = raibfduLRate+tmp8/innerLp; raibfddLRate = raibfddLRate+tmp9/innerLp;
    iraibfdSumRate = iraibfdSumRate+tmp1/innerLp; iraibfduLRate = iraibfduLRate+tmp2/innerLp; iraibfddLRate = iraibfddLRate+tmp3/innerLp;
    rafddSumRate = rafddSumRate+tmp4/innerLp; rafdduLRate = rafdduLRate+tmp5/innerLp; rafdddLRate = rafdddLRate+tmp6/innerLp;
    softNullSumRate = softNullSumRate+tmp10/innerLp; softNulluLRate=softNulluLRate+tmp11/innerLp; softNulldLRate=softNulldLRate+tmp12/innerLp;
end
raibfdSumRateAve=raibfdSumRate/nn;raibfduLRateAve=raibfduLRate/nn; raibfddLRateAve = raibfddLRate/nn;
iraibfdSumRateAve=iraibfdSumRate/nn;iraibfduLRateAve=iraibfduLRate/nn; iraibfddLRateAve = iraibfddLRate/nn;
rafddSumRateAve=rafddSumRate/nn;rafdduLRateAve=rafdduLRate/nn; rafdddLRateAve = rafdddLRate/nn;
softNullSumRateAve=softNullSumRate/nn;softNulluLRateAve=softNulluLRate/nn; softNulldLRateAve = softNulldLRate/nn;
%% plot
plt1 = {'k-o','k-x','k-*','k-^','k-v','k->','k-<'};
plt2 = 'r-';
plt3 = 'b-';
plt4 = 'm-';
%% UL
figure;
for ii = 1:bPsSetLen
    plot(bAdcSet,raibfduLRateAve(:,ii),plt1{ii},'lineWidth',2,'markersize',8);hold on
end
plot(bAdcSet,iraibfduLRateAve*ones(1,bAdcSetLen),plt2,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,rafdduLRateAve*ones(1,bAdcSetLen),plt3,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,softNulluLRateAve,plt4,'lineWidth',2,'markersize',8);

set(gca,'FontSize',12);
xlabel('ENOB (bit)','FontSize', 15);
ylabel('UL Rate (bps/Hz)','FontSize',15);
legend('RAIBFD,b=1','RAIBFD,b=2','RAIBFD,b=3','RAIBFD,b=4','RAIBFD,b=5','RAIBFD,b=6','RAIBFD,b=Inf',...
    'I-RAIBFD','RAFDD','SoftNull','FontSize',15);
grid on
%% DL
figure;
for ii = 1:bPsSetLen
    plot(bAdcSet,raibfddLRateAve(:,ii),plt1{ii},'lineWidth',2,'markersize',8);hold on
end
plot(bAdcSet,iraibfddLRateAve*ones(1,bAdcSetLen),plt2,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,rafdddLRateAve*ones(1,bAdcSetLen),plt3,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,softNulldLRateAve,plt4,'lineWidth',2,'markersize',8);
set(gca,'FontSize',12);
xlabel('ENOB (bit)','FontSize', 15);
ylabel('DL Rate (bps/Hz)','FontSize',15);
legend('RAIBFD,b=1','RAIBFD,b=2','RAIBFD,b=3','RAIBFD,b=4','RAIBFD,b=5','RAIBFD,b=6','RAIBFD,b=Inf',...
    'I-RAIBFD','RAFDD','SoftNull','FontSize',15);
grid on
%% Sum-rate
figure;
for ii = 1:bPsSetLen
    plot(bAdcSet,raibfdSumRateAve(:,ii),plt1{ii},'lineWidth',2,'markersize',8);hold on
end
plot(bAdcSet,iraibfdSumRateAve*ones(1,bAdcSetLen),plt2,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,rafddSumRateAve*ones(1,bAdcSetLen),plt3,'lineWidth',2,'markersize',8);hold on
plot(bAdcSet,softNullSumRateAve,plt4,'lineWidth',2,'markersize',8);
set(gca,'FontSize',12);
xlabel('ENOB (bit)','FontSize', 15);
ylabel('Sum-Rate (bps/Hz)','FontSize',15);
legend('RAIBFD,b=1','RAIBFD,b=2','RAIBFD,b=3','RAIBFD,b=4','RAIBFD,b=5','RAIBFD,b=6','RAIBFD,b=Inf',...
    'I-RAIBFD','RAFDD','SoftNull','FontSize',15);
grid on



