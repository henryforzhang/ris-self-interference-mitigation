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
% pathloss
f = 2.4e9;
c = 3e8;
lma = c/f;

%% channel
dRxTxIbfd = 6*lma/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRisIbfd = lma/2;
dRisShiftibfd = dRxTxIbfd/2;%lma/2*MrisSet(2)/2;
bPsSet = [1 2 3 4 5 6];
bPsSetLen = length(bPsSet);
%% Full Duplex 
Md = 8;
sicDbSet = zeros(MrisLen,mrMtLen,bPsSetLen);
nMonte = 20;
for kk = 1 : MrisLen
    MrisSize = MrisSet(kk,:);
    Mris = MrisSize(1)*MrisSize(2);
    for ii = 1 : mrMtLen
        MtSize = MtSet(ii,:);
        MrSize = MrSet(ii,:);
        Mr = MrSize(1)*MrSize(2);
        [HrBt,HbrR,HbrBt] = gen_nearFieldChan0516(MtSize,MrSize,MrisSize,dRxTxIbfd,dRxRisIbfd,dRisShiftibfd,f);
        for jj = 1 : bPsSetLen
            bPs = bPsSet(jj);
            sicDbTmp = 0;
            for nn = 1 : nMonte
                if nn == 1 || mod(nn, 5) == 0
                    fprintf(['\n',num2str(MrisSize(1)),'x',num2str(MrisSize(2)),' RIS, ', num2str(MrSize(1)),'x',num2str(MrSize(2)),...
                        ' Array, ', 'bPs = ',num2str(bPs),', nn = ',num2str(nn),', ',datestr(now,"HH:MM"), '\n']);
                end
                phi = rand(Mris,1)*2*pi;
                [D, P] = moBased_SIC_Algo(HbrR,HrBt,HbrBt,Md,phi,bPs);
                sic = norm((HbrR*D*HrBt+HbrBt)*P,'fro')^2/Mr;
                sicDbTmp = sicDbTmp-pow2db(sic); 
            end
            sicDbSet(kk,ii,jj) = sicDbTmp/nMonte;
        end
    end
end
%%
plt = {'k-o','k--o';'r-x','r--x';'b-*','b--*'};
figure;
for kk = 1 : MrisLen
    for ii = 1 : mrMtLen
        plot(bPsSet,squeeze(sicDbSet(kk,ii,:)),plt{kk,ii},'linewidth',2,'markersize',8)
        hold on
    end
end
grid on
xlabel('Phase Resolution (bit)')
ylabel('SIM amount: -10log_(10)(\kappa) (dB)');
legend('ULA, 4x4 RIS','URA, 4x4 RIS',...
    'ULA, 8x8 RIS','URA, 8x8 RIS',...
    'ULA, 16x16 RIS','URA, 16x16 RIS');
















%% full duplex communication system
%
%
clear;

% basic setting
%% RIS setting
MrisSet = [8,8];
Mris = MrisSet(1)*MrisSet(2);
%% Bs setting
MrSet = [8,1;4,2];
Mr = 8;
% MrSetMax = max(MrSet);
MtSet = [8,1;4,2];
Mt = 8;
% pathloss
f = 2.4e9;
c = 3e8;
lma = c/f;

%% channel
dRxTx = 6*lma/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRis = lma/2;
dRisShift = dRxTx/2;%lma/2*MrisSet(2)/2;


%% Full Duplex 
% load('./saveData0701/4x4RIS_0701_1306_22.mat');
% load('./saveData0701/8x8RIS_0701_1349_22.mat')
% load('./saveData0701/8x8RIS_0702_0110_22.mat')
% load('./saveData0701/16x16RIS_0701_1407_22.mat')
bSet = [2,3,4,5,6];
bSetLen = length(bSet);
Md = 8;
sicDbSet = zeros(2,bSetLen);
phi = rand(Mris,1)*2*pi;
for hh = 1 : 2
    [Hrat,Harr,Had,~] = gen_nearFieldChan0516(MtSet(hh,:),MrSet(hh,:),MrisSet,dRxTx,dRxRis,dRisShift,f);
    for mm = 1 : bSetLen
        b = bSet(mm);
        [D, P,test] = moBased_SIC_Algo(Harr,Hrat,Had,Md,phi,b);
        sic = norm((Harr*D*Hrat+Had)*P,'fro')^2/Mr;
        sicDbSet(hh,mm) = -pow2db(sic);
    end    
end
%%
plt1 = {'k-o','k--o'};
plt2 = {'r-x','r--x'};
plt3 = {'b-*','b--*'};
figure;
for ii = 1 : 2
    plot(bSet,sicDbSet(ii,:),plt2{ii},'linewidth',2,'markersize',8)
    hold on
end
grid on
xlabel('Phase Resolution (bit)');
ylabel('Self-Interference Cancellation (dB)');
title([num2str(MrisSet(1)),'x',num2str(MrisSet(2)),' RIS'])
% legend('ULA, 4x4 RIS','URA, 4x4 RIS');
legend('ULA, 8x8 RIS','URA, 8x8 RIS');
% legend('ULA, 16x16 RIS','URA, 16x16 RIS');
















