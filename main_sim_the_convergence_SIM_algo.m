%% full duplex communication system
%
%
clear;
% basic setting
%% RIS setting
MrisSet = [8,8];
Mris = MrisSet(1)*MrisSet(2);
%% Bs setting
MrSet = [8,1];
Mr = MrSet(1)*MrSet(2);
% MrSetMax = max(MrSet);
MtSet = [8,1];
Mt = MtSet(1)*MtSet(2);
% pathloss
fIBFD = 2.4e9;
c = 3e8;
lmaIBFD = c/fIBFD;
bPs = inf;
%% channel
dRxTxIBFD = 6*lmaIBFD/2; % 接收天线中心与发射天线阵列中心相隔距离
dRxRisIBFD = lmaIBFD/2;
dRisShift = dRxTxIBFD/2;%lma/2*MrisSet(2)/2;
%% Full Duplex （self interference initialization）
phi = rand(Mris,1)*2*pi;
MdSet = [1,2,3,4,5,6,7,8];
MdSetLen = length(MdSet);
DSet = zeros(MdSetLen,Mris,2);
len = sum(MdSet);
PSet = zeros(Mt,len,2);
plt = {'ULA','URA'};
for hh = 1 : 2
    figure;
    [HrBt,HbrR,HbrBt] = gen_nearFieldChan0516(MtSet(hh,:),MrSet(hh,:),MrisSet,dRxTxIBFD,dRxRisIBFD,dRisShift,fIBFD);
%     iraibfdSumRate = cal_IRAIBFD_sumRate(GmaD,HrBt,HdR,HdBt,dLTxPow,usrFdNoisePow,GmaU,HbrR,HrU,HbrU,txUsrPow,bsFdNoisePow,phi);
    for mm = 1 : MdSetLen
        Md = MdSet(mm);
        [D, P, testTmp] = moBased_SIC_Algo(HbrR,HrBt,HbrBt,Md,phi,bPs);
        DSet(mm,:,hh) = diag(D);
        PSet(:,Md*(Md-1)/2+(1:Md),hh) = P;
        semilogy(testTmp,'lineWidth',2);
        hold on
    end    
    grid on
    set(gca,'FontSize',12);
    xlabel('Iteration Times','FontSize',15);
    ylabel('Cost Function','FontSize',15);
    xlim([1 1e4])
    legend('M_d = 1','M_d = 2','M_d = 3','M_d = 4','M_d = 5','M_d = 6','M_d = 7','M_d = 8','FontSize',15);
    title([plt{hh},', ',num2str(MrisSet(1)),'x',num2str(MrisSet(2)),' RIS'],'FontSize',15)
end
save(['.\saveData0821\',num2str(MrisSet(1)),'x',num2str(MrisSet(2)),'RIS_',num2str(bPs),'-bitRIS_'...
    datestr(now,"mmDD_HHMM_YY")],'PSet','DSet')













