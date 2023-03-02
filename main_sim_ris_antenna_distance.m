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
tmpSet = [1,2,3,4,5,6,7,8];
dRxRisIbfdSet = lma/2*tmpSet;
dRxRisIbfdSetLen = length(dRxRisIbfdSet);
dRisShiftibfd = dRxTxIbfd/2;%lma/2*MrisSet(2)/2;
bPs = Inf;
%% Full Duplex 
Md = 8;
sicDbSet = zeros(MrisLen,mrMtLen,dRxRisIbfdSetLen);
nMonte = 20;
for kk = 1 : MrisLen
    MrisSize = MrisSet(kk,:);
    Mris = MrisSize(1)*MrisSize(2);
    for ii = 1 : mrMtLen
        MtSize = MtSet(ii,:);
        MrSize = MrSet(ii,:);
        Mr = MrSize(1)*MrSize(2);
        for jj = 1 : dRxRisIbfdSetLen
            dRxRisIbfd = dRxRisIbfdSet(jj);
            [HrBt,HbrR,HbrBt] = gen_nearFieldChan0516(MtSize,MrSize,MrisSize,dRxTxIbfd,dRxRisIbfd,dRisShiftibfd,f);
            sicDbTmp = 0;
            for nn = 1 : nMonte
                if nn == 1 || mod(nn, 5) == 0
                    fprintf(['\n',num2str(MrisSize(1)),'x',num2str(MrisSize(2)),' RIS, ', num2str(MrSize(1)),'x',num2str(MrSize(2)),...
                        ' Array, ', 'tmpIdx = ',num2str(tmpSet(jj)),', nn = ',num2str(nn),', ',datestr(now,"HH:MM"), '\n']);
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
        plot(tmpSet,squeeze(sicDbSet(kk,ii,:)),plt{kk,ii},'linewidth',2,'markersize',8)
        hold on
    end
end
grid on
xlabel('Distance between atennna and RIS (\frac{\lambda}{2})')
ylabel('SIM amount: -10log_(10)(\kappa) (dB)');
legend('ULA, 4x4 RIS','URA, 4x4 RIS',...
    'ULA, 8x8 RIS','URA, 8x8 RIS',...
    'ULA, 16x16 RIS','URA, 16x16 RIS');















