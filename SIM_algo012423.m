function [D,P] = SIM_algo012423(MrisSet,MrSet,MtSet,f,bPs,dRxTx,dRxRis,dRisShift,Md)
[HrBt,HbrR,HbrBt] = gen_nearFieldChan0516(MtSet,MrSet,MrisSet,dRxTx,dRxRis,dRisShift,f);
[D, P, testTmp] = moBased_SIC_Algo(HbrR,HrBt,HbrBt,Md,phi,bPs);
end