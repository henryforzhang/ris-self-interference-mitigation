function y = quantize_RIS_element(D,b)
%   Author  : Wei Zhang (wzhang@fudan.edu.cn)
%   Date    : 02/02/23

d = diag(D);
dAngle = angle(d);
dLen = length(d);
Nb = 2^b;
delta = 2*pi/Nb; % quantization interval
tmp = 2*pi-delta/2;
quantLevel = (0:Nb-1)/Nb*2*pi;

for dd = 1:dLen
    phase = dAngle(dd);
    if phase<0
        phase=phase+2*pi;
    end
    if phase> tmp
        phase=phase-2*pi;
    end
    [~,idx] = min(abs(phase-quantLevel));
    phaseQuant = quantLevel(idx);
    dAngle(dd) = phaseQuant;
end
y = diag(exp(1j*dAngle));
end