function y = gen_pathloss(dSet,f)
c = 3e8;
lma = c/f;
y = -20*log10(lma./(4*pi*dSet));
end