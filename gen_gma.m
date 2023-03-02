function y = gen_gma(betaDb)
len = length(betaDb);
y = zeros(len);
for ii = 1:len
    y(ii,ii) = db2mag(-betaDb(ii));
end
end