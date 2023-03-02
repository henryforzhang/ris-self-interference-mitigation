function AB = krb(A,B)
% KRB Khatri-Rao product

[I,F]=size(A);
[J,F1]=size(B);

if F~=F1
   szerror(' Error in krb.m - The matrices must have the same number of columns')
end

AB=zeros(I*J,F);
for f=1:F
   ab=B(:,f)*A(:,f).';
   AB(:,f)=ab(:);
end