function [ p ] = AssignProbability6copy20181030_pureHf6_Llimited(Hf6left,n,n1,addprobability,entropy_penalty,entropy_penalty2)
if n>24
    n=24;
end
dGn1 = (1-n/n1).*(-Hf6left+entropy_penalty2-log(nchoosek(24,n1))) + entropy_penalty - (log(nchoosek(24,n))-log(nchoosek(24,n1)));
p = exp(dGn1)*addprobability; 
end

