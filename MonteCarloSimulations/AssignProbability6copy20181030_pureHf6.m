function [ p ] = AssignProbability6copy20181030_pureHf6(H2BPDC,formic,formate_formic,water,n,n1,delta,x4,x5,addprobability,entropy_penalty)
DMF = log(1-exp(formic)-exp(water)); %DMF concentration
N=1+0.479*exp(formate_formic); %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
% dG10 = -2.15-log(0.479) ; %0.371 ligand and formic acid equilibrium as measured by experiment
% deltaS = 2.5+1.5+3.7+ delta;%excess amount of translational and rotational entropy gain for multiple connections
% dG1=dG10-deltaS/2;
% dG2=-7.7-deltaS;%7.7 is the experimental data from mass spectrum analysis with alpha1=0.5
% N=1+0.479*exp(formate_formic); %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
alpha1=exp(formic+formate_formic)./(exp(formic+formate_formic)+exp(-x4+2*water+formate_formic)+exp(-x5+2*DMF));
alpha2=exp(-x4+2*water+formate_formic)./(exp(formic+formate_formic)+exp(-x4+2*water+formate_formic)+exp(-x5+2*DMF));
formicstar=alpha1.*(formic)+alpha2.*(2*(water+2*log(1-exp(water)))-formate_formic)+(1-alpha1-alpha2).*(2*(DMF+2*log(1-exp(DMF)))-formate_formic);
X=formicstar-(H2BPDC-log(N));
% X=formicstar-1/2*H2BPDC;
dG1=-8.6015-delta;
dGn1 = (n-n1)*dG1+(n-n1)*X + entropy_penalty;
p = exp(dGn1)*addprobability; 
end

