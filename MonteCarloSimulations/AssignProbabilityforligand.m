function [ p ] = AssignProbabilityforligand(H2BPDC,formic,formate_formic,water,gama,headpenalty)
DMF = log(1-exp(formic)-exp(water)); %DMF concentration
N=1+0.479*exp(formate_formic); %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
K1=2.6;
K2=3.3;
alpha1=exp(formic+log(gama(2)))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
alpha2=K2*exp(2*(water+log(gama(1))))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
formicstar=alpha1.*(formic+log(gama(2)))+alpha2.*(2*(water+log(gama(1))))+(1-alpha1-alpha2).*((DMF+log(gama(3)))+(water+log(gama(1)))); % a parameter to describe all the terminating ligands on the SBUs
mixentropy=-alpha1.*log(alpha1)-alpha2.*log(alpha2)-(1-alpha1-alpha2).*log(1-alpha1-alpha2);
% mixentropy=0;
X=formicstar-(H2BPDC/2-log(N));% a parameter to account for all the concentration related terms
dG1=+1.0+alpha2*log(K2)+(1-alpha1-alpha2)*log(K1);%the Gibbs energy change of a single step of ligand replacement on the SBU,including the entropy of releasing a molecule. the dG(ligand replacement) = -RTlnK = -RTln(0.012)= +4.42 RT. 
% dG1=2.5;%the Gibbs energy change of a single step of ligand replacement on the SBU,including the entropy of releasing a molecule. the dG(ligand replacement) = -RTlnK = -RTln(0.012)= +4.42 RT. 

dST=12;%The entropy of releasing a molecule is dG(entropy)=-12RT (corresponding to 100 Jmol-1K-1).
dGn2 = 2*(dG1+mixentropy) -2*dST + 2*X + 2*headpenalty;
p = exp(-dGn2)./(1+exp(-dGn2)); 
end

