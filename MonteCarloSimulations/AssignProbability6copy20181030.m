function [ p ] = AssignProbability6copy20181030(H2BPDC,formic,formate_formic,water,n_head,n_lateral,delta,addprobability,gama,n1)
DMF = log(1-exp(formic)-exp(water)); %DMF concentration
N=1+0.479*exp(formate_formic); %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
n = (n_head + n_lateral);%number of connections to other SBUs
K1=2.6;
K2=3.3;
alpha1=exp(formic+log(gama(2)))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
alpha2=K2*exp(2*(water+log(gama(1))))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
formicstar=alpha1.*(formic+log(gama(2)))+alpha2.*(2*(water+log(gama(1))))+(1-alpha1-alpha2).*((DMF+log(gama(3)))+(water+log(gama(1)))); % a parameter to describe all the terminating ligands on the SBUs
mixentropy=0;
X=formicstar-(H2BPDC/2-log(N));% a parameter to account for all the concentration related terms
dG1=+1.0+alpha2*log(K2)+(1-alpha1-alpha2)*log(K1);%the Gibbs energy change of a single step of ligand replacement on the SBU,including the entropy of releasing a molecule. the dG(ligand replacement) = -RTlnK = -RTln(0.012)= +4.42 RT. 
dST=12;%The entropy of releasing a molecule is dG(entropy)=-9.02RT (corresponding to 75 Jmol-1K-1). the overal dG1=-9.02+4.42=-4.6 RT
entropy_penalty=4.2-delta;%consider the loss of entropy of Hf6 with five connections C(12,5)=792. Then log(C(12,5))=6.67 to the maximum; but the SBU in the MOFs are not that ordered either, we also have an entropy depending on the connection and symmetry of the SBU. In general, if there are more than three connections for the SBU, the Hf12 SBU can have only ln(6*2)= 2.5, while the Hf6 SBU can have ln(6*2) = 2.5. As a result the entropy penalty is roughly 6.67-2.5= 4.2 R.
head_penalty=0.5;%from DFT calculation
dGn1 = (n-n1)*(dG1+mixentropy-dST/2+X)+(n_head-1/2*n1)*head_penalty + entropy_penalty; %head_penalty is the difference in the enthalpy change of the head sites and the lateral sites  
p = exp(dGn1)*addprobability; %addprobability is the probability of adding an SBU in the simulation.
end

