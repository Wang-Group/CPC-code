%NRTL model to calculate activity coefficient for water(1)-formic
%acid(2)-DMF(3)
%system
%data from Fluid Phase Equilibria 275(2009) 27-32
%data from Fluid Phase Equilibria 429(2016) 254-265
%data from Fluid Phase Equilibria 415(2016) 51-57
%NRTL parameters in J/mol
A(1,2)=983.6;
A(2,1)=2491.1;
A(1,3)=872.3;
A(3,1)=425.9;
A(2,3)=-298.6;
A(3,2)=-2252.5;
Tau=A/(8.314*(298+150));
G=exp(-Tau*0.3);
x(:,1)=exp(water);
x(:,2)=exp(formic);
x(:,3)=1-x(:,1)-x(:,2);
Gama=zeros(length(water),3);
for i=1:length(water)
Gama(i,:)=exp(x(i,:)*(Tau.*G)./(x(i,:)*G) +  x(i,:)*(G./(ones(3,1)*x(i,:)*G).*(Tau-ones(3,1)*((x(i,:)*(Tau.*G))./(x(i,:)*G))))');
end