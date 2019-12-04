% close all
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data')
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes')
load('comparethenucleationstep20190909_delta0.mat')
load('comparegrowthrates20190909_delta=0.mat')
expdata=xlsread('sem learning-6.0-reduce.xlsx','E2:J85');
expthickness=expdata(:,1)';
expthicknessstd=expdata(:,2)';
explateralsize=expdata(:,3)';
explateralsizestd=expdata(:,4)';
expratio=expdata(:,5)';
expratiostd=expdata(:,6)';

thickness1=sort(crystalgrowth,2);
stdthickness1=std(thickness1(:,2:end-1),1,2);
thickness1=mean(thickness1(:,2:end-1),2);
thickness2=nuclearate*40;
select1=1:length(expthickness);

figure
plot(thickness1,thickness2,'*')

scale = max(thickness(select1));
thickness=thickness/scale;

select2=1:78;
select2=intersect(select2,find(exp(formic)/0.37+exp(water)/0.56>=1.07));

expdata=xlsread('sem learning-6.0-reduce.xlsx','AJ2:AL85');
Hf12_Hf6=log(expdata(:,1))';
ligandmetalratio=expdata(:,2)';
formate_formic=log(expdata(:,3))';
expdata=xlsread('sem learning-6.0-reduce.xlsx','L2:N85');
Hf=log(expdata(:,1))';
water=log(expdata(:,2))';
formic=log(expdata(:,3))';
formate=formic+formate_formic;
expdata=xlsread('sem learning-6.0-reduce.xlsx','AN2:AN85');
H2BPDC=log(expdata(:,1))';


%% figure of simulated thickness vs experimental thickness
figure
% p=polyfit(thickness(select2)',expthickness(select2),1);
plot(thickness1(select2),expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
errorbarxy(thickness1(select2),expthickness(select2),stdthickness1(select2)/2,expthicknessstd(select2)/2);
% hold on;
% f=polyval(p,0:0.01:1.2);
% plot(0:0.01:1.2,f,'k','LineWidth',1.5);
% hold off
set(gca,'FontSize',14);
xlabel('Simulated Growth Rate (SBU/step)');
ylabel('Thickness (Exp) / nm');
title('Simulated Growth Rate and Experimental Thickness')
% xlim([0 1.1])
ylim([0 110])
set(gca,'XTick',0:0.25:1,'YTick',0:50:100,'FontSize',14);


figure
% p=polyfit(thickness(select2)',expthickness(select2),1);
plot(thickness2(select2),expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
% hold on;
% f=polyval(p,0:0.01:1.2);
% plot(0:0.01:1.2,f,'k','LineWidth',1.5);
% hold off
set(gca,'FontSize',14);
xlabel('Simulated Nucl. Rate (SBU/step)');
ylabel('Thickness (Exp) / nm');
title('Simulated Nucl. Rate and Experimental Thickness')
% xlim([0 1.1])
ylim([0 110])
set(gca,'XTick',0:0.1:1,'YTick',0:50:100,'FontSize',14);


figure
% p=polyfit(thickness(select2)',expthickness(select2),1);
plot((thickness1(select2).*thickness2(select2)),expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
% hold on;
% f=polyval(p,0:0.01:1.2);
% plot(0:0.01:1.2,f,'k','LineWidth',1.5);
% hold off
set(gca,'FontSize',14);
xlabel('Simulated Nucl.* Growth Rate (SBU/step)');
ylabel('Thickness (Exp) / nm');
title('Simulated Nucl.* Growth Rate and Experimental Thickness')
% xlim([0 1.1])
ylim([0 110])
set(gca,'XTick',0:0.1:1,'YTick',0:50:100,'FontSize',14);
% data2save=[thickness(select3)' thicknessstd(select3) expthickness(select3)' expthicknessstd(select3)'];
% save('simulated_and_experimental_thickness.txt','data2save','-ascii')

figure
% p=polyfit(thickness(select2)',expthickness(select2),1);
plot(thickness1(select2),thickness2(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
% hold on;
% f=polyval(p,0:0.01:1.2);
% plot(0:0.01:1.2,f,'k','LineWidth',1.5);
% hold off
set(gca,'FontSize',14);
xlabel('Simulated Growth Rate (SBU/step)');
ylabel('Simulated Nucl. Rate (SBU/step)');
title('Simulated Nucl. vs. Growth Rates')
% xlim([0 1.1])
% ylim([0 110])
set(gca,'XTick',0:0.25:1,'YTick',0:0.1:1,'FontSize',14);


%figure of experimental thickness vs experimental lateral size
figure
% p=polyfit(explateralsize(select2)',expthickness(select2)',1);
plot(explateralsize,expthickness,'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy(explateralsize,expthickness,explateralsizestd/2,expthicknessstd/2);
hold on;
% f=polyval(p,0:2200);
% plot(0:2200,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',14);
% title('Lateral size and Thickness')
% xlabel('Lateral Size (Exp) / nm','FontSize',14);
% ylabel('Thickness (Exp) / nm','FontSize',14);
xlim([0 2200])
ylim([0 120])
set(gca,'XTick',0:1000:2000,'YTick',0:50:100,'FontSize',16);

data2save=[explateralsize' explateralsizestd' expthickness' expthicknessstd'];
save('lateralsize_and_thickness.txt','data2save','-ascii')


%figure of synthetic condition vs expthickness
figure
% p=polyfit(H2BPDC(select2)',expthickness(select2)',1);
plot(H2BPDC(select2)/2.303,expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy(H2BPDC(select2)/2.303,expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-10.5:0.01:-9.4);
% plot(-10.5:0.01:-9.4,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('H_{2}BPDC solubility and Thickness')
xlabel('log([H_{2}BPDC])');
ylabel('Thickness (Exp) / nm');
xlim([-4.5 -3.9])
ylim([0 120])
set(gca,'XTick',-5.5:0.2:0,'YTick',0:50:100,'FontSize',16);

% data2save=[H2BPDC' expthickness' expthicknessstd'];
% save('H2BPDC_and_thickness.txt','data2save','-ascii')



figure
% p=polyfit(water(select2)',expthickness(select2)',1);
plot(water(select2)/2.303,expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy(water(select2)/2.303,expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-4.5:0.01:0);
% plot(-4.5:0.01:0,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('Water Concentration and Thickness')
xlabel('log([H_{2}O])');
ylabel('Thickness (Exp) / nm');
xlim([-1.8 -0.3])
ylim([0 120])
set(gca,'XTick',-2:0.5:0,'YTick',0:50:100,'FontSize',16);

figure
% p=polyfit(formic(select2)',expthickness(select2)',1);
plot(formic(select2)/2.303,expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy(formic(select2)/2.303,expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-4.5:0.01:0);
% plot(-4.5:0.01:0,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('Formic Acid Concentration and Thickness')
xlabel('log([HCO_{2}H])');
ylabel('Thickness (Exp) / nm');
xlim([-1 -0.25])
ylim([0 120])
set(gca,'XTick',-2:0.25:0,'YTick',0:50:100,'FontSize',16);


figure
% p=polyfit(formate_formic(select2)',expthickness(select2)',1);
plot(formate_formic(select2)/2.303,expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy(formate_formic(select2)/2.303,expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-4.5:0.01:0);
% plot(-4.5:0.01:0,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('[HCO_{2}^{-}]/[HCO_{2}H] and Thickness')
xlabel('log([HCO_{2}^{-}]/[HCO_{2}H])');
ylabel('Thickness (Exp) / nm');
% xlim([-2.5 0.2])
ylim([0 120])
% set(gca,'XTick',-2:1:0,'YTick',0:50:100,'FontSize',16);


figure
% p=polyfit(formate_formic(select2)',expthickness(select2)',1);
plot((Hf12(select2)+log(12))/2.303,expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy((Hf12(select2)+log(12))/2.303,expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-4.5:0.01:0);
% plot(-4.5:0.01:0,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('[Hf] and Thickness')
xlabel('log([Hf])');
ylabel('Thickness (Exp) / nm');
xlim([-3.3 -2])
ylim([0 120])
set(gca,'XTick',-4.5:0.5:0,'YTick',0:50:100,'FontSize',16);



figure
% p=polyfit(formate_formic(select2)',expthickness(select2)',1);
plot((ligandmetalratio(select2)),expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy((ligandmetalratio(select2)),expthickness(select2),[],expthicknessstd(select2)/2);
hold on;
% f=polyval(p,-4.5:0.01:0);
% plot(-4.5:0.01:0,f,'k','LineWidth',1.5);
hold off
set(gca,'FontSize',16);
title('H_{2}BPDC/Hf and Thickness')
xlabel('H_{2}BPDC/Hf');
ylabel('Thickness (Exp) / nm');
xlim([-0.3 4.5])
ylim([0 120])
set(gca,'XTick',-1:1:10,'YTick',0:50:100,'FontSize',16);
hold on
plot(0.75*ones(1,121),0:1:120,'r','linewidth',1.5)
plot(1*ones(1,121),0:1:120,'b','linewidth',1.5)
plot(0.5*ones(1,121),0:1:120,'k','linewidth',1.5)
% 

save('simulatethickness.mat','-v7.3')

%%

range=0:0.01:1;
[qx,qy]=meshgrid(range,range);

F1 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',expthickness(select2)','linear','none');
qz1=F1(qx,qy);
qz1(qx+qy>1)=0;
qz1(qz1<0)=0;

figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
title('Nanoplate Thickness')
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'y','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
hold off

cH2BPDC = (exp(H2BPDC(select2))- min(exp(H2BPDC(select2))))/(max(exp(H2BPDC(select2)))-min(exp(H2BPDC(select2))));
F2 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cH2BPDC','linear','none');
qz2=F2(qx,qy);
qz2(qx+qy>1)=0;
qz2(qz2<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz2,9,'w','linewidth',2,'showtext','on');
contour(range,range,qz2,9,'w','linewidth',2);
title('[H_{2}BPDC] and Thickness')
hold off

cformate_formic = (exp(formate_formic(select2))- min(exp(formate_formic(select2))))/(max(exp(formate_formic(select2)))-min(exp(formate_formic(select2))));
F3 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cformate_formic','linear','none');
qz3=F3(qx,qy);
qz3(qx+qy>1)=0;
qz3(qz3<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz3,10,'w','linewidth',2,'showtext','on');
contour(range,range,qz3,10,'w','linewidth',2);
title('[HCO_{2}^{-}]/[HCO_{2}H] and Thickness')
hold off

cHf12_Hf6 = (exp(Hf12_Hf6(select2))- min(exp(Hf12_Hf6(select2))))/(max(exp(Hf12_Hf6(select2)))-min(exp(Hf12_Hf6(select2))));
F4 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cHf12_Hf6','linear','none');
qz4=F4(qx,qy);
qz4(qx+qy>1)=0;
qz4(qz4<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz4,9,'w','linewidth',2,'showtext','on');
contour(range,range,qz4,9,'w','linewidth',2);
title('[Hf_{12}]/[Hf_{6}] and Thickness')
hold off


BPDC_formic=exp(H2BPDC)./exp(formic);
cBPDC_formic = (BPDC_formic(select2)- min(BPDC_formic(select2)))/(max(BPDC_formic(select2))-min(BPDC_formic(select2)));
F5 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cBPDC_formic','linear','none');
qz5=F5(qx,qy);
qz5(qx+qy>1)=0;
qz5(qz5<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz5,9,'w','linewidth',2,'showtext','on');
contour(range,range,qz5,9,'w','linewidth',2);
title('[H_{2}BPDC]/[HCO_{2}H] and Thickness')
hold off

cHf = (exp(Hf(select2))- min(exp(Hf(select2))))/(max(exp(Hf(select2)))-min(exp(Hf(select2))));
F6 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cHf','linear','none');
qz6=F6(qx,qy);
qz6(qx+qy>1)=0;
qz6(qz6<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz6,9,'w','linewidth',2,'showtext','on');
contour(range,range,qz6,9,'w','linewidth',2);
title('[Hf] and Thickness')
hold off


cL2M = (ligandmetalratio(select2)- min(ligandmetalratio(select2)))/(max(ligandmetalratio(select2))-min(ligandmetalratio(select2)));
F7 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cL2M','linear','none');
qz7=F7(qx,qy);
qz7(qx+qy>1)=0;
qz7(qz7<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz7,9,'w','linewidth',2,'showtext','on');
contour(range,range,qz7,9,'w','linewidth',2);
title('L/M and Thickness')
hold off

addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes\MonteCarlosimulation')
Calculatingactivitycoefficient;
gama=Gama;
DMF = log(1-exp(formic)-exp(water)); %DMF concentration
N=1+0.479*exp(formate_formic); %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
K1=2.6;
K2=3.3;
alpha1=exp(formic+log(gama(2)))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
alpha2=K2*exp(2*(water+log(gama(1))))./(K1*exp(water+log(gama(1))).*exp(DMF+log(gama(3)))+K2*exp(2*(water+log(gama(1))))+exp(formic+log(gama(2))));
formicstar=alpha1.*(formic+log(gama(2)))+alpha2.*(2*(water+log(gama(1))))+(1-alpha1-alpha2).*((DMF+log(gama(3)))+(water+log(gama(1)))); % a parameter to describe all the terminating ligands on the SBUs
X=formicstar-(H2BPDC/2-log(N));
BPDC_formic3 =X.*(1.5.*exp(Hf12_Hf6)+1)./(exp(Hf12_Hf6)+1);
cBPDC_formic3=((BPDC_formic3)-min((BPDC_formic3)))/(max((BPDC_formic3))-min((BPDC_formic3)));
figure
plot((cBPDC_formic3(select2)),expthickness(select2),'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
errorbarxy((cBPDC_formic3(select2)),expthickness(select2),[],expthicknessstd(select2)/2);
axis([0 1 0 100]);
ylabel('Crystal Thickness / nm','fontsize',18);
xlabel('Normalized \DeltaG_{a}','fontsize',18);
set(gca,'xtick',0:0.5:1,'ytick',0:50:100,'fontsize',18)


F8 = scatteredInterpolant(exp(water(select2))',exp(formic(select2))',cBPDC_formic3(select2)','linear','none');
qz8=F8(qx,qy);
qz8(qx+qy>1)=0;
qz8(qz8<0)=0;
figure
imagesc(range,range,qz1)
axis([0 0.5 0.1 0.5]);
colorbar;
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
% contour(range,range,qz8,20,'w','linewidth',2,'showtext','on');
contour(range,range,qz8,9,'w','linewidth',2);
title('\DeltaG_{a} and Thickness','fontsize',16)
hold off


