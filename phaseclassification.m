%% load experimental data
clear
% close all
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data\xrd-UiO-67');
load('pxrdanalysis.mat');

%load the synthetic conditions
[expdata txt]=xlsread('UiO-67-BPDC-20190304-150-120-100- 20190308.xlsx','D2:I286');

H2BPDCadd=expdata(:,1);
Hf=expdata(:,2);
water=expdata(:,3);
formic=expdata(:,5);
T=expdata(:,4);
%change the order of synthetic conditions to match the pxrd pattern
H2BPDCadd=H2BPDCadd(samplenum);
Hf=Hf(samplenum);
water=water(samplenum);
formic=formic(samplenum);
T=T(samplenum);
txt=txt(samplenum);

Xvariables=[ water formic Hf T H2BPDCadd ];


%% classification
%1 - UiO; 2 - hcp; 3 - hxl; 4 - Hf-formate; 5 - H2BPDC; 6 - broadUiO

barlevel=[0.1 0.05 0.05 0.05 0.05];
for i=1:5
    phase(:,i) = coeff(:,i) >= barlevel(i);
end

phase([6 10:12 14:15 17 101:102 109 127 128 221 231 240 248 254],2)=0;
phase([35 51 56 136 191],2)=1;
phase([12 14:15 136 147 195 203 215 228 236 244 255],3)=0;
broadUiO=[9 11 13 16 138 177 212 214 217 220 224 229 231 233 237 239 241 263];
phase(broadUiO,3)=0;
phase([6 108 121 183:186 188 237],3)=1;
phase([236 240 263],1)=0;
phase([10 11 60 61 89 102 133 149:150 152 231 240],1)=1;
phase(broadUiO,1)=1;
phase([12 63 90 209 240],4)=0;

addpoint = [0.33 0.11 mean(Xvariables(:,3)) 150 mean(Xvariables(:,5));
    0.27 0.12 mean(Xvariables(:,3)) 150 mean(Xvariables(:,5));
    0.38 0.13 mean(Xvariables(:,3)) 150 mean(Xvariables(:,5)); 
    0.37 0.08 mean(Xvariables(:,3)) 150 mean(Xvariables(:,5));
    0.45 0.03 mean(Xvariables(:,3)) 150 mean(Xvariables(:,5));];
addphase = [0 0 0 0 0;
            1 0 0 0 0;
            0 0 0 0 0;
            1 0 0 0 0;
            1 0 0 0 0;];
Xvariables = cat(1,Xvariables,addpoint);
phase = cat(1,phase,addphase);
phase(sum(phase,2)==0,5)=1;


%% Random Forest
close all
% define the phase to study
studyphase = 5;
% phasedata = sum(phaseprevious(:,studyphase),2)>=1;
phasedata= sum(phase(:,studyphase),2)>=1;
% phasedata(errorassign)= phaseprevious(errorassign,studyphase);

%random forest to view the importance of all the variables
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes\SMOTE')
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes')
trainsize=1.8;
trainindex=randi(size(Xvariables,1),[round(trainsize*size(Xvariables,1)) 1]);
trainindex=unique(trainindex);
testindex=setdiff((1:size(Xvariables,1))',trainindex);

[XXvariables Xphasedata]=SMOTE(Xvariables(trainindex,:),phasedata(trainindex,:));

%create a random forest
bagMdl=TreeBagger(1000,XXvariables(trainindex,:),Xphasedata(trainindex,:),'MaxNumSplits',4,'MinLeafSize',5,...
'PredictorNames',{'Water','Formic acid','Halfnium','Temperature','Ligand'},'Method','classification',...
'OOBPrediction','Off','OOBPredictorImportance','on','InBagFraction',1);
view(bagMdl.Trees{1},'mode','graph')
figure;
oobErrorBaggedEnsemble = oobError(bagMdl);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';

%calculate the importance parameters
imp = bagMdl.OOBPermutedPredictorDeltaError;
imp(imp<0)=0;

figure(1);
bar(imp);
title('Importance of Variables');
ylabel('Predictor importance estimates');
xlabel('Predictors');
h = gca;
h.XTickLabel = bagMdl.PredictorNames;
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';
set(h,'fontsize',14)

%the confusion matrix
figure(2)
Ypredict=predict(bagMdl,Xvariables(trainindex,:));
Y1predict=zeros(length(Ypredict),1);
for i=1:length(Ypredict)
    Y1predict(i)=str2num(Ypredict{i});
end
Y1predict=logical(Y1predict);
C = confusionmat(Y1predict,logical(phasedata(trainindex,:)));
cm = confusionchart(gather(Y1predict),gather(logical(phasedata(trainindex,:))));
title('Confusion Matrix for the Training Set');
set(gca,'fontsize',14)
xlabel('True Class');
ylabel('Predicted Class');

%the confusion matrix
figure(3)
Ypredict=predict(bagMdl,Xvariables(testindex,:));
Y1predict=zeros(length(Ypredict),1);
for i=1:length(Ypredict)
    Y1predict(i)=str2num(Ypredict{i});
end
Y1predict=logical(Y1predict);
C = confusionmat(Y1predict,logical(phasedata(testindex,:)));
cm = confusionchart(gather(Y1predict),gather(logical(phasedata(testindex,:))));
title('Confusion Matrix for the Test Set');
set(gca,'fontsize',14)
xlabel('True Class');
ylabel('Predicted Class');


%% see the distribution
figure
F1 = scatteredInterpolant(Xvariables(:,1),Xvariables(:,2),double(real(phase(:,1))),'linear','none');
F2 = scatteredInterpolant(Xvariables(:,1),Xvariables(:,2),double(real(phase(:,2))),'linear','none');
F3 = scatteredInterpolant(Xvariables(:,1),Xvariables(:,2),double(real(phase(:,3))),'linear','none');
F4 = scatteredInterpolant(Xvariables(:,1),Xvariables(:,2),double(real(phase(:,4))),'linear','none');
F5 = scatteredInterpolant(Xvariables(:,1),Xvariables(:,2),double(real(phase(:,5))),'linear','none');

range=0:0.01:1;
[qx,qy]=meshgrid(range,range);
qz1=F1(qx,qy);
qz1(qx+qy>1)=0;
qz1(qz1<0)=0;
qz2=F2(qx,qy);
qz2(qx+qy>1)=0;
qz2(qz2<0)=0;
qz3=F3(qx,qy);
qz3(qx+qy>1)=0;
qz3(qz3<0)=0;
qz4=F4(qx,qy);
qz4(qx+qy>1)=0;
qz4(qz4<0)=0;
qz5=F5(qx,qy);
qz5(qx+qy>1)=0;
qz5(qz5<0)=0;
clearvars qz
qz(1,:,:)=qz1;
qz(2,:,:)=qz2;
qz(3,:,:)=qz3;
qz(4,:,:)=qz4;
qz(5,:,:)=qz5;
figure
imagesc(range,range,qz2)
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
plot(range,range(end:-1:1),'k','linewidth',1);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
hold off

%%
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes\Multivariate ADTree\adtree package\adtree\ADTree_AdaBoost')
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes\Multivariate ADTree\adtree package\adtree\ADTree_AdaBoost\ADTree_Model')

Numtrees=1000;
clearvars phasepredict

%amplify the number of data accoding to the density of the data
N1=hist(Xvariables(:,1))/size(phase,1);
N2=hist(Xvariables(:,1))/size(phase,2);
n1=floor(Xvariables(:,1)/0.1)+1;
n2=floor(Xvariables(:,2)/0.1)+1;
datadensity=N1(n1).*N2(n2);
maxdensity=max(datadensity);
amplifynum=floor(sqrt(maxdensity./datadensity)/2);
maxamplify=max(amplifynum);
phaseuse=phase;
Xuse=Xvariables;
while maxamplify>1
    phaseuse=cat(1,phaseuse,phase(amplifynum>1,:));
    Xuse=cat(1,Xuse,Xvariables(amplifynum>1,:));
    amplifynum=amplifynum-1;
    maxamplify=max(amplifynum);
end

for j=1:5
    phasedata=double((phaseuse(:,j)>=1));
    data150=find(Xuse(:,4)>0);

    if sum(phasedata(data150))/length(data150)<=0.3
        [simXvariables simphasedata]=SMOTE(Xuse(data150,[1 2]),phasedata(data150));
    elseif sum(phasedata(data150))/length(data150)>=0.7
        [simXvariables simphasedata]=SMOTE(Xuse(data150,[1 2]),1-phasedata(data150));
        simphasedata=1-simphasedata;
    else
        simXvariables=Xuse(data150,[1 2]);
        simphasedata=phasedata(data150,:);
    end
[Y Z]=meshgrid(range,range);
Waterscan=reshape(Y,[],1);
Formicscan=reshape(Z,[],1);
checkpoints=[Waterscan Formicscan];
clearvars answers

for i=1:Numtrees
dataselect=randi(length(simphasedata),[length(simphasedata) 1]);
Xvariablesuse=simXvariables(dataselect,:);
phasedatause=double(simphasedata(dataselect,:));
phasedatause(phasedatause==0)=-1;
mode = 5;
GorL = 2;
Tt = 5;
stop = -2;
[ADTreeModel] = ADTree(Xvariablesuse,phasedatause,Tt,mode,GorL,stop);
[class,score,finalresult(1)] = ADTree_Model_Evaluation(ADTreeModel,checkpoints,ones(size(checkpoints,1),1));
answers(:,i)=score;
i
end
answers=mean(answers,2);
minans=min(answers);
maxans=max(answers);
answers=(answers-minans)/(maxans-minans);
answers(Waterscan + Formicscan>1)=0;
phasepredict(j,:,:)=reshape(answers,size(Y));
end

range=0:0.01:1;
drawlevel(1,:)=[0.9 0.9];
drawlevel(2,:)=[0.9 0.9];
drawlevel(3,:)=[0.9 0.9];
drawlevel(4,:)=[0.9 0.9];
drawlevel(5,:)=[0.9 0.9];
for j=1:5
figure(j+5)
imagesc(range, range,squeeze(phasepredict(j,:,:)))
hold on
contour(range, range,squeeze(phasepredict(j,:,:)),drawlevel(j,:),'r','linewidth',1.5)
goodpoint=find(phase(:,j)==1);
badpoint=find(phase(:,j)==0);
plot(Xvariables(goodpoint,1),Xvariables(goodpoint,2),'ro')
plot(Xvariables(badpoint,1),Xvariables(badpoint,2),'go')
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
plot(range,range(end:-1:1),'k','linewidth',1);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
colorbar;
hold off
end


%% draw the phase diagram
drawlevel(1,:)=[0.25 0.25];
drawlevel(2,:)=[0.9 0.9];
drawlevel(3,:)=[0.9 0.9];
drawlevel(4,:)=[0.5 0.5];
drawlevel(5,:)=[0.25 0.25];

figure(11)
contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'r');
C1=contourc(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:));

hold on
contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'b');
C2=contourc(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:));

contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'k');
C3=contourc(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:));

contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'g');
C4=contourc(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:));

contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'p');
C5=contourc(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:));

hold off

figure(12)
transparency=0.9;
axis([0 1 0 1]); 
hold on;

middlepoint=find(C5(2,:)>1);
pe=fill([C5(1,2:middlepoint(2)-1) C5(1,middlepoint(2)+1:end)],[C5(2,2:middlepoint(2)-1) C5(2,middlepoint(2)+1:end)],'r','linewidth',1.5);
alpha(transparency)

pd=fill(C4(1,2:end),C4(2,2:end),'b','linewidth',1.5);
alpha(transparency)

middlepoint=find(C1(2,:)>1);
pa=fill([C1(1,2:end) 0], [C1(2,2:end) 0],'c','linewidth',1.5);
alpha(transparency)

middlepoint=find(C2(2,:)>1);
pb=fill([C2(1,2:end) 0],[C2(2,2:end) 0],'y','linewidth',1.5);
alpha(transparency)
middlepoint=find(C3(2,:)>1);
pc=fill([C3(1,2:end) 0],[C3(2,2:end) 0],'g','linewidth',1.5);
alpha(transparency)

xlabel('H_2O (mol%)');
ylabel('HCO_2H (mol%)');

legend([pa,pb,pc,pd,pe],'UiO-67','hcp-UiO-67','hxl-UiO-67','Hf-Formate','H_{2}BPDC');
box on
legend('boxoff')
set(gca,'fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);

hold off

%% plot the phase diagram and synthetic parameters
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data')

[expdata txt]=xlsread('UiO-67-BPDC-20190304-150-120-100- 20190308.xlsx','X2:AE286');


Hf12_Hf6=expdata(:,3);
Hf12_Hf6=Hf12_Hf6(samplenum);
cHf12_Hf6=(log10(Hf12_Hf6)-min(log10(Hf12_Hf6)))/(max(log10(Hf12_Hf6))-min(log10(Hf12_Hf6)));

H2BPDC=expdata(:,2);
H2BPDC=H2BPDC(samplenum);
cH2BPDC=(log10(H2BPDC)-min(log10(H2BPDC)))/(max(log10(H2BPDC))-min(log10(H2BPDC)));


formate_formic=expdata(:,end-1)./(expdata(:,end)+0.01);
formate_formic=formate_formic(samplenum);
logformate_formic=log10(formate_formic);
logformate_formic(isnan(logformate_formic))=0;
logformate_formic(isinf(logformate_formic))=0;
cformate_formic=(logformate_formic-min(logformate_formic))/(max(logformate_formic)-min(logformate_formic));

BPDC_formic=H2BPDC./(formic+0.001);
cBPDC_formic=(log10(BPDC_formic)-min(log10(BPDC_formic)))/(max(log10(BPDC_formic))-min(log10(BPDC_formic)));

BPDC_formic2=H2BPDC./(formic+0.001).^2;
cBPDC_formic2=(log10(BPDC_formic2)-min(log10(BPDC_formic2)))/(max(log10(BPDC_formic2))-min(log10(BPDC_formic2)));


F1 = scatteredInterpolant(water,formic,cHf12_Hf6,'linear','none');
range=0:0.01:1;
[qx,qy]=meshgrid(range,range);
qz1=F1(qx,qy);
qz1(qx+qy>1)=0;
qz1(qz1<0)=0;

figure
imagesc(range,range,qz1)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 1])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
contour(range,range,qz1,20,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('[Hf_{12}]/[Hf_{6}]','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('[Hf_{12}]/[Hf_{6}] and Phases','fontsize',16)
hold off
%[co1,co2,co3,co4,co5],


F2 = scatteredInterpolant(water,formic,cH2BPDC,'linear','none');
qz2=F2(qx,qy);
qz2(qx+qy>1)=0;
qz2(qz2<0)=0;

figure
imagesc(range,range,qz2)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 1])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
co=contour(range,range,qz2,20,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('[H_{2}BPDC]','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('[H_{2}BPDC] and Phases','fontsize',16)
hold off


F3 = scatteredInterpolant(water,formic,cformate_formic,'linear','none');
qz3=F3(qx,qy);
qz3(qx+qy>1)=0;
qz3(qz3<0)=0;

figure
imagesc(range,range,qz3)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 1])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
co=contour(range,range,qz3,5,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('[HCO_{2}^{-}]/[HCO_{2}H]','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('[HCO_{2}^{-}]/[HCO_{2}H] and Phases','fontsize',16)
hold off



F4 = scatteredInterpolant(water,formic,cBPDC_formic,'linear','none');
qz4=F4(qx,qy);
qz4(qx+qy>1)=0;
qz4(qz4<0)=0;

figure
imagesc(range,range,qz4)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 0.8])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
co=contour(range,range,qz4,20,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('Phases','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('[H_{2}BPDC]/[HCO_{2}H] distribution','fontsize',16)
hold off

F5 = scatteredInterpolant(water,formic,cBPDC_formic2,'linear','none');
qz5=F5(qx,qy);
qz5(qx+qy>1)=0;
qz5(qz5<0)=0;

figure
imagesc(range,range,qz5)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 0.8])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
co=contour(range,range,qz5,20,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('Phases','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('[H_{2}BPDC]/[HCO_{2}H]^{2} distribution','fontsize',16)
hold off



addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\codes\MonteCarlosimulation')
Calculatingactivitycoefficient;
gama=Gama;
DMF = (1-formic-water); %DMF concentration
N=1+0.479*formate_formic; %pKa of benzoic acid at 150 C is 2.37 and the pKa of formic acid at 150 C is 2.05, 10^(2.05-2.37)=0.479
K1=2.6;
K2=3.3;
alpha1=((formic+0.001)*gama(2))./(K1*((water+0.001)*gama(1)).*(DMF*gama(3))+K2*(water+0.001).^2*gama(1)+(formic+0.001)*gama(2));
alpha2=K2*(water+0.001).^2*gama(1)./(K1*((water+0.001)*gama(1)).*(DMF*gama(3))+K2*(water+0.001).^2*gama(1)+(formic+0.001)*gama(2));
formicstar=((formic+0.001)*gama(2)).^(alpha1).*(gama(1)^2*(water+0.001).^2).^(alpha2).*((water+0.001)*gama(1)).*(DMF*gama(3)).^(1-alpha1-alpha2); % a parameter to describe all the terminating ligands on the SBUs
X=log(formicstar)-(log(H2BPDC)/2-log(N));
BPDC_formic3 =-log((H2BPDC./N)./(formicstar+0.001).^2).*((0.75*Hf12_Hf6+1)./(Hf12_Hf6+1))+0.5*0.5*(Hf12_Hf6./(Hf12_Hf6+1));
cBPDC_formic3=(BPDC_formic3-min((BPDC_formic3)))/(max((BPDC_formic3))-min((BPDC_formic3)));

F6 = scatteredInterpolant(water,formic,cBPDC_formic3,'linear','none');
qz6=F6(qx,qy);
qz6(qx+qy>1)=0;
qz6(qz6<0)=0;

figure
imagesc(range,range,qz6)
axis([-0.01 1 -0.01 1]);
colorbar;
caxis([0 1])
set(gca,'ydir','normal','fontsize',18,'xtick',0:0.25:1,'ytick',0:0.25:1);
hold on
co=contour(range,range,qz6,20,'w','linewidth',0.5);
xlabel('water (fraction)','fontsize',18);
ylabel('formic acid (fraction)','fontsize',18);
xlim([0 1])
ylim([0 1])
co5=contour(range,range,squeeze(phasepredict(5,:,:)),drawlevel(5,:),'r','linewidth',2);
co4=contour(range,range,squeeze(phasepredict(4,:,:)),drawlevel(4,:),'b','linewidth',2);
co1=contour(range,range,squeeze(phasepredict(1,:,:)),drawlevel(1,:),'c','linewidth',2);
co2=contour(range,range,squeeze(phasepredict(2,:,:)),drawlevel(2,:),'k','linewidth',2);
co3=contour(range,range,squeeze(phasepredict(3,:,:)),drawlevel(3,:),'g','linewidth',2);
legend('Normalized \DeltaG','H_{2}BPDC','Hf-Formate','UiO-67','hcp-UiO-67','hxl-UiO-67','fontsize',14,'location','northeast');
title('Normalized \DeltaG and Phases','fontsize',16)
hold off


%calculate phase percentage
coeff=zeros(size(pxrd,1),size(standard,1));
coeffweights=[1 1 1 1 1];
coeff(:,5)= pxrdint(:,15)/stanint(5,15).*phase(1:end-5,5)*coeffweights(5);
coeff(:,4)= pxrdint(:,6)/stanint(4,6).*phase(1:end-5,4)*coeffweights(4);
coeff(:,2)= pxrdint(:,1)/stanint(2,1).*phase(1:end-5,2)*coeffweights(2);
coeff(:,1)= (pxrdint(:,4)-coeff(:,2)*stanint(2,4))/stanint(1,4).*phase(1:end-5,1)*coeffweights(1);
coeff(:,3)= (pxrdint(:,2)-coeff(:,2)*stanint(2,2)-coeff(:,1)*stanint(1,2))/stanint(3,2).*phase(1:end-5,3)*coeffweights(3);
for i=1:size(coeff,1)
coeff(i,:)=coeff(i,:)/sum(coeff(i,1:5));
end

ligand2metal=[1 0.75 0.5 0 0];
L2M=(phase(1:end-5,:).*coeff)*ligand2metal';



figure
plot(cBPDC_formic3,L2M,'s','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10);
xlim([0.4 1])
set(gca,'fontsize',18,'xtick',0:0.5:1,'ytick',0:0.5:1);
xlabel('Normalized \DeltaG','fontsize',18)
ylabel('Ligand/Metal Ratio in Sample','fontsize',18)



