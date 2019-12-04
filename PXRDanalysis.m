%% load the data
close all 
clear all
foldername ='C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data\xrd-UiO-67\pxrddata\';
cd(foldername)
thedir=dir([foldername '*.ASC']);
Numsample=length(thedir);
data = load(thedir(1).name);
twotheta=3:0.020:20;
stepsize=0.020;
pxrd=zeros(length(thedir),ceil((twotheta(end)-twotheta(1))/stepsize+1));
for i = 1:length(thedir)         
    data = load(thedir(i).name);
    if data(1,1)-twotheta(1)<-0.0001
       data(find(data(:,1)-twotheta(1)<-0.0001),:)=[];
    end
    if data(end,1)-twotheta(end)>0.0001
       data(find(data(:,1)-twotheta(end)>0.0001),:)=[];
    end
    checkstepsize=data(2,1)-data(1,1);
    data(:,2)=data(:,2)/max(data(:,2));
    if abs(checkstepsize-stepsize) < 0.0001
        startpoint=find(twotheta>=data(1,1),1,'first');
           pxrd(i,startpoint:startpoint+size(data,1)-1)=data(:,2);
    end
end


%% load standard spectra
cd 'C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data\xrd-UiO-67\standard\';
clearvars standard
data=load('uio67.xye');
databegin=find(data(:,1)>=3,1,'first');
standard(1,:)=data(databegin:end,2)'/max(data(:,2));
data=load('hcp.xye');
standard(2,:)=data(databegin:end,2)'/max(data(:,2));
data=load('hxl.xye');
standard(3,:)=data(databegin:end,2)'/max(data(:,2));
data=load('ZrFA.xye');
standard(4,:)=data(databegin:end,2)'/max(data(:,2));
data=load('BPDC.xye');
standard(5,:)=data(databegin:end,2)'/max(data(:,2));

figure
plot(twotheta,standard(1,:)');
hold on
plot(twotheta,standard(2,:)');
plot(twotheta,standard(3,:)');
plot(twotheta,standard(4,:)');
plot(twotheta,standard(5,:)');


% figure
% imagesc(twotheta,1:length(thedir),pxrd)
% set(gca,'ydir','normal')
% xlabel('2-theta')
% ylabel('Samples')
% hold on
% plot(twotheta,standard(1,:)'*Numsample,'r','linewidth',1);
% plot(twotheta,standard(2,:)'*Numsample,'b','linewidth',1);
% plot(twotheta,standard(3,:)'*Numsample,'g','linewidth',1);
% plot(twotheta,standard(4,:)'*Numsample,'k','linewidth',1);
% plot(twotheta,standard(5,:)'*Numsample,'y','linewidth',1);
% legend('UiO','hcp-UiO','hxl-UiO','Hf-formate','BPDC');

%% background removal

addpath('C:\Users\Wave\Desktop\desktop\crystalgrowth_simulation\codes_for_publication\xrd-UiO-67\bf')
testline=pxrd(1,:);
signalrange=mean(standard,1);
availablepts=find(signalrange<0.006);

figure
meanpxrd=mean(pxrd,1);
plot(meanpxrd,'b')
hold on
plot(availablepts,meanpxrd(availablepts),'bo')
toadd=[40 82 163:164 411:416 837 851];
toremove=[29:42 63:66 77:101 167:170 194:199 216:228 230:233 237:241 242:248 256:268 276:278 359:372 377:382 402:413 417:418 428:443 460:469 472:475 485:498 504:508 524:535 559:582 586:600 651:674 680:695 745:753 776:796 790:795 799:814 817:835];
availablepts=setdiff(availablepts,toremove);
availablepts=union(availablepts,toadd);


figure
meanpxrd=mean(pxrd,1);
plot(meanpxrd,'b')
hold on
plot(availablepts,meanpxrd(availablepts),'bo')
plot(availablepts,meanpxrd(availablepts),'r')
% hold on
% button=1;
% numpts=0;
% while button==1
% numpts=numpts+1;
% [xpts(numpts) ypts(numpts) button]=ginput(1);
% plot(round(x(numpts)),meanpxrd(round(x(numpts))),'o')
% end

for i=1:size(pxrd,1)
 i
 pxrd(i,:)=smooth(pxrd(i,:),8);
[ycorr(:,i),yfit(:,i)] = bf(pxrd(i,:)',availablepts(1:end),3,'linear'); 
end

potentialerror=[96 99 100 101 175 177 178 186 187 188 193 194 195 224 225];

correctpxrd=pxrd-yfit';
correctpxrd(correctpxrd<0)=0;
% for i=1:size(correctpxrd,1)
% correctpxrd(i,:)=correctpxrd(i,:)/max(correctpxrd(i,:));
% end


%update standard
data=load('UiO.ASC');
databegin=find(data(:,1)>=3,1,'first');
standard(1,:)=data(databegin:end,2)'/max(data(:,2));
data=load('hcp.ASC');
standard(2,:)=data(databegin:end,2)'/max(data(:,2));
data=load('hxl.ASC');
standard(3,:)=data(databegin:end,2)'/max(data(:,2));
data=load('ZrFA.xye');
standard(4,:)=data(databegin:end,2)'/max(data(:,2));
data=load('BPDC.ASC');
standard(5,:)=data(databegin:end,2)'/max(data(:,2));

for i=1:size(standard,1)
 i
 standard(i,:)=smooth(standard(i,:),8);
[yscorr(:,i),ysfit(:,i)] = bf(standard(i,:)',availablepts(1:end),3,'linear'); 
end

correctstand=standard-ysfit';
correctstand(correctstand<0)=0;
for i=1:size(correctstand,1)
correctstand(i,:)=correctstand(i,:)/max(correctstand(i,:));
end

figure
hold on
for i=1:size(correctstand,1)
plot(twotheta,correctstand(i,:))
end

figure
imagesc(twotheta,1:size(pxrd,1),correctpxrd)
set(gca,'ydir','normal');
hold on
% colormap('jet')
xlabel('2\theta')
ylabel('Samples')

plot(twotheta,correctstand(1,:)'*Numsample,'r','linewidth',0.5);
plot(twotheta,correctstand(2,:)'*Numsample,'b','linewidth',0.5);
plot(twotheta,correctstand(3,:)'*Numsample,'g','linewidth',0.5);
plot(twotheta,correctstand(4,:)'*Numsample,'k','linewidth',0.5);
plot(twotheta,correctstand(5,:)'*Numsample,'y','linewidth',0.5);
legend('UiO-67','hcp-UiO-67','hxl-UiO-67','Hf-Formate','H_{2}BPDC');
legend('boxoff')
legend('TextColor','green','Linewidth',2)
set(gca,'fontsize',14);



%% component decomposition
range(1,:)=[3.68 4.34]; % for hcp phase
range(2,:)=[5.00 5.46];% for hcp + hxl
range(3,:)=[5.48 6.1]; % for hcp+UiO
range(4,:)=[6.26 6.96]; % for hcp + UiO
range(5,:)=[7.26 7.58]; % for unknown phase
range(6,:)=[8.50 9.14];% for Hf-FA
range(7,:)=[9.16 9.58]; % for UiO + hxl + hcp
range(8,:)=[9.60 10.34]; % for Hf-FA
range(9,:)=[10.60 11.20]; % for UiO + hcp + hxl
range(10,:)=[11.28 11.70]; % for UiO + hcp + hxl
range(11,:)=[12 12.28]; % for hcp + hxl
range(12,:)=[12.42 12.88]; % for Hf-FA
range(13,:)=[13.18 13.62]; % for Hf-FA
range(14,:)=[14.26 14.6]; % for hxl
range(15,:)=[16.60 17.38]; % for BPDC
range(16,:)=[19.44 19.68]; % for UiO
range(17,:)=[19.76 20]; % for UiO

pxrdint=zeros(size(pxrd,1),size(range,1));
stanint=zeros(size(standard,1),size(range,1));
for i = 1: size(range,1)
    startpoint=find(twotheta>=range(i,1),1,'first');
    endpoint=find(twotheta<=range(i,2),1,'last');
    pxrdint(:,i)=sum(correctpxrd(:,startpoint:endpoint),2);
    stanint(:,i)=sum(correctstand(:,startpoint:endpoint),2);
end
% pxrdint(pxrdint<0)=0;
cutlim=0.03;
coeff=zeros(size(pxrd,1),size(standard,1));
% first round of guess values
coeff(:,5)= pxrdint(:,15)/stanint(5,15);
coeff(coeff(:,5)<cutlim,5)=0;
coeff(:,4)= pxrdint(:,6)/stanint(4,6);
coeff(coeff(:,4)<cutlim,4)=0;
coeff(:,2)= pxrdint(:,1)/stanint(2,1);
coeff(coeff(:,2)<2*cutlim,2)=0;

question1 = setdiff(find(coeff(:,2)<5*cutlim),find(coeff(:,2)==0));
test1 = pxrdint(question1,2)-coeff(question1,2)*stanint(2,2);
test2 = pxrdint(question1,3)-coeff(question1,2)*stanint(2,3);
test3 = pxrdint(question1,4)-coeff(question1,2)*stanint(2,4);
test4 = pxrdint(question1,7)-coeff(question1,2)*stanint(2,7);
test5 = pxrdint(question1,9)-coeff(question1,2)*stanint(2,9);
test6 = pxrdint(question1,10)-coeff(question1,2)*stanint(2,10);
test=(test1<0) + (test2<0) + (test3<0) + (test4<0) + (test5<0) + (test6<0);
coeff(question1(find(test>=2)),2)=0;
question1=question1(find(test<2));

coeff(:,1)= (pxrdint(:,4)-coeff(:,2)*stanint(2,4))/stanint(1,4);
coeff(coeff(:,1)<cutlim,1)=0;

question2 = setdiff(find(coeff(:,1)<0.1),find(coeff(:,1)==0));
getque21=question2(find(coeff(question2,2)>0.1));
coeff(getque21,1)=pxrdint(getque21,end)/stanint(1,end);
getque22=question2(find(coeff(question2,2)==0));
coeff(getque22,1)=(pxrdint(getque22,4)/stanint(1,4)+...
    pxrdint(getque22,10)/stanint(1,10)+...
    pxrdint(getque22,end)/stanint(1,end))/3;
question2=setdiff(question2,getque21);
question2=setdiff(question2,getque22);


coeff(:,3)= (pxrdint(:,2)-coeff(:,2)*stanint(2,2)-coeff(:,1)*stanint(1,2))/stanint(3,2);
coeff(coeff(:,3)<cutlim,3)=0;
question3 = setdiff(find(coeff(:,3)<0.1),find(coeff(:,3)==0));

getque32=question3(find(coeff(question3,2)==0));
getque33=question3(find(coeff(question3,1)==0));
getque34=intersect(getque32,getque33);
coeff(getque34,1)=(pxrdint(getque34,2)/stanint(3,2)+...
    pxrdint(getque34,7)/stanint(3,7)+...
    pxrdint(getque34,9)/stanint(3,9))/3;

getque33=setdiff(getque33,getque34);
coeff(getque34,1)=((pxrdint(getque34,7)-coeff(getque34,2)*stanint(2,7))/stanint(3,7)+...
    (pxrdint(getque34,9)-coeff(getque34,2)*stanint(2,9))/stanint(3,9))/2;

getque32=setdiff(getque32,getque34);
coeff(getque34,1)=(pxrdint(getque34,2)-coeff(getque34,1)*stanint(1,2))/stanint(3,2);

question3=setdiff(question3,getque32);
question3=setdiff(question3,getque33);
question3=setdiff(question3,getque34);

question=union(question2,question3);

%check
simpxrdint=coeff*stanint;
figure
imagesc(1:size(range,1),1:size(pxrd,1),pxrdint-simpxrdint)
set(gca,'ydir','normal');
xlabel('channels')
ylabel('Samples')

% figure
% for i=1:size(correctpxrd,1)
% plot(twotheta,correctpxrd(i,:))
% hold on
% plot(twotheta,pxrd(i,:))
% plot(twotheta,twotheta*0,'r')
% plot(range(4,1)*ones(1,12),-0.1:0.1:1,'r')
% plot(range(4,2)*ones(1,12),-0.1:0.1:1,'r')
% hold off
% title(strcat(num2str(i),'_sample',thedir(i).name))
% legend(num2str(coeff(i,1)),num2str(pxrdint(i,1)));
% legend('Location','southwest');
% pause(1.5)
% end

%% find the sample names
cd ..
[expdata, txt, raw]=xlsread('UiO-67-BPDC-20190304-150-120-100- 20190308.xlsx','A2:B286');

for i = 1:length(thedir)
    filename=thedir(i).name;
   samplename{i} = filename(1:end-4);
end

samplenum=zeros(length(thedir),1);
for i=1:length(thedir)
[LIA, LOCB]=ismember(samplename{i},txt);
if LIA 
    samplenum(i)=LOCB;
end
end

thedir(find(samplenum==0)).name

coeff(samplenum==0,:)=[];
pxrdint(samplenum==0,:)=[];
pxrd(samplenum==0,:)=[];
correctpxrd(samplenum==0,:)=[];
samplename(samplenum==0)=[];
samplenum(samplenum==0)=[];

save('pxrdanalysis.mat');
