%% load experimental data
clear
close all
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data');
xpdata=xlsread('syntheticconditions.xlsx','AJ2:AP94');
%load the concentration of different species in log scale
H2BPDC = log(expdata(:,1))';
formic= log(expdata(:,2))';
Hf12= log(expdata(:,3))';
Hf12_Hf6= log(expdata(:,4))';
formate_formic= log(expdata(:,5))';
water=log(expdata(:,6)');
ligandmetalratio=expdata(:,7)';
expdata=xlsread('syntheticconditions.xlsx','AT2:AV94');
expthickness=expdata(:,1)';
explateralsize=expdata(:,2)';
expratio=expdata(:,3)';
select1=find(exp(water)+1.847*exp(formic)>0.6131);


addsequence=[0 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
Numpoints=length(addsequence);

H2BPDC = H2BPDC(select1(1))+addsequence;
formic = formic(select1(1))+ 0*addsequence;
Hf12 = Hf12(select1(1))+0*addsequence;
Hf12_Hf6 = Hf12_Hf6(select1(1))+0*addsequence;
formate_formic = formate_formic(select1(1))+ 0*addsequence;
water = water(select1(1))+0*addsequence;
ligandmetalratio = ligandmetalratio(select1(1))+0*addsequence;
expthickness=expthickness(select1(1))+0*addsequence;
explateralsize=explateralsize(select1(1))+0*addsequence;
expratio=expratio(select1(1))+0*addsequence;

%% several parameters
delta=0;%a parameter to adjust the configurational entropy in the MOF per SBU, in the unit of R. The larger delta suggests a more disorderd MOF. delta = 0 for ideal MOFs.
n1=10;
Hf12overHf6rates=10;
addprobability=0.01;

%number of repeats to obtain statiscal measurements of each point
repeats=50;

%the Maximum steps for nucleation
Maxtry =50000;
%the maximum crystal size in this simulation
Maxcrystalsize=500;
%the maximum steps in the simulation
TimeSteps=10000;

nucleationsize = 20; %set the initial step to nucleation, in which SBU removal is suppressed to expedite nulceation speeds.




%store the direction vector
%Hf12_up
direction(1,1,:)=[0 -3^0.5/3 6^0.5/3 0 -3^0.5/3 -6^0.5/3];
direction(1,2,:)=[3^0.5/2*3^0.5/3 1/2*3^0.5/3 6^0.5/3 3^0.5/2*3^0.5/3 1/2*3^0.5/3 -6^0.5/3];
direction(1,3,:)=[-3^0.5/2*3^0.5/3 1/2*3^0.5/3 6^0.5/3 -3^0.5/2*3^0.5/3 1/2*3^0.5/3 -6^0.5/3];
direction(1,4,:)=[1 0 0 1 0 0];
direction(1,5,:)=[-1 0 0 -1 0 0];
direction(1,6,:)=[1/2 3^0.5/2 0 1/2 3^0.5/2 0];
direction(1,7,:)=[-1/2 3^0.5/2 0 -1/2 3^0.5/2 0];
direction(1,8,:)=[1/2 -3^0.5/2 0 1/2 -3^0.5/2 0];
direction(1,9,:)=[-1/2 -3^0.5/2 0 -1/2 -3^0.5/2 0];

%Hf12_down
direction(2,1,:)=[0 3^0.5/3 6^0.5/3 0 3^0.5/3 -6^0.5/3];
direction(2,2,:)=[-3^0.5/2*3^0.5/3 -1/2*3^0.5/3 6^0.5/3 -3^0.5/2*3^0.5/3 -1/2*3^0.5/3 -6^0.5/3];
direction(2,3,:)=[3^0.5/2*3^0.5/3 -1/2*3^0.5/3 6^0.5/3 3^0.5/2*3^0.5/3 -1/2*3^0.5/3 -6^0.5/3];
direction(2,4,:)=[1 0 0 1 0 0];
direction(2,5,:)=[-1 0 0 -1 0 0];
direction(2,6,:)=[1/2 3^0.5/2 0 1/2 3^0.5/2 0];
direction(2,7,:)=[-1/2 3^0.5/2 0 -1/2 3^0.5/2 0];
direction(2,8,:)=[1/2 -3^0.5/2 0 1/2 -3^0.5/2 0];
direction(2,9,:)=[-1/2 -3^0.5/2 0 -1/2 -3^0.5/2 0];


%Hf6_up
direction(3,1,:)=[0 -3^0.5/3 6^0.5/3 0 3^0.5/3 -6^0.5/3];
direction(3,2,:)=[3^0.5/2*3^0.5/3 1/2*3^0.5/3 6^0.5/3 -3^0.5/2*3^0.5/3 -1/2*3^0.5/3 -6^0.5/3];
direction(3,3,:)=[-3^0.5/2*3^0.5/3 1/2*3^0.5/3 6^0.5/3 3^0.5/2*3^0.5/3 -1/2*3^0.5/3 -6^0.5/3];
direction(3,4,:)=[1 0 0 1 0 0];
direction(3,5,:)=[-1 0 0 -1 0 0];
direction(3,6,:)=[1/2 3^0.5/2 0 1/2 3^0.5/2 0];
direction(3,7,:)=[-1/2 3^0.5/2 0 -1/2 3^0.5/2 0];
direction(3,8,:)=[1/2 -3^0.5/2 0 1/2 -3^0.5/2 0];
direction(3,9,:)=[-1/2 -3^0.5/2 0 -1/2 -3^0.5/2 0];

for i=4:9
    direction(3,i,4:6)=20000*[1 1 1];
end

%Hf6_down
direction(4,1,:)=[0 3^0.5/3 6^0.5/3 0 -3^0.5/3 -6^0.5/3];
direction(4,2,:)=[-3^0.5/2*3^0.5/3 -1/2*3^0.5/3 6^0.5/3 3^0.5/2*3^0.5/3 1/2*3^0.5/3 -6^0.5/3];
direction(4,3,:)=[3^0.5/2*3^0.5/3 -1/2*3^0.5/3 6^0.5/3 -3^0.5/2*3^0.5/3 1/2*3^0.5/3 -6^0.5/3];
direction(4,4,:)=[1 0 0 1 0 0];
direction(4,5,:)=[-1 0 0 -1 0 0];
direction(4,6,:)=[1/2 3^0.5/2 0 1/2 3^0.5/2 0];
direction(4,7,:)=[-1/2 3^0.5/2 0 -1/2 3^0.5/2 0];
direction(4,8,:)=[1/2 -3^0.5/2 0 1/2 -3^0.5/2 0];
direction(4,9,:)=[-1/2 -3^0.5/2 0 -1/2 -3^0.5/2 0];

for i=4:9
    direction(4,i,4:6)=20000*[1 1 1];
end

%initial HfSBUpositions
        HfSBUdis=[ 0     0;
                   0     0;
                   0.342/2     0;];%the distance between the two Hf6 in the Hf12 SBU is 0.342 measured in inter SBU distance between adjacent Hf6 linked by BPDC; the 0 is for Hf6 SBU;

        centropos=[0 0 0]';


%% 
%record the doping of Hf6 into the NMOF of Hf12 SBUs
   Zr6dopinglevel=zeros(1,Numpoints);
%record the simulated aspect ratio of the NMOF in the unit of inter-SBU
%distance
   aspect_ratio=zeros(1,Numpoints);
%record the simulated thickness of the NMOF in the unit of inter-SBU
%distance
   thickness=zeros(1,Numpoints);
%record the density of the NMOF in the units of Number of SBU/(inter-SBU
%distance)^3
   density=zeros(1,Numpoints);
%record the growthrate of the NMOF
   growthrate=zeros(1,Numpoints);
%record the NMOF size in the X direction
   Xsize=zeros(1,Numpoints);
%record the NMOF size in the Y direction
   Ysize=zeros(1,Numpoints);
%record the NMOF size in three directions
   crystallength1=zeros(Numpoints,repeats);
   crystallength2=zeros(Numpoints,repeats);
   crystallength3=zeros(Numpoints,repeats);
%record the nucleation shape of the NMOF
   crystalsizeintial=zeros(Numpoints,repeats);
   crystalsize1intial=zeros(Numpoints,repeats);
   crystalsize2intial=zeros(Numpoints,repeats);
   crystalsize3intial=zeros(Numpoints,repeats);
%record the NMOF growth rate in three directions
   crystalgrowth1=zeros(Numpoints,repeats);
   crystalgrowth2=zeros(Numpoints,repeats);
   crystalgrowth3=zeros(Numpoints,repeats);
%record the total growth rate of the NMOF in the unit of SBU/step
   crystalgrowth=zeros(Numpoints,repeats);
%record the size of the NMOF in the unit of number of SBUs
   crystalsize=zeros(Numpoints,repeats);
%record the number of steps of crystal growth
   t=zeros(Numpoints,repeats);
%record temperary thickness and aspectratios
    thicknesstrial=zeros(Numpoints,repeats);
    aspectratiotrial=zeros(Numpoints,repeats);
    
% the process of crystalgrowth
saveformic=formic;
saveformate_formic=formate_formic;

for k=1:length(addsequence) % k represents the k th point
    formic=saveformic+addsequence(k);
    formate_formic=saveformate_formic-addsequence(k);

for l=1:repeats
    trial =1;
Maxtime = TimeSteps;
    %initialize;
        pHf12=exp(Hf12_Hf6(k))*Hf12overHf6rates/(1+Hf12overHf6rates*exp(Hf12_Hf6(k)));%Hf12 possibility, the percentage of Hf12 in the total amounts of Hf12 and Hf6
        crystalsize(k,l)=1;
        %mark if the center is Zr6 or not
        Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
        %set SBU positions
        SBUpos=zeros(crystalsize(k,l),6);
        SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
        SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
        SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
        t(k,l)=0;
        Calculatingactivitycoefficient;%calculating the activity coefficient of the solvent components using NRTL model 
        %initialize the probability of including a Hf12 SBU into the crystal
        p=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
        for i=0:36
            for j=0:36
                                p(i+1,j+1)=AssignProbability6copy20181030(H2BPDC(k),formic(k),formate_formic(k),water(k),i,j,delta,addprobability,Gama(k,:),n1);

                if isnan(p(i+1,j+1))
                    p(i+1,j+1)=1;
                end
            end
        end
        for i=0:36
            for j=0:36
                if p(i+1,j+1)>1
                   p(i+1,j+1)=1;
                end
            end
        end
        
        
         p_Hf6=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
        for i=0:36
            for j=0:36
                                p_Hf6(i+1,j+1)=AssignProbability6copy20181030_Hf6(H2BPDC(k),formic(k),formate_formic(k),water(k),i,j,delta,addprobability,Gama(k,:),n1);

                if isnan(p_Hf6(i+1,j+1))
                    p_Hf6(i+1,j+1)=1;
                end
            end
        end
        for i=0:36
            for j=0:36
                if p_Hf6(i+1,j+1)>1
                   p_Hf6(i+1,j+1)=1;
                end
            end
        end       
% add the probability that the linkage between two SBUs in the solid MOF actually exists  
       pligand = AssignProbabilityforligand(H2BPDC(k),formic(k),formate_formic(k),water(k),Gama(k,:),0);
       pligandhead = AssignProbabilityforligand(H2BPDC(k),formic(k),formate_formic(k),water(k),Gama(k,:),0.5);
% add the kinetic terms
       pkineticsHf6 = 12/18*exp(H2BPDC(k))/exp(max(H2BPDC));
       pkineticsHf12 = exp(H2BPDC(k))/exp(max(H2BPDC));
       
%grow the crystal
trial=0;
while t(k,l) < Maxtime
    index6=[];
    if crystalsize(k,l)< 30
        t(k,l)=0;
    end

    %initialization
            %potential positions of the connected SBUs
            potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
            for i=1:9
                potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));

            %find the connecting SBUs that belong to SBUpos2. use index2 to store their
            %index in SBUpos2 and use internalconnect to record their index in potentialSBUconnect2  
              [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');

                index2(index2==0)=[];
                index2=mod(index2-1,crystalsize(k,l))+1;
            %convert the index in internalconnect to the connected SBU in the crystal
                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize(k,l))+1;% index of the connected SBU in crystal
            %to decide if the connection is to the upper side or the lower side
                internalupperlower=(index3<=size(potentialSBUconnect,1));
        
            %to decide if it is lateral connection or head connection
                internalhead=(mod(floor((index3-1)/crystalsize(k,l)),9)+1<=3);
                internallateral=(1-internalhead);

            %delete unmatched connections 1 for keep the connection, 0 for remove the
                 %connection
                headmatchness = (2-ceil(Hf12vsHf6(index1)/2)).*(2-ceil(Hf12vsHf6(index2)/2)).*(Hf12vsHf6(index1)~=Hf12vsHf6(index2)) ...
                     + (ceil(Hf12vsHf6(index1)/2)-1).*(ceil(Hf12vsHf6(index2)/2)-1).*(Hf12vsHf6(index1)==Hf12vsHf6(index2)) ...
                     + (ceil(Hf12vsHf6(index1)/2)-1).*(2-ceil(Hf12vsHf6(index2)/2)).*(internalupperlower.*( mod(Hf12vsHf6(index1),2) ~= mod(Hf12vsHf6(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6(index1),2) == mod(Hf12vsHf6(index2),2))) ...
                     + (2-ceil(Hf12vsHf6(index1)/2)).*(ceil(Hf12vsHf6(index2)/2)-1).*(internalupperlower.*( mod(Hf12vsHf6(index1),2) == mod(Hf12vsHf6(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6(index1),2) ~= mod(Hf12vsHf6(index2),2)));
                connectmatch = internalhead.*headmatchness + internallateral;

            %delete unmatched SBU connections from index2
                index2(connectmatch==0)=[];
    
     %reinitialize;           
     if crystalsize(k,l) <1          
            crystalsize(k,l)=1;
            %mark if the center is Zr6 or not
            Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
            SBUpos=zeros(crystalsize(k,l),6);
            SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
            SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
            potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
     end
     sumadd=0;
     roundnum = 0;
    while length(index6)~= crystalsize(k,l) || (sumadd~=0 && roundnum <100)
        roundnum = roundnum +1;
        %reinitialize;
        if (crystalsize(k,l)<=1)
            crystalsize(k,l)=1;
            Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
            SBUpos=zeros(crystalsize(k,l),6);
            SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
            SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
                index2=1;
                index1=1;
                index3=1;
        end
        
        %delete isolated ones
            deleteisolate = 1-ismember(1:crystalsize(k,l),index2);
            SBUpos(deleteisolate==1,:) = [];
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
            Hf12vsHf6(deleteisolate==1,:) = [];
            crystalsize(k,l)=crystalsize(k,l)-sum(deleteisolate);
           
        %determine the indexes again
            potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
            for i=1:9
                potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            
              [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');

                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize(k,l))+1;% index of the connected SBU in crystal
                
        %check if the upperlower match, remove those that do not match        
                internalhead=(mod(floor((index3-1)/crystalsize(k,l)),9)+1<=3);
                internallateral=(1-internalhead);
                internalupperlower=(index3<=size(potentialSBUconnect,1));
                index2(index2==0)=[];
                index8=index2;
                index2=mod(index2-1,crystalsize(k,l))+1;
                actualupperlower=(2-ceil(Hf12vsHf6(index2)/2)).*(index8>crystalsize(k,l))+(ceil(Hf12vsHf6(index2)/2)-1).*internalupperlower;
                connectmatch=internalhead.*(internalupperlower==actualupperlower)+internallateral;
                index2(connectmatch==0)=[];
                index1(connectmatch==0)=[];
                index3(connectmatch==0)=[];
                internalupperlower(connectmatch==0)=[];
                actualupperlower(connectmatch==0)=[];
        
        %reinitialize;   
             if (crystalsize(k,l)<=1)
                crystalsize(k,l)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k,l),6);
                SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
                index2=1;
                index1=1;
                index3=1;
                internalconnect = zeros(18,1);
             end
             
        %check if the headtype match, remove those that do not match        
            internalhead=(mod(floor((index3-1)/crystalsize(k,l)),9)+1<=3);
            internallateral=(1-internalhead);
            headmatchness = (2-ceil(Hf12vsHf6(index1)/2)).*(2-ceil(Hf12vsHf6(index2)/2)).*(Hf12vsHf6(index1)~=Hf12vsHf6(index2)) ...
                + (ceil(Hf12vsHf6(index1)/2)-1).*(ceil(Hf12vsHf6(index2)/2)-1).*(Hf12vsHf6(index1)==Hf12vsHf6(index2)) ...
                + (ceil(Hf12vsHf6(index1)/2)-1).*(2-ceil(Hf12vsHf6(index2)/2)).*(internalupperlower.*( mod(Hf12vsHf6(index1),2) ~= mod(Hf12vsHf6(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6(index1),2) == mod(Hf12vsHf6(index2),2))) ...
                + (2-ceil(Hf12vsHf6(index1)/2)).*(ceil(Hf12vsHf6(index2)/2)-1).*(internalupperlower.*( mod(Hf12vsHf6(index1),2) == mod(Hf12vsHf6(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6(index1),2) ~= mod(Hf12vsHf6(index2),2)));
            connectmatch = internalhead.*headmatchness + internallateral;
            index2(connectmatch==0)=[];
            index1(connectmatch==0)=[];
            internalhead(connectmatch==0)=[];
            internallateral(connectmatch==0)=[];
            internalupperlower(connectmatch==0)=[];
      
        %reinitialize;      
            if length(index2) <=1
                crystalsize(k,l)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k,l),6);
                SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
                index2 = 1;
                index1 = 1;
                internalhead =1;
                internallateral=0;
                internalupperlower = 1;
                internalconnect = zeros(18,1);
            end
     
    %find number of interal connections and store in index6;
            [index2, Ia] = sort(index2);
            index1 = index1(Ia);
            internalhead = internalhead(Ia);
            internallateral = internallateral(Ia);
            internalupperlower = internalupperlower(Ia);
            index4 = diff([index2; index2(end)+1]);
            index5 = find(index4>=1);
            index6 = diff([0; index5]);
            sumadd =  sum(abs(sort(index1)-sort(index2)));
    end

    %if the crystal falls apart, keep the larger piece
        index11=index1;
        index22=index2;
        %find one crystal
        if crystalsize(k,l) > 1
           crystal1=[];
           cc = tabulate(index11);
           [yy start]=max(cc(:,2));
           start=find(index11==cc(start,1));
           crystal1=union(crystal1,index11(start),'stable');
           crystal1=union(crystal1,index22(start),'stable');
           index11(start)=[];
           index22(start)=[];
           kk=1;
            while kk<length(crystal1)
                hh=find(index11==crystal1(kk));
                crystal1=union(crystal1,index11(hh),'stable');
                crystal1=union(crystal1,index22(hh),'stable');
                index11(hh)=[];
                index22(hh)=[];
                kk = kk + 1;
            end
            
            %find another crystal
            while ~isempty(index11)
                crystal2=[];
                cc = tabulate(index11);
                [yy start]=max(cc(:,2));
                start=find(index11==cc(start,1));
                crystal2=union(crystal2,index11(start),'stable');
                crystal2=union(crystal2,index22(start),'stable');
                index11(start)=[];
                index22(start)=[];
                kk=1;
                while kk<length(crystal2)
                    hh=find(index11==crystal2(kk));
                    crystal2=union(crystal2,index11(hh),'stable');
                    crystal2=union(crystal2,index22(hh),'stable');
                    index11(hh)=[];
                    index22(hh)=[];
                    kk = kk + 1;
                end

            end
            
            %decide if these two crystals can be linked

            if length(crystal1)>=length(crystal2)
                deleteother = 1-ismember(1:crystalsize(k,l),crystal1);
            else
                deleteother = 1- ismember(1:crystalsize(k,l),crystal2);
            end

            deleteother = deleteother';
        end        
           

    %determine the number of head connections for each node
        internal_n_head = zeros(crystalsize(k,l),1);
        internal_n_head(1) = sum(2*internalhead(1:index5(1)) + internallateral(1:index5(1)).*((Hf12vsHf6(index1(1:index5(1)))>=3)+(Hf12vsHf6(1)>=3)));
        for i=2:crystalsize(k,l)
            internal_n_head(i) = sum(2*internalhead(index5(i-1)+1:index5(i)) + internallateral(index5(i-1)+1:index5(i)).*((Hf12vsHf6(index1(index5(i-1)+1:index5(i)))>=3)+(Hf12vsHf6(i)>=3)));
        end
        internal_n_lateral = 2*index6 - internal_n_head;
     internal_n_head(internal_n_head>36)=36;
     internal_n_lateral(internal_n_lateral>36)=36;
     internal_n_headactual = zeros(crystalsize(k,l),1);
     internal_n_lateralactual = zeros(crystalsize(k,l),1);
     for i=1:crystalsize(k,l)
         internal_n_headactual(i) = sum(randsrc(1,ceil(internal_n_head(i)),[[1 0];[pligandhead 1-pligandhead]]));
         internal_n_lateralactual(i) = sum(randsrc(1,ceil(internal_n_lateral(i)),[[1 0];[pligand 1-pligand]]));
     end
     
    %determine the breaking linkages
        if crystalsize(k,l) > nucleationsize 
              p11= diag(p(round(internal_n_headactual)+1,round(internal_n_lateralactual)+1)).*(Hf12vsHf6<=2);
              p12= diag(p_Hf6(round(internal_n_headactual)+1,round(internal_n_lateralactual)+1)).*(Hf12vsHf6>2);
              p1=p11+p12;

            delete = zeros(crystalsize(k,l),1);
                for i=1:crystalsize(k,l)
                        delete(i) = randsrc(1,1,[[1 0];[p1(i) 1-p1(i)]]);
                end
            if sum(delete)>=crystalsize(k,l)
                    delete(1)=0;
            end
            delete=1-(1-delete).*(1-deleteother);
        else
            delete=zeros(crystalsize(k,l),1);
        end
  

    %get outside connections
    outsideconnect=1-internalconnect;
   
    %get outside connection whether Hf12 or Hf6
        outsideHf12vsHf6=outsideconnect.*randsrc(length(outsideconnect),1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]);
        IsoutsideHf6=(outsideHf12vsHf6>=3);
        outsideConnectHf12vsHf6=zeros(size(outsideHf12vsHf6));
        for i=1:18
            outsideConnectHf12vsHf6(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i)=Hf12vsHf6;
        end
        IsoutsideConnectHf6=(outsideConnectHf12vsHf6>=3);

    %outside connection head vs lateral
        outsideheadconnect=zeros(crystalsize(k,l)*18,1);
        outsideheadconnect([1:3*crystalsize(k,l) 9*crystalsize(k,l)+1:12*crystalsize(k,l)])=1;
        outsidelateralconnect=1-outsideheadconnect;
        upperlower=zeros(size(outsideheadconnect));
        upperlower(1:end/2)=1;%upper connection
        upperlower(1+end/2:end)=0; %lower connection

    %change some of the SBU type to make it logical
    %first term: lateral connection, no need to change the SBU type
    %second term: head connection, two Hf6 SBUs, the direction of the SBU is
    %the same as the one connected to
    %third term: head connection, two Hf12 SBUs, the direction of the SBU is
    %the opposite to the one connected to
    %fourth term: head connection, connected to Hf6, with Hf12
    %fifth term: head connection, connected to Hf12, with Hf6
  
        outsideHf12vsHf6 = outsidelateralconnect.*outsideHf12vsHf6 ... 
                 + outsideheadconnect.*IsoutsideConnectHf6.*IsoutsideHf6.*outsideConnectHf12vsHf6 ...
                 + outsideheadconnect.*(1-IsoutsideConnectHf6).*(1-IsoutsideHf6).*(mod(outsideConnectHf12vsHf6,2)+1) ...
                 + outsideheadconnect.*IsoutsideConnectHf6.*(1-IsoutsideHf6).*((1-upperlower).*(outsideConnectHf12vsHf6-2) + upperlower.*(5-outsideConnectHf12vsHf6)) ...
                 + outsideheadconnect.*(1-IsoutsideConnectHf6).*IsoutsideHf6.*(upperlower.*(outsideConnectHf12vsHf6+2) + (1-upperlower).*(5-outsideConnectHf12vsHf6));
             
    %add SBUs
        add=0;    
        addsum=sum(add);
        growthgoal=1;
     while addsum < growthgoal &&  t(k,l)<=Maxtime
         t(k,l) = t(k,l) + 1;
      %assign connection probability 
      %first term: connect on head site
      %second term: connect on lateral site and connect between a Hf6 and a Hf6
      %fourth term: connect with one internal node and the connection is to
      %lateral site and is connection between a Hf12 and a Hf6
      %fifth term: connect with one internal node and the connection is to
      %lateral site and is connection between a Hf6 and a Hf12
      %sixth term: connect with one internal node and the connection is to
      %lateral site and is connection between a Hf12 and a Hf12

        add = outsideconnect.*outsideheadconnect.*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]) ...
            + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]).*randsrc(length(outsideconnect),1,[[1 0];[pkineticsHf6 1-pkineticsHf6]]) ...
            + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]).*randsrc(length(outsideconnect),1,[[1 0];[pkineticsHf6 1-pkineticsHf6]]) ...
            + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]).*randsrc(length(outsideconnect),1,[[1 0];[pkineticsHf12 1-pkineticsHf12]]) ...
            + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]).*randsrc(length(outsideconnect),1,[[1 0];[pkineticsHf12 1-pkineticsHf12]]);

        %get rid of those connected to the deleted ones
        deletemark=zeros(size(add));
        for i=1:18
            deletemark(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i)=delete;
        end

        add = add - add.*deletemark;

        %get the SBU positions of the attached SBUs
        addSBUs=zeros(sum(add),6);
        select = find(add==1);
        if ~isempty(select)
            addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*(1-upperlower(select))*ones(1,3).*potentialSBUconnect2(select,:);
            addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*(1-upperlower(select))*ones(1,3).*(potentialSBUconnect2(select,:) - ones(length(select),1)*2*HfSBUdis(:,1)');
            addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*upperlower(select)*ones(1,3).*(potentialSBUconnect2(select,:) + ones(length(select),1)*2*HfSBUdis(:,1)');
            addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*upperlower(select)*ones(1,3).*potentialSBUconnect2(select,:);
            addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
            addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
            upordown = randsrc(length(select),1,[[1 0];[1/2 1/2]]);
            addSBUs(:,1:3) = addSBUs(:,1:3) + outsidelateralconnect(select).*(1-IsoutsideHf6(select))*ones(1,3).*(potentialSBUconnect2(select,:) + upordown*2*HfSBUdis(:,1)');
            addSBUs(:,4:6) = addSBUs(:,4:6) + outsidelateralconnect(select).*(1-IsoutsideHf6(select))*ones(1,3).*(potentialSBUconnect2(select,:) - (1-upordown)*2*HfSBUdis(:,1)');
            addSBUs(:,1:3) = addSBUs(:,1:3) + outsidelateralconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
            addSBUs(:,4:6) = addSBUs(:,4:6) + outsidelateralconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
        end
        addHf12vsHf6=outsideHf12vsHf6(select);
    
       % get rid of edge ones
        remove=[find(addSBUs(:,4)>10000); find(addSBUs(:,4)<-10000)];
        if isempty(remove)~=1
            addSBUs(remove,:)=[];  
            addHf12vsHf6(remove)=[];
        end
        addsum=sum(add)-length(remove);
     end
     
  %delete the detached SBUs
    [SBUpos Ia] =setdiff(SBUpos,SBUpos(delete==1,:),'rows','stable');
    Hf12vsHf6 = Hf12vsHf6(Ia);
        %reinitialize;      
            if size(SBUpos,1) <=1
                crystalsize(k,l)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k,l),6);
                SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            end
  %add the attached SBUs
    [SBUpos Ia Ib]=union(SBUpos,addSBUs,'rows','stable');
    if ~isempty(addSBUs)
      Hf12vsHf6 = cat(1,Hf12vsHf6(Ia),addHf12vsHf6(Ib));
    end
    crystalsize(k,l)=size(SBUpos,1);
    
  %find contradictionary
    SBUpos=round(SBUpos*10000)/10000;
    [SBUpos Ia]=sortrows(SBUpos);
    Hf12vsHf6=Hf12vsHf6(Ia);
        SBUjudge=diff(cat(1,SBUpos,SBUpos(end,:)),1);
    judge = (SBUjudge(:,1)==0).*(SBUjudge(:,2)==0).*(abs(SBUjudge(:,3)-2*HfSBUdis(3,1))<1);
    judge(end)=0;
    SBUpos(judge==1,:)=[];
    Hf12vsHf6(judge==1)=[];
        %reinitialize;      
            if size(SBUpos,1) <=1
                crystalsize(k,l)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k,l),6);
                SBUpos(1:crystalsize(k,l),1:3)=centropos(1:crystalsize(k,l)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k,l),4:6)=centropos(1:crystalsize(k,l)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k,l)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k,l)*(i-1)+1:crystalsize(k,l)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            end
            
  %update SBUpos2
    SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
  %update crystalsize(k,l)
    crystalsizeold=crystalsize(k,l);
    crystalsize(k,l)=size(SBUpos,1);
    crystallength1(k,l)=(3*std(SBUpos2(:,1)));
    crystallength2(k,l)=(3*std(SBUpos2(:,2)));
    crystallength3(k,l)=(3*std(SBUpos2(:,3)));
    
    if t(k,l)==1
        crystalsizeintial(k,l)= crystalsize(k,l);
        crystalsize1intial(k,l)= crystallength1(k,l);
        crystalsize2intial(k,l)= crystallength2(k,l);
        crystalsize3intial(k,l)= crystallength3(k,l);
    end

        if crystalsize(k,l) > Maxcrystalsize
            crystalgrowth1(k,l)=(crystallength1(k,l)-crystalsize1intial(k,l))/t(k,l);
            crystalgrowth2(k,l)=(crystallength2(k,l)-crystalsize2intial(k,l))/t(k,l);
            crystalgrowth3(k,l)=(crystallength3(k,l)-crystalsize3intial(k,l))/t(k,l);
            crystalgrowth(k,l)=(crystalsize(k,l)-crystalsizeintial(k,l))/t(k,l);
            t(k,l)=Maxtime;
        end

 fprintf('water=%f, expratio =%f, after %d steps, the crystal size becomes %.1f by %.1f by %.1f, at %d  and crystalsize(k,l) of %d\n',exp(water(k)),expratio(k),t(k,l),crystallength1(k,l),crystallength2(k,l),crystallength3(k,l),k,crystalsize(k,l));
 fprintf('delta, %f \n',delta);
 if crystalsize(k,l)<=nucleationsize
%      t(k,l)=0;
     trial = trial +1;
 end
 if trial >= Maxtry
    Maxtime =0;
 end
end

   Hf6=(Hf12vsHf6>=3);
   Zr6dopinglevel(k,l)=sum(Hf6)/crystalsize(k,l);
   aspectratiotrial(k,l)=crystallength3(k)/max(crystallength1(k,l),crystallength2(k,l));
   thicknesstrial(k,l)=crystallength3(k,l);
   fprintf('H2BPDC = %f, formic = %f, Hf12 = %f \n',H2BPDC(k),formic(k),Hf12(k));
   fprintf('Hf12_Hf6 = %f, formate_formic = %f, water = %f \n',Hf12_Hf6(k),formate_formic(k),water(k));

figure(3)
hold on;
plot(crystalgrowth(k,l),expthickness(k),'*','Color',[k/93 1-k/93 1]);

ylim([0 100]);
xlabel('growthrate');
ylabel('expthicikness');
title(strcat('delta= ',num2str(delta)));
pause(0.01);

end
end

figure
plot(H2BPDC/2.303,log(mean(crystalgrowth,2))/2.303,'*');
figure_FontSize = 14;
ylabel('log (crystal growth rate)','fontsize',figure_FontSize,'fontname','Arial')
xlabel('log (c[H_{2}BPDC])','fontsize',figure_FontSize,'fontname','Arial')
set(gca,'XTick',-0.5:0.2:0.2);
set(gca,'YTick',-0.3:0.1:0);
ylim([-0.4 -0.1])
set(gca,'fontsize',figure_FontSize,'fontname','Arial')

[p s] = polyfit((H2BPDC)'/2.303,log(mean(crystalgrowth,2))/2.303,1);
[y delta]=polyval(p,(formic)/2.303,s);
hold on
plot(formic/2.303,y,'r','linewidth',2.5);
legend('Average over 50 runs',strcat('Slope = ',num2str(p(2))))
legend('boxoff')
data2save=[formic' mean(crystalgrowth,2)];

save('comparegrowthrates_BPDC.txt','-ascii','data2save');
save('BPDCdependent_homogeneous.mat')

