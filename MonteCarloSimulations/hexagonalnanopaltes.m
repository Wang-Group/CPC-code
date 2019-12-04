%% load experimental data
clear
close all
addpath('C:\Users\Wave\Desktop\desktop\Papers\manuscript-CPC-20190621\crystalgrowth_simulation\data');
expdata=xlsread('20181005.xlsx','AJ2:AP97');
H2BPDC = log(expdata(:,1))';
formic= log(expdata(:,2))';
Hf12= log(expdata(:,3))';
Hf12_Hf6= log(expdata(:,4))';
formate_formic= log(expdata(:,5))';
water=log(expdata(:,6)');
ligandmetalratio=expdata(:,7)';
expdata=xlsread('20181005.xlsx','AT2:AV97');
expthickness=expdata(:,1)';
explateralsize=expdata(:,2)';
expratio=expdata(:,3)';

%% several parameters
% global HfSBUdis crystalsize(k) centropos Hf12vsHf6 pHf12 SBUpos SBUpos2 potentialSBUconnect potentialSBUconnect2 direction 
x5=10;%DMF coordination penalty
x4=2;%water coordination penalty
% n1=(formic-min(formic))/(max(formic)-min(formic))*(2-4);
n1(1:length(H2BPDC))=6;
n2(1:length(H2BPDC))=3;
delta=0.25;%0.75;%delta +3.5816+delta-X; max(X)=3.2926 min(X)=2.8263 current(X)=3.0605
head_penalty=0.5;
entropy_penalty=6.5;
addprobability=0.01;

Maxtry =10000;
Maxcrystalsize=2500;
TimeSteps=20000;
annealstep = 1000;
nucleationsize = 20; %set the initial step to nucleation
select1=[1:10:91 65:69 91:96];


% Maxcrystalsize = 1000; %set the number limit of the SBU in crystal

Numpoints=length(expthickness);
% Maxtime = 10^9;%set the maximum of time for growth of one crystal

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
   Zr6dopinglevel=zeros(1,Numpoints);
   aspect_ratio=zeros(1,Numpoints);
   thickness=zeros(1,Numpoints);
   density=zeros(1,Numpoints);
   growthrate=zeros(1,Numpoints);
   Xsize=zeros(1,Numpoints);
   Ysize=zeros(1,Numpoints);
   crystallength1=zeros(1,Numpoints);
   crystallength2=zeros(1,Numpoints);
   crystallength3=zeros(1,Numpoints);
   crystalsize=zeros(1,Numpoints);
   t=zeros(1,Numpoints);


for k=select1 % k represents the k th point
    trial =1;
Maxtime = TimeSteps;
    %initialize;
        pHf12=exp(Hf12_Hf6(k))/(1+exp(Hf12_Hf6(k)));%Hf12 possibility
        crystalsize(k)=1;
        %mark if the center is Zr6 or not
        Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down

        SBUpos=zeros(crystalsize(k),6);
        SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
        SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
        SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
        t(k)=0;
        %initialize the probability of including a Hf12 SBU into the crystal
        p=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
        for i=0:36
            for j=0:36
%                 p(i+1,j+1)=AssignProbability6copy(H2BPDC(k),formic(k),Hf12(k),formate_formic(k),water(k),i,j,delta);
                                p(i+1,j+1)=AssignProbability6copy20181030(H2BPDC(k),formic(k),formate_formic(k),water(k),Hf12_Hf6(k),i,j,n1(k),delta,head_penalty,x4,x5,addprobability,entropy_penalty);

                if isnan(p(i+1,j+1))
                    p(i+1,j+1)=1;
                end
            end
        end
        for i=0:36
            for j=0:36
%                 if p(i+1,j+1)<1
%                    p(i+1,j+1)=(1-p(i+1,j+1))*(p(i+1,j+1)/(1-p(i+1,j+1))+sum(sum(p(2:i,2:j)./(1-p(2:i,2:j)))));
%                 end
                if p(i+1,j+1)>1
                   p(i+1,j+1)=1;
                end
            end
        end
        
        
         p_Hf6=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
        for i=0:36
            for j=0:36
%                 p_Hf6(i+1,j+1)=AssignProbability6copy(H2BPDC(k),formic(k),Hf12(k),formate_formic(k),water(k),i,j,delta);
                                p_Hf6(i+1,j+1)=AssignProbability6copy20181030_Hf6(H2BPDC(k),formic(k),formate_formic(k),water(k),Hf12_Hf6(k),i,j,n2(k),delta,head_penalty,x4,x5,addprobability,entropy_penalty);

                if isnan(p_Hf6(i+1,j+1))
                    p_Hf6(i+1,j+1)=1;
                end
            end
        end
        for i=0:36
            for j=0:36
%                 if p(i+1,j+1)<1
%                    p(i+1,j+1)=(1-p(i+1,j+1))*(p(i+1,j+1)/(1-p(i+1,j+1))+sum(sum(p(2:i,2:j)./(1-p(2:i,2:j)))));
%                 end
                if p_Hf6(i+1,j+1)>1
                   p_Hf6(i+1,j+1)=1;
                end
            end
        end       
        
        
        
%grow the crystal

while t(k) < Maxtime
    index6=[];
    %initialization
            %potential positions of the connected SBUs
            potentialSBUconnect=zeros(crystalsize(k)*9,6);
            for i=1:9
                potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));

            %find the connecting SBUs that belong to SBUpos2. use index2 to store their
            %index in SBUpos2 and use internalconnect to record their index in potentialSBUconnect2  
            [internalconnect, index2] = ismember(round(round(10000*potentialSBUconnect2)/10),round(round(10000*SBUpos2)/10),'rows');
%         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);

                index2(index2==0)=[];
                index2=mod(index2-1,crystalsize(k))+1;
            %convert the index in internalconnect to the connected SBU in the crystal
                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal
            %to decide if the connection is to the upper side or the lower side
                internalupperlower=(index3<=size(potentialSBUconnect,1));
        
            %to decide if it is lateral connection or head connection
                internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
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
     if crystalsize(k) <1          
            crystalsize(k)=1;
            %mark if the center is Zr6 or not
            Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
            SBUpos=zeros(crystalsize(k),6);
            SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
            SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
            potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
     end
     sumadd=0;
     roundnum = 0;
    while length(index6)~= crystalsize(k) || (sumadd~=0 && roundnum <100)
        roundnum = roundnum +1;
        %reinitialize;
        if (crystalsize(k)<=1)
            crystalsize(k)=1;
            Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
            SBUpos=zeros(crystalsize(k),6);
            SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
            SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
                index2=1;
                index1=1;
                index3=1;
        end
        
        %delete isolated ones
            deleteisolate = 1-ismember(1:crystalsize(k),index2);
            SBUpos(deleteisolate==1,:) = [];
            SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
            Hf12vsHf6(deleteisolate==1,:) = [];
            crystalsize(k)=crystalsize(k)-sum(deleteisolate);
           
        %determine the indexes again
            potentialSBUconnect=zeros(crystalsize(k)*9,6);
            for i=1:9
                potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            
            [internalconnect, index2] = ismember(round(round(10000*potentialSBUconnect2)/10),round(round(10000*SBUpos2)/10),'rows');

%         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);

                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal
                
        %check if the upperlower match, remove those that do not match        
                internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
                internallateral=(1-internalhead);
                internalupperlower=(index3<=size(potentialSBUconnect,1));
                index2(index2==0)=[];
                index8=index2;
                index2=mod(index2-1,crystalsize(k))+1;
                actualupperlower=(2-ceil(Hf12vsHf6(index2)/2)).*(index8>crystalsize(k))+(ceil(Hf12vsHf6(index2)/2)-1).*internalupperlower;
                connectmatch=internalhead.*(internalupperlower==actualupperlower)+internallateral;
                index2(connectmatch==0)=[];
                index1(connectmatch==0)=[];
                index3(connectmatch==0)=[];
                internalupperlower(connectmatch==0)=[];
                actualupperlower(connectmatch==0)=[];
        
        %reinitialize;   
             if (crystalsize(k)<=1)
                crystalsize(k)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k),6);
                SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
                index2=1;
                index1=1;
                index3=1;
                internalconnect = zeros(18,1);
             end
             
        %check if the headtype match, remove those that do not match        
            internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
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
                crystalsize(k)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k),6);
                SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
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
        if crystalsize(k) > 1
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
%                 crystal1z=mean(SBUpos2([crystal1,crystal1+crystalsize(k)],3));
%                 crystal2z=mean(SBUpos2([crystal2,crystal2+crystalsize(k)],3));
%                 diffz=crystal2z-crystal1z;
%                 if abs(diffz)>1
%                     SBUpos(crystal2,[3 6])=SBUpos(crystal2,[3 6])-diffz* ones(length(crystal2),2);
%                 end
%                 crystal1=union(crystal1,crystal2);
            end
            
            %decide if these two crystals can be linked

            if length(crystal1)>=length(crystal2)
                deleteother = 1-ismember(1:crystalsize(k),crystal1);
            else
                deleteother = 1-ismember(1:crystalsize(k),crystal2);
            end
%             if length(crystal1)>=1/2*crystalsize(k)
%                 deleteother = 1-ismember(1:crystalsize(k),crystal1);
%             else
%                 deleteother = ismember(1:crystalsize(k),crystal1);
%             end
            deleteother = deleteother';
        end        
           

    %determine the number of head connections for each node
        internal_n_head = zeros(crystalsize(k),1);
        internal_n_head(1) = sum(2*internalhead(1:index5(1)) + internallateral(1:index5(1)).*((Hf12vsHf6(index1(1:index5(1)))>=3)+(Hf12vsHf6(1)>=3)));
        for i=2:crystalsize(k)
            internal_n_head(i) = sum(2*internalhead(index5(i-1)+1:index5(i)) + internallateral(index5(i-1)+1:index5(i)).*((Hf12vsHf6(index1(index5(i-1)+1:index5(i)))>=3)+(Hf12vsHf6(i)>=3)));
        end
        internal_n_lateral = 2*index6 - internal_n_head;
     internal_n_head(internal_n_head>36)=36;
     internal_n_lateral(internal_n_lateral>36)=36;
    %determine the breaking linkages
%     fprintf('internal_n_head = %d, internal_n_lateral = %d \n',internal_n_head,internal_n_lateral)
%     internal_n_head
%     internal_n_lateral
%     
        if crystalsize(k) > nucleationsize 
              p11= diag(p(round(internal_n_head)+1,round(internal_n_lateral)+1)).*(Hf12vsHf6<=2);
              p12= diag(p_Hf6(round(internal_n_head)+1,round(internal_n_lateral)+1)).*(Hf12vsHf6>2);
              p1=p11+p12;
            delete = zeros(crystalsize(k),1);
                for i=1:crystalsize(k)
                        delete(i) = randsrc(1,1,[[1 0];[p1(i) 1-p1(i)]]);
                end
            if sum(delete)>=crystalsize(k)
                    delete(1)=0;
            end
            delete=1-(1-delete).*(1-deleteother);
        else
            delete=zeros(crystalsize(k),1);
        end
  

    %get outside connections
    outsideconnect=1-internalconnect;
   
    %get outside connection whether Hf12 or Hf6
        outsideHf12vsHf6=outsideconnect.*randsrc(length(outsideconnect),1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]);
        IsoutsideHf6=(outsideHf12vsHf6>=3);
        outsideConnectHf12vsHf6=zeros(size(outsideHf12vsHf6));
        for i=1:18
            outsideConnectHf12vsHf6(crystalsize(k)*(i-1)+1:crystalsize(k)*i)=Hf12vsHf6;
        end
        IsoutsideConnectHf6=(outsideConnectHf12vsHf6>=3);

    %outside connection head vs lateral
        outsideheadconnect=zeros(crystalsize(k)*18,1);
        outsideheadconnect([1:3*crystalsize(k) 9*crystalsize(k)+1:12*crystalsize(k)])=1;
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
     while addsum < growthgoal &&  t(k)<=Maxtime
         t(k) = t(k) + 1;
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
            + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]) ...
            + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]) ...
            + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]) ...
            + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[addprobability 1-addprobability]]);

      %get rid of those connected to the deleted ones
        deletemark=zeros(size(add));
        for i=1:18
            deletemark(crystalsize(k)*(i-1)+1:crystalsize(k)*i)=delete;
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
                crystalsize(k)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k),6);
                SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            end
  %add the attached SBUs
    [SBUpos Ia Ib]=union(SBUpos,addSBUs,'rows','stable');
    if ~isempty(addSBUs)
      Hf12vsHf6 = cat(1,Hf12vsHf6(Ia),addHf12vsHf6(Ib));
    end
    crystalsize(k)=size(SBUpos,1);
    
  %find contradictionary
%     SBUpos=round(SBUpos*10000)/10000;
    [SBUpos Ia]=sortrows(SBUpos);
    Hf12vsHf6=Hf12vsHf6(Ia);
        SBUjudge=diff(cat(1,SBUpos,SBUpos(end,:)),1);
    judge = (SBUjudge(:,1)==0).*(SBUjudge(:,2)==0).*(abs(SBUjudge(:,3)-2*HfSBUdis(3,1))<1);
    judge(end)=0;
    SBUpos(judge==1,:)=[];
    Hf12vsHf6(judge==1)=[];
        %reinitialize;      
            if size(SBUpos,1) <=1
                crystalsize(k)=1;
                Hf12vsHf6=randsrc(1,1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]); % 1 2 3 4 for Hf12_up Hf12_down Hf6_up Hf6_down
                SBUpos=zeros(crystalsize(k),6);
                SBUpos(1:crystalsize(k),1:3)=centropos(1:crystalsize(k)) + HfSBUdis(:,ceil(Hf12vsHf6/2));
                SBUpos(1:crystalsize(k),4:6)=centropos(1:crystalsize(k)) - HfSBUdis(:,ceil(Hf12vsHf6/2)); 
                SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
                potentialSBUconnect=zeros(crystalsize(k)*9,6);
                for i=1:9
                    potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
                end
                potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
            end
            
  %update SBUpos2
    SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
  %update crystalsize(k)
    crystalsizeold=crystalsize(k);
    crystalsize(k)=size(SBUpos,1);
    crystallength1(k)=(3*std(SBUpos2(:,1)));
    crystallength2(k)=(3*std(SBUpos2(:,2)));
    crystallength3(k)=(3*std(SBUpos2(:,3)));

%         if crystalsize(k) > Maxcrystalsize
%             crystallength1(k)=crystallength1(k)*Maxtime/t(k);
%             crystallength2(k)=crystallength2(k)*Maxtime/t(k);
%             crystallength3(k)=(crystallength3(k)-0.7)*Maxtime/t(k)+0.7;
%             t(k)=Maxtime;
%         end

 fprintf('water=%f, expratio =%f, after %d steps, the crystal size becomes %.1f by %.1f by %.1f, at %d  and crystalsize(k) of %d\n',exp(water(k)),expratio(k),t(k),crystallength1(k),crystallength2(k),crystallength3(k),k,crystalsize(k));
 fprintf('delta, n1, entropy penalty = %f %f %f \n',delta,n1(k),entropy_penalty);
 if crystalsize(k)<=nucleationsize
%      t(k)=0;
     trial = trial +1;
 end
 if trial >= Maxtry
    Maxtime =0;
 end
  if crystalsize(k) >= Maxcrystalsize
    Maxtime =0;
 end
 if crystalsize(k) > 10 && mod(t(k),100)==1
    xsize = ceil((crystallength1(k)+0.001)/5)*5;
    xmean=round(mean(SBUpos2(:,1))/5)*5;
    ysize = ceil((crystallength2(k)+0.001)/5)*5;
    ymean=round(mean(SBUpos2(:,2))/5)*5;
    zsize = ceil((crystallength3(k)+0.001)/2)*2;
    zmean=round(mean(SBUpos2(:,3))/2)*2;
    figure(1)
        title(strcat('delta= ',num2str(delta),'head_penalty= ',num2str(head_penalty),'n1= ',num2str(n1(k))));
    subplot(2,2,1)
%         scatter3(SBUpos2(:,1),SBUpos2(:,2),SBUpos2(:,3),'fill');
        shp = alphaShape(SBUpos(:,1),SBUpos(:,2),SBUpos(:,3));
        plot(shp,'FaceColor','y','EdgeColor','none')
        camlight('right')
        lighting gouraud
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        xlim([xmean-xsize xmean+xsize]);
        ylim([ymean-ysize  ymean+ysize]);
        zlim([zmean-zsize zmean+zsize]);
    subplot(2,2,2)
        scatter(SBUpos2(:,1),SBUpos2(:,2),'fill');
        xlabel('X');
        ylabel('Y');
        xlim([xmean-xsize xmean+xsize]);
        ylim([ymean-ysize  ymean+ysize]);
    subplot(2,2,3)
        scatter(SBUpos2(:,1),SBUpos2(:,3),'fill');
        xlabel('X');
        ylabel('Z');
        xlim([xmean-xsize xmean+xsize]);
        ylim([zmean-zsize zmean+zsize]);
    subplot(2,2,4)
        scatter(SBUpos2(:,2),SBUpos2(:,3),'fill');
        xlabel('Y');
        ylabel('Z');
        xlim([ymean-ysize ymean+ysize]);
        ylim([zmean-zsize zmean+zsize]);
    pause(0.0001)

 end
end

   Hf6=(Hf12vsHf6>=3);
   Zr6dopinglevel(k)=sum(Hf6)/crystalsize(k);
   aspect_ratio(k)=crystallength3(k)/max(crystallength1(k),crystallength2(k));
   thickness(k)=crystallength3(k);
   Xsize(k)=crystallength1(k);
   Ysize(k)=crystallength2(k);
   volume=crystallength1(k)*crystallength2(k)*crystallength3(k);
   density(k)=crystalsize(k)/volume;
%    growthrate(k)=crystalsize(k)/t;
   fprintf('aspect_ratio is %f, density is %f\n',aspect_ratio(k),density(k));
   fprintf('H2BPDC = %f, formic = %f, Hf12 = %f \n',H2BPDC(k),formic(k),Hf12(k));
   fprintf('Hf12_Hf6 = %f, formate_formic = %f, water = %f \n',Hf12_Hf6(k),formate_formic(k),water(k));
% figure(2)
% hold on;
% plot(aspect_ratio(k),expratio(k),'*');
% % plot(expthickness,expthickness);
% xlabel('aspect_ratio');
% ylabel('expratio');
% title(strcat('delta= ',num2str(delta),'head_penalty= ',num2str(head_penalty),'n1= ',num2str(n1(k))));
% % title(strcat('delta = ',num2str(delta(mn))));
% pause(0.01);
% 
% figure(3)
% hold on;
% plot(crystallength3(k),expthickness(k),'*');
% % plot(0:1:100,0:1:100);
% % xlim([0 150]);
% ylim([0 100]);
% % plot(expthickness,expthickness);
% xlabel('thickness');
% ylabel('expthicikness');
% title(strcat('delta= ',num2str(delta),'head_penalty= ',num2str(head_penalty),'n1= ',num2str(n1(k))));
% % title(strcat('delta = ',num2str(delta(mn))));
% pause(0.01);
% % if crystalsize(k)<100
% %     aspect_ratio(k)=1000;
% %     break;
% % end

save(strcat('n1=',num2str(n1(k)),'k_',num2str(k),'_nanoplates.mat'));
end


% %% annealing
% % annealstep = 100;
% load('n1=8_nanoplates.mat')
% deleteratio=0.01;
% surfacepercentchange = 0.05;
% 
% for k=95 % k represents the k th point
% %     n1(k)= Zr6dopinglevel(k)*24+(1-Zr6dopinglevel(k))*36;
% n1(k) = 8;
% %     entropy_penalty=0;
% Maxtime = annealstep;
% 
%         %reinitialize the probability of including a Hf12 SBU into the crystal
%         p=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
%         for i=0:36
%             for j=0:36
% %                 p(i+1,j+1)=AssignProbability6copy(H2BPDC(k),formic(k),Hf12(k),formate_formic(k),water(k),i,j,delta(mn));
%                                 p(i+1,j+1)=AssignProbability6copy20181030(H2BPDC(k),formic(k),formate_formic(k),water(k),Hf12_Hf6(k),i,j,n1(k),delta,head_penalty,x4,x5,addprobability,entropy_penalty);
% 
%                 if isnan(p(i+1,j+1))
%                     p(i+1,j+1)=1;
%                 end
%             end
%         end
%         for i=0:36
%             for j=0:36
% %                 if p(i+1,j+1)<1
% %                    p(i+1,j+1)=(1-p(i+1,j+1))*(p(i+1,j+1)/(1-p(i+1,j+1))+sum(sum(p(2:i,2:j)./(1-p(2:i,2:j)))));
% %                 end
%                 if p(i+1,j+1)>1
%                    p(i+1,j+1)=1;
%                 end
%             end
%         end
%         
%         
%          p_Hf6=zeros(37,37);%the first number represents number of head connections %the second one represents number of lateral connections
%         for i=0:36
%             for j=0:36
% %                 p(i+1,j+1)=AssignProbability6copy(H2BPDC(k),formic(k),Hf12(k),formate_formic(k),water(k),i,j,delta(mn));
%                                 p_Hf6(i+1,j+1)=AssignProbability6copy20181030_Hf6(H2BPDC(k),formic(k),formate_formic(k),water(k),Hf12_Hf6(k),i,j,n1(k),delta,head_penalty,x4,x5,addprobability,entropy_penalty);
% 
%                 if isnan(p_Hf6(i+1,j+1))
%                     p_Hf6(i+1,j+1)=1;
%                 end
%             end
%         end
%         for i=0:36
%             for j=0:36
% %                 if p(i+1,j+1)<1
% %                    p(i+1,j+1)=(1-p(i+1,j+1))*(p(i+1,j+1)/(1-p(i+1,j+1))+sum(sum(p(2:i,2:j)./(1-p(2:i,2:j)))));
% %                 end
%                 if p_Hf6(i+1,j+1)>1
%                    p_Hf6(i+1,j+1)=1;
%                 end
%             end
%         end       
%         
%         
%         
% %grow the crystal
% 
% while t(k)-TimeSteps < Maxtime
%     oldcrystalsize=crystalsize(k);
%     t(k)=t(k)+1;
%     index6=[];
%     %initialization
%     %potential positions of the connected SBUs
%         potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
%           end
%             potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));
%             SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
%     %find the connecting SBUs that belong to SBUpos2. use index2 to store their
%         %index in SBUpos2 and use internalconnect to record their index in potentialSBUconnect2  
%             [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%             index2(index2==0)=[];
%             index2=mod(index2-1,crystalsize(k))+1;
%         %convert the index in internalconnect to the connected SBU in the crystal
%             index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%             index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal
%             
%     %get outside connections
%     outsideconnect=1-internalconnect; 
% 
%     %get outside connection whether Hf12 or Hf6
%         outsideHf12vsHf6=outsideconnect.*randsrc(length(outsideconnect),1,[[1 2 3 4];[pHf12/2 pHf12/2 1/2-pHf12/2 1/2-pHf12/2]]);
%         IsoutsideHf6=(outsideHf12vsHf6>=3);
%         outsideConnectHf12vsHf6=zeros(size(outsideHf12vsHf6));
%         for i=1:18
%             outsideConnectHf12vsHf6(crystalsize(k)*(i-1)+1:crystalsize(k)*i)=Hf12vsHf6;
%         end
%         IsoutsideConnectHf6=(outsideConnectHf12vsHf6>=3);
% 
%     %outside connection head vs lateral
%         outsideheadconnect=zeros(crystalsize(k)*18,1);
%         outsideheadconnect([1:3*crystalsize(k) 9*crystalsize(k)+1:12*crystalsize(k)])=1;
%         outsidelateralconnect=1-outsideheadconnect;
%         upperlower=zeros(size(outsideheadconnect));
%         upperlower(1:end/2)=1;%upper connection
%         upperlower(1+end/2:end)=0; %lower connection
% 
%     %change some of the SBU type to make it logical
%         %first term: lateral connection, no need to change the SBU type
%         %second term: head connection, two Hf6 SBUs, the direction of the SBU is
%         %the same as the one connected to
%         %third term: head connection, two Hf12 SBUs, the direction of the SBU is
%         %the opposite to the one connected to
%         %fourth term: head connection, connected to Hf6, with Hf12
%         %fifth term: head connection, connected to Hf12, with Hf6
%   
%         outsideHf12vsHf6 = outsidelateralconnect.*outsideHf12vsHf6 ... 
%                  + outsideheadconnect.*IsoutsideConnectHf6.*IsoutsideHf6.*outsideConnectHf12vsHf6 ...
%                  + outsideheadconnect.*(1-IsoutsideConnectHf6).*(1-IsoutsideHf6).*(mod(outsideConnectHf12vsHf6,2)+1) ...
%                  + outsideheadconnect.*IsoutsideConnectHf6.*(1-IsoutsideHf6).*((1-upperlower).*(outsideConnectHf12vsHf6-2) + upperlower.*(5-outsideConnectHf12vsHf6)) ...
%                  + outsideheadconnect.*(1-IsoutsideConnectHf6).*IsoutsideHf6.*(upperlower.*(outsideConnectHf12vsHf6+2) + (1-upperlower).*(5-outsideConnectHf12vsHf6));
%              
%     %add SBUs
%       %assign connection probability 
%       %first term: connect on head site
%       %second term: connect on lateral site and connect between a Hf6 and a Hf6
%       %fourth term: connect with one internal node and the connection is to
%       %lateral site and is connection between a Hf12 and a Hf6
%       %fifth term: connect with one internal node and the connection is to
%       %lateral site and is connection between a Hf6 and a Hf12
%       %sixth term: connect with one internal node and the connection is to
%       %lateral site and is connection between a Hf12 and a Hf12
% 
%         add = outsideconnect.*outsideheadconnect.*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]) ...
%             + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]) ...
%             + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*IsoutsideHf6.*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]) ...
%             + outsideconnect.*outsidelateralconnect.*IsoutsideConnectHf6.*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]) ...
%             + outsideconnect.*outsidelateralconnect.*(1-IsoutsideConnectHf6).*(1-IsoutsideHf6).*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]);
% 
%       %get the SBU positions of the attached SBUs
%         addSBUs=zeros(sum(add),6);
%         select = find(add==1);
%         if ~isempty(select)
%             addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*(1-upperlower(select))*ones(1,3).*potentialSBUconnect2(select,:);
%             addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*(1-upperlower(select))*ones(1,3).*(potentialSBUconnect2(select,:) - ones(length(select),1)*2*HfSBUdis(:,1)');
%             addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*upperlower(select)*ones(1,3).*(potentialSBUconnect2(select,:) + ones(length(select),1)*2*HfSBUdis(:,1)');
%             addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*(1-IsoutsideHf6(select)).*upperlower(select)*ones(1,3).*potentialSBUconnect2(select,:);
%             addSBUs(:,1:3) = addSBUs(:,1:3) + outsideheadconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
%             addSBUs(:,4:6) = addSBUs(:,4:6) + outsideheadconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
%             upordown = randsrc(length(select),1,[[1 0];[1/2 1/2]]);
%             addSBUs(:,1:3) = addSBUs(:,1:3) + outsidelateralconnect(select).*(1-IsoutsideHf6(select))*ones(1,3).*(potentialSBUconnect2(select,:) + upordown*2*HfSBUdis(:,1)');
%             addSBUs(:,4:6) = addSBUs(:,4:6) + outsidelateralconnect(select).*(1-IsoutsideHf6(select))*ones(1,3).*(potentialSBUconnect2(select,:) - (1-upordown)*2*HfSBUdis(:,1)');
%             addSBUs(:,1:3) = addSBUs(:,1:3) + outsidelateralconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
%             addSBUs(:,4:6) = addSBUs(:,4:6) + outsidelateralconnect(select).*IsoutsideHf6(select)*ones(1,3).*potentialSBUconnect2(select,:);
%         end
%         addHf12vsHf6=outsideHf12vsHf6(select);
%         
%         [addSBUs Ia]=uniquetol(addSBUs,'ByRows',true);
%         addHf12vsHf6=addHf12vsHf6(Ia);
%      
%       %get rid of edge ones
%         remove=[find(addSBUs(:,4)>10000); find(addSBUs(:,4)<-10000)];
%         if isempty(remove)~=1
%             addSBUs(remove,:)=[];  
%             addHf12vsHf6(remove)=[];
%         end
%         addsum=size(addSBUs,1); 
%  
%         
%     %propose to add the attached SBUs
%         [SBUpos_2b Ia Ib]=union(SBUpos,addSBUs,'rows','stable');
% %         if ~isempty(addSBUs)
%             Hf12vsHf6_2b = cat(1,Hf12vsHf6(Ia),addHf12vsHf6(Ib));
% %         end
%         crystalsize(k)=size(SBUpos_2b,1);
%     
%     %find contradictionary
% %         SBUpos_2b=round(SBUpos_2b*10000)/10000;
%         [SBUpos_2b Ia]=sortrows(SBUpos_2b);
%         Hf12vsHf6_2b=Hf12vsHf6_2b(Ia);
%         SBUjudge=diff(cat(1,SBUpos_2b,SBUpos_2b(end,:)),1);
%         judge = (SBUjudge(:,1)==0).*(SBUjudge(:,2)==0).*(abs(SBUjudge(:,3)-2*HfSBUdis(3,1))<1);
%         judge(end)=0;
%         
%         SBUpos_2b(judge==1,:)=[];
%         Hf12vsHf6_2b(judge==1)=[]; 
%         
%     %get rid of some of the added SBUs
%         keep=ismembertol(addSBUs,SBUpos_2b,'ByRows',true);
%         remove=find(keep==0);
%         if isempty(remove)~=1
%             addSBUs(remove,:)=[];  
%             addHf12vsHf6(remove)=[];
%         end
%         addsum=size(addSBUs,1); 
%         
%         
%         %propose to add the attached SBUs
%         [SBUpos_2b Ia Ib]=union(SBUpos,addSBUs,'rows','stable');
% %         if ~isempty(addSBUs)
%             Hf12vsHf6_2b = cat(1,Hf12vsHf6(Ia),addHf12vsHf6(Ib));
% %         end
%         crystalsize(k)=size(SBUpos_2b,1);
%               
%     %delete isolated ones
%     deleteisolate = 1;
%     while sum(deleteisolate>0)
%         potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos_2b + reshape(squeeze(direction(Hf12vsHf6_2b,i,:)),size(SBUpos_2b,1),size(SBUpos_2b,2));
%           end
%         potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));            
%         SBUpos2=cat(1,SBUpos_2b(:,1:3),SBUpos_2b(:,4:6));
%         [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%         index2(index2==0)=[];
%         index8=index2;        
%         index2=mod(index2-1,crystalsize(k))+1;
%         index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%         index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal 
%         
%         deleteisolate = 1-ismember(1:crystalsize(k),index2);
%         SBUpos_2b(deleteisolate==1,:) = [];
%         Hf12vsHf6_2b(deleteisolate==1,:) = [];   
%         crystalsize(k)=size(SBUpos_2b,1);
%     end 
%        
%         
%     %check if the upperlower match, remove those that do not match  
%         internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
%         internallateral=(1-internalhead);
%         internalupperlower=(index3<=size(potentialSBUconnect,1));
%         actualupperlower=(2-ceil(Hf12vsHf6_2b(index2)/2)).*(index8>crystalsize(k))+(ceil(Hf12vsHf6_2b(index2)/2)-1).*internalupperlower;
%         connectmatch=internalhead.*(internalupperlower==actualupperlower)+internallateral;
%         index2(connectmatch==0)=[];
%         index1(connectmatch==0)=[];
%         index3(connectmatch==0)=[];
%         internalupperlower(connectmatch==0)=[];
%         actualupperlower(connectmatch==0)=[];       
%         
%     %check if the headtype match, remove those that do not match    
%        
%         internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
%         internallateral=(1-internalhead);
%         headmatchness = (2-ceil(Hf12vsHf6_2b(index1)/2)).*(2-ceil(Hf12vsHf6_2b(index2)/2)).*(Hf12vsHf6_2b(index1)~=Hf12vsHf6_2b(index2)) ...
%                 + (ceil(Hf12vsHf6_2b(index1)/2)-1).*(ceil(Hf12vsHf6_2b(index2)/2)-1).*(Hf12vsHf6_2b(index1)==Hf12vsHf6_2b(index2)) ...
%                 + (ceil(Hf12vsHf6_2b(index1)/2)-1).*(2-ceil(Hf12vsHf6_2b(index2)/2)).*(internalupperlower.*( mod(Hf12vsHf6_2b(index1),2) ~= mod(Hf12vsHf6_2b(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6_2b(index1),2) == mod(Hf12vsHf6_2b(index2),2))) ...
%                 + (2-ceil(Hf12vsHf6_2b(index1)/2)).*(ceil(Hf12vsHf6_2b(index2)/2)-1).*(internalupperlower.*( mod(Hf12vsHf6_2b(index1),2) == mod(Hf12vsHf6_2b(index2),2))+(1-internalupperlower).*( mod(Hf12vsHf6_2b(index1),2) ~= mod(Hf12vsHf6_2b(index2),2)));
%         connectmatch = internalhead.*headmatchness + internallateral;
%         index2(connectmatch==0)=[];
%         index1(connectmatch==0)=[];
%         internalhead(connectmatch==0)=[];
%         internallateral(connectmatch==0)=[];
%         internalupperlower(connectmatch==0)=[]; 
%         index9=mod(index8-1,crystalsize(k))+1;
%         remove=~ismember(index9,index2);
%         if sum(remove)>0
%             SBUpos_2b(index9(remove==1),:)=[];
%             Hf12vsHf6_2b(index9(remove==1))=[];
%             crystalsize(k)=size(SBUpos_2b,1);
%         end
%         
%     %find number of interal connections and store in index6;
%     %delete isolated ones
%     deleteisolate = 1;
%     while sum(deleteisolate>0)
%          potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos_2b + reshape(squeeze(direction(Hf12vsHf6_2b,i,:)),size(SBUpos_2b,1),size(SBUpos_2b,2));
%           end
%         potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));            
%         SBUpos2=cat(1,SBUpos_2b(:,1:3),SBUpos_2b(:,4:6));
%         [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%         index2(index2==0)=[];
%         index8=index2;        
%         index2=mod(index2-1,crystalsize(k))+1;
%         index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%         index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal 
%         
%         deleteisolate = 1-ismember(1:crystalsize(k),index2);
%         SBUpos_2b(deleteisolate==1,:) = [];
%         Hf12vsHf6_2b(deleteisolate==1,:) = [];   
%         crystalsize(k)=size(SBUpos_2b,1);
%     end       
% 
%      %get rid of some of the added SBUs
%         keep=ismembertol(addSBUs,SBUpos_2b,'ByRows',true);
%         remove=find(keep==0);
%         if isempty(remove)~=1
%             addSBUs(remove,:)=[];  
%             addHf12vsHf6(remove)=[];
%         end
%         addsum=size(addSBUs,1);    
%     
%       %add the attached SBUs
%         [SBUpos Ia Ib]=union(SBUpos,addSBUs,'rows','stable');
% %         if ~isempty(addSBUs)
%             Hf12vsHf6 = cat(1,Hf12vsHf6(Ia),addHf12vsHf6(Ib));
% %         end
%         crystalsize(k)=size(SBUpos,1);
% 
%     %delete isolated ones
%     deleteisolate = 1;      
%     while sum(deleteisolate>0)
%          potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
%           end
%         potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));            
%         SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
%         [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%         index2(index2==0)=[];
%         index8=index2;        
%         index2=mod(index2-1,crystalsize(k))+1;
%         index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%         index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal 
%         
%         deleteisolate = 1-ismember(1:crystalsize(k),index2);
%         SBUpos(deleteisolate==1,:) = [];
%         Hf12vsHf6(deleteisolate==1,:) = [];   
%         crystalsize(k)=size(SBUpos,1);
%     end             
%         
%       %determine the connection again
%         potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
%           end
%         potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));            
%         SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
%         [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%         index2(index2==0)=[];
%         index8=index2;        
%         index2=mod(index2-1,crystalsize(k))+1;
%         index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%         index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal 
%    
%       %obtain internal connections
%         internalhead=(mod(floor((index3-1)/crystalsize(k)),9)+1<=3);
%         internallateral=(1-internalhead);
%         internalupperlower=(index3<=size(potentialSBUconnect,1));
%         
%         [index2, Ia] = sort(index2);
%         index1 = index1(Ia);
%         internalhead = internalhead(Ia);
%         internallateral = internallateral(Ia);
%         internalupperlower = internalupperlower(Ia);
%         index4 = diff([index2; index2(end)+1]);
%         index5 = find(index4>=1);
%         index6 = diff([0; index5]);  
%         
%     %determine the number of head connections for each node
%         internal_n_head = zeros(crystalsize(k),1);
%         internal_n_head(1) = sum(2*internalhead(1:index5(1)) + internallateral(1:index5(1)).*((Hf12vsHf6(index1(1:index5(1)))>=3)+(Hf12vsHf6(1)>=3)));
%         for i=2:crystalsize(k)
%             internal_n_head(i) = sum(2*internalhead(index5(i-1)+1:index5(i)) + internallateral(index5(i-1)+1:index5(i)).*((Hf12vsHf6(index1(index5(i-1)+1:index5(i)))>=3)+(Hf12vsHf6(i)>=3)));
%         end
%         internal_n_lateral = 2*index6 - internal_n_head;
%         internal_n_head(internal_n_head>36)=36;
%         internal_n_lateral(internal_n_lateral>36)=36;
%      
%     %determine the breaking linkages
%          addsum=crystalsize(k)-oldcrystalsize;
%          if addsum <0
%              addsum = 0;
%          end
%          p11= diag(p(round(internal_n_head)+1,round(internal_n_lateral)+1)).*(Hf12vsHf6<=2);
%          p12= diag(p_Hf6(round(internal_n_head)+1,round(internal_n_lateral)+1)).*(Hf12vsHf6>2);
%          p1=p11+p12;
%          p1=p1/sum(p1)*deleteratio*addsum;
%          p1(p1>1)=1;
%          a=1;
%          while sum(p1)<deleteratio*addsum
%             a=a*1.001;
%             p1=p11+p12;
%             p1=p1/sum(p1)*deleteratio*addsum*a;
%             p1(p1>1)=1;
%          end
%          delete = zeros(crystalsize(k),1);
%            for i=1:crystalsize(k)
%               delete(i) = randsrc(1,1,[[1 0];[p1(i) 1-p1(i)]]);
%            end
%      
%   %delete the detached SBUs
%         [SBUpos Ia] =setdiff(SBUpos,SBUpos(delete==1,:),'rows','stable');
%         Hf12vsHf6 = Hf12vsHf6(Ia);
%         crystalsize(k)=size(SBUpos,1);
%         
%     %delete isolated ones
%     deleteisolate = 1;
%     while sum(deleteisolate>0)
%          potentialSBUconnect=zeros(crystalsize(k)*9,6);
%           for i=1:9
%             potentialSBUconnect(crystalsize(k)*(i-1)+1:crystalsize(k)*i,:) = SBUpos + reshape(squeeze(direction(Hf12vsHf6,i,:)),size(SBUpos,1),size(SBUpos,2));
%           end
%         potentialSBUconnect2 = cat(1,potentialSBUconnect(:,1:3),potentialSBUconnect(:,4:6));            
%         SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
%         [internalconnect, index2] = ismember(round(round(10^3*potentialSBUconnect2)/10),round(round(10^3*SBUpos2)/10),'rows');
% %         [internalconnect, index2] = ismembertol(potentialSBUconnect2,SBUpos2,'ByRows',true);
%         index2(index2==0)=[];
%         index8=index2;        
%         index2=mod(index2-1,crystalsize(k))+1;
%         index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%         index1=mod(index3-1,crystalsize(k))+1;% index of the connected SBU in crystal 
%         
%         deleteisolate = 1-ismember(1:crystalsize(k),index2);
%         SBUpos(deleteisolate==1,:) = [];
%         Hf12vsHf6(deleteisolate==1,:) = [];   
%         crystalsize(k)=size(SBUpos,1);
%     end 
%        
%   %update SBUpos2
%     SBUpos2=cat(1,SBUpos(:,1:3),SBUpos(:,4:6));
%   %update crystalsize(k)
%     crystalsize(k)=size(SBUpos,1);
%     crystallength1(k)=(3*std(SBUpos2(:,1)));
%     crystallength2(k)=(3*std(SBUpos2(:,2)));
%     crystallength3(k)=(3*std(SBUpos2(:,3)));
% 
% %         if crystalsize(k) > Maxcrystalsize
% %             crystallength1(k)=crystallength1(k)*Maxtime/t(k);
% %             crystallength2(k)=crystallength2(k)*Maxtime/t(k);
% %             crystallength3(k)=(crystallength3(k)-0.7)*Maxtime/t(k)+0.7;
% %             t(k)=Maxtime;
% %         end
% 
%  fprintf('water=%f, expratio =%f, after %d steps, the crystal size becomes %.1f by %.1f by %.1f, at %d  and crystalsize(k) of %d\n',exp(water(k)),expratio(k),t(k),crystallength1(k),crystallength2(k),crystallength3(k),k,crystalsize(k));
%  fprintf('delta, n1, entropy penalty = %f %f %f \n',delta,n1(k),entropy_penalty);
% 
%  if crystalsize(k) > 10
%     xsize = ceil((crystallength1(k)+0.001)/5)*5;
%     xmean=round(mean(SBUpos2(:,1))/5)*5;
%     ysize = ceil((crystallength2(k)+0.001)/5)*5;
%     ymean=round(mean(SBUpos2(:,2))/5)*5;
%     zsize = ceil((crystallength3(k)+0.001)/2)*2;
%     zmean=round(mean(SBUpos2(:,3))/2)*2;
%     figure(1)
% 
%     subplot(2,2,1)
% %         scatter3(SBUpos2(:,1),SBUpos2(:,2),SBUpos2(:,3),'fill');
%         shp = alphaShape(SBUpos(:,1),SBUpos(:,2),SBUpos(:,3));
%         plot(shp,'FaceColor','y','EdgeColor','none')
% 
%         camlight('right')
%         lighting gouraud
% 
%         xlabel('X');
%         ylabel('Y');
%         zlabel('Z');
%         xlim([xmean-xsize xmean+xsize]);
%         ylim([ymean-ysize  ymean+ysize]);
%         zlim([zmean-zsize zmean+zsize]);
%     subplot(2,2,2)
%         scatter(SBUpos2(:,1),SBUpos2(:,2),'fill');
%         xlabel('X');
%         ylabel('Y');
%         xlim([xmean-xsize xmean+xsize]);
%         ylim([ymean-ysize  ymean+ysize]);
%     subplot(2,2,3)
%         scatter(SBUpos2(:,1),SBUpos2(:,3),'fill');
%         xlabel('X');
%         ylabel('Z');
%         xlim([xmean-xsize xmean+xsize]);
%         ylim([zmean-zsize zmean+zsize]);
%     subplot(2,2,4)
%         scatter(SBUpos2(:,2),SBUpos2(:,3),'fill');
%         xlabel('Y');
%         ylabel('Z');
%         xlim([ymean-ysize ymean+ysize]);
%         ylim([zmean-zsize zmean+zsize]);
%     pause(0.0001)
% 
%  end
% end
% 
% 
% end
% save(strcat('n1=',num2str(n1(k)),'k_',num2str(k),'_refinednanoplates.mat'));
% 
