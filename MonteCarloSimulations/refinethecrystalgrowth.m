load('UiO67_startpoint2.mat');
%% annealing
% for n1k=8:2:18
load('UiO67_startpoint20181218.mat');
n1=n1k;
annealstep = 300;
entropy_penalty2=0;
entropy_penalty=0;
deletepercent = 0.9;
surfacepercentchange = 0.1;
insidepenalty=1;
Hf6left=log(exp(Hf6)-exp(ligand)/6);
% H2BPDC = -9.3833;
H2BPDC=-10.5;
%     n1= 2*mean(index6(index7));  

%     entropy_penalty=0;
Maxtime = annealstep;

%reinitialize the probability of including a Hf12 SBU into the crysta      
        
         p=zeros(25,1);
        for i=0:24
                p(i+1)=AssignProbability6copy20181030_pureHf6_Llimited(Hf6left,i,n1,addprobability,entropy_penalty,entropy_penalty2);
%                 p(i+1)=AssignProbability6copy20181030_pureHf6(H2BPDC,formic,formate_formic,water,i,n1,delta,x4,x5,addprobability,entropy_penalty);
                if isnan(p(i+1))
                    p(i+1)=1;
                end
        end
        for i=0:24
                if p(i+1)>1
                   p(i+1)=1;
                end
        end       
        
% grow the crystal
crystalsizeold=size(SBUpos,1);
while t - TimeSteps < Maxtime
    if crystalsizeold>size(SBUpos,1)+50
        deletepercent=deletepercent - 0.01;
    elseif crystalsizeold < size(SBUpos,1)- 50
        deletepercent=deletepercent + 0.01;
    end  
    t=t+1;
    index6=[];
    crystalsizeold=size(SBUpos,1);
    %initialization
    %determine the indexes again     
        potentialSBUconnect=zeros(crystalsize*12,3);
            for i=1:12
                potentialSBUconnect(crystalsize*(i-1)+1:crystalsize*i,:) = SBUpos + reshape(ones(size(SBUpos,1),1)*squeeze(direction(i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            [internalconnect, ~] = ismember(round(round(1000*potentialSBUconnect)/10),round(round(1000*SBUpos)/10),'rows');
   
    %get outside connections
    outsideconnect=1-internalconnect;
     %add SBUs
      %assign connection probability 
        add = outsideconnect.*randsrc(length(outsideconnect),1,[[1 0];[surfacepercentchange 1-surfacepercentchange]]);
  %add the attached SBUs
        %get the SBU positions of the attached SBUs
        addSBUs=zeros(sum(add),3);
        select = find(add==1);
        if ~isempty(select)
            addSBUs = potentialSBUconnect(select,:);
        end
%         addSBUs = addSBUs(1:end-1,:);
        addsum=size(addSBUs,1);   

       %consider spatial restriction for growth
        shp = alphaShape(SBUpos(:,1),SBUpos(:,2),SBUpos(:,3));
        spatialresitriction = inShape(shp,addSBUs(:,1),addSBUs(:,2),addSBUs(:,3));
        for kk= find(spatialresitriction == 1)
            spatialresitriction(kk) = randsrc(1,1,[[1 0];[insidepenalty 1-insidepenalty]]);
        end
        Ib=find(spatialresitriction == 1);
        addSBUs(Ib,:)=[];
        addsum=size(addSBUs,1);          
        
    SBUpos=union(SBUpos,addSBUs,'rows','stable'); 
    SBUpos=uniquetol(SBUpos,'ByRows',true);
     
    
    crystalsize=size(SBUpos,1);
        
    %determine the indexes again     
        potentialSBUconnect=zeros(crystalsize*12,3);
            for i=1:12
                potentialSBUconnect(crystalsize*(i-1)+1:crystalsize*i,:) = SBUpos + reshape(ones(size(SBUpos,1),1)*squeeze(direction(i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            [internalconnect, index2] = ismembertol(potentialSBUconnect,SBUpos,'ByRows',true);
                index2(internalconnect==0)=[];
                index2=mod(index2-1,crystalsize)+1;
                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize)+1;% index of the connected SBU in crystal

     if ~isempty(index2)
    %find number of interal connections and store in index6;
            [index2, Ia] = sort(index2);
            index1 = index1(Ia);
            index4 = diff([index2; index2(end)+1]);
            index5 = find(index4>=1);
            index6 = diff([0; index5]);
%             index6(index6>12)=12;
     end        
        addsum=crystalsize-crystalsizeold;
     %determine the breaking linkages
            p1=p(2*index6+1);
%               p1=p1/sum(p1)*sum(outsideconnect)*surfacepercentchange; 
              p1=p1/sum(p1)*deletepercent*addsum;
              p1(p1>1)=1;
              a=1;
              while sum(p1)<deletepercent*addsum
                  a=a*1.001;
                  p1=p(2*index6+1);
                  p1=p1/sum(p1)*addsum*a;
                  p1(p1>1)=1;
              end
            delete = zeros(crystalsize,1);
                for i=1:crystalsize
                        delete(i) = randsrc(1,1,[[1 0];[p1(i) 1-p1(i)]]);
                end
   
  %delete the detached SBUs
    [SBUpos Ia] =setdiff(SBUpos,SBUpos(delete==1,:),'rows','stable');
    crystalsize=size(SBUpos,1);

    %potential positions of the connected SBUs
            potentialSBUconnect=zeros(crystalsize*12,3);
            for i=1:12
                potentialSBUconnect(crystalsize*(i-1)+1:crystalsize*i,:) = SBUpos + reshape(ones(size(SBUpos,1),1)*squeeze(direction(i,:)),size(SBUpos,1),size(SBUpos,2));
            end
            %find the connecting SBUs that belong to SBUpos. use index2 to store their
            %index in SBUpos and use internalconnect to record their index in potentialSBUconnect2  
%             [internalconnect, index2] = ismember(round(round(10000*potentialSBUconnect2)/10),round(round(10000*SBUpos)/10),'rows');
            [internalconnect, index2] = ismember(round(round(1000*potentialSBUconnect)/10),round(round(1000*SBUpos)/10),'rows');
                index2(index2==0)=[];
                index2=mod(index2-1,crystalsize)+1;
            %convert the index in internalconnect to the connected SBU in the crystal
                index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
                index1=mod(index3-1,crystalsize)+1;% index of the connected SBU in crystal
        %delete isolated ones
            deleteisolate = 1-ismember(1:crystalsize,index2);
            SBUpos(deleteisolate==1,:) = [];
            crystalsize=crystalsize-sum(deleteisolate);   
    
  %update crystalsize
    crystallength1=(3*std(SBUpos(:,1)));
    crystallength2=(3*std(SBUpos(:,2)));
    crystallength3=(3*std(SBUpos(:,3)));

 fprintf('water=%f, after %d steps , crystalsize of %d\n',exp(water),t,crystalsize);
 fprintf('delta, n1, entropy penalty = %f %f %f \n',delta,n1,entropy_penalty);
    
 if crystalsize > 10
    xsize = ceil((crystallength1+0.001)/5)*5;
    xmean=round(mean(SBUpos(:,1))/5)*5;
    ysize = ceil((crystallength2+0.001)/5)*5;
    ymean=round(mean(SBUpos(:,2))/5)*5;
    zsize = ceil((crystallength3+0.001)/2)*2;
    zmean=round(mean(SBUpos(:,3))/2)*2;
    figure(1)

    subplot(2,2,1)
%         scatter3(SBUpos(:,1),SBUpos(:,2),SBUpos(:,3),'fill');
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
        axis equal;
    subplot(2,2,2)
        scatter(SBUpos(:,1),SBUpos(:,2),'fill');
        xlabel('X');
        ylabel('Y');
        xlim([xmean-xsize xmean+xsize]);
        ylim([ymean-ysize  ymean+ysize]);
        axis equal;
    subplot(2,2,3)
        scatter(SBUpos(:,1),SBUpos(:,3),'fill');
        xlabel('X');
        ylabel('Z');
        xlim([xmean-xsize xmean+xsize]);
        ylim([zmean-zsize zmean+zsize]);
        axis equal;
    subplot(2,2,4)
        scatter(SBUpos(:,2),SBUpos(:,3),'fill');
        xlabel('Y');
        ylabel('Z');
        xlim([ymean-ysize ymean+ysize]);
        ylim([zmean-zsize zmean+zsize]);
        axis equal;
    pause(0.0001)

 end
end

% 
%    %draw a shell figure
%       ZrSBUpositions = SBUpos;
%    
%             %potential positions of the connected SBUs
%             potentialSBUconnect=zeros(crystalsize*12,3);
%             for i=1:12
%                 potentialSBUconnect(crystalsize*(i-1)+1:crystalsize*i,:) = SBUpos + reshape(ones(size(SBUpos,1),1)*squeeze(direction(i,:)),size(SBUpos,1),size(SBUpos,2));
%             end
% 
%             %find the connecting SBUs that belong to SBUpos. use index2 to store their
%             %index in SBUpos and use internalconnect to record their index in potentialSBUconnect  
%             [internalconnect, index2] = ismember(round(round(1000*potentialSBUconnect)/10),round(round(1000*SBUpos)/10),'rows');
%                 index2(index2==0)=[];
%                 index2=mod(index2-1,crystalsize)+1;
%             %convert the index in internalconnect to the connected SBU in the crystal
%                 index3=find(internalconnect==1); % the corresponding index in potentialSBUconnect2
%                 index1=mod(index3-1,crystalsize)+1;% index of the connected SBU in crystal   
%             %find number of interal connections and store in index6;
%             [index2, Ia] = sort(index2);
%             index1 = index1(Ia);
%             index4 = diff([index2; index2(end)+1]);
%             index5 = find(index4>=1);
%             index6 = diff([0; index5]);
%             index7=find(index6<12);
%             outershell=SBUpos(index7,:);
%             outlineSBU=SBUpos(index7,:);
%             index8=find(index6<6);
%             outershelledge=SBUpos(index7,:);
%             edgeSBU=SBUpos(index8,:);
%             figure(2)
%                 plot3(SBUpos(:,1),SBUpos(:,2),SBUpos(:,3));
%                 xlabel('X');
%                 ylabel('Y');
%                 zlabel('Z');
%                 xlim([xmean-xsize xmean+xsize]);
%                 ylim([ymean-ysize  ymean+ysize]);
%                 zlim([zmean-zsize zmean+zsize]);
%              axis equal;


save(strcat('n1=',num2str(n1),'_20181218.mat'));
% save('20181218_1.mat');


% draw a surface figure
% coarsedraw=30;
%                 x=linspace(min(SBUpos(:,1))-2,max(SBUpos(:,1))+2,coarsedraw);
%                 y=linspace(min(SBUpos(:,2))-2,max(SBUpos(:,2))+2,coarsedraw);
%                 z=linspace(min(SBUpos(:,3))-2,max(SBUpos(:,3))+2,coarsedraw);                
%                 data = zeros(length(x),length(y),length(z));
%              ZrSBUpositions = sortrows(SBUpos,[1 2 3]);
%                 for m1=1:size(data,1)-1
%                     for m2=1:size(data,2)-1
%                         for m3=1:size(data,3)-1
%                             selectSBUs1=intersect(find(ZrSBUpositions(:,1)>=x(m1)),find(ZrSBUpositions(:,1)<x(m1+1)));
%                             selectSBUs2=intersect(find(ZrSBUpositions(:,2)>=y(m2)),find(ZrSBUpositions(:,2)<y(m2+1)));
%                             selectSBUs3=intersect(find(ZrSBUpositions(:,3)>=z(m3)),find(ZrSBUpositions(:,3)<z(m3+1)));
%                             selectSBUs=intersect(selectSBUs1,selectSBUs2);
%                             selectSBUs=intersect(selectSBUs,selectSBUs3);
%                             data(m1,m2,m3)=length(selectSBUs);
%                         end
%                     end
%                 end
%               
% data = interp3(data,3,'linear');
% figure
% fv = isosurface(data,0.1);
% % p2 = patch(fv,'FaceColor','red','EdgeColor','none');
% % isonormals(data,p2)
% p1 = patch(fv,'FaceColor','red','EdgeColor','none');
% view(3)
% daspect([1,1,1])
% axis tight
% camlight
% camlight(-80,-10)
% lighting gouraud
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
figure
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
title(strcat('n1 = ',num2str(n1)));
% end
