%{
Copyright (c) 2015, Sok Hong Kuan, Kuang Ye Chow
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}

function [] = ADTreeModelDraw(ADTreeModel,std_name)

a0    = ADTreeModel.a0;
stump = ADTreeModel.stump;
Node = 0;
% Node ( put parent node number for each node)
Ptemp = 1;
for i = 1 : length(stump)
    
    D     = stump{i}.P_entry;  %obtain the entry node number for the new stump
    entry = Ptemp(D);
    Loc   = length(Node)+1;    %node number of new stump itself
    Node  = [Node,entry,Loc,Loc];
    Ptemp = [Ptemp, length(Node)-1,length(Node)];
    
end
treeplot(Node,'k.');
xstring = strcat('T stumps :',num2str(length(stump)));
xlabel(xstring);

[x,y] = treelayout(Node);

text(x(1),y(1),num2str(a0),'BackgroundColor',[1,1,0]);

off = 1;
for i = 1 : length(stump)
    Dec    = stump{i}.decision;
    DecStr = strcat('S',num2str(i));
    text(x(1+off),y(1+off),DecStr,'BackgroundColor',[1,.6,1]);   
    text(x(2+off),y(2+off),num2str( stump{i}.a),'BackgroundColor',[1,1,0]);
    text(x(3+off),y(3+off),num2str( stump{i}.b),'BackgroundColor',[1,1,0]);
    off = off + 3;

end
matlabFig  = strcat(std_name,'.fig'); 
saveas(gcf,matlabFig);

end

