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

function [] = LADTreeModelWrite(LADTreeModel,std_name)
rule = LADTreeModel.rule;
%==========================================================================
% open file to write
%==========================================================================
filename  = strcat(std_name,'_rule');
fid       = fopen(filename,'w');

fprintf(fid,'---------------------------------------------------------------------------\n');

for i = 1 : length(rule)
    
    Dec   = rule{i}.beta;
    str = strcat('S',num2str(i),'(num):');
    fprintf(fid,str);
    str = [];
    for p = 1 : length(Dec) -1
        str = strcat(str,num2str(Dec(p)),'*[',num2str(p),']+');
    end
    str = strcat(str,num2str(Dec(end)),'*[',num2str(p+1),'] > 0');
    fprintf(fid,str);
    fprintf(fid,'\n');
    fprintf(fid,'---------------------------------------------------------------------------\n');
end

fclose(fid);
end
