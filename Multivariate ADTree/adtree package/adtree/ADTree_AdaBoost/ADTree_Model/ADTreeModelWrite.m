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

function [] = ADTreeModelWrite(ADTreeModel,std_name)
stump = ADTreeModel.stump;
%==========================================================================
% open file to write
%==========================================================================
filename  = strcat(std_name,'_rule'); 
fid       = fopen(filename,'w');

     fprintf(fid,'---------------------------------------------------------------------------\n');
     
for i = 1 : length(stump)
    
    Dec   = stump{i}.decision;
    str = strcat('S',num2str(i));
    switch(Dec.att_type)
        case 1
            %==============================================================
            str = strcat(str,'(cat):[',num2str(Dec.att),']=',num2str(Dec.param));
            fprintf(fid,str);  
            %==============================================================
        case 2
            %==============================================================
            str = strcat(str,'(num):');
            if length(Dec.att) == 1
                if Dec.param(end) == 1
                    str = strcat(str,'[',num2str(Dec.att),']-', num2str(Dec.param(3)),'> 0');
                    fprintf(fid,str);
                else
                    str = strcat(str,'[',num2str(Dec.att),']-', num2str(Dec.param(3)),'< 0');
                    fprintf(fid,str);
                end
            else
                if Dec.param(end) == 1
                    %==========================================================
                    for j = 1 : length(Dec.param)-2;
                        str = strcat(str,num2str(Dec.param(j)),'*[',num2str(Dec.att(j)),']+');
                    end
                    
                    str = strcat(str,num2str(Dec.param(end-1)),' > 0');
                    fprintf(fid,str);
                    %==========================================================
                else
                    %==========================================================
                    B       = length(Dec.att);
                    axislen = Dec.param(1:B);
                    axiscen = Dec.param(B+1:2*B);
                    
                    if Dec.param(1) > 0
                        %[(att - axiscen)/axislen]^2 -1 > 0
                        %enclose neg
                        for j = 1 : length(axislen)
                            str = strcat(str,'[([',num2str(Dec.att(j)),']-',num2str(axiscen(j)),')/',num2str(axislen(j)),']2 + ');
                        end
                        str = strcat(str,'-1 > 0');
                        
                    else
                        %[(att - axiscen)/axislen]^2 -1 < 0
                        %enclose pos
                        axislen = -axislen;
                        for j = 1 : length(axislen)
                            str = strcat(str,'[([',num2str(Dec.att(j)),']-',num2str(axiscen(j)),')/',num2str(axislen(j)),']2 + ');
                        end
                        str = strcat(str,'-1 < 0');
                    end
                    fprintf(fid,str);
                    %==========================================================
                end
            end
            %==============================================================
        case 3
            %==============================================================
            str = strcat(str,'(com):');
            if Dec.param(end) == 1
                %==========================================================
                str = strcat(str,num2str(Dec.param(1)),'*r[',num2str(Dec.att),']+',num2str(Dec.param(2)),'*i[',num2str(Dec.att),']+');
                str = strcat(str,num2str(Dec.param(end-1)),' > 0');
                fprintf(fid,str);
                %==========================================================
            else
                %==========================================================
                B       = 2;
                axislen = Dec.param(1:B);
                axiscen = Dec.param(B+1:2*B);
                
                if Dec.param(1) > 0
                    %[(att - axiscen)/axislen]^2 -1 > 0
                    %enclose neg
                    str = strcat(str,'[(r[',num2str(Dec.att),']-',num2str(axiscen(1)),')/',num2str(axislen(1)),']2 + [(i[',num2str(Dec.att),']-',num2str(axiscen(2)),')/',num2str(axislen(2)),']2');
                    str = strcat(str,'-1 > 0');
                    
                else
                    %[(att - axiscen)/axislen]^2 -1 < 0
                    %enclose pos
                    axislen = -axislen;
                    str = strcat(str,'[(r[',num2str(Dec.att),']-',num2str(axiscen(1)),')/',num2str(axislen(1)),']2 + [(i[',num2str(Dec.att),']-',num2str(axiscen(2)),')/',num2str(axislen(2)),']2');
                    str = strcat(str,'-1 < 0');
                end
                fprintf(fid,str);
                %==========================================================
            end
            %==============================================================
    end
    fprintf(fid,'\n');
    fprintf(fid,'---------------------------------------------------------------------------\n');
end
%==========================================================================
% close file
fclose(fid);
%==========================================================================
end
