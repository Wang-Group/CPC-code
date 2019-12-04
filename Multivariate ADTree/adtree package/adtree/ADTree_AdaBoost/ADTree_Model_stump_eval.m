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

%==========================================================================
% stump evaluation
%==========================================================================
function c = ADTree_Model_stump_eval(att_type,att_label,att,param)
att = att(:,att_label);
switch(att_type)
    case 1
        c = (att == param); 
    case 2
        %------------------------------------------------------------------
        [A,B] = size(att);
        if B == 1
            d = bsxfun(@minus,att,param(3));
            if param(4) == 1
                c = d > 0;
            else
                c = d < 0;
            end
        else
            %--------------------------------------------------------------
            if param(end) == 1
                d = [att,ones(A,1)]*param(1:end-1)';
            else
                R = reshape(param(2*B+1:end-1),B,B)';
                q = (att-repmat(zeta(B+1:2*B),A,1)) *R';
                N = param(1:B);
                if param(1) > 0
                    for i = 1 : size(q,2)
                        q(:,i) = q(:,i)/N(i);
                    end
                    d   = sum(q.^2,2) - 1;
                else
                    for i = 1 : size(q,2)
                        q(:,i) = q(:,i)/N(i);
                    end
                    d   = -(sum(q.^2,2) - 1);
                end
            end
            %--------------------------------------------------------------
            c = d > 0;        
        end
        
        %------------------------------------------------------------------
    case 3
        A = size(att,1);
        B = 2;
        if param(end) == 1
            d = [att,ones(A,1)]*param(1:end-1)';
        else
            R = reshape(param(2*B+1:end-1),B,B)';
            q = (att-repmat(zeta(B+1:2*B),A,1)) *R';
            N = param(1:B);
            if param(1) > 0
                for i = 1 : size(q,2)
                    q(:,i) = q(:,i)/N(i);
                end
                d   = sum(q.^2,2) - 1;
            else
                for i = 1 : size(q,2)
                    q(:,i) = q(:,i)/N(i);
                end
                d   = -(sum(q.^2,2) - 1);
            end
        end
        c = d > 0;
end
end