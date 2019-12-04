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
%  Test            : ADTree Model
%==========================================================================
function [class,score,result] = ADTree_Model_Evaluation(ADTreeModel,te_data,te_label)
result            = struct('tp',{},'tn',{},'fp',{},'fn',{},'acc',{});

a0    = ADTreeModel.a0;
stump = ADTreeModel.stump;

score = zeros(size(te_label,1),1);
score(:) = a0;
myP = true(size(te_label,1),1);


for i = 1 : length(stump) 
    P_entry  = stump{i}.P_entry;
    decision = stump{i}.decision;
    a        = stump{i}.a;
    b        = stump{i}.b;
        
    c = ADTree_Model_stump_eval(decision.att_type,decision.att,te_data,decision.param);
    
    score(myP(:,P_entry) & c) = score(myP(:,P_entry) & c) + a;
    score(myP(:,P_entry) &~c) = score(myP(:,P_entry) &~c) + b;
    myP = [myP, myP(:,P_entry) & c,myP(:,P_entry) & ~c];
end


class = sign(score);
%==========================================================================
% compare with the true label : 
% confusion matrix
%==========================================================================
tp = sum((te_label ==  1) == (class ==  1));
tn = sum((te_label == -1) == (class == -1));
fp = sum((te_label == -1) == (class ==  1));
fn = sum((te_label ==  1) == (class == -1));

result(1).tp  = tp;
result(1).tn  = tn;
result(1).fp  = fp;
result(1).fn  = fn;
result(1).acc = (tp+tn)/(tp+tn+fp+fn);

end