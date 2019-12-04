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

function varargout = LADTree_Model_Evaluation(LADTreeModel,te_data,te_label)
result            = struct('tp',{},'tn',{},'fp',{},'fn',{},'acc',{});

rule = LADTreeModel.rule;

score = zeros(size(te_data,1),1);
score(:) = LADTreeModel.a0;

myP  = true(size(te_data,1),1);


for i = 1 : length(rule)
    P    = rule{i}.precondition;
    beta = rule{i}.beta;
    a    = rule{i}.a;
    b    = rule{i}.b;
    

    c = te_data*beta > 0;
    
    score(myP(:,P) &  c) = score(myP(:,P) &  c) + a;
    score(myP(:,P) & ~c) = score(myP(:,P) & ~c) + b;
    
    myP = [myP, myP(:,P) & c, myP(:,P) & ~c];
end

class = sign(score);
%==========================================================================
% compare with the true label : 
% confusion matrite_data
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

if nargout > 1
    varargout{1} = class;
    varargout{2} = score;
    varargout{3} = result;
else
    varargout{1} = score;
end
end 