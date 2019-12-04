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

% *************************************************************************
%  LADTree (2 classes):
%  LogitBoost: An adaptive Newton algorithm for fitting an additive logistic
%              regression model
%  LADTree   : An extension to LogitBoost where locality is introduced to
%              the additive linear model
% *************************************************************************

% y : vector of yi values
% X : N x p matrix of xi values
% p : N x 1 vector of fitted probabilities with i-th element p(xi;Bold)
% W : N x N diagonal matrix of weights with i-th element p(xi;Bold)(1-p(xi;Bold))

% probability p(c=1|x) denoted as p = exp(F)./(exp(F)+exp(-F))

function [LADTreeModel] = LADTree(X,label,T,stop)
% Number of training samples and number of features
[N,d] = size (X);
% label yi need to be either {0,1} instead of -1,+1
y  = (label+1)/2;

Epsilon = 1;
% LADTree Model setup
LADTreeModel      = struct('a0',{},'rule',{});
rule              = struct('precondition',{},'beta',{},'a',{},'b',{});

% 1 Initialization
% 1.a weight
w = ones(N,1)/N;
% 1.b probability estimates p(x)
p  = 0.5*ones(N,1);
% 1.c precondition set
P{1} = true(N,1);
% 1.d root decision rule
LADTreeModel(1).a0 = 0.5*log(sum(w(y == 1))/sum(w(y == 0)));

Ft = zeros(N,1);
Emp = zeros(T,1);
for t = 1 : T
    
    % 2.1 compute working response and weights
    Z         =  sqrt((1-p)./p);
    temp      = -sqrt(p./(1-p));
    Z(y == 0) = temp(y == 0);
    
    Z(Z> 4)   =  4;
    Z(Z<-4)   = -4;
    
    
    w = p.*(1-p);
    % 2.2 fit the function f{i}
    
    % IRLS-LARS regularization path
    [b info] = lasso(repmat(w.^0.5,1,d).*X,Z, stop, true, false);
%     [b info] = elasticnet(repmat(w.^0.5,1,d).*X, Z, 1e-3, stop, true, false);
    % Information based model selection
%     [~,in] = min(info.GCV(:,2:end));
    
    
    % 2.3 Compute Z-score
    C         = mat2cell(b(:,2:end),size(b,1),ones(1, size(b,2)-1));
    if false, mytemp = C; clear C; C{1} = mytemp{end}; clear mytemp; end
    NumP      = length(P);
    NumC      = length(C);
    
    wP_Tc1Tc2 = zeros(NumP,NumC);
    wN_Tc1Tc2 = zeros(NumP,NumC);
    wP_Tc1Fc2 = zeros(NumP,NumC);
    wN_Tc1Fc2 = zeros(NumP,NumC);
    Zt        = zeros(NumP,NumC);
    
    for iP    = 1 : NumP
        c1    = P{iP};
        w_Fc1 = sum(w(~c1));
        
        for iC = 1 : NumC
            bvec     = C{iC};
            c2       = X*bvec > 0; 
            Tc1Tc2   = c1    &  c2;
            Tc1Fc2   = c1    & ~c2;
            
            P_Tc1Tc2 = (label>0) & Tc1Tc2;
            N_Tc1Tc2 = (label<0) & Tc1Tc2;
            P_Tc1Fc2 = (label>0) & Tc1Fc2;
            N_Tc1Fc2 = (label<0) & Tc1Fc2;
            
            wP_Tc1Tc2(iP,iC) =  sum(w(P_Tc1Tc2));
            wN_Tc1Tc2(iP,iC) =  sum(w(N_Tc1Tc2));
            wP_Tc1Fc2(iP,iC) =  sum(w(P_Tc1Fc2));
            wN_Tc1Fc2(iP,iC) =  sum(w(N_Tc1Fc2));
        end
        
        Zt(iP,:) = 2*sum(sqrt(wP_Tc1Tc2(iP,:).*wN_Tc1Tc2(iP,:))+...
        sqrt(wP_Tc1Fc2(iP,:).*wN_Tc1Fc2(iP,:)),1) + w_Fc1;
    end

    % 2.4 Decision rule updating
    [MinZt,bestP] = min(Zt,[],1);
    [~    ,bestC] = min(MinZt);
    bestP         = bestP(bestC);
    
    alphapos      = 0.5*log((wP_Tc1Tc2(bestP,bestC)+Epsilon)/(wN_Tc1Tc2(bestP,bestC)+Epsilon));
    alphaneg      = 0.5*log((wP_Tc1Fc2(bestP,bestC)+Epsilon)/(wN_Tc1Fc2(bestP,bestC)+Epsilon));
                
    LADTreeModel(1).rule{t}.precondition = bestP;
    LADTreeModel(1).rule{t}.beta         = b(:,bestC+1);
    LADTreeModel(1).rule{t}.a            = alphapos;
    LADTreeModel(1).rule{t}.b            = alphaneg;
    
    % 2.5 update F(x) = Fm-1(x) + 1/2*fm(x) and p(x) = exp(F)./(exp(F+exp(-F))
    Ft = LADTreeModel_eval(LADTreeModel,X);
    p = exp(Ft)./(exp(Ft)+exp(-Ft));
    p(find(isnan(p) == 1)) = 1;
    p( p>1 ) = 1;
    p( p<0 ) = 0;
    
    % 2.6 Update precondition set
    P{NumP+1} = P{bestP} &  c2;
    P{NumP+2} = P{bestP} & ~c2;
    
    % plot empirical loss at each boosting step
    Emp(t,1) = sum(log(1+exp(-label.*Ft)));
    
end
% plot(1:T,Emp);
end


