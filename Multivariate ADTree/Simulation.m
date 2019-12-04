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

% Preliminary results with simulated data sets which might give better
% insight to the structural advantage

% linearly separable
% Univariate
% Multivariate + linearly separable

close all;clear;clc;
% -------------------------------------------------------------------------
% Dataset [X,y]
%
n   = 1000;
p   = 2;
sig = eye(p);
cor = 0.8;
sig = sig*(1-cor)+cor;
% mu1 = [1,2]; mu2 = [2,1];
mu1 = [2,1]; mu2 = [5,1];
Xpos = mvnrnd(mu1,sig,n/2);
Xneg = mvnrnd(mu2,sig,n/2);


X = [Xpos;Xneg]; y = [ones(n/2,1);-ones(n/2,1)];
[X,m,d] = normalize(X);
Xpos = X(1:n/2,:);
Xneg = X(n/2+1:end,:);
X1 = X(:,1);
X2 = X(:,2);
%}
% -------------------------------------------------------------------------
% XOR Dataset 
% n   = 1000;
% p   = 2;
% sig = eye(p);
% cor = 0;
% sig = sig*(1-cor)+cor;
% % mu1 = [1,2]; mu2 = [2,1];
% mu1 = [2,2]; mu2 = [-2,-2];mu3 = [2,-2]; mu4 = [-2,2];
% Xpos1 = mvnrnd(mu1,sig,n/4);
% Xpos2 = mvnrnd(mu2,sig,n/4);
% Xneg1 = mvnrnd(mu3,sig,n/4);
% Xneg2 = mvnrnd(mu4,sig,n/4);
% 
% X = [Xpos1;Xpos2;Xneg1;Xneg2]; y = [ones(n/2,1);-ones(n/2,1)];
% [X,m,d] = normalize(X);
% Xpos = X(1:n/2,:);
% Xneg = X(n/2+1:end,:);
% X1 = X(:,1);
% X2 = X(:,2);


% -------------------------------------------------------------------------
% Regularized LADTree
[LADTreeModel] = LADTree(X,y,3,0);
LADTreeModel_draw(LADTreeModel)

figure;hold on
plot(Xpos(:,1),Xpos(:,2),'r.');
plot(Xneg(:,1),Xneg(:,2),'b.');
numrule = length(LADTreeModel.rule);
LargeMarginClassifier = false;
beta = [0;0];
if LargeMarginClassifier
    for i = 1 : numrule
        beta = beta + LADTreeModel.rule{i}.beta;
    end
    m    = -beta(1)/beta(2);
    if beta(1) == 0
        plot([min(X1),max(X1)],[0,0],'k')
    elseif beta(2) == 0
        plot([0,0],[min(X2),max(X2)],'k')
    else
        plot([min(X1),max(X1)],m*[min(X1),max(X1)],'k')
    end
else
    off = 0;
    for i = 1 : numrule
        beta = LADTreeModel.rule{i}.beta;
        m    = -beta(1)/beta(2);
        if LADTreeModel.rule{i}.precondition~=1
            if beta(1) == 0
                plot([min(X1),max(X1)],[0,0],'g')
                text(max(X1)+off,0+off,num2str(i));
            elseif beta(2) == 0
                plot([0,0],[min(X2),max(X2)],'g')
                text(0+off,max(X2)+off,num2str(i));
            else
                plot([min(X1),max(X1)],m*[min(X1),max(X1)],'g')
                text(max(X1)+off,m*max(X1)+off,num2str(i));
            end
        else
            if beta(1) == 0
                plot([min(X1),max(X1)],[0,0],'k')
                text(max(X1)+off,0+off,num2str(i));
            elseif beta(2) == 0
                plot([0,0],[min(X2),max(X2)],'k')
                text(0+off,max(X2)+off,num2str(i));
            else
                plot([min(X1),max(X1)],m*[min(X1),max(X1)],'k')
                text(max(X1)+off,m*max(X1)+off,num2str(i));
            end
        end
        off = off + 0.005;
    end
end
hold off
% -------------------------------------------------------------------------
%{
[ADTreeModel] = ADTree_2015([],X,[],y,[],5,1,1);
figure;hold on
plot(Xpos(:,1),Xpos(:,2),'r.');
plot(Xneg(:,1),Xneg(:,2),'b.');
numrule = length(ADTreeModel.stump);
for i = 1 : numrule
    pause
    att = ADTreeModel.stump{i}.decision.att;
    param = ADTreeModel.stump{i}.decision.param;
    if att == 1
        plot([param(3),param(3)],[min(X2),max(X2)],'k');
    else
        plot([min(X1),max(X1)],[param(3),param(3)],'k');
    end
end
hold off
%}



