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

close all;clear all;clc;
%==========================================================================
% load Dataset
% cur_dir = cd;
% cd 'Dataset'
dataset = 'Glomnewfeat';
% load(dataset);
% label(label == 2) = -1;
% cd(cur_dir);
load('C:\Users\wktey\Dropbox\Research - Archive\Literature\Texture\Code\newfeatvec.mat')
load('C:\Users\wktey\Dropbox\Research - Archive\Literature\Texture\Code\class.mat')
feat(isnan(feat))=0;
% load('impfeat.mat')
% impfeat = [4,5,6,381:383,252:269,50,99,148,274,427,51:94,100:143,149:192,275:318,428:471];
% [Allfeatures mu d] = normalize(Allfeatures);
%==========================================================================
% parameters: 
%      = positive value (the number of active features)
% T     (number of boosting cycles)
% mode = 1 (lasso)
%        2 (elasticnet)

% If STOP is negative, STOP is an integer
%   that determines the desired number of non-zero variables. If STOP is
%   positive, it corresponds to an upper bound on the L1-norm of the BETA
%   coefficients. Setting STOP to zero (default) yields the entire
%   regularization path.

T    = 5;
mode = 2;
stop = -40; 
[LADTreeModel] = LADTree(feat,class,T,mode,stop);

%==========================================================================
cur_dir = cd;
cd 'LADTree_Model';
save([dataset,'_LADTree_model'],'LADTreeModel');
LADTreeModelDraw(LADTreeModel,dataset);
LADTreeModelWrite(LADTreeModel,dataset);
cd(cur_dir);
load('C:\Users\wktey\Dropbox\Research - Archive\Literature\Texture\Code\testingfeatvec.mat')
load('C:\Users\wktey\Dropbox\Research - Archive\Literature\Texture\Code\testclass.mat')
feat_vec(isnan(feat_vec))=0;
% % [feat_vec mu d] = normalize(features);
[rclass,score,finalresult(1)] = LADTree_Model_Evaluation(LADTreeModel,feat,class);