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

This folder contains two subfolders:
  1. ADTree_AdaBoost
  2. ADTree_LogitBoost
==============================================================================
Folder: ADTree_AdaBoost
------------------------------------------------------------------------------
- contains two folders:
  1. Dataset (contains the training dataset)
  2. ADTree_Model (the learned model will be saved here)

Inside ADTree_Model, there are two function files:
  1. ADTreeModelDraw 
     - visualization of ADTree model
  2. ADTreeModelWrite
     - describes all the decision nodes of the ADTree model

ADTree.m
- based on AdaBoost, the first adaptive boosting in the literature
- there are three different base learner options for ADTree:
  1. univariate ADTree (Weak_Classifier_Univariate_Threshold)
  2. multivariate ADTree (Weak_Classifier_LDA) 
  3. sparse ADTree (Weak_Classifier_SLDA > sldaWeighted > elasticnetWeighted
                    > larsenWeighted)

Both choldelete and cholinsert are used in larsenWeighted for fast computation

ADTree_Model_Evaluation.m
- return the score for test sample using ADTree model

ADTree_Model_stump_eval.m
- evaluate each decision stump which consists of a decision node and two 
  child prediction nodes

demoADTree.m
- this is an example file on how to build an ADTree based on dataset, the 
  ADTree model learned will be stored in the ADTree_Model

==============================================================================
Folder: ADTree_LogitBoost
------------------------------------------------------------------------------
- contains two folders:
  1. Dataset (contains the training dataset)
  2. LADTree_Model (the learned model will be saved here)

Inside LADTree_Model foler, there are two function files:
  1. LADTreeModelDraw
     - visualization of LADTree model
  2. LADTreeModelWrite
     - describes all the decision nodes of the LADTree model

LADTree.m
- based on LogitBoost, an adaptive Newton algorithm for fitting additive 
  logistic regression model 
- LADTree extends LogitBoost to form a decision tree by introducing locality 
  into the additive logistic regression model
- the proposed regularized LADTree allows any standard linear regularization
  seamlessly without any modification to the available solvers
- this m-file can be modified easily to use any standard solvers. For example,
  the written m-file uses the SpaSM, Matlab toolbox for sparse statistical 
  modeling to perform
  1. lasso 
  2. elasticnet

Both lasso and elasticnet problems are solved using larsen (together with
choldelete and cholinsert for fast computation)

LADTree_Model_Evaluation.m
- return the score for test sample using the LADTree model


