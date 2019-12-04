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

function [ADTreeModel] = ADTree(data,label,T,mode,GorL,stop)
%==========================================================================
% INITIALIZATION
%--------------------------------------------------------------------------
% m, total number of training instances
m = length(label);
%--------------------------------------------------------------------------
% Initialize Precondition set
P{1}              = true(m,1);
%--------------------------------------------------------------------------
% Initialize ADTree model
ADTreeModel       = struct('a0',{},'stump',{});
stump             = struct('P_entry',{},'decision',{},'a',{},'b',{});
decision          = struct('att_type',{},'att',{},'param',{});
%--------------------------------------------------------------------------
% step 1 : set each training instance's weight
weight            = ones(m,1)./m;
%--------------------------------------------------------------------------
% step 2 : set Precondition set = {true} and obtain the a0
W_Pos             = sum(weight(label ==  1));
W_Neg             = sum(weight(label == -1));
a0                = 0.5*log(W_Pos/W_Neg);
ADTreeModel(1).a0 = a0;
%--------------------------------------------------------------------------
% step 3 : reweight the training instances before start boosting process
weight            = weight.*exp(-bsxfun(@times,a0,label));
%==========================================================================

% INDUCTION


numD   = size(data,2);


for t = 1 : T
    %----------------------------------------------------------------------
    % 1. Weak Classifier Generation : BC (base condition set)
    clear BC
    BC        = struct('att_type',{},'att',{},'param',{});
    
    NumP      = length(P);
    
    switch(mode)
        case 1
            %--------------------------------------------------------------
            % numerical attribute
            for i = 1 : numD
                BC{i}.att_type = 2;
                BC{i}.att      = i;
                BC{i}.param    = ...
                    Weak_Classifier_Univariate_Threshold(data(:,i),label,weight);
            end
            %--------------------------------------------------------------
        case 5
            if GorL == 1
                BC{1}.att_type = 2;
                BC{1}.att      = 1:size(data,2);
                [param,~,~] = Weak_Classifier_LDA([data,label],weight);
                BC{1}.param = [param,1];
            else
                for iP = 1 : NumP
                    c1 = P{iP};
                    BC{iP}.att_type = 2;
                    if length(unique(label(c1,1))) == 1
                        param = zeros(1,size(data,2)+1);
                    else
                        [param,~,~] = Weak_Classifier_LDA([data(c1,:),label(c1,1)],weight(c1,1));
                    end
                    BC{iP}.att      = 1:size(data,2);
                    BC{iP}.param = [param,1];
                end
            end
        case 6
            if GorL == 1
                BC{1}.att_type = 2;
                param = Weak_Classifier_SLDA(data,label,weight,stop);
                BC{1}.att      = 1:size(data,2);
                BC{1}.param = [param,1];
                
            else
                for iP    = 1 : NumP
                    c1    = P{iP};
                    BC{iP}.att_type = 2;
                    if length(unique(label(c1,1))) == 1
                        param = zeros(1,size(data,2)+1);
                    else
                        param = Weak_Classifier_SLDA(data(c1,:),label(c1,1),weight(c1,1),stop);
                    end
                    BC{iP}.att      = 1:size(data,2);
                    BC{iP}.param = [param,1];
                end
            end
    end
    %----------------------------------------------------------------------
    % 2. ADTree Model Combination
    % 2.1 Z-value cost function evaluation
    % 2.2 Identify the best combination
    % 2.3 Expanding two new prediction nodes
    %----------------------------------------------------------------------
    Epsilon = 1;
    %----------------------------------------------------------------------
    % step 2.1 : Z-value cost function
    %----------------------------------------------------------------------
    
    NumBC     = length(BC);
    
    wP_Tc1Tc2 = zeros(NumP,NumBC);
    wN_Tc1Tc2 = zeros(NumP,NumBC);
    wP_Tc1Fc2 = zeros(NumP,NumBC);
    wN_Tc1Fc2 = zeros(NumP,NumBC);
    Zt        = zeros(NumP,NumBC);
    aAll      = zeros(NumP,NumBC);
    bAll      = zeros(NumP,NumBC);
    %------------------------------
    for iP    = 1 : NumP
        c1    = P{iP};
        w_Fc1 = sum(weight(~c1));
        
        if ((mode == 5 || mode == 8) && GorL == 2)|| (mode == 6 && GorL == 2)
            iC = iP;
            bc       = BC{iC};
            switch(bc.att_type)
                case 1
                    c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_cat_att,bc.param);
                case 2
                    c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,data,bc.param);
                case 3
                    c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_com_att,bc.param);
            end
            
            Tc1Tc2   = c1    &  c2;
            Tc1Fc2   = c1    & ~c2;
            
            P_Tc1Tc2 = (label>0) & Tc1Tc2;
            N_Tc1Tc2 = (label<0) & Tc1Tc2;
            P_Tc1Fc2 = (label>0) & Tc1Fc2;
            N_Tc1Fc2 = (label<0) & Tc1Fc2;
            
            wP_Tc1Tc2(iP,iC) =  sum(weight(P_Tc1Tc2));
            wN_Tc1Tc2(iP,iC) =  sum(weight(N_Tc1Tc2));
            wP_Tc1Fc2(iP,iC) =  sum(weight(P_Tc1Fc2));
            wN_Tc1Fc2(iP,iC) =  sum(weight(N_Tc1Fc2));
            aAll(iP,iC)      =  0.5*log((wP_Tc1Tc2(iP,iC)+Epsilon)/(wN_Tc1Tc2(iP,iC)+Epsilon));
            bAll(iP,iC)      =  0.5*log((wP_Tc1Fc2(iP,iC)+Epsilon)/(wN_Tc1Fc2(iP,iC)+Epsilon));
        else
            for iC = 1 : NumBC
                bc       = BC{iC};
                switch(bc.att_type)
                    case 1
                        c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_cat_att,bc.param);
                    case 2
                        c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,data,bc.param);
                    case 3
                        c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_com_att,bc.param);
                end
                
                Tc1Tc2   = c1    &  c2;
                Tc1Fc2   = c1    & ~c2;
                
                P_Tc1Tc2 = (label>0) & Tc1Tc2;
                N_Tc1Tc2 = (label<0) & Tc1Tc2;
                P_Tc1Fc2 = (label>0) & Tc1Fc2;
                N_Tc1Fc2 = (label<0) & Tc1Fc2;
                
                wP_Tc1Tc2(iP,iC) =  sum(weight(P_Tc1Tc2));
                wN_Tc1Tc2(iP,iC) =  sum(weight(N_Tc1Tc2));
                wP_Tc1Fc2(iP,iC) =  sum(weight(P_Tc1Fc2));
                wN_Tc1Fc2(iP,iC) =  sum(weight(N_Tc1Fc2));
                %----------------------------------------------------------
                aAll(iP,iC)      =  0.5*log((wP_Tc1Tc2(iP,iC)+Epsilon)/(wN_Tc1Tc2(iP,iC)+Epsilon));
                bAll(iP,iC)      =  0.5*log((wP_Tc1Fc2(iP,iC)+Epsilon)/(wN_Tc1Fc2(iP,iC)+Epsilon));
                %----------------------------------------------------------
            end
        end
        Zt(iP,:) = 2*sum(sqrt(wP_Tc1Tc2(iP,:).*wN_Tc1Tc2(iP,:))+...
            sqrt(wP_Tc1Fc2(iP,:).*wN_Tc1Fc2(iP,:)),1) + w_Fc1;
        
        %------------------------------------------------------------------
        % get rid of those < 50% accuracy
        %------------------------------------------------------------------
        Ztemp = Zt(iP,:);
        Ztemp((sign(aAll(iP,:))~= sign(bAll(iP,:)))==0) = inf;
        
        aAlltemp = aAll(iP,:);
        bAlltemp = bAll(iP,:);
        Ztemp((aAlltemp == 0) | (bAlltemp == 0)) = inf;
        Zt(iP,:) = Ztemp;
        
        
    end
    %----------------------------------------------------------------------
    % step 2 : Identify the best combination
    %----------------------------------------------------------------------
    [MinZt,bestP] = min(Zt,[],1);
    [~    ,bestC] = min(MinZt);
    bestP         = bestP(bestC);
    %----------------------------------------------------------------------
    % step 3 : Expanding two new prediction nodes
    %----------------------------------------------------------------------
    a  = aAll(bestP,bestC);
    b  = bAll(bestP,bestC);
    
    %======================================================================
    % 3. ADTree Model Updating
    %======================================================================
    
    %----------------------------------------------------------------------
    % Updating P
    %----------------------------------------------------------------------
    bc       = BC{bestC};
    switch(bc.att_type)
        case 1
            c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_cat_att,bc.param);
        case 2
            c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,data,bc.param);
        case 3
            c2       = ADTree_Model_stump_eval(bc.att_type,bc.att,tr_com_att,bc.param);
    end
    P{NumP+1} = P{bestP} &  c2;
    P{NumP+2} = P{bestP} & ~c2;
    
    %----------------------------------------------------------------------
    % Updating ADTree Model
    %----------------------------------------------------------------------
    ADTreeModel(1).stump{t}.P_entry      = bestP   ;
    ADTreeModel(1).stump{t}.decision     = bc      ;
    ADTreeModel(1).stump{t}.a            = a       ;
    ADTreeModel(1).stump{t}.b            = b       ;
    
    %======================================================================
    % 4. Weightage Updating
    %======================================================================
    rt = zeros(size(P{bestP},1),1);
    rt(P{bestP} & P{NumP+1}) = a;
    rt(P{bestP} & P{NumP+2}) = b;
    
    weight = weight.*exp(-rt.*label);
end
end