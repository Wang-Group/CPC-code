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
% n dimensional hyperplane
%==========================================================================
function varargout = Weak_Classifier_SLDA(att,label,weight,stop) 

[m1,n1]     = size(att);

index       = find(label ==  1);
jndex       = find(label == -1);

PosWeight   = sum(weight(index));
NegWeight   = sum(weight(jndex));

weight      = weight';
%==========================================================================
% weighted centroid
%==========================================================================
PosMean     = weight(index)*att(index,:)/PosWeight; 
NegMean     = weight(jndex)*att(jndex,:)/NegWeight; 

Centroid    = PosWeight*PosMean + NegWeight*NegMean;

PosDev      = att-repmat(PosMean,[m1,1]); 
NegDev      = att-repmat(NegMean,[m1,1]); 

Dev         = att-repmat(Centroid,[m1,1]);

% %==========================================================================
% % Scatter Matrix
% %==========================================================================
% Sw = ( PosDev(index,:)'*(PosDev(index,:).*repmat(weight(index)',[1,n1])) + ...
%        NegDev(jndex,:)'*(NegDev(jndex,:).*repmat(weight(jndex)',[1,n1])) ); 
% 
% St          = Dev'*(Dev.*repmat(weight',[1,n1]));
% 
% Sb          = St - Sw;
% 
% %==========================================================================
% % Generalized eigendecomposition
% %==========================================================================
% [V,D]       = eig(Sb,Sw);
% 
% [~,kndex]= max(abs(diag(D)));
% 
% param       = V(:,kndex)';
% 
% if param*PosMean' < param*NegMean'
%     param   = -param;
% end
l1 = label;
l2 = label;
l1(l1 == -1) = 0;
l2(l2 ==  1) = 0;
l2(l2 == -1) = 1;
slda_label = [l1,l2];
[param,theta,info] = sldaWeighted(att,slda_label,1e-6,stop,1,100,1e-6,false,weight');

param = param';

param       = param/realsqrt(param*param');

PosDist     = param*att(index,:)';
NegDist     = param*att(jndex,:)';

Mu1         = -mean(PosDist);
Mu2         = -mean(NegDist);

if Mu1 > Mu2
    temp    = Mu2;
    Mu2     = Mu1;
    Mu1     = temp;
    Sigma1  = std(NegDist);
    Sigma2  = std(PosDist);
else
    Sigma1  = std(PosDist);
    Sigma2  = std(NegDist);
end

n = size(param,2);
param(n+1) = (Sigma1*Mu2 + Sigma2*Mu1)/(Sigma1 + Sigma2);

cost     = size(find(label.*(param*[att';ones(1,m1)])'<0),1);

%==========================================================================
if nargout > 1
    varargout{1} = param;
    varargout{2} = cost;
else
    varargout{1} = param;
end
%==========================================================================
end