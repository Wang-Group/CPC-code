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
function varargout = Weak_Classifier_LDA(Dataset,dist)
att         = Dataset(:,1:end-1);
label       = Dataset(:,end);

[m,n]     = size(att);

index       = find(label ==  1);
jndex       = find(label == -1);

Posdist   = sum(dist(index));
Negdist   = sum(dist(jndex));

dist      = dist';
%==========================================================================
% disted centroid
%==========================================================================
PosMean     = dist(index)*att(index,:)/Posdist; 
NegMean     = dist(jndex)*att(jndex,:)/Negdist; 

Centroid    = Posdist*PosMean + Negdist*NegMean;

PosDev      = att-repmat(PosMean,[m,1]); 
NegDev      = att-repmat(NegMean,[m,1]); 

Dev         = att-repmat(Centroid,[m,1]);

%==========================================================================
% Scatter Matrix
%==========================================================================
Sw = ( PosDev(index,:)'*(PosDev(index,:).*repmat(dist(index)',[1,n])) + ...
       NegDev(jndex,:)'*(NegDev(jndex,:).*repmat(dist(jndex)',[1,n])) ); 

St          = Dev'*(Dev.*repmat(dist',[1,n]));

Sb          = St - Sw;

%==========================================================================
% Generalized eigendecomposition
%==========================================================================
[V,D]       = eig(Sb,Sw);

[~,kndex]= max(abs(diag(D)));

param       = V(:,kndex)';

if param*PosMean' < param*NegMean'
    param   = -param;
end

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

cost     = size(find(label.*(param*[att';ones(1,m)])'<0),1);

%==========================================================================
varargout{1} = param;
if nargout > 1 
    varargout{2} = cost;
    if nargout > 2
        varargout{3} = @(x)sign([x,ones(size(x,1),1)]*param');
    end
end
%==========================================================================
end