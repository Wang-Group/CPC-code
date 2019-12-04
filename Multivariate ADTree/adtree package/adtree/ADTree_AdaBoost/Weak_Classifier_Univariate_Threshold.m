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

function varargout = Weak_Classifier_Univariate_Threshold(att,label,weight)
param(1) = 1;
param(2) = 0;
param(3) = 0;

index  = find(label >  0);
jndex  = find(label <= 0);

PosWeight = sum(weight(index));
NegWeight = sum(weight(jndex));

weight      = weight';

PosMean    = weight(index)*att(index,:)/PosWeight;
NegMean    = weight(jndex)*att(jndex,:)/NegWeight;

PosCentroid  = PosMean;
NegCentroid  = NegMean;

PosDev       = att-PosCentroid;
NegDev       = att-NegCentroid;

if PosCentroid > NegCentroid
    Mu1      = NegCentroid;
    Mu2      = PosCentroid;
    Sigma1   = std(NegDev);
    Sigma2   = std(PosDev);
    param(4) = 1;
else
    Mu1      = PosCentroid;
    Mu2      = NegCentroid;
    Sigma1   = std(PosDev);
    Sigma2   = std(NegDev);
    param(4) = 2;
end

param(3) = (Sigma1*Mu2 + Sigma2*Mu1)/(Sigma1 + Sigma2);
cost = size(find(label.*(att > param(3))<0),1);
%==========================================================================
if nargout > 1
    varargout{1} = param;
    varargout{2} = cost;
else
    varargout{1} = param;
end
%==========================================================================
end