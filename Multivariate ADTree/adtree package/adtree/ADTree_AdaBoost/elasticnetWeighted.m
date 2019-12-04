%{
    MCopyright (c) 2015, Sok Hong Kuan, Kuang Ye Chow
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

function [b steps info] = elasticnetWeighted(X, y, delta, stop, Gram, storepath, verbose, weight)
[n,p]=size(X);
Gram = [];
% if n is approximately a factor 10 bigger than p it is faster to use a
% precomputed Gram matrix rather than Cholesky factorization when solving
% the partial OLS problem. Make sure the resulting Gram matrix is not
% prohibitively large.
if (n/p) > 10 && p < 1000
  Gram = X'*X + delta*eye(p);
end

[b steps] = larsenWeighted(X, y, delta, stop, Gram, storepath, verbose, weight);

b = (1 + delta)*b;

if nargout == 3 % only compute if asked for
  info.steps = steps;
  if (delta < eps)
    % delta is zero, use minimum-norm solution
    b0 = pinv(X)*y;
  else
    % delta is non-zero, use ridge regression solution
    [U D V] = svd(X, 'econ');
    b0 = V*diag(1./(diag(D).^2 + delta))*D*U'*y;
  end
  penalty0 = sum(abs(b0)); % L1 constraint size of low-bias model
  indices = (1:p)';
  
  if storepath % for entire path
    q = info.steps + 1;
    info.lambda = zeros(1,q);
    info.df = zeros(1,q);
    info.Cp = zeros(1,q);
    info.AIC = zeros(1,q);
    info.BIC = zeros(1,q);
    info.GCV = zeros(1,q);
    info.s = zeros(1,q);
    sigma2e = sum((y - X*b0).^2)/n; % Mean Square Error of low-bias model
    for step = 1:q
      A = indices(b(:,step) ~= 0); % active set
      [U D] = svd(X(:,A), 'econ');
      d2 = diag(D).^2;
      info.df(step) = sum(sum(U.*(U*diag(d2./(d2 + delta)))));
      % compute godness of fit measurements Cp, AIC and BIC
      r = y - X(:,A)*b(A,step); % residuals
      rss = sum(r.^2); % residual sum-of-squares
      info.Cp(step) = rss/sigma2e - n + 2*info.df(step);
      info.AIC(step) = rss + 2*sigma2e*info.df(step);
%       info.AIC(step) = n*log(rss) + 2*info.df(step);
      info.GCV(step) = rss/((n-info.df(step)).^2);
      info.BIC(step) = rss + log(n)*sigma2e*info.df(step);
      % compute L1 penalty constraints s and lambda
      info.s(step) = sum(abs(b(A,step)))/penalty0;
      if (step == 1)
        info.lambda(step) = max(2*abs(X'*y));
      else
        r2 = y - X(:,A)*b(A,step)/(1 + delta);
        info.lambda(step) = median(2*abs(X(:,A)'*r2 - delta/(1 + delta)*b(A,step)));
      end
    end
    
  else % for single solution
    info.steps = steps;
    A = indices(b ~= 0); % active set
    % compute L1 penalty constraints s and lambda at solution
    info.s = sum(abs(b))/penalty0;
    [U D] = svd(X(:,A), 'econ');
    d2 = diag(D).^2;
    info.df = sum(sum(U.*(U*diag(d2./(d2 + delta)))));
    if isempty(A)
      info.lambda = max(2*abs(X'*y));
    else
      r2 = y - X(:,A)*b(A)/(1 + delta);
      info.lambda = median(2*abs(X(:,A)'*r2 - delta/(1 + delta)*b(A)));
    end
  end

end