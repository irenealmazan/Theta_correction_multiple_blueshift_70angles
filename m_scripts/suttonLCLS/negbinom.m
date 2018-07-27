function y = negbinom(k,kbar,m)
%UNTITLED Summary of this function goes here
%   Detailed explanati-on goes here

%y=exp(gammaln(k+m) ... 
%    -gammaln(k+1) ... 
%    -gammaln(m)) ...
%    .*(1+m/kbar).^(-k) ... 
%    .*(1+kbar/m).^(-m);

%same answer, different algebraic form
y= (m./(m+kbar)).^m .* ...
    gamma(m+k)./(gamma(m).*gamma(k+1)) .* ...
    (kbar./(m+kbar)).^k;

end