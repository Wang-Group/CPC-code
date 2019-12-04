function [y] = Relufunc(a,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

y=a(1)*(x<=a(2))+(a(1)+(x-a(2))*a(3)).*(x>a(2));

end

