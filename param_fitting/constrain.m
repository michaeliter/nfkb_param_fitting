function p = constrain(params,a,b)
%CONSTRAIN Transform values in [-inf, inf] to [a,b]
% a(1,1:5) = log10(a(1,1:5));
% b(1,1:5) = log10(b(1,1:5));
p = (b-a).*exp(params)./(exp(params)+1) + a;
end

