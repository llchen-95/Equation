function [e,f] = cqglTridiagLU(coeff)

n = length(coeff.a);
e = zeros(n, 1); f = e;

e(1) = coeff.b(1);
f(1) = coeff.c(1) / coeff.b(1);

for i = 2:n
   e(i) = coeff.b(i) - coeff.a(i) * f(i-1);
   f(i) = coeff.c(i) / e(i);
end

end

