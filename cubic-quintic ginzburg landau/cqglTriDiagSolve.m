function v = cqglTriDiagSolve(d,a,e,f,v)

n = length(d);

% --- forward substitution
v(1) = d(1) / e(1);
for i = 2:n
    v(i) = (d(i) - a(i) * v(i-1))/e(i);
end

for i = n-1:-1:1
   v(i) = v(i) - f(i) * v(i+1); 
end

end

