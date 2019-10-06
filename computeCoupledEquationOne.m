function A = computeCoupledEquationOne(A, B, chi, gamma, beta, delta, xi)
    % obtain the second derivative
    % w.r.t. x for A
    [Ax, Ay] = gradient(A);
    [Axx, Axy] = gradient(Ax);
    
    A = chi * A + gamma * Axx - beta * (abs(A)^2) * A - delta * (abs(A)^4) * A - xi * (abs(B)^2) * A;
end

