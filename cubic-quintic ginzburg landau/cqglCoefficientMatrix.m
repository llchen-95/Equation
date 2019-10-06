function coeff = cqglCoefficientMatrix(params, grid)

% --- coefficients of the tridiagonal system
a = (1 - ((grid.dt * -params.v)/(grid.dx)) + grid.dt * ...
    (params.gamma_1/2) * (2 / grid.dx^2)) * ones(grid.nx-2, 1);
b = (- grid.dt * (params.gamma_1/2) * (1/grid.dx^2)) * ones(grid.nx-3, 1);
c = ((grid.dt * -params.v)/(grid.dx)) - (grid.dt * (params.gamma_1/2) ...
    * (1/grid.dx^2)) * ones(grid.nx-3, 1);

% --- fix coefficient boundary nodes
% b(1) = 1; b(end) = 1;
% c(1) = 0; a(end) = 0;

coeff.a = a;
coeff.b = b;
coeff.c = c;
AA = diag(a) + diag(b, 1) + diag(c, -1);
coeff.AA = AA;
end

