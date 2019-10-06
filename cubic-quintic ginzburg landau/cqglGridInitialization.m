function grid = cqglGridInitialization(nx, nt, L, tmax)

% --- Mesh Spacing & time step
grid.dx = L / (nx - 1);
grid.dt = tmax / (nt - 1);
grid.nx = nx;
grid.nt = nt;

% --- Create arrays to save output
grid.x = linspace(-L, L, nx); 
grid.t = linspace(0, tmax, nt);
grid.U = zeros(nx, nt);

% set boundary conditions
grid.U(:, 1) = 0;
grid.u0 = 0; grid.uL = 0;

end

