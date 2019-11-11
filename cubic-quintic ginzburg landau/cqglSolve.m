function grid = cqglSolve(caseType, caseLetter)
% cqgl: Computes
% slow variying amplitude (A & B)
% quintic ginzburg landau equation by
% caseType and caseLetter
%
% Synopsis: cqglSolveSLowVaryingAmplitude
%           cqglSolveSLowVaryingAmplitude(caseType, caseLetter):
%           
%
% Input:    caseType = the type of case
%           of interest. Case 1 (bright
%           -bright, front-front, and train solitrary waves),
%           Case 2 (Similar to case 1),
%           Case 3 (progressive bright-bright, train solitary waves,
%           front-front, dark-dark, W-dark-W dark waves)
%           Default: 1
%
%           caseLetter = The type of
%           wave of interest
%           Default: "A"
%
%
% Output:   value = amplitude at a point

if(nargin <= 1)
    caseType = 1;
    caseLetter = "B";
end

% --- Obtain required parameters and
% constraints
param = cqglComputeParams(caseType, caseLetter);
constraints = cqglComputeConstraints(caseType, caseLetter);

% --- Obtain the amplitude for the given
% casetype
amp = cqglSlowVaryingAmplitude(caseType);

% --- actually compute the value of the
% constraints
[alpha_1, alpha_2] = constraints.alpha(param.delta_1, param.delta_2);
alpha = [alpha_1, alpha_2];
nsq = constraints.nsq(param.b, alpha(1,1), param.beta_1, param.chi_1, param.xi_1);
musq = constraints.musq(param.b, alpha(1,1), param.beta_1, param.chi_1, param.xi_1);

L = constraints.L(param.b, param.delta_1, param.xi_1, param.chi_1, alpha(1,1), param.beta_1);
xi_2r = constraints.xi_2r(alpha(1,1), alpha(1,2), param.beta_2, param.beta_1, param.gamma_1, param.gamma_2, param.xi_1 ...
    ,param.delta_1, param.delta_2);
xi_2i = constraints.xi_2i(alpha(1,1), alpha(1,2), param.beta_1, param.beta_2, param.gamma_1, param.gamma_2, param.xi_1);

param.xi_2 = xi_2r + i * xi_2i;
gamma_1i = constraints.gamma_1i(alpha(1,1), alpha(1,2), param.gamma_2, param.beta_1, param.xi_1, param.beta_2, param.xi_2);
param.gamma_1 = param.gamma_1 + i * gamma_1i;
R = constraints.R(param.chi_1, alpha(1,1), param.gamma_1);
w = constraints.w(R, param.v, param.k1,param.gamma_1, param.chi_1);
chi2 = constraints.chi2(alpha(1,1), alpha(1,2), param.gamma_1, param.gamma_2, param.chi_1);
k2 = constraints.k2(R, param.v, param.k1, param.gamma_1, param.chi_1, chi2, param.gamma_2);
omega2 = constraints.omega_2(param.v, k2, R, param.gamma_2);
omega1 = constraints.omega_1(param.v, param.k1, R, param.gamma_1);

%param.xi_2 = xi_2r + i*xi_2i;
param.omega_1 = omega1;
param.omega_2 = omega2;
param.chi_2 = chi2;
param.w = w;
param.R = R;
param.musq = musq;
param.nsq = nsq;
param.k2 = k2;
param.alpha_1 = alpha(1,1);
param.alpha_2 = alpha(1,2);
param.L = L;

% --- Initialize the grid
nx = 100; % number of positional entried
nt = 100; % number of time entries
L = 20; % maximal positional entry
tmax = 40; % maximal time entry
grid = cqglGridInitialization(nx, nt, L, tmax);

% -- obtain coefficient matrix
coeff = cqglCoefficientMatrix(param, grid);

% --- compute LU Factorization of
% coefficient matrix
% [e,f] = cqglTridiagLU(coeff);

AA = grid.U;

for c = 1:grid.nx
    for j = 1:grid.nt
       x = grid.x(1, c);
       t = grid.t(1, j); 
       if grid.t(1, j) == 0 ||  grid.x(1, c) == 0
           grid.U(c, j) = amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k1, param.omega_1 ...
               ,param.alpha_1, param.b, param.L) +  amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k1, param.omega_1 ...
               ,param.alpha_1, param.b, param.L) * (0.1* rand(1));
       end
    end
end



% % loop over time steps and compute
% % ginzburg landau
% for m = 2:grid.nt
%      x = grid.x(1, m);
%      t = grid.t(1, m);
%      ampA = amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k1, param.omega_1 ...
%      ,param.alpha_1, param.b, param.L);
%      ampB = amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k2, param.omega_2 ...
%      ,param.alpha_2, param.b, param.L);
%      linear = grid.dt * param.chi_1 * ampA +  ...
%           - grid.dt * param.gamma_1 * abs(ampA)^2* ampA - grid.dt * param.delta_1 * abs(ampA)^4 * ampA ...
%           - param.delta_1 * param.xi_1 * abs(ampB)^2 * ampA;
%      d = grid.U(:, m-1) + [0; (grid.dt * ( param.gamma_1/(2 * grid.dx^2))).* grid.U(1:end-2, m-1); 0] ...
%          - [0; (grid.dt * ((2 * param.gamma_1)/(2 * grid.dx^2))).* grid.U(2:end-1, m-1); 0] ...
%          + [0; (grid.dt * (param.gamma_1/(grid.dx^2))).* grid.U(3:end, m-1); 0];
%      d = d + linear;
%      d(1) = grid.u0; d(end) = grid.uL;  % overwrite boundary condition values
%      grid.U(:,m) = cqglTriDiagSolve(d,coeff.a,e,f,grid.U(:, m-1)); % solve the system
% end

amplitude = grid.U;

Uo(1) = grid.U(1);
Uo(2:grid.nx-1) = 0;
Uo(grid.nx) = grid.U(grid.nx);
Un(1) = grid.U(1); Un(grid.nx) = grid.U(grid.nx);
for m = 2:grid.nt
   for ii = 1:grid.nx-2
     x = grid.x(1, ii);
     t = grid.t(1, m);
     ampA = amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k1, param.omega_1 ...
     ,param.alpha_1, param.b, param.L);
     ampB = amp.A(x, t, sqrt(param.nsq), param.R, param.w, param.k2, param.omega_2 ...
     ,param.alpha_2, param.b, param.L);
     d(ii) = (ampA/grid.dt) + param.chi_1 * ampA - ...
          param.beta_1 * abs(ampA)^2* ampA - param.delta_1 * abs(ampA)^4 * ampA ...
          - param.xi_1 * abs(ampB)^2 * ampA;
      amplitude(ii, nt) = ampA;
   end
   
   UU = coeff.AA\d';
   Un = [Un(1), UU', Un(nx)];
   grid.U(m, :) = Un;
   Uo = Un;
end

grid.Amp = amplitude;
disp("yoan");
