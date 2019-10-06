function constraints = constraintsCQCGLE(case_const, case_subref)
% Provides a set of parameter
% constraints specific to the cubic
% quintic complex ginzburg landau
% equation based on input case param and
% case subreference provided in the
% paper `Dissipative Solitons in a
% Generalize Coupled Cubic-Quintic
% Ginzburg Landau Eqautions`
% Input:
%   case_const: int valued experimental case
%   case_subref: string valued case sub
%   reference
% Output:
%   constraints: struct comprised of
%   provided case specific cqcgle
%   constraints

% sanitize inputs
if (isa(case_const, "double") ~= 1|| isa(case_const, "single") ~= 1) && isa(case_subref, "char") ~= 1   
   error('Error. \Invalid input params')
   return 
end

switch case_const
    case 1
        % call the case 1 subroutine
    case 2
        % call the case 2 subroutine
    case 3 
        % call the case 3 subroutine
    otherwise
        error('Error. \Please provide a case number in the range of 1 to 3')
        return

end


function constraints = case1CQCGLEConstraints(case_subref)
% Provides a set of parameter
% constraints specific to the cubic
% quintic complex ginzburg landau
% equation for case 1 provided in the
% paper `Dissipative Solitons in a
% Generalize Coupled Cubic-Quintic
% Ginzburg Landau Eqautions`
% Input:
%   case_subref: string valued case sub
%   reference
% Output:
%   constraints: struct comprised of
%   provided case specific cqcgle
%   constraints
% sanitize inputs
if isa(case_subref, "char") ~= 1   
   error('Error. \Invalid input params')
   return 
end


% define usefule anonymous functions
P = @(a2, y2) 3 * a2 * imag(y2) - real(y2) + 2 * (a2^2) * real(y2);
N = @(a2, y2) (-3+4*(a2^2))*imag(y2) - 8 * a2 * real(y2);
S1 = @(a1, delta1) -8*a1*imag(delta1) - 3*real(delta1) + 4 * (a1^2) * real(delta);
S2 = @(a1, y1) 8 * a1 * imag(y1) - 3 * real(y1) + 4 * (a1^2) * real(y1);
S3 = @(a1, y2) 4 * a1 * imag(y1) - real(y1) + 4 * (a1^2) * real(y1);
S4 = @(a1) 9 + 40 * (a1^2) + 16 * (a1^4);
S5 = @(a1, y1) (-3 + 4 *(a1^2))*imag(y1) - 8 * a1 * real(y1);
S6 = @(a1,y1) (-1 + 4*(a1^2))*imag(y1) - 4 * a1 * real(y1);

% solve the quadratic relation for a1
% and a2
p_a1 = @(delta1)[-4 * real(delta), 8 * imag(delta), 3*real(delta)];
p_a2 = @(delta2, y2)[4 * (-real(y2)*imag(delta2) + (imag(y2) * real(delta2))), ...
    -8 * (imag(y2)*imag(delta2) + (real(y2) * real(delta2))), ...
    3*real(y2)*imag(delta2) - 3 * imag(y2) * real(delta2)];

% TODO: implement logic deciphering
% wether to take the root value less
% than 0 or greater than zero
a1 = routs(p_a1(delta1));
a2 = roots(p_a2(delta2, y2));

y1_r = 0; y1_i = imag(y1);
y2_i = imag(y2); y2_r = real(y2);

beta_1r = real(beta1); beta_1i = imag(beta1);
beta_2r = real(beta2); beta_2i = imag(beta2);


x2 = ((4 * a2 * y2_i -y2_r + 4 * (a2^2)*y2_r)* x1)/(4 * a1 * y1_i);
e2_i = (-3*a1* beta_1i * beta_2i * y1_i + 2*(a1^2)*beta_2i * beta_1r * y1_i + ...
    beta_1r * (-beta_2i * y1_i + (y2_i - 2 * (a2^2)*y2_i + 3*a2*y2_r)*e1_i) + ...
    beta_1i * ((-1 + 2*(a2^2)) * y2_i - 3*a2*y2_r)*e1_r)/ ...
    (y1_i* (-3*a1*e1_i - e1_r + 2*(a1^2)*e1_r));

beta_2r = -N(a2, y2) * real(delta1) * real(e1)^3 * real(e2) + 8 * N(a2, y2) * a1^6 * real(delta1) * real(delta1)^3 * real(e2) + ...
    a1 * (8 * P(a2, y2) * real(b1)^2*imag(delta2) * (real(b1)*iamag(e1) - imag(b1)*real(e1))) - 9 * N(a2, y2)*real(delta1) * imag(e1)*real(e1)^2 * real(e2) - ...
    4*(a1^5)*(8 * P(a2, y2) * real(b1)^2 * imag(delta2) * (-real(b1)*imag(e1)+imag(b1)*real(e1)) + 9 * N(a2, y2) * real(delta1)*imag(e1)*(real(e1)^2)*real(e2)) + ...
    


end


end

