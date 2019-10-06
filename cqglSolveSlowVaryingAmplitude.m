function amp = cqglSlowVaryingAmplitude(caseType)
% cqglComputeParams: Computes
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
end

if caseType ==  1
    amp.A = @(x, t, n, r, w, k1, omega_1, alpha_1, b) ...
    (n * exp(0.5*(r*x+w*t) + 1i*(k1 * x - omega_1 * t)))/...
    ((1+b*exp((r*x+w*t)) + L* exp(r*x+w*t))^(0.5 + i*alpha_1));
    amp.B = @(x, t, mu, r, w, k2, omega_2, alpha_2, b) ...
    (n * exp(0.5*(r*x+w*t) + 1i*(k1 * x - omega_2 * t)))/...
    ((1+b*exp((r*x+w*t)) + L* exp(r*x+w*t))^(0.5 + i*alpha_2));
elseif caseType == 2
    amp.A = @(x, t, n, r, w, k1, omega_1, alpha_1, b) ...
    (n * exp((r*x+w*t) + 1i*(k1 * x - omega_1 * t)))/...
    ((1+b*exp(2*(r*x+w*t)))^(0.5 + i*alpha_1));

    amp.B = @(x, t, mu, r, w, k2, omega_2, alpha_2, b) ...
    (mu * exp((r*x+w*t) + 1i*(k2 * x - omega_2 * t)))/...
    ((1+b*exp(2*(r*x+w*t)))^(0.5 + i*alpha_2));
end

end

