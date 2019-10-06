function param = cqglComputeParams(caseType, caseLetter)
% cqglComputeParams: Computes
% necessary parameters for the 2D cubic
% quintic ginzburg landau equation by
% caseType and caseLetter
%
% Synopsis: cqglComputeParams
%           cqglComputeParams(caseType, caseLetter):
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
% Output:   result = structure comprise
%           parameters

if(nargin <= 1)
    caseType = 1;
    caseLetter = "A";
end

% --- Define Parameters
if(caseType == 1)
   if(caseLetter == "B" || caseLetter == "b")
       param.gamma_1 = 0;
       param.beta_1 = 1.2 - (0.7 * i);
       param.delta_1 = 2.4+(0.6 * i);
       param.xi_1 = 1.6 + (0.5 * i);
       param.chi_1 = 0.1;
       param.gamma_2 = 0 - (1.1*i);
       param.delta_2 = 2.75 - i;
       param.beta_2 = 0.4 + i;
       param.b = 3;
       param.v = 2;
       param.k1 = 5;
   end
end

end

