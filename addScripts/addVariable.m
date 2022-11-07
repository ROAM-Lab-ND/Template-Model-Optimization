%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% addVariable %%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to create CasADi variables and track them in the
% appropriate indices. This function stores the generated variables and
% indices in the NLP (problem) struct.

% INPUTS:
%   problem - Structure used to pass necessary information to NLP solver
%   name - Variable name to add to problem struct
%   len - Size of variable to add to problem struct
%   first - Initial guess of variable value to seed NLP solver
%   lb - Lower bound value permissible for variable
%   ub - Upper bound value permissible for variable

% OUTPUTS:
%   problem - Structure used to pass necessary information to NLP solver
%   var - CasADi variable generated based off input arguments
%   varInds - Updated indices for tracking instances of variable creation

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [problem, var, varInds] =...
    addVariable(problem, name, len, first, lb, ub)
    
    % Have to import CasADi to generate CasADi variables
    import casadi.*

    % Generate symbolic based on inputted name and length
    var = MX.sym(name, len);

    % Store variable in problem structure with initial value and bounds
    problem.vars = {problem.vars{:}, var};
    problem.varsInit = [problem.varsInit; first];
    problem.varsLB = [problem.varsLB; lb];
    problem.varsUB = [problem.varsUB; ub];

    % Return updated indices
    indEnd = length(problem.varsLB);
    indBegin = indEnd - len + 1;
    varInds = indBegin:indEnd;

end