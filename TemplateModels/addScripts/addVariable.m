%% %%%%%%%%%%%%%%%%%%%%%%%%% structInitBSLIP %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to create CasADi variables and tracking them in the
% appropriate indices. This function stores the generated variables and
% indices in the NLP (problem) struct.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [problem,var,varInds] = addVariable(problem, name, len, first, lb, ub)

    import casadi.*

    % Generate symbolic
    var = MX.sym(name,len);

    % Store variable in problem structure with initial value and bounds
    problem.vars = {problem.vars{:},var};
    problem.varsInit = [problem.varsInit;first];
    problem.varsLB = [problem.varsLB;lb];
    problem.varsUB = [problem.varsUB;ub];

    % If index is included store the index
    indEnd = length(problem.varsLB);
    indBegin = indEnd-len+1;
    varInds = indBegin:indEnd;

end