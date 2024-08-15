%% %%%%%%%%%%%%%%%%%%%%%%%%% structInitBSLIP %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 6 July 2021
% Last Updated: 7 July 2021

% This function is used to create constraints for recognition by the NLP.
% The constraints are stored in the NLP (problem) struct for future passing
% to the NLP solver.

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function [problem] = addConstraint(problem,constr,lb,ub)
    
    % Store the constraint and bounds in the problem structure
    problem.constraints = {problem.constraints{:},constr};
    problem.constraintsLB = [problem.constraintsLB;lb];
    problem.constraintsUB = [problem.constraintsUB;ub];

end