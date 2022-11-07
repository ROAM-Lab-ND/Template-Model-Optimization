%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% footUpdate %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 10 December 2021
% Last Updated: 10 December 2021

% This function is used to update the optimized foot vector for graphing

% INPUTS:
%   varS - Structure used for variable and post-optimization storage

% OUTPUTS:
%   varS - Structure used for variable and post-optimization storage

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function varS = footUpdate(varS)
    
    % Create dummy counter for animated plotting
    count = [1,1];
    
    % For each step taken by template model
    for i = 1:varS.params.steps
        
        % Append to dummy counter
        count = [count,(count(end-1:end)-1)];
        
    end

    % Initialize foot position update vector
    fFix = [];
    
    % For each foot position stored
    for i = 1:(length(varS.optims.fOpt) - 1)
        
        % Check if current and next foot position are duplicates
        if abs((varS.optims.fOpt(i) - varS.optims.fOpt(i+1))) < 1e-5
            
            % Append foot position index
            fFix = [fFix, i];
            
        end
        
    end
    
    % Clear all foot positions marked as duplicates
    varS.optims.fOpt(fFix(:))=[];
    
end