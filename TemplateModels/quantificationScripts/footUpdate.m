%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% footUpdate %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 10 December 2021
% Last Updated: 10 December 2021

% This function is used to update the optimized foot vector for graphing

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function vS = footUpdate(vS)
    
    % Create dummy counter for animated plotting (needed in order to choose the
    % correct foot value)
    count = [1,1];
    
    for i = 1:vS.params.steps
        
        count = [count,(count(end-1:end)-1)];
        
    end

    % Remove duplicates from foot variable for animating plot (needed in order
    % to choose double support and single support foot locations correctly)
    fFix = [];
    
    for i=1:length(vS.optims.fOpt)-1
        
        if abs((vS.optims.fOpt(i)-vS.optims.fOpt(i+1)))<1e-5
            
            fFix=[fFix,i];
            
        end
        
    end
    
    vS.optims.fOpt(fFix(:))=[];
    
end