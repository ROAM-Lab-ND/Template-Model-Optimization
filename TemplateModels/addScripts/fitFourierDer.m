%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% fitFourierDer %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 12 April 2022
% Last Updated: 12 April 2022

% This function is a callable function for evaluating the derivative of the
% fourier fit to the human COM velocity data

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function varOut = fitFourierDer(cVals, varIn, strVal)

    import casadi.* 
    
    % Determine Number of Harmonics to Evaluate
    len = length(cVals);
    
    % Evaluate Fourier Fit for Vert, Fourier Fit w/ Polyfit for Hor
    if strcmp(strVal, 'Vert')
        
        varOut = 0;
        
        for i = 1:(len-2)/2
            
            varOut = varOut - i*cVals(end)*...
                (cVals(2*i)*sin(i*cVals(end)*varIn) -...
                cVals(2*i+1)*cos(i*cVals(end)*varIn));
            
        end
        
    else
        
        varOut = cVals(1);
        
        for i = 1:(len-4)/2
            
            varOut = varOut - i*cVals(end)*...
                (cVals(2*i+2)*sin(i*cVals(end)*varIn) -...
                cVals(2*i+3)*cos(i*cVals(end)*varIn));
            
        end
        
    end

end