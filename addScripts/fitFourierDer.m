%% %%%%%%%%%%%%%%%%%%%%%%%%%%% fitFourierDer %%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 12 April 2022
% Last Updated: 12 April 2022

% This function is a callable function for evaluating the derivative of the
% fourier fit to the human COM velocity data

% INPUTS:
%   cVals - Coefficient terms for the Fourier fit or Fourier w/ Poly fit
%   varIn - Independent variable to evaluate the Fourier fit at
%   strVal - String used for denoting between Vertical and Horizontal CoM

% OUTPUTS:
%   varOut - Fourier fit evaluation at the varIn value

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function varOut = fitFourierDer(cVals, varIn, strVal)

    import casadi.* 
    
    % Determine Number of Harmonics to Evaluate
    len = length(cVals);
    
    % Evaluate Fourier Fit for Vert, Fourier Fit w/ Polyfit for Hor
    if strcmp(strVal, 'Vert')
        
        % Initialize varOut with constant term (0 since derivative)
        varOut = 0;
        
        % Append each harmonic to varOut
        for i = 1:(len-2)/2
            
            varOut = varOut - i*cVals(end)*...
                (cVals(2*i)*sin(i*cVals(end)*varIn) -...
                cVals(2*i+1)*cos(i*cVals(end)*varIn));
            
        end
        
    else
        
        % Initialize varOut with constant and linear terms (constant = 0)
        varOut = cVals(1);

        % Append each harmonic to varOut
        for i = 1:(len-4)/2
            
            varOut = varOut - i*cVals(end)*...
                (cVals(2*i+2)*sin(i*cVals(end)*varIn) -...
                cVals(2*i+3)*cos(i*cVals(end)*varIn));
            
        end
        
    end

end