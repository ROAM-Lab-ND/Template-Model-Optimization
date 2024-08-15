%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% colPoints %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 December 2021
% Last Updated: 7 December 2021

% This function is used to build the collocation points

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function dynS = colPoints(method, polyDeg)

    import casadi.*

    % Get collocation points
    dynS.colPoints =...
            [0 collocation_points(polyDeg, method)];

    % Coefficients of the collocation equation
    dynS.colEqCoeff = zeros(polyDeg + 1, polyDeg + 1);

    % Coefficients of the continuity equation
    dynS.contEqCoeff = zeros(polyDeg + 1, 1);

    % Coefficients of the quadrature function
    dynS.quadFuncCoeff = zeros(polyDeg + 1, 1);

    % Construct polynomial basis
    for i = 1:(polyDeg + 1)
        
      % Construct Lagrange polynomials to get the polynomial basis at the collocation point
      coeffLagrange = 1;
      
      for j = 1:(polyDeg + 1)
          
        if j ~= i
            
          coeffLagrange = conv(coeffLagrange,...
              [1, -dynS.colPoints(j)]);
          coeffLagrange = coeffLagrange/...
              (dynS.colPoints(i) - dynS.colPoints(j));
          
        end
        
      end
      
      % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
      dynS.contEqCoeff(i) = polyval(coeffLagrange, 1.0);

      % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
      polyDer = polyder(coeffLagrange);
      
      for j=1:polyDeg+1
          
        dynS.colEqCoeff(i,j) = polyval(polyDer,...
            dynS.colPoints(j));
        
      end

      % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
      polyInt = polyint(coeffLagrange);
      dynS.quadFuncCoeff(i) = polyval(polyInt, 1.0);
      
    end
    
end