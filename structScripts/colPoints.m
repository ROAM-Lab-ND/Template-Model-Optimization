%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% colPoints %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% Author: David Kelly
% Created: 7 December 2021
% Last Updated: 7 December 2021

% This function is used to build the collocation points

% INPUTS:
%   method - method to be used for collocation point generation
%   polyDeg - degree of collocation point generation

% OUTPUTS:
%   dynS - Structure used for dynamics storage

%% %%%%%%%%%%%%%%%%%%%%%%%% FUNCTION CREATION %%%%%%%%%%%%%%%%%%%%%%%%%% %%

function dynS = colPoints(method, polyDeg)

    % Import CasADi toolbox
    import casadi.*

    % Create collocation points
    dynS.colPoints =...
            [0 collocation_points(polyDeg, method)];

    % Initialize coefficients of the collocation equation
    dynS.colEqCoeff = zeros(polyDeg + 1, polyDeg + 1);

    % Initialize coefficients of the continuity equation
    dynS.contEqCoeff = zeros(polyDeg + 1, 1);

    % Initialize coefficients of the quadrature function
    dynS.quadFuncCoeff = zeros(polyDeg + 1, 1);

    % For desired degree of polynomials
    for i = 1:(polyDeg + 1)
        
      % Intiliaze Lagrange polynomials coefficients
      coeffLagrange = 1;
      
      % For desired degree of polynomials
      for j = 1:(polyDeg + 1)
          
          % Check if i and j are equivalent
          if j ~= i
              
              % Append Lagrange polynomial coefficients
              coeffLagrange = conv(coeffLagrange,...
                  [1, -dynS.colPoints(j)]);
              coeffLagrange = coeffLagrange/...
                  (dynS.colPoints(i) - dynS.colPoints(j));
              
          end
        
      end
      
      % Calculate continuity equation coefficient
      dynS.contEqCoeff(i) = polyval(coeffLagrange, 1.0);

      % Calculate derivative of the Lagrange polynomial coefficients
      polyDer = polyder(coeffLagrange);
      
      % For desired degree of polynomials
      for j=1:polyDeg+1
          
          % Calculate collocation equation coefficients
          dynS.colEqCoeff(i,j) = polyval(polyDer,...
              dynS.colPoints(j));
          
      end

      % Calculate integral of the Lagrange polynomial coefficients
      polyInt = polyint(coeffLagrange);
      
      % Calculate quadrature function coefficients
      dynS.quadFuncCoeff(i) = polyval(polyInt, 1.0);
      
    end
    
end