% osc2mean
%   Conversion from osculating to mean orbital elements
%   mean_elem = osc2mean_iterative(osc_elem,J2_flag) computes 
%   J2-perturbed mean orbital elements using the iterative approach given
%   in Alfriend's text, with a Newton-Raphson solver using the Jacobian of
%   the nonlinear transformation.  
% 
%   INPUTS:
%   osc_elem = [a;e;i;O;w;M]: osculating Keplerian orbital element vector
% 
%   J2_flag: Flags whether J2 should be considered. 
%       J2_flag = 1 -> J2 is enabled, calculate according to algorithm
%       J2_flag = 0 -> J2 is disabled, osc elements = mean elements
%       DEFAULT J2_flag = 1
% 
%   OUTPUTS:
%   mean_elem: Mean Keplerian orbital element vector

function mean_elem = osc2mean(osc_elem,J2_flag)

    % Check inputs
    if (nargin < 2) || isempty(J2_flag)
        J2_flag = 1;
    end
    if (nargin < 1) || isempty(osc_elem)
        error('Must input mean elements set');
    end
    
    % Format input to column vector and set tolerance
    osc_elem = osc_elem(:);
    tol = 1e-08;
    
    % With J2, run iterative method
    if J2_flag == 1

        % Convert to osculating equinoctial elements
        osc_equioe = koe2equioe(osc_elem);
        
        % Convert to mean equinoctial elements
        equi_c_mean = osc2mean_NRiterator(osc_equioe, tol);
        
        % Convert to mean keplerian elements
        mean_elem = equioe2koe(equi_c_mean);

        % Reduce numerical issues
        %mean_elem(3:6) = wrapTo2Pi(mean_elem(3:6));
        mean_elem(3:5) = wrapTo2Pi(mean_elem(3:5));
        
         
        % If an angle is close to 2pi, make it 0
        for i=3:6
            if 2*pi - abs(mean_elem(i)) < 1e-12
                mean_elem(i) = 0;
            end
        end

    % Without J2, elements are equal
    else
        mean_elem = osc_elem;

        % Reduce numerical issues
        %mean_elem(3:6) = wrapTo2Pi(mean_elem(3:6));
        mean_elem(3:5) = wrapTo2Pi(mean_elem(3:5));
         
        % If an angle is close to 2pi, make it 0
        for i=3:6
            if 2*pi - abs(mean_elem(i)) < 1e-12
                mean_elem(i) = 0;
            end
        end
    end