% mean2osc
%   Conversion from mean to osculating orbital elements
%   osc_elem = mean2osc(mean_elem,J2_flag) computes J2-perturbed
%   osculating orbital elements from a set of mean orbital elements. 
% 
%   INPUTS:
%   mean_elem = [a;e;i;O;w;M]: mean Keplerian orbital element vector
% 
%   J2_flag: Flags whether J2 should be considered. 
%       J2_flag = 1 -> J2 is enabled, calculate according to algorithm
%       J2_flag = 0 -> J2 is disabled, osc elements = mean elements
%       DEFAULT J2_flag = 1
% 
%   OUTPUTS:
%   osc_elem: osculating Keplerian orbital element vector

function osc_elem = mean2osc(mean_elem,J2_flag)

    % Check inputs
    if (nargin < 2) || isempty(J2_flag)
        J2_flag = 1;
    end
    if (nargin < 1) || isempty(mean_elem)
        error('Must input mean elements set');
    end
    
    % Format input to column vector and set tolerance
    mean_elem = mean_elem(:);

    % With J2, run method
    if J2_flag == 1
        
        % Convert to mean equinoctial elements
        mean_equi_elem = koe2equioe(mean_elem);

        % Convert to osculating equinoctial elements
        [~, osc_equi_elem] = mean_osc_closed_equi(mean_equi_elem,J2_flag);

        % Convert to keplerian elements
        osc_elem = equioe2koe(osc_equi_elem)';
        
        % Format output
        osc_elem = osc_elem(:);
        
    % Without J2, elements are equal
    else
        osc_elem = mean_elem;
    end
end