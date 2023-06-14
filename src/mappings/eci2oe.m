function [oe,eVec,hVec,mechEnergy] = eci2oe ( mu, cartesianState )
%{
ECI2OE converts position , r , and velocity , v , vectors in the ECI frame
to
orbital elements .
Inputs :
r_ijk - 3 x1 vector of radius in ECI frame [ m ]
v_ijk - 3 x1 vector of velocity in ECI frame [ m ]
mu - gravitational parameter of central body [ m ^ -3/ s ^2]
Outputs :
oe = struct with the following fields :
a - semi - major axis of orbit [ m ]
e - eccentricity of orbit
i - inclination of orbit [ rad ]
Om - right ascension of the ascending node [ rad ]
w - argument of periapsis [ rad ]
anom - true anomaly [ rad ]
ang - extra placeholder angle for special cases [ deg ]
This function is able to handle any orbit type , including equatorial and
circular orbits . The output ang will output the following values for
each orbit type :
1) Elliptical equatorial
Longitude of periapsis , Pi = Om + w
2) Circular inclined
Argument of latitude , u = w + anom
3) Circular equatorial
True latitude , lambda = Om + w + anom
Otherwise , ang will be NaN
%}
r_ijk = cartesianState(1:3);
v_ijk = cartesianState(4:6);

v = norm( v_ijk ) ;
r = norm( r_ijk ) ;
% Create angular momentum , mean motion , and eccentricity vectors :
hVec = cross( r_ijk , v_ijk ) ;
h = norm( hVec ) ;
nVec = cross([0 , 0 , 1] , hVec ) ;
n = norm( nVec ) ;
eVec = (1/mu)*(( v^2 - mu/r ) * r_ijk - dot( r_ijk , v_ijk ) * v_ijk ) ;
e = norm( eVec ) ;
% Compute orbit shape / size :
mechEnergy = 0.5* v^2 - mu/r ;
if e ~= 1 % Elliptical , Circular , Hyperbolic orbits
a = - mu/(2* mechEnergy ) ;
else % Parabolic orbit
a = inf ;
end
if hVec(3)/ h < 1 && hVec(3)/ h > -1
    % Compute the orientation of the orbit
    i = acos( hVec(3)/ h ) ;
    Om = acos( nVec(1)/ n ) ;
    w = acos( dot( nVec , eVec )/( n * e ) ) ;
    anom = acos( dot( eVec , r_ijk )/( e * r ) ) ;
    % Place angles in the correct range :
    if nVec(2) < 0
        Om = 2*pi - Om ;
    end
    if eVec(3) < 0
        w = 2*pi - w ;
    end
    if dot( r_ijk , v_ijk ) < 0
        anom = 2*pi - anom ;
    end
    % Account for special cases
    % if i == 0 && e ~= 0 % Elliptical equatorial
    % % Provide the longitude of periapsis ( PI = Om + w )
    % ang = acos( eVec(1)/ e ) ;
    % if eVec(2) < 0
    % ang = 360 - ang ;
    % end
    % elseif i ~= 0 && e == 0 % Circular inclined
    % % Provide the argument of latitude ( u = w + anom )
    % ang = acos( dot( nVec , r_ijk )/( n * r ) ) ;
    % if r_ijk(3) < 0
    % ang = 360 - ang ;
    % end
    % elseif i == 0 && e == 0 % Circular equatorial
    % % Provide the true latitude ( lambda = Om + w + anom )
    % ang = acos( r_ijk(1)/ r ) ;
    % if r_ijk(2) < 0
    % ang = 360 - ang ;
    % end
    % else
    % % Default
    % ang = NaN ;
    % end
    % Assign everything to struct :
    oe = zeros(6,1);
    oe(1,1) = a ;
    oe(2,1) = e ;
    oe(3,1) = i ;
    oe(4,1) = Om ;
    oe(5,1) = w ;
    
    % Mean anomaly
    anom = true2mean(anom, e);
    oe(6,1) = anom ;
else
    oe = cartesianState;
    return
end

end