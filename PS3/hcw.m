%% Hill-Clohessy-Wilshire equations
% Inputs: state = [x,y,z,xdot,ydot,zdot]
% Outputs: stateDot = [xdot,ydot,zdot,xddot,yddot,zddot]

function stateDot = hcw(t,state,n)
    x = state(1);
    z = state(3);
    xdot = state(4);
    ydot = state(5);
    zdot = state(6);

    xddot = 2*n*ydot + 3*n^2*x;
    yddot = -2*n*xdot;
    zddot = -n^2*z;

    stateDot = [xdot;ydot;zdot;xddot;yddot;zddot];

end