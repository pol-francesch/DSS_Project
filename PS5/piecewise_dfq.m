function [state_eci_c,state_eci_d] = piecewise_dfq(rv_init_c, rv_init_d,man_dVs, t_man, options, tspan, stepSize, mu,rE,J2)
    % Initialize state arrays
    state_eci_c = [];
    state_eci_d = []; % zeros(steps,6);
    
    % Manually propagate to each maneuver time
    t0 = 0;
    r0_c = rv_init_c;
    r0_d = rv_init_d;
    true_tspan = [];
    
    for j=1:length(t_man)
        % Propagate to maneuver
        tspan_to_man = t0:stepSize:t_man(j);
        [~, chief]  = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_man,r0_c,options);
        [~, deputy] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_man,r0_d,options);

        oe_c = rv2singular_oe(mu, chief(end,1:3), chief(end,4:6));
        qns_oe_c = singular2qns(oe_c);

        sma0 = oe_c(1); ecc0 = oe_c(2); u0 = qns_oe_c(6);
        inc0 = oe_c(3); raan0 = oe_c(4);
        
        % Convert maneuver to ECI
        dV_RTN = man_dVs(:,j);
        thetadot0_0 = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(u0))^2; % RTN coords / ECI frame
        w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame
        T_RTN2ECI = rtn2eci_rot(inc0,raan0,u0);

        dV_ECI = T_RTN2ECI * (dV_RTN + cross(w_RTNinECI, T_RTN2ECI'*[0;0;0])); %deputy(end,1:3)'));

        % Perform maneuver
        post_man = deputy(end,:) + [0,0,0,dV_ECI'];
    
        % Add to common state
        if tspan_to_man(end) == t_man(j)
            chief = chief(1:end-1,:);
            deputy = deputy(1:end-1,:);
            tspan_to_man = tspan_to_man(1:end-1);
        end
        % state_eci_d = [state_eci_d; deputy; post_man];
        % true_tspan  = [true_tspan, tspan_to_man, t_man(j)];
        state_eci_c = [state_eci_c; chief];
        state_eci_d = [state_eci_d; deputy];
        true_tspan  = [true_tspan, tspan_to_man];
        
        r0_c = chief(end,:);
        r0_d = post_man;
        t0   = t_man(j);
    end
    
    % Finish the propagation 
    if t_man(end) < tspan(end)
        tspan_to_end = t0:stepSize:tspan(end);
        [~, chief]  = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_end,r0_c,options);
        [~, deputy] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_end,r0_d,options);
        
        state_eci_c = [state_eci_c; chief];
        state_eci_d = [state_eci_d; deputy];
        true_tspan  = [true_tspan, tspan_to_end];
    end
    
    tspan = true_tspan;
    % state_eci_d = interp1(true_tspan, state_eci_d, tspan);
end

