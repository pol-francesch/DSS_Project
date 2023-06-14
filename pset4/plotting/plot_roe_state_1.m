% Plot dex-dey, dix-diy, dlambda-dsma
function plot_roe_state_1(roe_set)
    dims = size(roe_set);
    sets = dims(2) / 6;
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,3),u_i(:,4));
        hold on
    end
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,5),u_i(:,6),'.');
        hold on
    end
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,2),u_i(:,1));
        hold on
    end
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end