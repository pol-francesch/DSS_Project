function plot_roe_state_w_init_and_control_bounds(roe_set, roe_init, ecc_max)
    % Eccentricity controller bounds
    th = 0:pi/50:2*pi;
    x_circ = ecc_max*cos(th) + roe_init(3);
    y_circ = ecc_max*sin(th) + roe_init(4);

    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    plot(roe_set(:,3),roe_set(:,4),'.','MarkerSize',10); hold on
    plot(roe_init(3), roe_init(4), '.', 'MarkerSize', 20);
    plot(x_circ, y_circ, 'r--');
    hold on
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    plot(roe_set(:,5),roe_set(:,6),'.','MarkerSize',5); hold on
    plot(roe_init(5), roe_init(6), '.', 'MarkerSize', 10);
    hold on
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    plot(roe_set(:,2),roe_set(:,1)); hold on
    plot(roe_init(2), roe_init(1), '.', 'MarkerSize', 20);
    hold on
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
    legend('Time history', 'Initial')
end

