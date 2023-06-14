% Plot dex-dey, dix-diy, dlambda-dsma
function plot_roe_state(roe_set)
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    plot(roe_set(:,3),roe_set(:,4),'.','MarkerSize',10); hold on
    hold on
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    plot(roe_set(:,5),roe_set(:,6),'.','MarkerSize',10); hold on
    hold on
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    plot(roe_set(:,2),roe_set(:,1),'.','MarkerSize',10); hold on
    hold on
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end