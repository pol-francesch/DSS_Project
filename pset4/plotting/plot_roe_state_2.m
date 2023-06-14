% Plot overlaid unperturbed and J2 for dex-dey, dix-diy, dlambda-dsma
function plot_roe_state_2(unptb_roe,j2_roe)
    dims = size(unptb_roe);
    sets = dims(2) / 6;
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,3),unptb_roe(:,4)) %,'.','MarkerSize',10); hold on
        plot(j2_roe(:,3),j2_roe(:,4),'-') %,'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,5),unptb_roe(:,6),'.','MarkerSize',10); hold on
        plot(j2_roe(:,5),j2_roe(:,6),'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,2),unptb_roe(:,1)); hold on
        plot(j2_roe(:,2),j2_roe(:,1),'-');
        hold on
    end
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end