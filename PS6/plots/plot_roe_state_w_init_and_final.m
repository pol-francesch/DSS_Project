function plot_roe_state_w_init_and_final(roe_set, roe_init, roe_final)
    % path colormap
    n = length(roe_set(:,3));
    left_color = [.85 .88 .93]; % lighter blue
    right_color = [0 0.25 0.7]; % darker blue
    colors = interp1([0, 1], [left_color; right_color], linspace(0, 1, n));
    cd = [uint8(colors*255) uint8(ones(n,1))].';
    blue = "#0072BD";

    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    plot(roe_final(3), roe_final(4), 'x','Color','red', 'MarkerSize', 7,'LineWidth',2);
    plot(roe_init(3), roe_init(4), 'o', 'Color',blue,'MarkerSize', 6,'LineWidth',.7);
    path = plot(roe_set(:,3),roe_set(:,4),'MarkerSize',20);
    drawnow
    set(path.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    
    plot(roe_set(end,3),roe_set(end,4),'.','Color',blue,'MarkerSize',20);
    hold on
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    plot(roe_final(5), roe_final(6), 'x','Color','red', 'MarkerSize', 7,'LineWidth',2);
    plot(roe_init(5), roe_init(6), 'o', 'Color',blue,'MarkerSize', 6,'LineWidth',.7);
    path = plot(roe_set(:,5),roe_set(:,6),'MarkerSize',20);
    drawnow
    set(path.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    
    plot(roe_set(end,5),roe_set(end,6),'.','Color',blue,'MarkerSize',20);
    hold on
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    p_tgt=plot(roe_final(2), roe_final(1), 'x','Color','red', 'MarkerSize', 7,'LineWidth',2);
    p_init=plot(roe_init(2), roe_init(1), 'o', 'Color',blue,'MarkerSize', 6,'LineWidth',.7);
    path = plot(roe_set(:,2),roe_set(:,1),'MarkerSize',20);
    drawnow
    set(path.Edge, 'ColorBinding','interpolated', 'ColorData',cd)
    
    p_final=plot(roe_set(end,2),roe_set(end,1),'.','Color',blue,'MarkerSize',20);
    path = plot([NaN NaN], [NaN NaN],'Color', blue, 'DisplayName', 'Path'); % for legend color mismatch
    hold on
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')

    legend([p_tgt,p_init,path,p_final],'Target','Initial','Time history','End state','Interpreter','latex');
end

