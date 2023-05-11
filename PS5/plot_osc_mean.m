function plot_osc_mean(t, mean_elems, osc_elems, ylabels)
    % Plots osculating and mean elements
    % Works for either absolute or relative
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    % Outputs:
    
    figure();
    for i = 1:1:6
        subplot(2,3,i); hold on;
        plot(t, osc_elems(:,i));
        plot(t, mean_elems(:,i),'--');
        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;
    end
    legend('Osculating','Mean')
end

% Plot overlaid unperturbed and J2 for dex-dey, dix-diy, dlambda-dsma
function plot_roe_state_2(unptb_roe,j2_roe)
    dims = size(unptb_roe);
    sets = dims(2) / 6;
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    for i=1:sets
        u_i = unptb_roe(:,1+6*(i-1):6+6*(i-1));
        j2_i = j2_roe(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,3),u_i(:,4)) %,'.','MarkerSize',10); hold on
        plot(j2_i(:,3),j2_i(:,4),'--') %,'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    for i=1:sets
        u_i = unptb_roe(:,1+6*(i-1):6+6*(i-1));
        j2_i = j2_roe(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,5),u_i(:,6)) %,'.','MarkerSize',10); hold on
        plot(j2_i(:,5),j2_i(:,6)) %'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    for i=1:sets
        u_i = unptb_roe(:,1+6*(i-1):6+6*(i-1));
        j2_i = j2_roe(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,2),u_i(:,1)); hold on
        plot(j2_i(:,2),j2_i(:,1),'--');
        hold on
    end
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end