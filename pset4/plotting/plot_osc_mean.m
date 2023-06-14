function plot_osc_mean(t, mean_elems, osc_elems, ylabels)
    % Plots osculating and mean elements
    % Works for either absolute or relative
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    % Outputs:
    
    for i = 1:1:6
        subplot(2,3,i); hold on;
        plot(t, osc_elems(:,i));
        plot(t, mean_elems(:,i),'--');
        if i > 3
            xlabel('Orbital Periods')
        end

        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;
    end
    legend('Osculating','Mean')
end

