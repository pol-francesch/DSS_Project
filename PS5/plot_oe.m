function plot_oe(t, elems, ylabels)
    % Plots orbital elements
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    % Outputs:
    dims = size(elems);
    
    for i = 1:1:dims(2)
        subplot(2,3,i); hold on;
        plot(t, elems(:,i));
        xlabel('Orbit Periods')
        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;
    end
end