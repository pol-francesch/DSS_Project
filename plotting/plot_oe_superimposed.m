function plot_oe_superimposed(t, elems1, elems2, ylabels)
    % Plots orbital elements of 2 sets, assuming the data is the same size
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    % Outputs:
    dims1 = size(elems1);
    num_orbits = round(max(t));
    
    for i = 1:1:dims1(2)
        subplot(2,3,i); hold on;
        plot(t, elems1(:,i));
        plot(t, elems2(:,i),'');
        xticks(1:1:num_orbits)
        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;

        if i > 3
            xlabel('Orbit Periods')
        end
    end
end