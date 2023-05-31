function plot_oe_error(t, elems, err, ylabels)
    % Plots orbital elements with error bounds
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    %   err         - error (standard deviation/covariance)
    %   ylabels     - labels for variables
    % Outputs:

    sigma_color = "#3da3ff";

    dims = size(elems);
    num_orbits = round(max(t));
    
    for i = 1:1:dims(2)
        subplot(2,3,i); hold on;
        % Plot the error bounds
        plus_sigma = elems(:,i) + err(:,i);
        minus_sigma = elems(:,i) - err(:,i);
%         plot(t, plus_sigma,'b')
%         plot(t, minus_sigma,'b')
        % Fill between
        patch([t,fliplr(t)], [minus_sigma',fliplr(plus_sigma')],'b','FaceColor',sigma_color,'EdgeColor','none','FaceAlpha',0.7)
        plot(t, elems(:,i),'Color','#194063');
        xticks([0 .1 .2 .3])
%         xticks(1:1:num_orbits)
        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;

        if i > 3
            xlabel('Orbit Periods')
        end
    end
end