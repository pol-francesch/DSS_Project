function plot_rtn_time_series(t,data,ylabels)
% Plots time series data for 3 components in a vertical stack of plots
    for j=1:size(data,2)
        subplot(3,1,j); hold on
        plot(t,data(:,j))
        ylabel(ylabels(j))
        grid on
    end
    xlabel('Orbit Periods')
end