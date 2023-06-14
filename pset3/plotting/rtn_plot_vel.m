function rtn_plot_vel(xdot,ydot,zdot)
    subplot(2,2,1)
    plot(ydot,xdot); grid on;
    xlabel('T (m/s)'); ylabel('R (m/s)')
    subplot(2,2,2)
    plot(zdot,xdot); grid on;
    xlabel('N (m/s)'); ylabel('R (m/s)')
    subplot(2,2,3)
    plot(ydot,zdot); grid on;
    xlabel('T (m/s)'); ylabel('N (m/s)');
    subplot(2,2,4)
    plot3(xdot,ydot,zdot); grid on; xlabel('R (m/s)'); ylabel('T (m/s)'); zlabel('N (m/s)');

end

