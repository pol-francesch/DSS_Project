function rtn_plot_pos(x,y,z)
    subplot(2,2,1)
    plot(y,x)
    xlabel('T (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,2)
    plot(z,x)
    xlabel('N (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,3)
    plot(y,z)
    xlabel('T (m)'); ylabel('N (m)'); axis equal; grid on
    subplot(2,2,4)
    plot3(x,y,z); xlabel('R (m)'); ylabel('T (m)'); zlabel('N (m)'); grid on
end

