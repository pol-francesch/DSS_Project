function plot_rtn(rho1)
    subplot(2,2,1); hold on;
    plot(rho1(:,2),rho1(:,1), 'LineWidth',2); 
    xlabel('T (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,2); hold on;
    plot(rho1(:,3),rho1(:,1), 'LineWidth',2); 
    xlabel('N (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,3); hold on;
    plot(rho1(:,2),rho1(:,3), 'LineWidth',2); 
    xlabel('T (m)'); ylabel('N (m)'); axis equal; grid on
    subplot(2,2,4); hold on; view(3);
    plot3(rho1(:,1),rho1(:,2),rho1(:,3), 'LineWidth',2); 
    xlabel('R (m)'); ylabel('T (m)'); zlabel('N (m)'); grid on
end