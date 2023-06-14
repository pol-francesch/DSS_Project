function plot_rtn_pos(rho,marker,linewidth)
    subplot(2,2,1); hold on;
    plot(rho(:,2),rho(:,1),marker,'LineWidth',linewidth); 
    xlabel('T (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,2); hold on;
    plot(rho(:,3),rho(:,1),marker,'LineWidth',linewidth); 
    xlabel('N (m)'); ylabel('R (m)'); axis equal; grid on
    subplot(2,2,3); hold on;
    plot(rho(:,2),rho(:,3),marker,'LineWidth',linewidth); 
    xlabel('T (m)'); ylabel('N (m)'); axis equal; grid on
    subplot(2,2,4); hold on; view(3);
    plot3(rho(:,1),rho(:,2),rho(:,3),marker,'LineWidth',linewidth); 
    xlabel('R (m)'); ylabel('T (m)'); zlabel('N (m)'); grid on
end