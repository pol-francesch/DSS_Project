function plot_single_orbit_earth(rE,z)
    [xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);

    hold on; axis equal; grid on;
    set(gca,'DefaultLineLineWidth',1)
    plot3(z(:,1),z(:,2), z(:,3), 'm');
    surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
    xlabel('Position along I (m)');
    ylabel('Position along J (m)');
    zlabel('Position along K (m)')
    legend('Orbit','Earth');
    view(3);
    hold off;
end

