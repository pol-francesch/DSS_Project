function earth_plot(r_body, chief, deputy)
    % Inputs: r_body - radius of body
    %         chief - position and velocity [x,y,z,xvel,yvel,zvel] of chief in body-fixed frame
    %         deputy - position and velocity [x,y,z,xvel,yvel,zvel] of deputy in body-fixed frame
    % Sanity check plot to verify orbit propagation
    [xE,yE,zE] = ellipsoid(0,0,0,r_body,r_body,r_body,20);

    figure(); hold on; axis equal; grid on;
    set(gca,'DefaultLineLineWidth',1)
    plot3(chief(:,1),chief(:,2), chief(:,3), 'm');
    plot3(deputy(:,1),deputy(:,2), deputy(:,3), 'g');
    surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
    title('Orbit around the Earth in ECI');
    xlabel('Position along I (m)');
    ylabel('Position along J (m)');
    zlabel('Position along K (m)')
    legend('Orbit 1', 'Orbit 2','Earth');
    view(3);
    hold off;
end