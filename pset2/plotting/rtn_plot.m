function rtn_plot(t,x,y,z,vel,linewidth)
    % Plots the RTN time history
    % Inputs:
    %   t       - time [orbit periods]
    %   x       - relative position in R normalized by sma_0 []
    %   y       - relative position in T normalized by sma_0 []
    %   z       - relative position in N normalized by sma_0 []
    %   vel     - whether this is a velocity
%     figure();
    set(0,'defaultTextInterpreter','latex');
    tiledlayout(6,4);
    linestyles = ["-", "--", "-."];
    linecolors = ["#0072BD", "red", "green"];
    fontSize = 20;

    % Position vs time 
    nexttile(1,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t, x(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('R (m/s)','FontSize',fontSize), else ylabel('R (m)','FontSize',fontSize), end
        grid on;
    end

    nexttile(9,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t, y(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(17,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t, z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end        
        grid on;
    end
    xlabel('Orbital Periods','FontSize',fontSize);

    % Position vs Position
    nexttile(3,[3,1]); hold on;
    for i=1:size(x,2)
        plot(x(:,i), y(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));% axis equal;
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(4,[3,1]); hold on;
    for i=1:size(x,2)
        plot(x(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));% axis equal;
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(15,[3,1]); hold on;
    for i=1:size(x,2)
        plot(y(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));% axis equal;
        if vel==1 xlabel('T (m/s)','FontSize',fontSize), else xlabel('T (m)','FontSize',fontSize), end
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end       
        grid on;
    end

    % Position vs Position vs Position
    nexttile(16,[3,1]); hold on; view(3);
    for i=1:size(x,2)
        plot3(x(:,i), y(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i)); %axis equal;
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end 
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end
        if vel==1 zlabel('N (m/s)','FontSize',fontSize), else zlabel('N (m)','FontSize',fontSize), end
        grid on;
    end
end

