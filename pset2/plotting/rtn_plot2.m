function rtn_plot2(t,x,y,z,vel, linewidth)
    % Plots the RTN time history with two columns in t
    % Inputs:
    %   t       - time [orbit periods]
    %   x       - relative position in R normalized by sma_0 []
    %   y       - relative position in T normalized by sma_0 []
    %   z       - relative position in N normalized by sma_0 []
    %   vel     - whether this is a velocity
    set(0,'defaultTextInterpreter','latex');
    tiledlayout(6,4);
    linestyles = ["-", "-"];
    linecolors = ["#0072BD", "red"];
    fontSize = 20;

    % Position vs time 
    nexttile(1,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t(:,i), x(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('R (m/s)','FontSize',fontSize), else ylabel('R (m)','FontSize',fontSize), end
        grid on;
    end

    nexttile(9,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t(:,i), y(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(17,[2 2]); hold on;
    for i=1:size(x,2)
        plot(t(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end        
        grid on;
    end
    xlabel('Orbital Periods','FontSize',fontSize);

    % Position vs Position
    nexttile(3,[3,1]); hold on;
    for i=1:size(x,2)
        plot(x(:,i), y(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(4,[3,1]); hold on;
    for i=1:size(x,2)
        plot(x(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end        
        grid on;
    end

    nexttile(15,[3,1]); hold on;
    for i=1:size(x,2)
        plot(y(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 xlabel('T (m/s)','FontSize',fontSize), else xlabel('T (m)','FontSize',fontSize), end
        if vel==1 ylabel('N (m/s)','FontSize',fontSize), else ylabel('N (m)','FontSize',fontSize), end       
        grid on;
    end

    % Position vs Position vs Position
    nexttile(16,[3,1]); hold on; view(3);
    for i=1:size(x,2)
        plot3(x(:,i), y(:,i), z(:,i),"LineStyle",linestyles(i), "LineWidth", linewidth, "Color", linecolors(i));
        if vel==1 xlabel('R (m/s)','FontSize',fontSize), else xlabel('R (m)','FontSize',fontSize), end
        if vel==1 ylabel('T (m/s)','FontSize',fontSize), else ylabel('T (m)','FontSize',fontSize), end
        if vel==1 zlabel('N (m/s)','FontSize',fontSize), else zlabel('N (m)','FontSize',fontSize), end
        grid on;
    end
end

