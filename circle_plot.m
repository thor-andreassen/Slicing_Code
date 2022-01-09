function h = circle_plot(x,y,r)
% function to plot a basic circle, and return the handle to the plot
% containing the circle. Inputs are the x, y location of the center of the
% circle, and the radius of the circle.
        hold on
        th = 0:pi/50:2*pi;
        xunit = r * cos(th) + x;
        yunit = r * sin(th) + y;
        h = plot(xunit, yunit);
        hold off
end