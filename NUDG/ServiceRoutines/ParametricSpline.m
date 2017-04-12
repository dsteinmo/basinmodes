function [ x_t, y_t, t ] = ParametricSpline(x , y)
%ParametricSpline Returns arclength-parameterized spline
%parameter is 't', x_t & y_t may be evaluated with 
% ppval(x_t,tt),ppval(y_t,tt) where tt is in [min(t),max(t)]

    n = length(x);
    t = zeros(n, 1);

    for i=2:n
        arc_length = sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2);
        t(i) = t(i-1) + arc_length;
    end

    x_t = spline(t, x);
    y_t = spline(t, y);

end
