function [a, p_options] = lamberts(type,s,c,r1,r2,mu,TOF)
    syms a_sym
    p_options = [0,0];

    % application of lamberts theorem to find a
    if (type == "1A") || (type == "2A")
        %% ELLIPSE
        alpha = 2*asin(sqrt(s/(2*a_sym)));
        beta = 2*asin(sqrt((s-c)/(2*a_sym)));

        if type == "1A"
            eqn = sqrt(mu)*TOF == (a_sym^(3/2))*...
                ((alpha - sin(alpha)) - (beta - sin(beta)));
            a = double(vpasolve(eqn,a_sym));
        else
            % type 2A
            a = double(vpasolve(...
                sqrt(mu)*TOF == (a_sym^(3/2))*...
                ((alpha - sin(alpha)) + (beta - sin(beta))),a_sym));
        end

        alpha = vpa(subs(alpha,a_sym,a));
        beta = vpa(subs(beta,a_sym,a));

        p_options(1) = 4*a*(s-r1)*(s-r2)*(sin((alpha+beta)/2)^2)/c^2;
        p_options(2) = 4*a*(s-r1)*(s-r2)*(sin((alpha-beta)/2)^2)/c^2;
    elseif (type == "1H") || (type == "2H")
        %% HYPERBOLA
        alpha = 2*asinh(sqrt(s/(2*abs(a_sym))));
        beta = 2*asinh(sqrt((s-c)/(2*abs(a_sym))));

        if type == "1H"
            a = double(vpasolve(...
                sqrt(mu)*TOF == (abs(a_sym)^(3/2))*...
                ((sinh(alpha) - alpha) - (sinh(beta) - beta)),a_sym));
        else
            % type 2H
            a = double(vpasolve(...
                sqrt(mu)*TOF == (abs(a_sym)^(3/2))*...
                ((sinh(alpha) - alpha) + (sinh(beta) - beta)),a_sym));
        end

        alpha = vpa(subs(alpha,a_sym,a));
        beta = vpa(subs(beta,a_sym,a));

        p_options(1) = 4*abs(a)*(s-r1)*(s-r2)*(sinh((alpha+beta)/2)^2)/c^2;
        p_options(2) = 4*abs(a)*(s-r1)*(s-r2)*(sinh((alpha-beta)/2)^2)/c^2;
    else
        a = 0;
        disp('ERROR: Orbit type invalid')
    end
end

