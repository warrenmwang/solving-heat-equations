classdef Solver
    methods(Static)
        function approx = gaussianQuadrature(a, b, f, c, nodes, n)
            % numerically approximate integral of f from a to b
            % works for any range [a,b] due to change of vars
            % Inputs:
            %    a (int) - lower bound
            %    b (int) - upper bound
            %    f (function handle) - function to integrate
            %    c (list of floats) - special gaussian coefficients
            %    nodes (list of floats) - special gaussian node points
            %    n (int) - order of the gaussian quadrature
            % Outputs:
            %    approx (float) - integral approximation value
            g = @(a,b,x) ((b-a)*x+(a+b))/2;
            approx = 0;
            for i = 1:n+1
                approx = approx + c(i)*f(g(a,b,nodes(i)));
            end
            approx = approx * (b-a)/2;
        end
        
        % function approx = compositeGaussianQuadrature(a,b,f,c,nodes,n,numSubIntervals)
        %     % approximates the integral for function f using gaussian quadratures
        %     % Inputs:
        %     %    a (int) - lower bound
        %     %    b (int) - upper bound
        %     %    f (function handle) - function to integrate
        %     %    c (list of floats) - special gaussian coefficients
        %     %    nodes (list of floats) - special gaussian node points
        %     %    n (int) - order of the gaussian quadrature
        %     %    numSubIntervals (int) - number of subintervals to divide (b-a) into
        %     % Outputs:
        %     %    approx (float) - integral approximation value
        %     approx = 0;
        %     integrationRange = (b-a)/numSubIntervals;
        
        %     % note start with a = a
        %     b = integrationRange;
        %     for i = 1:numSubIntervals
        %         approx = approx + Solver.gaussianQuadrature(a, b, f, c, nodes, n);
        %         a = b;
        %         b = b + integrationRange;
        %     end
        % end

        function y = L_0(i, x, h)
            % L^{i}_{0}
            % lagrange interpolation for the 0th (from top left to bottom right) piece of the tent/phi functions
            % Inputs:
            %    i (int) - tent function idx
            %    x (float) - input point to interpolate at
            %    h (float) - step size on interval
            % Outputs:
            %    y (float) - interpolated value
            y = (x-i*h)/(-h);
        end

        function y = L_1(i, x, h)
            % L^{i}_{1}
            % lagrange interpolation for the 1th (from bottom left to top right) piece of the tent/phi functions
            % Inputs:
            %    i (int) - tent function idx
            %    x (float) - input point to interpolate at
            %    h (float) - step size on interval
            % Outputs:
            %    y (float) - interpolated value
            y = (x-(i-1)*h)/h;
        end

        function y = L_0_prime(i, h)
            % d/dx [L^{i}_{0}]
            % derivative of lagrange interpolation for the 0th (from top left to bottom right) piece of the tent/phi functions
            % Inputs:
            %    i (int) - tent function idx
            %    x (float) - input point to interpolate at
            %    h (float) - step size on interval
            % Outputs:
            %    y (float) - interpolated value
            y = (1/h) - i;
        end

        function y = L_1_prime(i, h)
            % d/dx [L^{i}_{1}]
            % derivative of lagrange interpolation for the 1th (from bottom left to top right) piece of the tent/phi functions
            % Inputs:
            %    i (int) - tent function idx
            %    x (float) - input point to interpolate at
            %    h (float) - step size on interval
            % Outputs:
            %    y (float) - interpolated value
            y = (1/h) - (i-1);
        end

        % function y = tentFunction(i, a, h, x)
        %     % special tent / phi function
        %     % Inputs:
        %     %    i (int) - idx of the tent function
        %     %    a (int) - lower bound of interval
        %     %    h (float) - step size
        %     %    x (float) - input to evaluate at
        %     % Outputs:
        %     %    y (float) - tent func output
        %     if x >= a + (i*h) && x <= a + ((i+1)*h)
        %         % in bounds, use left lagrange interpolation
        %         y = Solver.L_1(i+1,x,h);
        %         return;
        %     elseif x > a + ((i+1)*h) && x <= a + ((i+2)*h)
        %         % in bounds, use right lagrange interpolation
        %         y = Solver.L_0(i+2,x,h);
        %         return;
        %     else
        %         % out of bounds
        %         y = 0.0;
        %         return;
        %     end
        % end

        function old_main()
            % setup the matrix for solving u coeffs
            a11 = gaussianQuadrature(0.00, 0.25, @(x) L_1_prime(1,h)^2, c, nodes, n);
            a12 = gaussianQuadrature(0.25, 0.50, @(x) L_1_prime(2,h) * L_0_prime(2,h), c, nodes, n);
            a13 = 0;
            
            a21 = gaussianQuadrature(0.25, 0.50, @(x) L_0_prime(2,h)*L_1_prime(2,h), c, nodes, n);
            a22 = gaussianQuadrature(0.25, 0.50, @(x) L_1_prime(2,h)^2, c, nodes, n) + gaussianQuadrature(0.5, 0.75, @(x) L_1_prime(3,h)*L_0_prime(3,h), c, nodes, n);
            a23 = gaussianQuadrature(0.50, 0.75, @(x) L_1_prime(3,h)*L_0_prime(3,h), c, nodes, n);

            a31 = 0;
            a32 = gaussianQuadrature(0.50, 0.75, @(x) L_0_prime(3,h)*L_1_prime(3,h), c, nodes, n);
            a33 = gaussianQuadrature(0.50, 0.75, @(x) L_1_prime(3,h)^2, c, nodes, n) + gaussianQuadrature(0.75, 1, @(x) L_0_prime(4,h)^2, c, nodes, n);

            b1 = -gaussianQuadrature(0.00, 0.25, @(x) L_1(1,x,h)*f(x), c, nodes, n) - gaussianQuadrature(0.25, 0.50, @(x) L_0(2,x,h)*f(x), c, nodes, n);
            b2 = -gaussianQuadrature(0.25, 0.50, @(x) L_1(2,x,h)*f(x), c, nodes, n) - gaussianQuadrature(0.50, 0.75, @(x) L_0(3,x,h)*f(x), c, nodes, n);
            b3 = -gaussianQuadrature(0.50, 0.75, @(x) L_1(3,x,h)*f(x), c, nodes, n) - gaussianQuadrature(0.75, 1.00, @(x) L_0(4,x,h)*f(x), c, nodes, n);

            A = [a11 a12 a13; ...
                 a21 a22 a23; ...
                 a31 a32 a33];
            b = [b1 b2 b3]';

            % solve
            sol = A\b;
            disp("Computed Coefficients [u_1, u_2, u_3]: ")
            disp(sol')

            x = 0:h:1;
            y = [0 sol(1) sol(2) sol(3) 0];
            fig = figure;
            plot(x,y)

            actual_sol = @(x) (x.^4 / 12) - (1/2).*x.^3 + (5/12).*x;
            x = linspace(0,1);
            fig = figure;
            plot(x,actual_sol(x),'r')
        end

        function main()
            clear; close all;
            format long;
            gaussianQuadrature = @(a,b,f,c,nodes,n) Solver.gaussianQuadrature(a,b,f,c,nodes,n);
            L_0_prime = @(i,h) Solver.L_0_prime(i,h);
            L_1_prime = @(i,h) Solver.L_1_prime(i,h);
            L_0 = @(i,x,h) Solver.L_0(i,x,h);
            L_1 = @(i,x,h) Solver.L_1(i,x,h);

            % solving
            % \frac{\partial^2 u}{\partial x^2} + (-x^2 + 3x) = 0

            % function f in diff eq
            f = @(x) -x^2 + 3*x;

            numTentFunctions = 3;
            a = 0;
            b = 1;
            h = (b-a)/(numTentFunctions+1);
            fprintf("h=%f\n",h);

            % gaussian quadrature values
            c = [5/9, 8/9, 5/9];
            nodes = [-sqrt(3/5), 0, sqrt(3/5)];
            n = 2;
            
            %TODO: assemble connectivity table

        end
    end
end
