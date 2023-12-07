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
        
        function y = L_0(i, x, h)
            % L^{i}_{0}
            % lagrange interpolation for the 0th (from top left to bottom right) piece of the tent/phi functions
            % Inputs:
            %    i (int) - tent function idx
            %    x (float) - input point to interpolate at
            %    h (float) - step size on interval
            % Outputs:
            %    y (float) - interpolated value
            if x >= (i-1)*h && x <= i*h
                y = (i*h - x)/h;
            else
                y = 0;
            end
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
            if x >= (i-1)*h && x <= i*h
                y = (x - (i-1)*h)/h;
            else
                y = 0;
            end
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
            y = -1/h;
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
            y = 1/h;
        end

        function [y,b] = createRow(i, h, f, numTentFunctions)
            % i (int) - ith idx
            % h (float) - step size
            % f (function handle)
            y = zeros(1,numTentFunctions);

            gaussianQuadrature = @(a,b,f,c,nodes,n) Solver.gaussianQuadrature(a,b,f,c,nodes,n);
            L_0_prime = @(i,h) Solver.L_0_prime(i,h);
            L_1_prime = @(i,h) Solver.L_1_prime(i,h);
            L_0 = @(i,x,h) Solver.L_0(i,x,h);
            L_1 = @(i,x,h) Solver.L_1(i,x,h);

            % gaussian quadrature values
            c = [5/9, 8/9, 5/9];
            nodes = [-sqrt(3/5), 0, sqrt(3/5)];
            n = 2;

            % if numTentFunctions == 3, then we only allow coefficients a1, a2, a3 -- therefore a0, a4 should be zero
            % therefore it's the bounds that we need to check if we need to ignore.
            % if i == 1, then the i_minus_1 value therefore needs to be ignored
            % if i == numTentFunctions, then the i_plus_1 value therefore needs to be ignored
            if i == 1
                a_i_minus_1 = 0;
            else
                a_i_minus_1 = gaussianQuadrature((i-1)*h, i*h, @(x) L_0_prime(i,h) * L_1_prime(i,h), c,nodes,n);
            end

            a_i = gaussianQuadrature((i-1)*h, i*h, @(x) L_1_prime(i,h) * L_1_prime(i,h), c,nodes,n) + gaussianQuadrature(i*h, (i+1)*h, @(x) L_0_prime(i+1,h) * L_0_prime(i+1,h), c,nodes,n);

            if i == numTentFunctions
                a_i_plus_1 = 0;
            else
                a_i_plus_1 = gaussianQuadrature(i*h, (i+1)*h, @(x) L_1_prime(i+1,h) * L_0_prime(i+1,h), c,nodes,n);
            end

            b = -gaussianQuadrature((i-1)*h, i*h, @(x) L_1(i,x,h) * f(x), c,nodes,n) - gaussianQuadrature(i*h, (i+1)*h, @(x) L_0(i+1,x,h) * f(x), c,nodes,n);

            y = [a_i_minus_1, a_i, a_i_plus_1];
        end

        function sol = solveF(a,b,f,n)
            % solves function f
            % with boundary conditions
            % u(a) = u(b) = 0
            % with n number of tent functions
            % returns the solution vector sol, u's coefficients

            h = (b-a)/(n+1);

            A = zeros(n,n+2);
            bs = zeros(1,n);

            % each row corresponds to a tentFunction
            for i = 1:n
                [y,b] = Solver.createRow(i, h, f, n);
                A(i,i:i+2) = y;
                bs(i) = b;
            end
            A = A(:, 2:end-1); % remove the first and last columns (extra padding)
            bs = bs'; % make b a column vector

            sol = A\bs;
        end

        function main()
            clear; close all;
            format long;

            % solving
            % \frac{\partial^2 u}{\partial x^2} = f(x)
            % \frac{\partial^2 u}{\partial x^2} = (x^2 - 3x)

            % function f in diff eq
            f = @(x) x^2 - 3*x;

            % show one good approximation
            numTentFunctions = 20;
            a = 0;
            b = 1;
            h = (b-a)/(numTentFunctions+1);
            
            sol = Solver.solveF(a,b,f,numTentFunctions);
            fig = figure;
            x = 0:h:1;
            y = zeros(numTentFunctions+2,1);
            for i = 1:numTentFunctions
                y(i+1) = sol(i); % coeffs are the values
            end
            plot(x,y)
            hold on;
            actual_sol = @(x) (x.^4 / 12) - (1/2).*x.^3 + (5/12).*x;
            x = linspace(a,b);
            plot(x,actual_sol(x),'r');
            legend('Approximation', 'Actual');
            title(sprintf('Approximation (n=%d) vs Actual', numTentFunctions));
            xlabel('x');
            ylabel('y');
            hold off;

            % show convergence with increasing number of tent functions
            % from little (bad approximation) to more (better
            % approximation)
            numTentFunctionsIter = 1:5;
            fig = figure;
            plot(x,actual_sol(x),'r');
            hold on;
            for i = 1:length(numTentFunctionsIter)
               n = numTentFunctionsIter(i);
               sol = Solver.solveF(a,b,f,n);
               h = (b-a)/(n+1); 
               x = a:h:b;
               y = zeros(n+2,1);
               for j = 1:n
                   y(j+1) = sol(j); % coeffs are the values
               end
               plot(x,y);
            end
            legend('Actual', 'n=1', 'n=2', 'n=3', 'n=4', 'n=5');
            title('Approximations vs Actual');
            xlabel('x');
            ylabel('y');
            hold off;
        end
    end
end
