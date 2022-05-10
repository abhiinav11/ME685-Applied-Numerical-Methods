%% Code by Abhinav Maheshwari
%% Roll no 190028
%% ME685 Programming Assignment 1

clc ;
clear all ;

%% mention the initial guess
Q = [5; 5; 5; 5; 5; 5] ;

%% Define the function
F = @(Q) [Q(1)^2 + Q(2)^2 + Q(3)^2 - 14;
    2*(Q(2)^2) + Q(3)^2 + 2*(Q(4)^2) - 35;
    Q(3)^2 + Q(4)^2 - 2*(Q(5)^2) - 10;
    Q(1) + Q(5) - 3;
    Q(4) + Q(6) - 4;
    Q(2) + Q(6) - 3] ;

max_iter = 6 ; % we have to do max 5 iterations, first iteration is for quessed value, so 6
function_value = F(Q);
l = length(Q) ;
E = ( 1e-5*max(1, norm(Q)) ) ; % determining the value of Epsilon, given relaxation factor is unity
Eps = E*eye(l) ; % define a diagonal matrix with all diagonal element equal to Epsilon
count = 0 ;

%% Define the columns for our table
Iteration_Index = [] ;
Guessed_values = {} ;
Corrections = {} ;
Updates = {} ;

%% Use the Newton-Raphson method for iteration
while ( (norm(function_value)*max(1, norm(Q)) > 1e-10) && (count < max_iter) )
    count = count + 1 ;
    for i=1:l
        Diff(:,i) = (F(Q + Eps(:, i)) - F(Q - Eps(:, i))) / (2*E) ; %% Numerical differentiation
    end
    x = [Diff, function_value] ;

    %% Using the Gauss Jordan Method
    for j = 1:(length(x)-1)
        A = x(j, :) ;
        A = A/A(j) ;
        x(j, :) = A ;     
        for k = 1:(length(x)-1)
            if (j ~= k)
                x(k,:) = A*(-1 * x(k, j)) + x(k, :) ;
            end
        end
    end
    dell_Q = -x(:, length(x)) ;

    %% Input values in table
    Iteration_Index(count) = count ;
    Guessed_values{count} = mat2str(Q, 5) ;    %% converting values from matrix to string
    Corrections{count} = mat2str(dell_Q, 5) ;  %% converting values from matrix to string
    Updates{count} = mat2str((Q+dell_Q), 5) ;  %% converting values from matrix to string

    %% Improving the initial guess
    Q = Q + dell_Q ;
    function_value = F(Q) ;
end

%% Fill values in table using the right format
Iteration_Index = Iteration_Index' ;
Guessed_values = Guessed_values' ;
Corrections = Corrections' ;
Updates = Updates' ;

%% Display output in form of a table
Tab = table(Iteration_Index, Guessed_values, Corrections, Updates) ;
disp(Tab) ;

