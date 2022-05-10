%% Abhinav Maheshwari
%% 190028
%% ME685 Programming Assignment 4

clear all ;
clc ;
format long ;

%% Specifying the values
theta = [63.3 52.2 38.1 27.9 19.2 11.4 6.3 2.9] ;
t = [10, 30, 60, 90, 130, 180, 250, 300] ;
Z = zeros(8, 2) ;

%% initial guess
a = 85 ;
b = 0.03 ;
delA = 1 ;
delB = 1 ;

%% Solving iterations till convergence
while (delA >= 1E-5 || delB >= 1E-5)
    Z(:, 1) = exp(-b*t) ;
    Z(:, 2) = (-a*t).*(exp(-b*t)) ;
    B = (theta - (a*exp(-b*t)))' ;
    Zt_Z = (Z')*Z ;                     %% Z(transpose)*Z
    Zt_B = (Z')*B ;                     %% Z(transpose)*B
    delta = Zt_Z\Zt_B ;                 %% calculating the delta matrix
    a = a + delta(1) ;
    b = b + delta(2) ;

    delA = abs(delta(1)) ;
    delB = abs(delta(2)) ;
end

%% regression coefficient 
theta_avg = mean(theta) ;
S_t = 0 ;
S_r = 0 ;
for i = 1:8
    S_t = S_t + (theta(i) - theta_avg)^2 ;
    S_r = S_r + (B(i))^2 ;
end
r = 100*sqrt((S_t - S_r) / S_t) ;

%% Printing results
fprintf('a = %f\n', a) ;
fprintf('b = %f\n', b) ;
fprintf('Regression coefficient r = %f%%\n', r) ;
