%% Abhinav Maheshwari
%% 190028
%% ME685 Programming Assignment III

clear all ;
clc ;
format long ;

%% Specifying the parameters

N = [101;201]; %% for mid of slab

dx1 = 1/(N(1)-1);
T_num = num_T(N(1));
T_numreq = T_num((N(1)+1)/2,:);
DT_Dx = grad_T_numerical(T_num,N(1));
dT_dx = grad_T_ana();
T_ana = analyt_T(dx1*((N(1)+1)/2));
S_1 = [dT_dx,DT_Dx,T_ana,T_numreq'];
        
dx2 = 1/(N(2)-1);
T_num = num_T(N(2));
T_numreq = T_num((N(2)+1)/2,:);
DT_Dx = grad_T_numerical(T_num,N(2));
dT_dx = grad_T_ana();
T_ana = analyt_T(dx2*((N(2)+1)/2));
S_2 = [dT_dx,DT_Dx,T_ana,T_numreq'];

disp('N = 101');
disp(S_1);
disp('N = 201');
disp(S_2);


%% functions

function T_ana = analyt_T(inputs)
    x = inputs;
    T_ana = (1-x)*ones(21,1);
    lambda = 9;
    dt = 5.188184292273184e-04;
    for t = 1:21
    for n = 1:100
        T_ana(t) = T_ana(t) - ((2*sin(n*pi*x))*(lambda + n*n*pi*pi*exp(-n*n*pi*pi*(dt*(t-1)) - lambda*dt*(t-1)))/(n*pi*(lambda + n*n*pi*pi)));
    end
    end
end

function dT_dx = grad_T_ana()
lambda = 9;
dT_dx = -ones(21,1);
x = 0;
dt = 5.188184292273184e-04;
	for t = 1:21
		for n = 1:100
		dT_dx(t) = dT_dx(t) - 2*cos(n*pi*x)*(lambda + n*n*pi*pi*exp(-((t-1)*dt)*(n*n*pi*pi+lambda)))/(n*n*pi*pi+lambda);
		end
	end
end

function DT_Dx = grad_T_numerical(T,N)
    dx = 1/(N-1);
    DT_Dx = zeros(21,1);
    for t = 1:21
		DT_Dx(t) = (T(2,t)-T(1,t))/dx;
    end
end

function T_num = num_T(inputs)
    N = inputs;
    lambda = 9;
        dx = 1/(N-1);
        dt = 5.188184292273184e-04;
        T = zeros(N,1);
        T(1) = 1; r = dt/(dx^2);
        
        for n = 1:20
            stop_loop = 1; [lr,lc] = size(T);
            Tn = T(:,lc);
            Tnp1_old = Tn; Tnp1_new = Tnp1_old;
            while stop_loop
                for i = 2:(N-1)
                    Tnp1_new(i) = (Tn(i) + r*Tnp1_new(i-1) + r*Tnp1_old(i+1))/(1+2*r+lambda*dt);
                end
                maxe = 0;
                for i = 1:101
                    e = abs(Tnp1_new(i)-Tnp1_old(i));
                    if e>maxe
                        maxe = e;
                    end
                end
                if (maxe)<=1e-16
                    T = [T,Tnp1_new]; stop_loop = 0;
                else
                    Tnp1_old = Tnp1_new;
                end
            end
        end
    T_num = T;
end

