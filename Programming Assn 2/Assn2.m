%% Abhinav Maheshwari
%% 190028
%% ME685  

%% domain
length = 1 ;
N = 1001 ;
dx = length / (N-1) ;

%% Problemm Intializa

T=zeros(1,N)
T(1)=1 ;
T(N)=0 ;
K=zeros(1,N)
K(1)=0.9 ;
K(N)=1 ;
for i=2:N
    T(i)=0.5 ;
end
eps = 1e-6 %% for convergence
iter = 0 ;
max_iter = 5000 %% maximum iterations given
for i= 2:N
    K(i)=1;
end

error = 1e5;

%% Gauss Seidel
while iter < max_iter && error > eps
    error = 0;
    for x=2:1000
        Ke = 2/((1/K(x))+(1/K(x+1)));
        Kw = 2/((1/K(x)+(1/K(x-1))));
        T_n = (Ke*T(x+1)+Kw*T(x-1)-((1-(x-1)*0.001)*dx^2))/(Ke+Kw);
        error = max(error,abs(T(x)-T_n));
        T(x)=T_n;
        K(x)=1-0.1*T(x);
    end
    iter = iter + 1;
end
inter = 1:99:1001
T(inter)
plot(T(inter))
