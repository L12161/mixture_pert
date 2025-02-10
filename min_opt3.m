function [X,FVAL] = min_opt3(W,fun)
N_loc = length(W);

c = @(x) nonlcon(x,N_loc,W);         

% probably upper bound and lower bound
LB = [0.5*ones(N_loc,1); 0.0001*ones(N_loc,1); 0.0001*ones(N_loc,1)];    
UB = [0.9999*ones(N_loc,1); 0.5*ones(N_loc,1); 0.9999*ones(N_loc,1)];

% x0 refers to the initial guess for the parameters. we start with all a
% values to be 0.5 and all b values to be 1/(1+exp(min epsilon))
x0 = [0.5*ones(N_loc,1); 1/(1+exp(min(W)))*ones(N_loc,1) ; 0.5*ones(N_loc,1)];


options = optimoptions('fmincon','Algorithm','sqp');
[X,FVAL,EXITFLAG] = fmincon(fun,x0,[],[],[],[],LB,UB,c,options);


end



function [c, ceq] = nonlcon(x,N_loc,W)
% non linear constraint function 

row = @(i,j) N_loc*(i-1) + j;

a = x(1:N_loc); 
b = x(N_loc+1:2*N_loc);
Beta = x(2*N_loc : end)

c = zeros(N_loc^2,1);
for i = 1:N_loc
    for j = 1:N_loc
        p = exp(W(i))/(W(i) + 100 -1);
        q = 1 / (W(i) + 100 -1);
        c(row(i,j)) = Beta(i)*(p/q)*(100-1)*b(i)*(1-a(j)) + (1-Beta(i))*a(i)*(1-b(j)) - exp(min(W(i),W(j)))*b(i)*(1-a(j)); %  equation 9 constraints, with the form LHS <= 0 
    end
end

ceq = [];
end