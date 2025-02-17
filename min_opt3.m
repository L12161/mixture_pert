function [X,FVAL] = min_opt3(W,fun, d)
N_loc = length(W);

    % paper notation -> MATLAB Variable name
    %       ai ->  x(1 : N_lev) 
    % bi ->  x(N_lev+1 : 2*N_lev)
    % p -> x(2*N_lev+1 : 3*N_lev)
    % q -> x(3*N_lev+1 : 4*N_lev)
    % Alpha -> x(4*N_lev+1)
    % mi -> alpha 

c = @(x) nonlcon(x,N_loc,W,d);         

% probably upper bound and lower bound
LB = [0.3*ones(N_loc,1); zeros(N_loc,1); zeros(N_loc,1); zeros(N_loc,1); zeros(1,1)];    
UB = [ones(N_loc,1); 0.3*ones(N_loc,1); ones(N_loc,1); ones(N_loc,1); ones(1,1)];

% x0 refers to the initial guess for the parameters. we start with all a
% values to be 0.5 and all b values to be 1/(1+exp(min epsilon))
%x0 = [0.5*ones(N_loc,1); 1/(1+exp(min(W)))*ones(N_loc,1) ; ones(N_loc,1)];
x0 = [0.5*ones(N_loc,1); zeros(N_loc,1) ; 0.7*ones(N_loc,1); 0.3*ones(N_loc,1); 0.5*ones(1,1)];




options = optimoptions('fmincon','Algorithm','interior-point');
[X,FVAL,EXITFLAG] = fmincon(fun,x0,[],[],[],[],LB,UB,c,options);


end



function [c, ceq] = nonlcon(x,N_loc,W,d)
% non linear constraint function 

row = @(i,j) N_loc*(i-1) + j;

a = x(1:N_loc); 
b = x(N_loc+1:2*N_loc);
Alpha = x(4*N_loc+1 : end);
p_all = x(2*N_loc+1 : 3*N_loc);
q_all = x(3*N_loc+1 : 4*N_loc);
c = zeros(N_loc^2+3,1);
ceq = zeros(N_loc,1);
for i = 1:N_loc

    p = p_all(i);
    q = q_all(i);
    ceq(i) = p_all(i) + q_all(i)*(d-1) - 1;
    for j = 1:N_loc
        %p = exp(W(i))/(exp(W(i)) + d -1);
        %q = 1 / (exp(W(i)) + d -1);
        
        c(row(i,j)) =               Alpha*p + (1-Alpha)*a(i)*(1-b(j)) - ...
            exp( min(W(i),W(j)) )*( Alpha*q + (1-Alpha)*b(i)*(1-a(j)))  ; % NO CONSTRAINTS FROM Alpha or BETA
        %c(row(i,j)) = (Alpha(i)/Alpha(j))*(p/q)*(d-1)*b(i)*(1-a(j)) + a(i)*(1-b(j)) - exp(min(W(i),W(j)))*b(i)*(1-a(j)); % trying out the A1/A2 format. 
        %c(row(i,j)) =  (1-Beta(i))*a(i)*(1-b(j)) - exp(min(W(i),W(j)))*b(i)*(1-a(j)); 
        %c(row(i,j)) = Alpha(i)*(p/q)*(d-1)*b(i)*(1-a(j)) + (1-Alpha(i))*a(i)*(1-b(j)) - exp(min(W(i),W(j)))*b(i)*(1-a(j)); %  equation 9 constraints, with the form LHS <= 0 
    end
    c(row(N_loc,N_loc)+i,1) = p - q*exp( min(W(i),W(j)) );

end 


end