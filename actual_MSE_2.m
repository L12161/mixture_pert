function [MSE,RE] = actual_MSE(a_list,b_list,data,N_loc,W_list,p_list,q_list,alpha_list)
%w list will contain the priv lvl for the users. There should be 90% lvl3 s
%5% lvl 2s and 5% lvl 1s 
repeat = 10;
N_user = length(data);

tab_true = tabulate(data);
% creates a table where the percentage of every entry showing up is
% presented 
count_true = tab_true(:,2);
% count_true holds the true frequency of occurance of the different classes
MSE = zeros(1,repeat);
RE = zeros(N_loc,repeat);

parfor i = 1:repeat
    count_est = est(W_list,N_loc,N_user,a_list,b_list,count_true,p_list,q_list,alpha_list);
    MSE(i) = norm(count_est-count_true,2)^2/N_user;
    RE(:,i) = abs(count_est-count_true)./count_true;
end

MSE = mean(MSE);
RE = mean(RE,2);

end


function est_count = est(W_list,N_loc,N_user,a_list,b_list,count_true,p_list,q_list,alpha_list)

est_count = zeros(N_loc,1);

parfor k = 1:N_loc

    n1 = count_true(k);
    n2 = N_user - n1;
    a = a_list(W_list(k)); 
    b = b_list(W_list(k)); 
    p = p_list(W_list(k));
    q = q_list(W_list(k));
    Alpha = alpha_list;
    y1 = randsrc(n1,1,[[0 1]; [1-a a]]);
    y2 = randsrc(n2,1,[[0 1]; [1-b b]]);
    raw_count = sum([y1;y2]);

    est_count(k) = (raw_count-Alpha*N_user*b-(1-Alpha)*(1-q) ) / (Alpha*a - Alpha*b + (1-Alpha)*q - (1-Alpha)*(1-q) );

end
    
end