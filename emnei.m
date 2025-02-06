%% Plot setting
clc; clear; close all;
set(0,'defaulttextinterpreter','latex'); 
set(0,'defaultlinelinewidth',1.3); 
set(0,'DefaultLineMarkerSize',6); 
set(0,'DefaultTextFontSize', 14); 
set(0,'DefaultAxesFontSize',14); 

figure;
h = plot([1 2],rand(6,2));
c = get(h,'Color');
Color = cell2mat(c);
Color = Color([1 2 5 3 4],:);
Color(3,:) = [0.13,0.55,0.13];
close all;


%%

N_lev = 3;  % the different stages of elsilon values
W = [1 1.2 2]; % these terms are? 
N_loc = 100;
N_user = 100000; 
% no of synthetic users 
Prob = [0.05 0.05 0.9];
alpha = N_loc*Prob; 
% this was the alpha for power law, used to generate the synthetic dataset
% here alpha comes out to be 5 5 90



W_list = [];
for j = 1:N_lev
    W_list = [W_list; j*ones(alpha(j),1)];
end     %  W_list contains 5 1's 5 2's and 90 3's after this loop 
W_list = W_list(randperm(N_loc));   % randomly taken 100 samples from the previous W  
W_tab = tabulate(W_list);

W_index = cell(N_lev,1);
for j = 1:N_lev
    W_index{j} = find(W_list==j);
end


data = generate_powlaw(N_user,N_loc);
thamenn