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
N_lev = 3; W = [1 1.2 2];
Prob = [0.05 0.05 0.9];      % same prob distribution for every attribute

% ---- attributes: same Prob, different domain size N_loc ----

N_loc_list = [74 9 28523 16 16 7 15 6 5 2 123 99 96 42 4];   % domain size per attribute (edit as needed)
% domain set to [74 9 28523 16 16 7 15 6 5 2 123 99 96 42 4] for UCI Adult
% dataset bassed on the dataset domain distribution. 
N_attr = length(N_loc_list);

N_user = 100000;   % users sampled to estimate the effective domain size each run

myEpsilon = [0.5:0.5:6];          % epsilon used PER attribute (x-axis)
N_epsilon = length(myEpsilon);

% ---- Moments-accountant knobs (now applied per-mechanism, not just the mixture) ----
% Each mechanism's OWN moment-accountant function (paper Prop E.1 for GRR,
% Prop E.2 for BUE/unary-type, Theorem 5.2 for the mixture) replaces the
% earlier Gaussian-surrogate formula -- see helper functions at the bottom.
MA_delta      = 1e-5;             % target delta for the (eps,delta) budget
MA_alpha_grid = 1.01:0.01:200;    % RDP orders (lambda) to minimize over

% ---- Repeated-run settings ----
% Only GRR and the opt3 mixture depend on domain_size (hence on
% generate_powlaw's draw); RAPPOR, OUE, and IDUE-opt0/1/2 do not, so their
% across-run deviation will be exactly zero -- that's expected, not a bug.
N_runs = 20;   % number of independent generate_powlaw draws to average over

% Per-run accumulators (3rd dimension = run index)
% 1=RAPPOR, 2=OUE, 3=GRR
result_MA_runs      = zeros(3, N_epsilon, N_runs);
result_min_MA_runs  = zeros(4, N_epsilon, N_runs);   % opt0, opt1, opt2, opt3
Pi_runs              = zeros(N_attr, N_epsilon, N_runs);
domain_size_runs     = zeros(N_attr, N_runs);

% ---- Per-attribute eps for RAPPOR and OUE, via their OWN moment-accountant ----
% (paper Prop E.2). These don't depend on domain_size, so compute once,
% outside the run/attribute loops, instead of redoing it 20x per epsilon.
% eps_naive_attr_list is kept too since it's still needed for GRR/mix's
% domain-dependent solve, and for the summary table.
eps_naive_attr_list  = zeros(1, N_epsilon);
eps_actual_RAPPOR_list = zeros(1, N_epsilon);
eps_actual_OUE_list    = zeros(1, N_epsilon);
for i = 1:N_epsilon
    eps_total = myEpsilon(i);
    eps_naive_attr_list(i)    = eps_total / N_attr;
    eps_actual_RAPPOR_list(i) = solve_eps_loose_mechanism(@alpha_rappor_query, eps_total, N_attr, MA_delta, MA_alpha_grid);
    eps_actual_OUE_list(i)    = solve_eps_loose_mechanism(@alpha_oue_query,    eps_total, N_attr, MA_delta, MA_alpha_grid);
end

% GRR and the mixture DO depend on domain_size (redrawn each run), so their
% per-attribute eps is tracked per (epsilon, run) -- k=1 shown/used since
% N_attr=1 in the current N_loc_list; generalizes if N_loc_list grows.
eps_actual_GRR_runs = zeros(N_epsilon, N_runs);
eps_actual_mix_runs = zeros(N_epsilon, N_runs);

%% ================= REPEATED RUNS =================
for run = 1:N_runs
    fprintf('=== Run %d / %d ===\n', run, N_runs);

    % ---- Effective domain size (matches PG3_synth_MSE.m), redrawn each run ----
    % Requires generate_powlaw.m on the path.
    domain_size_list = N_loc_list;
    domain_size_runs(:,run) = domain_size_list(:);
    fprintf('    domain_size = %s\n', mat2str(domain_size_list));

    result_sum_MA     = zeros(3, N_epsilon);
    result_min_sum_MA = zeros(4, N_epsilon);   % opt0, opt1, opt2, opt3
    Pi_sum             = zeros(N_attr, N_epsilon);

    % ---- LOOP OVER PER-ATTRIBUTE EPSILON (x-axis) ----
    for i = 1:N_epsilon
        eps_total = myEpsilon(i);   % TOTAL per-user budget for this x-axis point

        eps_actual_RAPPOR = eps_actual_RAPPOR_list(i);   % precomputed above (domain-independent)
        eps_actual_OUE    = eps_actual_OUE_list(i);      % precomputed above (domain-independent)

        % ---------- accumulate over attributes ----------
        for k = 1:N_attr
            N_loc = N_loc_list(k);
            alpha = N_loc*Prob;
            domain_size = domain_size_list(k);   % from this run's generate_powlaw sample

            % GRR and the mixture DO depend on domain_size, so solve per
            % (epsilon, attribute, run) using their own moment-accountant
            % functions (paper Prop E.1 for GRR, Theorem 5.2 for the mixture).
            eps_actual_GRR = solve_eps_loose_mechanism(@(e,l) alpha_grr_query(e,domain_size,l), eps_total, N_attr, MA_delta, MA_alpha_grid);
            eps_actual_mix = solve_eps_loose_mechanism(@(e,l) alpha_mix_query(e,domain_size,l), eps_total, N_attr, MA_delta, MA_alpha_grid);
            if k == 1
                eps_actual_GRR_runs(i,run) = eps_actual_GRR;
                eps_actual_mix_runs(i,run) = eps_actual_mix;
            end

            % IDUE-opt0/1/2 objectives
            fun0 = @(x) alpha*((x(N_lev+1:end)-x(N_lev+1:end).^2)./((x(1:N_lev)-x(N_lev+1:end)).^2))...
                + max((1-x(1:N_lev)-x(N_lev+1:end))./(x(1:N_lev)-x(N_lev+1:end)));
            fun1 = @(x) alpha*(exp(x)./((exp(x)-1).^2));
            fun2 = @(b) alpha*((b-b.^2)./((0.5-b).^2)) + 1;

            % mixture objective (min_opt3)
            fun3 = @(x) x(4*N_lev+1)*sum( (x(2*N_lev+1:3*N_lev).*(1-x(2*N_lev+1:3*N_lev))) ...
                            ./ (x(2*N_lev+1:3*N_lev)-x(3*N_lev+1:4*N_lev)).^2 ) ...
                + sum( alpha(:) .* (1-x(4*N_lev+1)) .* x(N_lev+1:2*N_lev).*(1-x(N_lev+1:2*N_lev)) ...
                            ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)).^2 ) ...
                + max( (1-x(4*N_lev+1)) .* (1-x(1:N_lev)-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)) );

            % ---- Which mechanism-specific eps feeds each optimizer ----
            % opt0: general BUE (free a,b) -- approximated at OUE's a=0.5
            %       (the known variance-minimizer for unary encoding), so it
            %       gets OUE's moment-accountant eps.
            % opt1: matches symmetric UE/RAPPOR's exact parametrization
            %       (constraint x(i)+x(j)<=min(W(i),W(j)), fun1's exp(x) form).
            % opt2: matches OUE EXACTLY (x0, fun2 both = OUE's (0.5, q) form).
            % opt3: mixture -- uses the beta-independent max(GRR,BUE) bound.
            temp_opt0 = W*eps_actual_OUE;
            temp_opt1 = W*eps_actual_RAPPOR;
            temp_opt2 = W*eps_actual_OUE;
            temp_opt3 = W*eps_actual_mix;

            % ---- Closed-form pure mechanisms, each at its OWN actual eps ----
            % RAPPOR
            result_sum_MA(1,i) = result_sum_MA(1,i) + N_loc*exp(eps_actual_RAPPOR/2)./((exp(eps_actual_RAPPOR/2)-1).^2);
            % OUE
            result_sum_MA(2,i) = result_sum_MA(2,i) + N_loc*4*exp(eps_actual_OUE)./((exp(eps_actual_OUE)-1).^2)+1;
            % GRR (perturbation): summed MSE = N_loc*(domain_size-2+e^eps)/(e^eps-1)^2
            result_sum_MA(3,i) = result_sum_MA(3,i) + N_loc*(domain_size-2+exp(eps_actual_GRR))./((exp(eps_actual_GRR)-1).^2);

            % ---- IDUE-opt0/1/2/3, each at its own mechanism-specific eps ----
            [~, f0m] = min_opt0(temp_opt0,fun0);
            [~, f1m] = min_opt1(temp_opt1,fun1);
            [~, f2m] = min_opt2(temp_opt2,fun2);
            [X3, f3m] = min_opt3(temp_opt3,fun3,domain_size);
            result_min_sum_MA(1,i) = result_min_sum_MA(1,i) + f0m;
            result_min_sum_MA(2,i) = result_min_sum_MA(2,i) + f1m;
            result_min_sum_MA(3,i) = result_min_sum_MA(3,i) + f2m;
            result_min_sum_MA(4,i) = result_min_sum_MA(4,i) + f3m;

            Pi_sum(k,i) = X3(4*N_lev+1);
        end
    end

    result_MA_runs(:,:,run)     = result_sum_MA;
    result_min_MA_runs(:,:,run) = result_min_sum_MA;
    Pi_runs(:,:,run)            = Pi_sum;
end

% ---- Aggregate across the N_runs draws: mean and std ----
result_MA_mean      = mean(result_MA_runs, 3);
result_MA_std       = std(result_MA_runs, 0, 3);
result_min_MA_mean  = mean(result_min_MA_runs, 3);
result_min_MA_std   = std(result_min_MA_runs, 0, 3);

result_MA        = result_MA_mean;      % kept for naming compatibility below
result_min_MA    = result_min_MA_mean;

%% ===================== MSE PLOT (PG2-style, no shading) =====================
figure('Position',[100,100,1000,1000]);
ah1 = TightPlots(1, 1, 800,[10 7],[80,80],[50,50,60,40],[70,10],'pixels');
axes(ah1(1));

xlim_min = min(myEpsilon); xlim_max = max(myEpsilon);

all_mse_data = [result_MA(1,:), result_MA(2,:), result_MA(3,:), ...
                result_min_MA(1,:), result_min_MA(2,:), result_min_MA(3,:), result_min_MA(4,:)];
mse_min_y = min(all_mse_data);
mse_max_y = max(all_mse_data);
mse_y_range = mse_max_y - mse_min_y;
mse_ylim_min = max(1, mse_min_y - 0.1*mse_y_range);
mse_ylim_max = mse_max_y + 0.1*mse_y_range;

% GRR gets its own color (Color has 5 rows)
Color_GRR = [0.85 0.33 0.10];   % orange-red for GRR

% Line widths: opt3 (mixture) bolder/thicker than the rest
LW_norm = 1.3;   % all mechanisms except opt3
LW_opt3 = 3.0;   % opt3 mixture: thicker

hold on;

% ---------- MA set only: SOLID lines ----------
hm1 = plot(myEpsilon, result_MA(1,:),      '-', 'Color', Color(1,:),  'LineWidth', LW_norm);
hm2 = plot(myEpsilon, result_MA(2,:),      '-', 'Color', Color(2,:),  'LineWidth', LW_norm);
hg  = plot(myEpsilon, result_MA(3,:),      '-', 'Color', Color_GRR,   'LineWidth', LW_norm);
hm3 = plot(myEpsilon, result_min_MA(1,:),  '-', 'Color', Color(3,:),  'LineWidth', LW_norm);
hm4 = plot(myEpsilon, result_min_MA(2,:),  '-', 'Color', Color(4,:),  'LineWidth', LW_norm);
hm5 = plot(myEpsilon, result_min_MA(3,:),  '-', 'Color', Color(5,:),  'LineWidth', LW_norm);
hm6 = plot(myEpsilon, result_min_MA(4,:),  '-', 'Color', [0.5 0 0.5], 'LineWidth', LW_opt3);

% ---------- Markers (one per mechanism, filled) ----------
MS_norm = 6; MS_opt3 = 9;
plot(myEpsilon, result_MA(1,:),     'o', 'Color', Color(1,:),  'MarkerSize', MS_norm, 'MarkerFaceColor', Color(1,:),  'LineStyle', 'none');
plot(myEpsilon, result_MA(2,:),     's', 'Color', Color(2,:),  'MarkerSize', MS_norm, 'MarkerFaceColor', Color(2,:),  'LineStyle', 'none');
plot(myEpsilon, result_MA(3,:),     'v', 'Color', Color_GRR,   'MarkerSize', MS_norm, 'MarkerFaceColor', Color_GRR,   'LineStyle', 'none');
plot(myEpsilon, result_min_MA(1,:), '^', 'Color', Color(3,:),  'MarkerSize', MS_norm, 'MarkerFaceColor', Color(3,:),  'LineStyle', 'none');
plot(myEpsilon, result_min_MA(2,:), '*', 'Color', Color(4,:),  'MarkerSize', MS_norm, 'MarkerFaceColor', Color(4,:),  'LineStyle', 'none');
plot(myEpsilon, result_min_MA(3,:), 'd', 'Color', Color(5,:),  'MarkerSize', MS_norm, 'MarkerFaceColor', Color(5,:),  'LineStyle', 'none');
plot(myEpsilon, result_min_MA(4,:), 'p', 'Color', [0.5 0 0.5], 'MarkerSize', MS_opt3, 'MarkerFaceColor', [0.5 0 0.5], 'LineStyle', 'none');

% ---------- Deviation markers: mean +/- std across the N_runs draws ----------
% RAPPOR/OUE/opt0/opt1/opt2 don't depend on domain_size, so their std is 0
% and these whiskers will simply not appear for those series.
EB_CapSize = 4;
errorbar(myEpsilon, result_MA(1,:),     result_MA_std(1,:),     'LineStyle','none', 'Color', Color(1,:),  'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_MA(2,:),     result_MA_std(2,:),     'LineStyle','none', 'Color', Color(2,:),  'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_MA(3,:),     result_MA_std(3,:),     'LineStyle','none', 'Color', Color_GRR,   'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_min_MA(1,:), result_min_MA_std(1,:), 'LineStyle','none', 'Color', Color(3,:),  'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_min_MA(2,:), result_min_MA_std(2,:), 'LineStyle','none', 'Color', Color(4,:),  'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_min_MA(3,:), result_min_MA_std(3,:), 'LineStyle','none', 'Color', Color(5,:),  'CapSize', EB_CapSize, 'HandleVisibility','off');
errorbar(myEpsilon, result_min_MA(4,:), result_min_MA_std(4,:), 'LineStyle','none', 'Color', [0.5 0 0.5], 'CapSize', EB_CapSize, 'LineWidth', 1.5, 'HandleVisibility','off');

xlim([xlim_min, xlim_max]);
ylim([mse_ylim_min, mse_ylim_max]);
set(gca, 'YScale', 'log');

set(gca, 'XTick', myEpsilon);
set(gca, 'XTickLabel', arrayfun(@(v) sprintf('%g', v), myEpsilon, 'UniformOutput', false));

yt_lo = floor(log10(mse_ylim_min));
yt_hi = ceil(log10(mse_ylim_max));
decade_ticks = 10.^(yt_lo:yt_hi);
half_decade_ticks = 5 * 10.^(yt_lo:yt_hi-1);
yticks_all = sort([decade_ticks, half_decade_ticks]);
yticks_all = yticks_all(yticks_all >= mse_ylim_min & yticks_all <= mse_ylim_max);
set(gca, 'YTick', yticks_all);
set(gca, 'YTickLabel', arrayfun(@(v) sprintf('%g', v), yticks_all, 'UniformOutput', false));
set(gca, 'YMinorTick', 'on');
set(gca, 'TickDir', 'out');
box on;

xlabel('$\epsilon_{\mathrm{total}}$', 'FontSize', 30, 'FontWeight', 'bold');
ylabel('MSE', 'FontSize', 24, 'FontWeight', 'bold');
legend([hm1 hm2 hg hm3 hm4 hm5 hm6], ...
    {'RAPPOR','OUE','GRR','IDUE-opt0','IDUE-opt1','IDUE-opt2','IDUE-opt3 mix (MA)'}, ...
     'FontSize', 15, 'Location', 'northeast');
grid on;

annotation('textbox',...
    [0.386 0.839909089283511 0.10950000231266 0.0823863652619449],...
    'String',{'High','Privacy'}, 'EdgeColor','none', 'FontSize', 16);
annotation('textbox',...
    [0.549 0.217045452919874 0.106500002241135 0.0823863652619449],...
    'String',{'Low','Privacy'}, 'EdgeColor','none', 'FontSize', 16);

print(gcf, sprintf('MSE_Comparison_IDLDP_multiattr_MA_mean%druns.png', N_runs), '-dpng', '-r600');

%% ===================== MSE TABLE (mean +/- std across N_runs) =====================
fprintf('\n');
fprintf('Mean +/- std across %d runs (each run redraws domain_size via generate_powlaw)\n', N_runs);
fprintf('%8s %16s %16s %16s %16s %16s %16s %16s\n', ...
    'Epsilon', 'RAPPOR', 'OUE', 'GRR', 'IDUE-opt0', 'IDUE-opt1', 'IDUE-opt2', 'Mix(opt3)');
for i = 1:N_epsilon
    fprintf('%8.2f %7.2f+-%6.2f %7.2f+-%6.2f %7.2f+-%6.2f %7.2f+-%6.2f %7.2f+-%6.2f %7.2f+-%6.2f %7.2f+-%6.2f\n', ...
        myEpsilon(i), ...
        result_MA(1,i),      result_MA_std(1,i), ...      % RAPPOR
        result_MA(2,i),      result_MA_std(2,i), ...      % OUE
        result_MA(3,i),      result_MA_std(3,i), ...      % GRR
        result_min_MA(1,i),  result_min_MA_std(1,i), ...  % IDUE-opt0
        result_min_MA(2,i),  result_min_MA_std(2,i), ...  % IDUE-opt1
        result_min_MA(3,i),  result_min_MA_std(3,i), ...  % IDUE-opt2
        result_min_MA(4,i),  result_min_MA_std(4,i));     % Mix (opt3)
end
fprintf('\n');

% report how the TOTAL budget splits per attribute under each mechanism's
% own moment accountant (GRR/mix averaged across the N_runs domain_size draws)
eps_actual_GRR_mean = mean(eps_actual_GRR_runs, 2);
eps_actual_mix_mean = mean(eps_actual_mix_runs, 2);
fprintf('\nEps_total\tNaive(=/N_attr)\tRAPPOR\t\tOUE\t\tGRR(mean)\tMix(mean)\n');
for i = 1:N_epsilon
    fprintf('%.2f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n', myEpsilon(i), ...
        eps_naive_attr_list(i), eps_actual_RAPPOR_list(i), eps_actual_OUE_list(i), ...
        eps_actual_GRR_mean(i), eps_actual_mix_mean(i));
end

%% ===================== Mechanism-specific moment-accountant functions =====================
% These replace the earlier Gaussian-surrogate ma_compose_eps/solve_sigma_for_eps.
% They implement the paper's OWN moment-accountant results directly:
%   - alpha_grr_query  : Proposition E.1 (GRR)
%   - alpha_bue_query  : Proposition E.2 (BUE / unary-type, homogeneous case)
%   - alpha_rappor_query, alpha_oue_query: BUE specialized to symmetric UE and OUE
%   - alpha_mix_query  : Theorem 5.2 upper bound, using the fact that
%       log(beta*e^u+(1-beta)*e^v) is convex in beta on [0,1], so its max
%       over beta in [0,1] is attained at an endpoint. This gives a bound
%       that's valid for WHATEVER mixing weight the optimizer picks,
%       without needing to know it in advance (avoids the circularity of
%       beta being an optimization output, not an input).
% Composition (Eq 5) and conversion (Eq 7) are exactly as before, just
% applied per-mechanism instead of only for the mixture.
%
% IMPORTANT: MA/RDP composition is NOT guaranteed to beat naive linear
% composition for every mechanism and (eps_total, N_attr) combination --
% confirmed empirically here: GRR benefits hugely at large domain sizes
% (its per-query divergence shrinks as d grows), while OUE/RAPPOR's
% per-query divergence does NOT shrink with d and can end up needing a
% SMALLER per-attribute eps than naive composition would allow at moderate
% eps_total. solve_eps_loose_mechanism() therefore always returns
% max(naive, MA-derived), so results never regress below naive composition.

function out = logsumexp_rows(M)
% Numerically stable log-sum-exp down each COLUMN of M (M is rows x N).
% Needed because the raw p^lambda*q^(1-lambda) terms overflow/underflow
% for large lambda (up to 200) combined with small q (large domain sizes).
    m = max(M, [], 1);
    out = m + log(sum(exp(M - m), 1));
end

function a = alpha_grr_query(eps, d, lambda)
% Paper Proposition E.1 (homogeneous case, p_i=p_j=p, q_i=q_j=q): moment
% accountant of ONE GRR query, vectorized over RDP order(s) lambda.
    logp = eps - log(exp(eps) + d - 1);
    logq = -log(exp(eps) + d - 1);
    t1 = lambda.*logp + (1-lambda).*logq;
    t2 = lambda.*logq + (1-lambda).*logp;
    t3 = (log(d-2) + logq) * ones(size(lambda));
    a = logsumexp_rows([t1; t2; t3]);
end

function a = alpha_bue_query(a_prob, b_prob, lambda)
% Paper Proposition E.2 (homogeneous case): moment accountant of ONE
% BUE/unary query with bit-flip probabilities (a_prob, b_prob), vectorized
% over lambda. Covers RAPPOR (symmetric a+b=1), OUE (a=0.5 fixed), and the
% general BUE used by opt0.
    loga = log(a_prob);   log1ma = log(1-a_prob);
    logb = log(b_prob);   log1mb = log(1-b_prob);
    T1 = logsumexp_rows([lambda.*loga+(1-lambda).*logb;  lambda.*log1ma+(1-lambda).*log1mb]);
    T2 = logsumexp_rows([lambda.*logb+(1-lambda).*loga;  lambda.*log1mb+(1-lambda).*log1ma]);
    a = T1 + T2;
end

function a = alpha_rappor_query(eps, lambda)
% Symmetric unary encoding (RAPPOR/SUE, paper's "UE"): a=p, b=q, p+q=1.
    p = exp(eps/2)/(exp(eps/2)+1);
    q = 1/(exp(eps/2)+1);
    a = alpha_bue_query(p, q, lambda);
end

function a = alpha_oue_query(eps, lambda)
% OUE: a=1/2 fixed, b=q=1/(e^eps+1). Also used for opt0 (approximation,
% since a=0.5 is the known variance-minimizer for free-a,b BUE) and opt2
% (exact match -- opt2's x0/fun2 ARE OUE's own parametrization).
    q = 1/(exp(eps)+1);
    a = alpha_bue_query(0.5, q, lambda);
end

function a = alpha_mix_query(eps, d, lambda)
% Mixture (Theorem 5.2): alpha_mix(lambda) <= log[beta*e^{a_GRR}+(1-beta)*e^{a_BUE}].
% beta-independent worst-case bound: max(a_GRR, a_BUE) (see header comment).
    aGRR = alpha_grr_query(eps, d, lambda);
    aBUE = alpha_oue_query(eps, lambda);
    a = max(aGRR, aBUE);
end

function eps_actual = solve_eps_loose_mechanism(alpha_query_fn, eps_total, N_attr, delta, alpha_grid)
% Finds the per-attribute eps whose N_attr-fold moment-accountant
% composition (Eq 5: alpha_total = N_attr*alpha_per_query; Eq 7:
% eps = min_lambda[(alpha_total+log(1/delta))/(lambda-1)]) hits eps_total,
% for the given mechanism's alpha_query_fn(eps_per_attr, lambda_grid).
% Returns max(naive, MA-derived) -- see header comment on why this matters.
    eps_naive = eps_total / N_attr;
    logdelta_term = log(1/delta) ./ (alpha_grid - 1);
    total_eps_at = @(eps_pa) min( N_attr .* alpha_query_fn(eps_pa, alpha_grid) + logdelta_term );

    lo = 1e-8; hi = max(eps_total, 1e-6);
    for it = 1:60
        if total_eps_at(hi) >= eps_total
            break;
        end
        hi = hi * 2;
    end
    for it = 1:80
        mid = 0.5*(lo+hi);
        if total_eps_at(mid) <= eps_total
            lo = mid;
        else
            hi = mid;
        end
        if (hi-lo) < 1e-10*max(1,hi)
            break;
        end
    end
    eps_actual = max(eps_naive, lo);
end