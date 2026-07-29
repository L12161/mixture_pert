% Single Attribute Synthetic Data Distribution 
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

N_loc_list = [100, 100, 100];   % domain size per attribute (edit as needed)
N_attr = length(N_loc_list);

N_user = 100000;   % users sampled to estimate the effective domain size each run

myEpsilon = [0.5:0.5:6];          % epsilon used PER attribute (x-axis)
N_epsilon = length(myEpsilon);

% ---- Moments-accountant knobs (mixture mechanism ONLY) ----
MA_delta       = 1e-5;   % target delta for the (eps,delta) budget
MA_alpha_grid  = 1.01:0.01:200;   % RDP orders to minimize over (used only if MA_alpha_fixed = [])
MA_alpha_fixed = [];     % [] = optimize alpha (tightest eps).  Set e.g. 10 to PIN a single RDP order.

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
eps_loose_list       = zeros(1, N_epsilon);       % per-attr eps for mixture (MA); run-independent
eps_naive_attr_list  = zeros(1, N_epsilon);       % per-attr eps for naive mechanisms; run-independent

%% ================= REPEATED RUNS =================
for run = 1:N_runs
    fprintf('=== Run %d / %d ===\n', run, N_runs);

    % ---- Effective domain size (matches PG3_synth_MSE.m), redrawn each run ----
    % Requires generate_powlaw.m on the path.
    domain_size_list = zeros(1, N_attr);
    for k = 1:N_attr
        data = generate_powlaw(N_user, N_loc_list(k));  % can be chaanged to implement different synthetic distributions 
        domain_size_list(k) = numel(unique(data));
    end
    domain_size_runs(:,run) = domain_size_list(:);
    fprintf('    domain_size = %s\n', mat2str(domain_size_list));

    result_sum_MA     = zeros(3, N_epsilon);
    result_min_sum_MA = zeros(4, N_epsilon);   % opt0, opt1, opt2, opt3
    Pi_sum             = zeros(N_attr, N_epsilon);

    % ---- LOOP OVER PER-ATTRIBUTE EPSILON (x-axis) ----
    for i = 1:N_epsilon
        epsilon = myEpsilon(i);     % epsilon per attribute (spent by NON-mixture mechs)

        % ---------- BUDGET SPLIT: x-axis epsilon is the TOTAL user budget ----------
        % The user has a fixed TOTAL privacy budget eps_total = epsilon, spent across
        % all N_attr attributes (sequential queries on the same user).
        eps_total = epsilon;

        % (A) NAIVE / sequential composition (RAPPOR, OUE, GRR, opt0/1/2):
        %     budget adds linearly  =>  each attribute gets eps_total / N_attr.
        eps_naive_attr = eps_total / N_attr;

        % (B) MOMENTS ACCOUNTANT (opt3 mixture ONLY):
        %     RDP composition is TIGHTER than linear, so for the SAME eps_total the
        %     per-attribute budget it permits is LARGER than eps_total/N_attr.
        %     Back-solve the per-attribute eps whose N_attr-fold Gaussian RDP
        %     composition equals (eps_total, delta):
        %       1) find sigma s.t. N_attr-fold composition hits eps_total,
        %       2) map that sigma to the SINGLE-release (per-attribute) eps.
        ma_eps = @(sigma) ma_compose_eps(sigma, N_attr, MA_delta, MA_alpha_grid, MA_alpha_fixed);
        sigma_star = solve_sigma_for_eps(ma_eps, eps_total);
        eps_loose  = ma_compose_eps(sigma_star, 1, MA_delta, MA_alpha_grid, MA_alpha_fixed);

        % eps_loose/eps_naive_attr don't depend on domain_size or run, but
        % it's cheap enough to just recompute every run rather than branch.
        eps_loose_list(i)      = eps_loose;       % per-attribute eps for the mixture (MA)
        eps_naive_attr_list(i) = eps_naive_attr;  % per-attribute eps for everyone else
        % --------------------------------------------------------------------------

        % ---------- accumulate over attributes ----------
        for k = 1:N_attr
            N_loc = N_loc_list(k);
            alpha = N_loc*Prob;
            domain_size = domain_size_list(k);   % from this run's generate_powlaw sample

            % IDUE-opt0/1/2 objectives (evaluated at the LOOSER eps via temp_MA below)
            fun0 = @(x) alpha*((x(N_lev+1:end)-x(N_lev+1:end).^2)./((x(1:N_lev)-x(N_lev+1:end)).^2))...
                + max((1-x(1:N_lev)-x(N_lev+1:end))./(x(1:N_lev)-x(N_lev+1:end)));
            fun1 = @(x) alpha*(exp(x)./((exp(x)-1).^2));
            fun2 = @(b) alpha*((b-b.^2)./((0.5-b).^2)) + 1;

            % mixture objective (min_opt3) also uses the LOOSER epsilon
            fun3 = @(x) x(4*N_lev+1)*sum( (x(2*N_lev+1:3*N_lev).*(1-x(2*N_lev+1:3*N_lev))) ...
                            ./ (x(2*N_lev+1:3*N_lev)-x(3*N_lev+1:4*N_lev)).^2 ) ...
                + sum( alpha(:) .* (1-x(4*N_lev+1)) .* x(N_lev+1:2*N_lev).*(1-x(N_lev+1:2*N_lev)) ...
                            ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)).^2 ) ...
                + max( (1-x(4*N_lev+1)) .* (1-x(1:N_lev)-x(N_lev+1:2*N_lev)) ./ (x(1:N_lev)-x(N_lev+1:2*N_lev)) );

            temp_naive = W*eps_naive_attr;   % naive per-attribute eps = eps_total/N_attr
            temp_MA    = W*eps_loose;         % MA per-attribute eps (opt3 mixture ONLY)

            % ---- Closed-form pure mechanisms at NAIVE per-attribute epsilon ----
            % RAPPOR
            result_sum_MA(1,i) = result_sum_MA(1,i) + N_loc*exp(eps_naive_attr/2)./((exp(eps_naive_attr/2)-1).^2);
            % OUE
            result_sum_MA(2,i) = result_sum_MA(2,i) + N_loc*4*exp(eps_naive_attr)./((exp(eps_naive_attr)-1).^2)+1;
            % GRR (perturbation): summed MSE = N_loc*(domain_size-2+e^eps)/(e^eps-1)^2
            % (outer N_loc = population size, inner domain_size = sampled effective
            % domain -- same convention as PG3_synth_MSE.m)
            result_sum_MA(3,i) = result_sum_MA(3,i) + N_loc*(domain_size-2+exp(eps_naive_attr))./((exp(eps_naive_attr)-1).^2);

            % ---- IDUE-opt0/1/2 at NAIVE epsilon; opt3 mixture at eps_loose (MA) ----
            [~, f0m] = min_opt0(temp_naive,fun0);
            [~, f1m] = min_opt1(temp_naive,fun1);
            [~, f2m] = min_opt2(temp_naive,fun2);
            [X3, f3m] = min_opt3(temp_MA,fun3,domain_size);   % <-- ONLY the mixture uses MA eps_loose
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

% report how the TOTAL budget splits per attribute under each accounting
fprintf('\nEps_total\tEps/attr(naive=/N_attr)\tEps/attr(MA mixture)\n');
for i = 1:N_epsilon
    fprintf('%.2f\t\t%.4f\t\t\t%.4f\n', myEpsilon(i), eps_naive_attr_list(i), eps_loose_list(i));
end

%% ===================== Moments-accountant helper functions =====================
function eps = ma_compose_eps(sigma, N, delta, alpha_grid, alpha_fixed)
% (eps,delta) for N-fold composition of a unit-sensitivity Gaussian mechanism
% at noise sigma, via RDP.  rho(a) = N * a / (2 sigma^2).
%
% RDP->(eps,delta):  eps(a) = N*a/(2 sigma^2) + log(1/delta)/(a-1),   a > 1.
%
% Mode A (alpha_fixed non-empty): use the SINGLE pinned order alpha_fixed.
%   Still a valid (eps,delta) bound, but generally looser than the optimum.
% Mode B (alpha_fixed = []): optimize alpha.  eps(a) is convex in a with
%   closed-form minimizer a* = 1 + sqrt( 2 sigma^2 ln(1/delta) / N ).
%   We evaluate at a* and also scan alpha_grid as a guard, taking the min.
    L = log(1/delta);

    if nargin >= 5 && ~isempty(alpha_fixed)
        % ---- Mode A: pinned alpha ----
        a = alpha_fixed;
        if a <= 1
            error('MA_alpha_fixed must be > 1 (RDP order).');
        end
        eps = N*a/(2*sigma^2) + L/(a - 1);
        return;
    end

    % ---- Mode B: optimized alpha ----
    a_star = 1 + sqrt( 2*sigma^2 * L / N );
    a_star = max(a_star, 1 + 1e-12);          % enforce a > 1
    eps_star = N*a_star/(2*sigma^2) + L/(a_star - 1);

    a = alpha_grid(:);
    a = a(a > 1);
    eps_grid = min( N.*a./(2*sigma^2) + L./(a - 1) );

    eps = min(eps_star, eps_grid);
end

function sigma = solve_sigma_for_eps(ma_eps_fun, eps_target)
% Solve ma_eps_fun(sigma) = eps_target.  ma_eps is decreasing in sigma.
% Bisection on sigma.
    lo = 1e-3; hi = 1e3;
    % ensure bracketing: eps(lo) large, eps(hi) small
    for it = 1:200
        mid = 0.5*(lo+hi);
        if ma_eps_fun(mid) > eps_target
            lo = mid;   % need more noise
        else
            hi = mid;   % too much noise
        end
        if (hi-lo) < 1e-9*max(1,hi)
            break;
        end
    end
    sigma = 0.5*(lo+hi);
end