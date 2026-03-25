function [concordance, minusconcordance, p_values, HR, log_HR, Conf_interv, Conf_interv_log] = Analiza(path, variable)


data = path;
% Wybór M1 lub M2
if strcmp(variable, 'struct_ness_1')
    struct_ness = cell2mat(data.struct_ness_1);
elseif strcmp(variable, 'struct_ness_2')
    struct_ness = cell2mat(data.struct_ness_2);
end

wiek = data.age;
wiek = cell2mat(wiek);
nan_indices = isnan(wiek);
mean_age = mean(wiek, 'omitnan');
wiek(nan_indices) = mean_age;

X = [struct_ness, wiek];

OS = 1 - cell2mat(data.os);
OS_time = cell2mat(data.os_time);
Censored = true(size(OS));
nan_indices = isnan(OS);
Censored(~nan_indices) = logical(OS(~nan_indices));


Time = double(OS_time);


[b, logl, H, stats] = coxphfit(X, Time, 'Censoring', Censored,'Baseline',0);

% Hazard ratio
HR = exp(b(1,:));
log_HR = log(HR);

% Granice ufności w logarytmie
se = stats.se;
log_hr = stats.beta;
log_ci_lower = log_hr - se * 1.96;
log_ci_upper = log_hr + se * 1.96;
Conf_interva_log = [log_ci_lower, log_ci_upper];
Conf_interv_log = Conf_interva_log(1,:);

% Przekształcenie logarytmiczne na liniową skalę
ci_lower = exp(log_ci_lower);
ci_upper = exp(log_ci_upper);
Conf_interva = [ci_lower, ci_upper];
Conf_interv = Conf_interva(1,:);

% Harrell's c-index
concordance = concordanceIndex(Time, Censored, struct_ness);
minusconcordance = concordanceIndex(Time, Censored, -(struct_ness));

% p-values
p_values = stats.p(1,:);

end
