clear; close all; clc;
% Script para ajustar PDFs a planetas rocosos (R < 2 R_earth)
clear; clc; close all;

% Cargar base de datos
data = readtable('1.rocky_planets.xlsx');
R = data.radius_earth;

% Limpiar NaNs, infinitos y valores negativos o cero
R = R(~isnan(R) & isfinite(R) & R > 0);

% Filtro de la poblacion: 0 <= R < 2.0
R_pob = R(R >= 0.0 & R < 2.0);
N_pob = length(R_pob);

if N_pob == 0
    error('No hay datos en este rango, revisar la tabla.');
end

% Bins usando la regla de la raiz de N
numBins = round(sqrt(N_pob));
if numBins < 3
    numBins = 3;
end
edges = linspace(min(R_pob), max(R_pob), numBins + 1);

fprintf('Poblacion: Rocky planets. N = %d. Bins = %d\n', N_pob, numBins);

% Ejecutar ajustes
fila_weibull = ajustar_weibull(R_pob, edges, numBins);
fila_beta = ajustar_beta(R_pob, numBins);
fila_lognormal = ajustar_lognormal(R_pob, edges, numBins);
fila_gamma = ajustar_gamma(R_pob, edges, numBins);
fila_pareto = ajustar_pareto(R_pob, edges, numBins);
fila_powerlaw = ajustar_powerlaw(R_pob, edges, numBins);

% Juntar todo en una tabla
tabla = [fila_weibull; fila_beta; fila_lognormal; fila_gamma; fila_pareto; fila_powerlaw];
disp('Resultados finales:');
disp(tabla);

% Guardar en Excel
nombre_excel = 'resultados_pdfs_AD_AIC_BIC_rocky_planets.xlsx';
writetable(tabla, nombre_excel);
fprintf('Listo! Tabla guardada en %s\n', nombre_excel);


%% Funciones para cada ajuste

function fila = ajustar_weibull(R_pob, edges, numBins)
    N = length(R_pob);
    pd = fitdist(R_pob, 'Weibull');
    
    A = pd.A; B = pd.B;
    fprintf('\n--- Weibull --- \nA=%.4f, B=%.4f\n', A, B);

    [h, p, D] = kstest(R_pob, 'CDF', pd);
    dec_KS = decision(h);

    obs = histcounts(R_pob, edges);
    esp = N * diff(cdf(pd, edges));
    [chi2, gl, p_chi2, dec_chi2] = calc_chi2(obs, esp, 2);

    F_func = @(x) cdf(pd, x);
    AD = calc_AD(R_pob, F_func);
    
    % Intentar sacar p-value de AD (a veces falla en MATLAB)
    try
        [h_ad, p_ad, ad_stat, cv_ad] = adtest(R_pob, 'Distribution', pd);
        dec_AD = decision(h_ad);
    catch
        p_ad = NaN; cv_ad = NaN; dec_AD = "Solo relativo";
    end

    [logL, AIC, BIC] = calc_AIC(R_pob, @(x) pdf(pd, x), 2);
    
    print_res("WEIBULL", h, p, D, dec_KS, chi2, gl, p_chi2, dec_chi2, AD, p_ad, dec_AD, logL, AIC, BIC);

    % Graficas
    x = linspace(min(R_pob), max(R_pob), 1000);
    plot_pdf(R_pob, edges, x, pdf(pd, x), 'Weibull', sprintf('N=%d\nA=%.2f B=%.2f', N, A, B));
    plot_cdf(R_pob, x, cdf(pd, x), 'Weibull', sprintf('KS = %.3f\np = %.3f', D, p));

    fila = crear_fila("Rocky planets", "Weibull", "A_escala", A, "B_forma", B, ...
        D, p, dec_KS, chi2, gl, p_chi2, dec_chi2, AD, p_ad, cv_ad, dec_AD, logL, AIC, BIC);
end

function fila = ajustar_beta(R_pob, numBins)
    N = length(R_pob);
    Rmin = min(R_pob); Rmax = max(R_pob);
    
    % Normalizar a [0,1] para la Beta
    R_norm = (R_pob - Rmin) ./ (Rmax - Rmin);
    R_norm(R_norm <= 0) = eps; 
    R_norm(R_norm >= 1) = 1 - eps;

    pd = fitdist(R_norm, 'Beta');
    alpha = pd.a; beta = pd.b;
    fprintf('\n--- Beta --- \nAlpha=%.4f, Beta=%.4f\n', alpha, beta);

    [h, p, D] = kstest(R_norm, 'CDF', pd);
    dec_KS = decision(h);

    edges_norm = linspace(0, 1, numBins + 1);
    obs = histcounts(R_norm, edges_norm);
    esp = N * diff(cdf(pd, edges_norm));
    [chi2, gl, p_chi2, dec_chi2] = calc_chi2(obs, esp, 2);

    F_func = @(x) cdf(pd, x);
    AD = calc_AD(R_norm, F_func);
    
    try
        [h_ad, p_ad, ad_stat, cv_ad] = adtest(R_norm, 'Distribution', pd);
        dec_AD = decision(h_ad);
    catch
        p_ad = NaN; cv_ad = NaN; dec_AD = "Solo relativo";
    end

    % AIC en escala original
    pdf_orig = @(x) pdf(pd, (x - Rmin) ./ (Rmax - Rmin)) ./ (Rmax - Rmin);
    [logL, AIC, BIC] = calc_AIC(R_pob, pdf_orig, 2);
    
    print_res("BETA", h, p, D, dec_KS, chi2, gl, p_chi2, dec_chi2, AD, p_ad, dec_AD, logL, AIC, BIC);

    % Graficas
    x_norm = linspace(0.001, 0.999, 1000);
    x = Rmin + x_norm .* (Rmax - Rmin);
    y = pdf(pd, x_norm) ./ (Rmax - Rmin);
    edges_orig = linspace(Rmin, Rmax, numBins + 1);
    
    plot_pdf(R_pob, edges_orig, x, y, 'Beta', sprintf('N=%d\nalpha=%.2f beta=%.2f', N, alpha, beta));
    plot_cdf(R_pob, x, cdf(pd, x_norm), 'Beta', sprintf('KS = %.3f\np = %.3f', D, p));

    fila = crear_fila("Rocky planets", "Beta", "alpha", alpha, "beta", beta, ...
        D, p, dec_KS, chi2, gl, p_chi2, dec_chi2, AD, p_ad, cv_ad, dec_AD, logL, AIC, BIC);
end

% (Aqui irian las funciones para Lognormal, Gamma, Pareto y Power Law, con el mismo estilo)

%% Funciones auxiliares

function dec = decision(h)
    if h == 0
        dec = "No se rechaza H0";
    else
        dec = "Se rechaza H0";
    end
end

function [chi2, gl, p, dec] = calc_chi2(obs, esp, n_params)
    % Filtrar donde lo esperado sea >= 5 (requisito de chi2)
    validos = isfinite(esp) & esp >= 5;
    obs = obs(validos); esp = esp(validos);
    
    gl = length(obs) - 1 - n_params;
    
    if gl <= 0 || isempty(obs)
        chi2 = NaN; p = NaN; dec = "No calculable";
        return;
    end
    
    chi2 = sum((obs - esp).^2 ./ esp);
    p = 1 - chi2cdf(chi2, gl);
    dec = decision(p < 0.05);
end

function AD = calc_AD(datos, F_func)
    datos = sort(datos(:));
    n = length(datos);
    F = F_func(datos);
    
    % Evitar problemas con log(0)
    F(F <= 0) = eps; F(F >= 1) = 1 - eps;
    
    i = (1:n)';
    AD = -n - mean((2*i - 1) .* (log(F) + log(1 - flipud(F))));
end

function [logL, AIC, BIC] = calc_AIC(datos, pdf_func, k)
    n = length(datos);
    val = pdf_func(datos);
    val(val <= 0 | ~isfinite(val)) = realmin; % Evitar log(0)
    
    logL = sum(log(val));
    AIC = 2*k - 2*logL;
    BIC = k*log(n) - 2*logL;
end

function print_res(nombre, h, p, D, dec_KS, chi2, gl, p_chi2, dec_chi2, AD, p_ad, dec_AD, logL, AIC, BIC)
    fprintf('\nPruebas para %s:\n', nombre);
    fprintf('KS: D=%.4f, p=%.4f -> %s\n', D, p, dec_KS);
    fprintf('Chi2: X2=%.2f, gl=%d, p=%.4f -> %s\n', chi2, gl, p_chi2, dec_chi2);
    fprintf('AD: stat=%.4f, p=%.4f -> %s\n', AD, p_ad, dec_AD);
    fprintf('AIC: %.2f, BIC: %.2f (LogL: %.2f)\n', AIC, BIC, logL);
end

function plot_pdf(R, edges, x, y, nombre, caja)
    figure('Color', 'w');
    histogram(R, 'BinEdges', edges, 'Normalization', 'pdf', 'FaceAlpha', 0.5);
    hold on; plot(x, y, 'LineWidth', 2);
    grid on; xlabel('Radio (R_\oplus)'); ylabel('Densidad');
    title(['Ajuste ' nombre]);
    legend('Datos', ['PDF ' nombre]);
    annotation('textbox', [0.15 0.6 0.2 0.2], 'String', caja, 'BackgroundColor', 'w', 'FitBoxToText', 'on');
end

function plot_cdf(R, x, y, nombre, caja)
    figure('Color', 'w');
    [f, xi] = ecdf(R);
    stairs(xi, f, 'LineWidth', 2); hold on;
    plot(x, y, 'LineWidth', 2);
    grid on; xlabel('Radio (R_\oplus)'); ylabel('Prob. Acumulada');
    title(['CDF ' nombre]);
    legend('Empirica', 'Teorica');
    annotation('textbox', [0.6 0.2 0.2 0.2], 'String', caja, 'BackgroundColor', 'w', 'FitBoxToText', 'on');
end

function fila = crear_fila(pob, dist, p1_nom, p1_val, p2_nom, p2_val, KS_stat, KS_p, KS_dec, Chi_stat, Chi_gl, Chi_p, Chi_dec, AD_stat, AD_p, AD_cv, AD_dec, logL, AIC, BIC)
    fila = table(string(pob), string(dist), string(p1_nom), p1_val, string(p2_nom), p2_val, ...
        KS_stat, KS_p, string(KS_dec), Chi_stat, Chi_gl, Chi_p, string(Chi_dec), ...
        AD_stat, AD_p, AD_cv, string(AD_dec), logL, AIC, BIC, ...
        'VariableNames', {'Poblacion','Distribucion','Parametro1','ValorParametro1','Parametro2','ValorParametro2',...
        'KS_Estadistico','KS_pvalue','KS_Decision','Chi2_Estadistico','Chi2_gl','Chi2_pvalue','Chi2_Decision',...
        'AD_Estadistico','AD_pvalue','AD_ValorCritico','AD_Decision','LogVerosimilitud','AIC','BIC'});
end