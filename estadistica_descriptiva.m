% Estadistica descriptiva de las poblaciones planetarias
% (Calcula media, mediana, desviacion, cuartiles, etc.)
clear; clc; close all;

%% 1. Configuracion
% Cambiar estas variables segun la poblacion que se vaya a analizar

archivo = '5.1 not_so_hot_jupiters_outliers.xlsx';
nombre_pob = "Not-so-hot Jupiters (con outliers)";
rango_pob = "R >= 6.0; P > 10 dias";
sufijo = "not_so_hot_jupiters_outliers";

% Otros ejemplos (descomentar segun sea necesario):
% archivo = '1.rocky_planets.xlsx';
% nombre_pob = "Rocky planets"; rango_pob = "0 <= R < 2.0"; sufijo = "rocky";

fprintf('\nLeyendo archivo: %s...\n', archivo);

%% 2. Cargar y limpiar datos
data = readtable(archivo, 'VariableNamingRule', 'preserve');
data.Properties.VariableNames = matlab.lang.makeValidName(data.Properties.VariableNames);

if ~ismember('radius_earth', data.Properties.VariableNames)
    error('Ojo: no se encontro la columna radius_earth en el archivo.');
end

R = data.radius_earth;

% Si la columna se leyo como texto, convertir a numero
if ~isnumeric(R)
    R = string(R);
    R = strrep(R, ',', '.');
    R = str2double(R);
end

% Filtrar NaNs y valores invalidos
R = R(~isnan(R) & isfinite(R) & R > 0);

if isempty(R)
    error('La poblacion no tiene datos validos.');
end

%% 3. Calcular estadisticos
N = length(R);

min_R = min(R);
max_R = max(R);
media = mean(R);
mediana = median(R);
desv_std = std(R);
varianza = var(R);

q1 = prctile(R, 25);
q3 = prctile(R, 75);
iqr = q3 - q1;
rango = max_R - min_R;

% Numero de bins con la regla de la raiz de N
num_bins = round(sqrt(N));
if num_bins < 3
    num_bins = 3;
end

[cuentas, bordes] = histcounts(R, num_bins);
[frec_max, idx_max] = max(cuentas);
intervalo_frecuente = sprintf('[%.3f , %.3f)', bordes(idx_max), bordes(idx_max+1));

%% 4. Crear tabla y mostrar
tabla = table( ...
    string(nombre_pob), string(rango_pob), N, ...
    min_R, max_R, media, mediana, desv_std, varianza, ...
    q1, q3, iqr, rango, num_bins, string(intervalo_frecuente), frec_max, ...
    'VariableNames', { ...
        'Poblacion', 'Intervalo', 'N', ...
        'Minimo', 'Maximo', 'Media', 'Mediana', 'DesvStd', 'Varianza', ...
        'Q1', 'Q3', 'IQR', 'Rango', 'NumBins', 'IntervaloMasFrecuente', 'FrecMaxima'});

fprintf('\n--- Resultados ---\n');
disp(tabla);

%% 5. Guardar en Excel
nombre_salida = sprintf('resumen_estadistico_%s.xlsx', sufijo);
writetable(tabla, nombre_salida);
fprintf('Tabla guardada en: %s\n', nombre_salida);