clear all;
clc;
rng(1234);     % Seed to control pseudo-random number generation.

cd '/Users/gabriel/Documents/CATO/2023-2/Laboratorio de Matlab'

format bank;

%% Importing Data

datos = xlsread('Datos_base.xlsx','Base_vis');

cpi_per         = datos(:,1);   % inflacion peru (var interanual)
cpi_usa         = datos(:,2);   % Inflacion US (variacion interanual)
ffr             = datos(:,3);   % Federal funds rate
gdp_per         = datos(:,4);   % PBI peru (niveles trimestrales)
ti_per          = datos(:,5);   % Terminos de intercambio Peru (niveles)
pm              = datos(:,6);   % Tasa de interes interbancaria promedio
dates           = datos(:,7);   % Fechas

y1      = [cpi_per cpi_usa ffr gdp_per ti_per pm]; % Abreviation of variables. Concatemos las series
names1  = ['     Inflation Peru     '
           '     Inflation US       ';
           '   Federal funds rate   ';
           '        GDP Peru        ';
           '     Terms of Trade     ';
           '  Interbank rate mean   '];
order   = [3 2 5 4 1 6]; % Put the order of the variables for VAR. The number represents the variable.

%% Comportamiento de las variables

% Grafico de las series
figure('Name','Evolucion de las series')
% Definir el primer grafico
subplot(3, 2, 1);
plot(dates, gdp_per, 'LineWidth', 1.5);
title('PBI de Peru');
xlabel('Fecha');
ylabel('PBI interanual');

% Definir el segundo gráfico 
subplot(3, 2, 2);
plot(dates, ti_per, 'LineWidth', 1.5); 
title('Terminos de intercambio');
xlabel('Fecha');
ylabel('TI en niveles');

% Definir el tercer grafico 
subplot(3, 2, 3);
plot(dates, cpi_per, 'LineWidth', 1.5); 
title('Inflacion de Peru');
xlabel('Fecha');
ylabel('Inflacion interanual');

% Definir el cuarto gráfico 
subplot(3, 2, 4);
plot(dates, pm, 'LineWidth', 1.5); 
title('Tasa de interes interbancaria');
xlabel('Fecha');
ylabel('Tasas en niveles');

subplot(3, 2, 5);
plot(dates, cpi_usa, 'LineWidth', 1.5); 
title('Inflacion de US');
xlabel('Fecha');
ylabel('Inflacion interanual');

subplot(3, 2, 6);
plot(dates, ffr, 'LineWidth', 1.5);
title('Tasa de los fondos federales de US');
xlabel('Fecha');
ylabel('tasa interanual');

%% Importamos series ya estacionarias

clear all;
clc;
rng(1234);

cd '/Users/gabriel/Documents/CATO/2023-2/Laboratorio de Matlab'

datos = xlsread('Datos_base.xlsx','Base_est');

cpi_per         = datos(:,1);   % inflacion peru (var interanual)
cpi_usa         = datos(:,2);   % Inflacion US (variacion interanual)
ffr             = datos(:,3);   % Federal funds rate
gdp_per         = datos(:,4);   % PBI peru (niveles trimestrales)
ti_per          = datos(:,5);   % Terminos de intercambio Peru (niveles)
pm              = datos(:,6);   % Tasa de interes interbancaria promedio
dates           = datos(:,7);   % Fechas

y1      = [cpi_per cpi_usa ffr gdp_per ti_per pm]; % Abreviation of variables. Concatemos las series
names1  = ['   Inflacion de Peru    '
           '    Inflacion de US     ';
           'Tasa de fondos federales';
           '      PBI de Peru       ';
           '         T.I.           ';
           '   Tasa interbancaria   '];
order   = [3 2 5 4 1 6]; % Put the order of the variables for VAR. The number represents the variable.

start_year      = 1996;             % Initial year
start_quarter   = 1;                % Initial quarter
end_year        = 2019;             % Last year
end_quarter     = 4;                % Last quarter

y1 = y1(find(dates==start_year+start_quarter/4):find(dates==end_year+end_quarter/4),:); 
% creamos secuencia de numeros para ordenar acorde a la fecha, pero como ya lo puse, no sirve

lag             = 4;                % Lags of VAR model
c               = 1;                % 0: no constant, 1: constant, 2:linear trend
T               = size(y1,1);       % Size of data
n               = size(y1,2);       % Number of variables
k               = size(y1,2)*lag+c; % Number of regressors
S               = 20;    

%% Order of variables in VAR model 
% esta parte tambien relacionada a SVAR
y       = [];
names   = ['                        '];

for i = 1:n
    
    y           = [y y1(:,order(i))];
    names(i,:)  = string(names1(order(i),:));

end

%% VAR Estimation

Y = y(lag+1:T,:); % Dropping first "lags" observation of y vector.
% en modelos VAR por cada rezago que pongo pierdo obs, por lo que las obs
% empiezan a partir del numero de rezagos
X = [];

switch(c) % determinamos el componente dterministico
    case 0        % No deterministic components.
        X = [];    
    case 1        % Include a constant.
        X = ones(T-lag,1);        
    case 2        % Include a linear trend.
        for t = 1:T-lag            
            X(t,:) = [1 t];            
        end
end

% Lets generate X matrix for lags of vector y. % creamos los regresores de
% la matriz como rezagos

for i = 1:lag
    x{i} = y(lag+1-i:T-i,:);
    X = [X x{i}];
end

B_ols = (X'*X)^(-1)*(X'*Y); % OLS estimator of the model.

e = Y-X*B_ols; % Reduced errors.

sigmae = e'*e/(T-lag-k); % Variance matrix of reduced errors. % formula de la varianza de la regresion SUR

%% Getting companion matrix F

AA = B_ols(c+1:end,:); % Coefficients of lags. Excluding coefficients of deterministic components.

F=[]; % Definimos matriz compañera

for i = 1:lag
    A{i}=AA(1+n*(i-1):n*i,:); % Generating A1,...,Ap matrices.
    F = [F A{i}];  % Generating first matrix row of companion matrix F.
end
 
F = [F;eye(n*(lag-1)),zeros(n*(lag-1),n)]; % Generating others matrix rows of companion matrix F.
J = [eye(n) zeros(n,n*(lag-1))]; % definimos matriz seleccionadora que permite volver del VAR(1) al VAR(p)

%% Analizamos estacionariedad

lambda = eig(F); % Calculating eigenvalues of companion matrix F
modulo = (real(lambda).^2 + imag(lambda).^2 ).^(1/2); % Calculating modulus of each eigenvalue.

disp('Modulus of eigenvalues');
disp(modulo);
if max(abs(modulo))<1
    disp('The model satisfy stationary conditions.');
elseif max(abs(modulo))>=1
    disp('The model does not satisfy stationary conditions.');
end

% Graph of unit circle and inverse roots

tempx = [-1:0.01:1];
tempy = (1-tempx.^2).^(1/2);
rootx = real(lambda)';
rooty = imag(lambda)';

figure('Name','Unit Circle')
h = scatter(rootx,rooty);
set(h,'MarkerEdgeColor','Blue')
set(h,'MarkerFaceColor','Blue')
hold on
h = plot(tempx,tempy,tempx,-tempy);
set(h(1),'Color','Black')
set(h(2),'Color','Black')
xlim([(min(tempx)-0.5) (max(tempx)+0.5)])
ylim([(min(tempy)-1.2) (max(tempy)+0.2)])
grid on
hold off
title('Inverse Roots of Characteristic Polynominal')

clear tempx tempy rootx rooty;

%% ACF graphs

figure('Name','ACF ');
subplot(3, 2, 1);
autocorr(gdp_per);
title('ACF del PBI de Peru');

subplot(3, 2, 2);
autocorr(ti_per);
title('ACF de los Terminos de intercambio');

subplot(3, 2, 3);
autocorr(cpi_per);
title('ACF de la Inflacion de Peru');

subplot(3, 2, 4);
autocorr(pm);
title('ACF de la Tasa de interes interbancaria');

subplot(3, 2, 5);
autocorr(cpi_usa);
title('ACF de la Inflacion de US');

subplot(3, 2, 6);
autocorr(ffr)
title('ACF de la Tasa de los fondos federales');

%% Shock of Structural Form - Cholesky

% El codigo permite realizar dos identificaciones, como voy a hallar los
% parametros de la forma estructural en base a los de la forma reducida

type_ident = 'bq'; % oir: orthogonal impulse response; bq: blanchard-quah
units_structural_shock = 1; % 0: one unit shock, 1: one standard deviation shock

if strcmp(type_ident,'oir')
    
    P = chol(sigmae)';
    u = inv(P)*e';
    
elseif strcmp(type_ident,'bq')
    
    FF = inv(eye(length(F))-F)';
    Psi1 = J*FF*J';
    Sigma_LR = Psi1*sigmae*Psi1';
    cSigma_LR = chol(Sigma_LR)';
    P = Psi1\cSigma_LR;
    u = inv(P)*e';
    
end

switch(units_structural_shock)
    case 0
        shocks = eye(n).*(diag(P).^(-1)); % one unit shock
    case 1
        shocks = eye(n); % one std deviation shock
end

for i = 1:n
    
    name_sh_struc(i,:) = strcat(' Structural Shock ','',num2str(i));
    
end

IRF_s  = zeros(n,S+1,n);
MSE_s  = zeros(n,S+1);
Ome_s  = zeros(n,S+1,n);
FEVD_s = zeros(n,S+1,n);
omega  = zeros(n,n);

% Generating IRFs and FEVD

for j = 0:S
    
    Phi = J*(F'^j)*J';
    Theta = Phi*P*shocks;
    IRF_s(:,j+1,:) = Theta;
    Ome_s(:,j+1,:) = Theta.^2;
    
    if j == 0
        MSE_s(:,j+1) = diag(Phi*sigmae*Phi');
    elseif j>=1
        MSE_s(:,j+1) = MSE_s(:,j) + diag(Phi*sigmae*Phi');
    end
   
end

Ome_s = cumsum(Ome_s,2);
MSE_s2 = (MSE_s).^(1/2);

for s = 0:S
    
    FEVD_s(:,s+1,:) = Ome_s(:,s+1,:)./MSE_s(:,s+1)*100;
    
end

% Generating Historical Decomposition

Thetas_HD = NaN(n,T-lag,n);

for j = 0:T-lag-1
    
    Phi = J*(F'^j)*J';
    Theta = Phi*P*shocks;
    Thetas_HD(:,j+1,:) = Theta;
   
end

for i = 1:T-lag
    
    y1hat1(i,:) = dot(Thetas_HD(1,1:i,1),u(1,i:-1:1));
    y1hat2(i,:) = dot(Thetas_HD(1,1:i,2),u(2,i:-1:1));
    y1hat3(i,:) = dot(Thetas_HD(1,1:i,3),u(3,i:-1:1));
    y1hat4(i,:) = dot(Thetas_HD(1,1:i,4),u(4,i:-1:1));
    y1hat5(i,:) = dot(Thetas_HD(1,1:i,5),u(5,i:-1:1));
    y1hat6(i,:) = dot(Thetas_HD(1,1:i,6),u(6,i:-1:1));

    y2hat1(i,:) = dot(Thetas_HD(2,1:i,1),u(1,i:-1:1));
    y2hat2(i,:) = dot(Thetas_HD(2,1:i,2),u(2,i:-1:1));
    y2hat3(i,:) = dot(Thetas_HD(2,1:i,3),u(3,i:-1:1));
    y2hat4(i,:) = dot(Thetas_HD(2,1:i,4),u(4,i:-1:1));
    y2hat5(i,:) = dot(Thetas_HD(2,1:i,5),u(5,i:-1:1));
    y2hat6(i,:) = dot(Thetas_HD(2,1:i,6),u(6,i:-1:1));

    y3hat1(i,:) = dot(Thetas_HD(3,1:i,1),u(1,i:-1:1));
    y3hat2(i,:) = dot(Thetas_HD(3,1:i,2),u(2,i:-1:1));
    y3hat3(i,:) = dot(Thetas_HD(3,1:i,3),u(3,i:-1:1));
    y3hat4(i,:) = dot(Thetas_HD(3,1:i,4),u(4,i:-1:1));
    y3hat5(i,:) = dot(Thetas_HD(3,1:i,5),u(5,i:-1:1));
    y3hat6(i,:) = dot(Thetas_HD(3,1:i,6),u(6,i:-1:1));

    y4hat1(i,:) = dot(Thetas_HD(4,1:i,1),u(1,i:-1:1));
    y4hat2(i,:) = dot(Thetas_HD(4,1:i,2),u(2,i:-1:1));
    y4hat3(i,:) = dot(Thetas_HD(4,1:i,3),u(3,i:-1:1));
    y4hat4(i,:) = dot(Thetas_HD(4,1:i,4),u(4,i:-1:1));
    y4hat5(i,:) = dot(Thetas_HD(4,1:i,5),u(5,i:-1:1));
    y4hat6(i,:) = dot(Thetas_HD(4,1:i,6),u(6,i:-1:1));

    y5hat1(i,:) = dot(Thetas_HD(5,1:i,1),u(1,i:-1:1));
    y5hat2(i,:) = dot(Thetas_HD(5,1:i,2),u(2,i:-1:1));
    y5hat3(i,:) = dot(Thetas_HD(5,1:i,3),u(3,i:-1:1));
    y5hat4(i,:) = dot(Thetas_HD(5,1:i,4),u(4,i:-1:1));
    y5hat5(i,:) = dot(Thetas_HD(5,1:i,5),u(5,i:-1:1));
    y5hat6(i,:) = dot(Thetas_HD(5,1:i,6),u(6,i:-1:1));

    y6hat1(i,:) = dot(Thetas_HD(6,1:i,1),u(1,i:-1:1));
    y6hat2(i,:) = dot(Thetas_HD(6,1:i,2),u(2,i:-1:1));
    y6hat3(i,:) = dot(Thetas_HD(6,1:i,3),u(3,i:-1:1));
    y6hat4(i,:) = dot(Thetas_HD(6,1:i,4),u(4,i:-1:1));
    y6hat5(i,:) = dot(Thetas_HD(6,1:i,5),u(5,i:-1:1));
    y6hat6(i,:) = dot(Thetas_HD(6,1:i,6),u(6,i:-1:1));

end

y1hat  = [y1hat1 y1hat2 y1hat3 y1hat4 y1hat5 y1hat6];
y2hat  = [y2hat1 y2hat2 y2hat3 y2hat4 y2hat5 y2hat6];
y3hat  = [y3hat1 y3hat2 y3hat3 y3hat4 y3hat5 y3hat6];
y4hat  = [y4hat1 y4hat2 y4hat3 y4hat4 y4hat5 y4hat6];
y5hat  = [y5hat1 y5hat2 y5hat3 y5hat4 y5hat5 y5hat6];
y6hat  = [y6hat1 y6hat2 y6hat3 y6hat4 y6hat5 y6hat6];

yhat     = [y1hat y2hat y3hat y4hat y5hat y6hat];
ysum     = [sum(y1hat,2) sum(y2hat,2) sum(y3hat,2) sum(y4hat,2) sum(y5hat,2) sum(y6hat,2)]; 
ydemean  = [detrend(Y(:,1),c-1) detrend(Y(:,2),c-1) detrend(Y(:,3),c-1) detrend(Y(:,4),c-1) detrend(Y(:,5),c-1) detrend(Y(:,6),c-1)];

clear y1hat y1hat1 y1hat2 y1hat3 y1hat4 y1hat5 y1hat6...
      y2hat y2hat1 y2hat2 y2hat3 y2hat4 y2hat5 y2hat6...
      y3hat y3hat1 y3hat2 y3hat3 y3hat4 y3hat5 y3hat6...
      y4hat y4hat1 y4hat2 y4hat3 y4hat4 y4hat5 y4hat6...
      y5hat y5hat1 y5hat2 y5hat3 y5hat4 y5hat5 y5hat6...
      y6hat y6hat1 y6hat2 y6hat3 y6hat4 y6hat5 y6hat6

%% Graph of IRF of Structural Shocks

zeroline = zeros(1,S+1);
irf_axis = 0:S;

switch(units_structural_shock)
    case 0
        figure('Name','IRF de los choques estructurales - One Unit Shock')
    case 1
        figure('Name','IRF de los choques estructurales - One Std. Dev. Shock')
end

% cont = 1;
% 
% for var = 1:n
% 
%     for sh = 1:n
% 
%         subplot(n,n,cont)
%         h = plot(irf_axis,zeroline,irf_axis,IRF_s(var,:,sh));
%         set(h(1),'Color','Black','LineWidth',0.5)
%         set(h(2),'Color','Black','LineWidth',1.5)
%         title(strcat('Response of ',' ',names(var,:),' to ',name_sh_struc(sh,:)))
%         axis tight;grid on;
%         cont = cont+1;
% 
%     end
% 
% end

% Respuesta de los Terminos de intercambio a un shock de pm de US
subplot(4, 2, 1);
plot(irf_axis, zeroline, irf_axis, IRF_s(3, :, 1), 'LineWidth',1.5);
title('Respuesta de los T.I. a un choque de p.m. de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta de los Terminos de intercambio a un shock de costos de US
subplot(4, 2, 2);
plot(irf_axis, zeroline, irf_axis, IRF_s(3, :, 2), 'LineWidth',1.5);
title('Respuesta de los T.I. a un choque inflacionario de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta del PBI peruano a un shock de pm de US
subplot(4, 2, 3);
plot(irf_axis, zeroline, irf_axis, IRF_s(4, :, 1), 'LineWidth',1.5);
title('Respuesta del PBI a un choque de p.m. de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta del PBI peruano a un shock de costos de US
subplot(4, 2, 4);
plot(irf_axis, zeroline, irf_axis, IRF_s(4, :, 2), 'LineWidth',1.5);
title('Respuesta del PBI a un choque inflacionario de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta de la inflacion peruana a un shock de pm de US
subplot(4, 2, 5);
plot(irf_axis, zeroline, irf_axis, IRF_s(5, :, 1), 'LineWidth',1.5);
title('Respuesta de la inflacion peruana a un choque de p.m. de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta de la inflacion peruana a un shock de costos de US
subplot(4, 2, 6);
plot(irf_axis, zeroline, irf_axis, IRF_s(5, :, 2), 'LineWidth',1.5);
title('Respuesta de la inflacion peruana a un choque inflacionario de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta de la PM peruana a un shock de pm de US
subplot(4, 2, 7);
plot(irf_axis, zeroline, irf_axis, IRF_s(6, :, 1), 'LineWidth',1.5);
title('Respuesta de la p.m. del BCR a un choque de p.m. de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

% Respuesta de la PM peruana a un shock de costos de US
subplot(4, 2, 8);
plot(irf_axis, zeroline, irf_axis, IRF_s(6, :, 2), 'LineWidth',1.5);
title('Respuesta de la p.m. del BCR a un choque inflacionario de US');
xlabel('Time');
ylabel('Response');
axis tight;
grid on;

%% Graph of FEVD

% colors = [
%     0.114 0.388 0.639; % Azul ligeramente oscuro
%     0.855 0.439 0.125; % Naranja ligeramente oscuro
%     0.647 0.165 0.165; % Rojo ligeramente oscuro
%     0.929 0.694 0.125; % Amarillo ligeramente oscuro
%     0.125 0.694 0.929; % Celeste ligeramente oscuro
%     0.071 0.392 0.071  % Verde ligeramente más oscuro
% ];
% 
% colors = [
%     0.647 0.165 0.165; % Rojo
%     0.855 0.439 0.125; % Naranja claro
%     0.85 0.85 0.85;       % Gris claro
%     0.70 0.70 0.70;       % Gris medio
%     0.5 0.5 0.5;       % Gris oscuro
%     0.35 0.35 0.35        % Gris más oscuro
% ];
% 
% switch(units_structural_shock)
%     case 0
%         figure('Name','FEVD de los choques estructuralees - One Unit Shock')
%     case 1
%         figure('Name','FEVD de los choques estructuralees - One Std. Dev. Shock')
% end
% 
%     subplot(1,4,1)
%     indx = zeros(S+1,n);
%     indx(:,:) = FEVD_s(3,:,:);
%     h = bar(indx,'stacked');
%     for i = 1:n
%         set(h(i),'FaceColor',colors(i,:))
%     end
%     axis tight
%     title(strcat('Descomposicion de',' ',names(3,:)))
%     legend(h,name_sh_struc,'Location','southoutside')
% 
%     subplot(1,4,2)
%     indx = zeros(S+1,n);
%     indx(:,:) = FEVD_s(4,:,:);
%     h = bar(indx,'stacked');
%     for i = 1:n
%         set(h(i),'FaceColor',colors(i,:))
%     end
%     axis tight
%     title(strcat('Descomposicion de',' ',names(4,:)))
%     legend(h,name_sh_struc,'Location','southoutside')
% 
%     subplot(1,4,3)
%     indx = zeros(S+1,n);
%     indx(:,:) = FEVD_s(5,:,:);
%     h = bar(indx,'stacked');
%     for i = 1:n
%         set(h(i),'FaceColor',colors(i,:))
%     end
%     axis tight
%     title(strcat('Descompisicion de',' ',names(5,:)))
%     legend(h,name_sh_struc,'Location','southoutside')
% 
%     subplot(1,4,4)
%     indx = zeros(S+1,n);
%     indx(:,:) = FEVD_s(6,:,:);
%     h = bar(indx,'stacked');
%     for i = 1:n
%         set(h(i),'FaceColor',colors(i,:))
%     end
%     axis tight
%     title(strcat('Descomposicion de',' ',names(6,:)))
%     legend(h,name_sh_struc,'Location','southoutside')
% 

colors = [
    0.647 0.165 0.165; % Rojo
    0.855 0.439 0.125; % Naranja claro
    0.85 0.85 0.85;    % Gris claro
    0.70 0.70 0.70;    % Gris medio
    0.5 0.5 0.5;       % Gris oscuro
    0.35 0.35 0.35     % Gris más oscuro
];

switch(units_structural_shock)
    case 0
        figure('Name','FEVD de los choques estructurales - One Unit Shock')
    case 1
        figure('Name','FEVD de los choques estructurales - One Std. Dev. Shock')
end

for i = 1:4
    subplot(1,4,i)
    indx = zeros(S+1,n);
    indx(:,:) = FEVD_s(i+2,:,:); 
    h = bar(indx,'stacked');
    for j = 1:n
        set(h(j),'FaceColor',colors(j,:))
    end
    axis tight
    title(strcat('Descomposicion de',' ',names(i+2,:))) 
end

%% Grafico del performance de la descomposicion historica de las variables de interes 

zeroline = zeros(1,T-lag);

switch(units_structural_shock)
    case 0
        figure('Name','Performance de la descomposicion historica - One Unit Shock')
    case 1
        figure('Name','Performance de la descomposicion historica - One Std. Dev. Shock')
end

subplot(4, 1, 1);
h1 = plot(dates(1 + lag:end, 1), zeroline, 'k-');
hold on;
h2 = plot(dates(1 + lag:end, 1), ysum(:, 3), 'k-', ...
    dates(1 + lag:end, 1), ydemean(:, 3), 'k-.');
set(h1(1), 'LineWidth', 0.5);
set(h2(1), 'LineWidth', 1.5);
set(h2(2), 'LineWidth', 1.5);
title(strcat(names(3, :)));
axis tight;
grid on;
hold off;
% legend('Sum of Contribution of Shocks', 'Demeaned Series', ...
%     'Orientation', 'Horizontal', 'Location', 'southoutside');

%(Variable 2)
subplot(4, 1, 2);
h1 = plot(dates(1 + lag:end, 1), zeroline, 'k-');
hold on;
h2 = plot(dates(1 + lag:end, 1), ysum(:, 4), 'k-', ...
    dates(1 + lag:end, 1), ydemean(:, 4), 'k-.');
set(h1(1), 'LineWidth', 0.5);
set(h2(1), 'LineWidth', 1.5);
set(h2(2), 'LineWidth', 1.5);
title(strcat(names(4, :)));
axis tight;
grid on;
hold off;
% legend('Sum of Contribution of Shocks', 'Demeaned Series', ...
%     'Orientation', 'Horizontal', 'Location', 'southoutside');

% Variable 3 
subplot(4, 1, 3);
h1 = plot(dates(1 + lag:end, 1), zeroline, 'k-');
hold on;
h2 = plot(dates(1 + lag:end, 1), ysum(:, 5), 'k-', ...
    dates(1 + lag:end, 1), ydemean(:, 5), 'k-.');
set(h1(1), 'LineWidth', 0.5);
set(h2(1), 'LineWidth', 1.5);
set(h2(2), 'LineWidth', 1.5);
title(strcat(names(5, :)));
axis tight;
grid on;
hold off;
% legend('Sum of Contribution of Shocks', 'Demeaned Series', ...
%     'Orientation', 'Horizontal', 'Location', 'southoutside');

% Variable 4
subplot(4, 1, 4);
h1 = plot(dates(1 + lag:end, 1), zeroline, 'k-');
hold on;
h2 = plot(dates(1 + lag:end, 1), ysum(:, 6), 'k-', ...
    dates(1 + lag:end, 1), ydemean(:, 6), 'k-.');
set(h1(1), 'LineWidth', 0.5);
set(h2(1), 'LineWidth', 1.5);
set(h2(2), 'LineWidth', 1.5);
title(strcat(names(6, :)));
axis tight;
grid on;
hold off;
% legend('Sum of Contribution of Shocks', 'Demeaned Series', ...
%     'Orientation', 'Horizontal', 'Location', 'southoutside');

%% Graph of Historical Decomposition of the variable

zeroline = zeros(1,T-lag);

colors = [
    0.647 0.165 0.165; % Rojo
    0.855 0.439 0.125; % Naranja claro
    0.85 0.85 0.85;       % Gris claro
    0.70 0.70 0.70;       % Gris medio
    0.5 0.5 0.5;       % Gris oscuro
    0.35 0.35 0.35        % Gris más oscuro
];

% colors = [
%     0.647 0.165 0.165; % Rojo
%     1.0 0.8 0.2; % Amarillo claro
%     0.85 0.85 0.85;       % Gris claro
%     0.70 0.70 0.70;       % Gris medio
%     0.5 0.5 0.5;       % Gris oscuro
%     0.35 0.35 0.35        % Gris más oscuro
% ];

switch(units_structural_shock)
    case 0
        figure('Name','Descomposicion historica - One Unit Shock')
    case 1
        figure('Name','Descomposicion historica - One Std. Dev. Shock')
end

    indx_pos = yhat(:,1+(3-1)*n:3*n);
    indx_pos(yhat(:,1+(3-1)*n:3*n)<0) = 0;
    indx_neg = yhat(:,1+(3-1)*n:3*n);
    indx_neg(yhat(:,1+(3-1)*n:3*n)>0) = 0;

    subplot(4,1,1)
    h1 = plot(dates(1+lag:end,1),zeroline,'k-');
    axis tight; grid on;
    hold on;
    h2 = bar(dates(1+lag:end,1),indx_pos,'stacked');
    h3 = bar(dates(1+lag:end,1),indx_neg,'stacked');
    for i = 1:n
        set(h2(i),'FaceColor',colors(i,:))
        set(h3(i),'FaceColor',colors(i,:))
    end
    title(strcat(names(3,:)))
    hold off;

    indx_pos = yhat(:,1+(4-1)*n:4*n);
    indx_pos(yhat(:,1+(4-1)*n:4*n)<0) = 0;
    indx_neg = yhat(:,1+(4-1)*n:4*n);
    indx_neg(yhat(:,1+(4-1)*n:4*n)>0) = 0;

    subplot(4,1,2)
    h1 = plot(dates(1+lag:end,1),zeroline,'k-');
    axis tight; grid on;
    hold on;
    h2 = bar(dates(1+lag:end,1),indx_pos,'stacked');
    h3 = bar(dates(1+lag:end,1),indx_neg,'stacked');
    for i = 1:n
        set(h2(i),'FaceColor',colors(i,:))
        set(h3(i),'FaceColor',colors(i,:))
    end
    title(strcat(names(4,:)))
    hold off;

    indx_pos = yhat(:,1+(5-1)*n:5*n);
    indx_pos(yhat(:,1+(5-1)*n:5*n)<0) = 0;
    indx_neg = yhat(:,1+(5-1)*n:5*n);
    indx_neg(yhat(:,1+(5-1)*n:5*n)>0) = 0;

    subplot(4,1,3)
    h1 = plot(dates(1+lag:end,1),zeroline,'k-');
    axis tight; grid on;
    hold on;
    h2 = bar(dates(1+lag:end,1),indx_pos,'stacked');
    h3 = bar(dates(1+lag:end,1),indx_neg,'stacked');
    for i = 1:n
        set(h2(i),'FaceColor',colors(i,:))
        set(h3(i),'FaceColor',colors(i,:))
    end
    title(strcat(names(5,:)))
    hold off;

    indx_pos = yhat(:,1+(6-1)*n:6*n);
    indx_pos(yhat(:,1+(6-1)*n:6*n)<0) = 0;
    indx_neg = yhat(:,1+(6-1)*n:6*n);
    indx_neg(yhat(:,1+(6-1)*n:6*n)>0) = 0;
        
    subplot(4,1,4)
    h1 = plot(dates(1+lag:end,1),zeroline,'k-');
    axis tight; grid on;
    hold on;
    h2 = bar(dates(1+lag:end,1),indx_pos,'stacked');
    h3 = bar(dates(1+lag:end,1),indx_neg,'stacked');
    for i = 1:n
        set(h2(i),'FaceColor',colors(i,:))
        set(h3(i),'FaceColor',colors(i,:))
    end
    title(strcat(names(6,:)))
    hold off;

legend(h2, name_sh_struc, 'Orientation', 'Horizontal', 'Location', 'southoutside')   


% %% Sign Restrictions
% 
% IRF_sr      = [];                   % Matrix where IRFs will be reshaped.
% ho          = 4;                    % Number of periods for sign restrictions.
% maxdraw     = 1000;                 % Number of IRFs which must satisfy sign restrictions.
% dr          = 0;                    
% count       = 0;                    
% S_IRF_all   = NaN(n,n,S+1,maxdraw); % Matrix where IRFs will be saved.
% S_IRF_all2  = NaN(n,S+1,n,maxdraw); % IRF's reshaped for graphic purposes.
% QQ          = NaN(n,n,maxdraw);     % Matrix where q matrices will be saved.
% u_shocks    = NaN(n,T-lag,maxdraw);
% 
% for i = 1:S+1
% 
%     IRF_sr(:,:,i) = IRF_s(:,i,:);
% 
% end
% 
% while dr<maxdraw
% 
%     M     = 0;
%     it    = 0;
%     count = count+1;
% 
%     while M == 0
% 
%         [q r]   = qr(normrnd(0,1,n,n)); % definimos matriz Q nxn
% 
%         for ii=1:n;
%             if r(ii,ii)<0
%                 q(:,ii)=-q(:,ii);
%             end
%         end
% 
%         for i = 1: ho
% 
%             Q(:,:,i) = IRF_sr(:,:,i)*q;
% 
%         end
% 
%         % Choque de la tasa de fondos federales
%         z1 = [Q(2,1,1:ho)<=0 Q(3,1,1:ho)<=0 Q(5,1,1:ho)<=0];
%         % Q(I,J) EL EFECTO DEL CHOQUE DEL CHOQUE J SOBRE LA VARIABLE I es
%         zz1 = all(z1);    
% 
%         if zz1 == 0
%             mz1 = [Q(2,1,1:ho)>=0 Q(3,1,1:ho)>=0 Q(5,1,1:ho)>=0]; % para verificar con la restriccion inversa
%             mzz1 = all(mz1);
%             if mzz1 == 1 
%                 q(:,1) = -q(:,1); % si cumple con restriccion inversa, se multiplica por -1 y nos quedamos con esa matriz
%             elseif mzz1 == 0
%                 break
%             end
%         end
% 
%         % si cumple con las restricciones se queda con esa matriz, si no
%         % cumple la descarto y continuo a la sigueinte
% 
%         % Choque de inflacion de US
%         z2 = [Q(1,2,1:ho)>=0 Q(3,2,1:ho)>=0 Q(4,2,1:ho)<=0 Q(5,2,1:ho)>=0 Q(6,2,1:ho)>=0];
%         zz2 = all(z2);      
% 
%         if zz2 == 0
%             mz2 = [Q(1,2,1:ho)<=0 Q(3,2,1:ho)<=0 Q(4,2,1:ho)>=0 Q(5,2,1:ho)<=0 Q(6,2,1:ho)<=0];
%             mzz2 = all(mz2);
%             if mzz2 == 1
%                 q(:,2) = -q(:,2);
%             elseif mzz2 == 0
%                 break
%             end
%         end
% 
%         % Choque de terminos de intercambio
%         z3 = [Q(1,3,1:ho)>=0];
%         zz3 = all(z3); 
% 
%         if zz3 == 0
%             mz3 = [Q(1,3,1:ho)<=0];
%             mzz3 = all(mz3);
%             if mzz3 == 1
%                 q(:,3) = -q(:,3);
%             elseif mzz3 == 0
%                 break
%             end
%         end
% 
%         % Choque de demanda de Peru
%         z4 = [Q(2,4,1:ho)>=0];
%         zz4 = all(z4); 
% 
%         if zz4 == 0
%             mz4 = [Q(2,4,1:ho)<=0];
%             mzz4 = all(mz4);
%             if mzz4 == 1
%                 q(:,4) = -q(:,4);
%             elseif mzz4 == 0
%                 break
%             end
%         end
% 
%         % Choque inflacion de Peru 
%         z5 = [Q(2,5,1:ho)<=0 Q(3,5,1:ho)>=0 Q(4,5,1:ho)<=0];
%         zz5 = all(z5); 
% 
%         if zz5 == 0
%             mz5 = [Q(2,5,1:ho)>=0 Q(3,5,1:ho)<=0 Q(4,5,1:ho)>=0];
%             mzz5 = all(mz5);
%             if mzz5 == 1
%                 q(:,5) = -q(:,5);
%             elseif mzz5 == 0
%                 break
%             end
%         end
% 
%         % Choque tasa de interes interbancaria Peru 
%         z6 = [Q(3,6,1:ho)<=0 Q(4,6,1:ho)<=0];
%         zz6 = all(z6); 
% 
%         if zz6 == 0
%             mz6 = [Q(3,6,1:ho)>=0 Q(4,6,1:ho)>=0];
%             mzz6 = all(mz6);
%             if mzz6 == 1
%                 q(:,6) = -q(:,6);
%             elseif mzz6 == 0
%                 break
%             end
%         end
% 
%         M = 1;
%         dr = dr+1;
%         %disp('Found it...')
%         %disp(dr)
%         for i = 1:S+1
%             S_IRF_all(:,:,i,dr) = IRF_sr(:,:,i)*q;
%         end
%         QQ(:,:,dr) = q;
%         V = P*q;
%         u_shocks(:,:,dr) = inv(V)*e';
% 
%     end 
% end
% 
% %disp('Hit ratio...')
% %disp(dr/count*100)
% 
% for i = 1:S+1
% 
%     S_IRF_all2(:,i,:,:) = S_IRF_all(:,:,i,:);
% 
% end
% 
% % Graph of IRFs
% 
% zeroline = zeros(1,S+1);
% irf_axis = 0:S;
% 
% figure('Name','Sign Restrictions IRF of Structural Shocks - One Std. Deviation Shock')
% cont=1;
% 
% for var = 1:n
% 
%     for sh = 1:n
% 
%         subplot(n,n,cont)
%         SS = sort(S_IRF_all2(var,:,sh,:),4);
%         SS_84 = SS(1,:,1,0.84*maxdraw);
%         SS_16 = SS(1,:,1,0.16*maxdraw);
%         h = plot(irf_axis,zeroline,irf_axis,median(SS,4),...,
%             irf_axis,SS_84,irf_axis,SS_16);
%         set(h(3),'Color','Red','LineWidth',0.5)
%         set(h(4),'Color','Red','LineWidth',0.5)
%         set(h(1),'Color','Black','LineWidth',0.5)
%         set(h(2),'Color','Black','LineWidth',1.5)
%         title(strcat('Response of ',' ',names(var,:),' to ',name_sh_struc(sh,:)))
%         axis tight;grid on;
%         cont = cont+1;
% 
%     end
% 
% end
% 
% 



