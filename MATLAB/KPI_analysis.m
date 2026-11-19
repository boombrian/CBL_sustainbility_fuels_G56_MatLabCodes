clear; clc; close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);          % so cylinder_volume.m (same folder) is found from any current folder

%% Units and constants
mm    = 1e-3;
dm    = 0.1;                 
bara  = 1e5;               
kWhr  = 1000*3600;          
volperc = 0.01;             
ppm   = 1e-6;                
g     = 1e-3;              

%molar masses with units [kg/mol]
M_NO2 = 46.0055e-3;          
M_CO2 = 44.0095e-3;         
M_exh = 28.8e-3;             % mean molar mass of diesel exhaust, assuming the value to be close to air

dataFile = fullfile(scriptDir, '..', 'Example Data Set', 'Example Data Set', 'Data', 'ExampleDataSet.txt');
CaResolution = 0.2;          % crank-angle resolution [deg] of the input data

% --- Engine geometry---
Cyl.Bore             = 104*mm;
Cyl.Stroke           = 85*mm;
Cyl.ConRod           = 136.5*mm;
Cyl.CompressionRatio = 21.5;
Cyl.TDCangle         = 0;    

% --- Operating point ---
N   = 1500;                  % engine speed [rpm] 
n_R = 2;                     % revolutions per work cycle

% --- Fuel property ---
LHV = 43.0e6;                % [J/kg] This is a value I found on the internet because I don't have the real values of our viewers yet.
%This one will be used when We get the real experiment data: 
% fuelFile = fullfile(scriptDir, '..', 'Experimental Data', 'FuelProperties.csv');
% fuelTab  = readtable(fuelFile, 'VariableNamingRule', 'preserve');
% fuelName = 'B7';                                                               
% LHV = fuelTab.('LHV_MJ_kg')(strcmp(fuelTab.('Fuel'), fuelName)) * 1e6;           

% --- Slow signals ---
meas.mdot_f   = 0.24*g;      % fuel mass flow [kg/s]  
meas.mdot_air = 9.5*g;       % air mass flow  [kg/s]  
meas.X_CO2    = 5.0*volperc; % CO2 mole fraction  
meas.X_NOx    = 600*ppm;     % NOx mole fraction

%This one will be used when We get the real experiment data: 

% slowFile = fullfile(scriptDir, '..', 'Experimental Data', 'SlowSignals.csv'); 
% slow     = readtable(slowFile, 'VariableNamingRule', 'preserve');
% disp(slow.Properties.VariableNames);   % prints the column names
%
% % Average each signal over the measurement period.
% meas.mdot_f   = mean(slow.('FuelFlow_g_s'), 'omitnan') * g;       
% meas.mdot_air = mean(slow.('AirFlow_g_s'),  'omitnan') * g;       
% meas.X_CO2    = mean(slow.('CO2_vol%'),     'omitnan') * volperc; 
% meas.X_NOx    = mean(slow.('NOx_ppm'),      'omitnan') * ppm;    


%% Load and reshape the data
dataIn  = readmatrix(dataFile);
NptsPerCycle = round(720/CaResolution);
Ncycles = size(dataIn,1)/NptsPerCycle;
if Ncycles ~= round(Ncycles)
    error('Number of rows (%d) is not a multiple of %d points per cycle.', size(dataIn,1), NptsPerCycle);
end

Ca = reshape(dataIn(:,1), NptsPerCycle, Ncycles);          % [deg], one column per cycle
p  = reshape(dataIn(:,2), NptsPerCycle, Ncycles) * bara;   % [Pa],  one column per cycle
CaVec = Ca(:,1);                                           % same crank-angle grid for every cycle

%% Pre-processing
% these will be included in future SSAs, which includes noise filtering and pegging

%% Cylinder volume (slider-crank)
% Uses cylinder_volume.m wrote by Rivan in the last SSA.
V  = cylinder_volume(CaVec, Cyl);                          % [m^3]
Vd = pi/4 * Cyl.Bore^2 * Cyl.Stroke;                       % displacement volume [m^3]
Vc = Vd / (Cyl.CompressionRatio - 1);                      % clearance volume [m^3]


%% Actual work per cycle

Vloop = [V; V(1)];
W_i = zeros(1, Ncycles);
for j = 1:Ncycles
    pLoop  = [p(:,j); p(1,j)];
    W_i(j) = trapz(Vloop, pLoop);                          % [J]
end

%% IMEP
IMEP = W_i / Vd;                                           % [Pa]

%% Indicated power and efficiency
P_i     = W_i * N / (60*n_R);                              % [W]
m_f_cyc = meas.mdot_f * 60*n_R / N;                        % fuel mass per cycle [kg]
eta_i   = W_i ./ (m_f_cyc * LHV);                          % [-] 

%% Indicated specific fuel consumption
isfc = meas.mdot_f ./ P_i;                                 % [kg/J]

%% isNOx, isCO2
mdot_exh = meas.mdot_air + meas.mdot_f;                    % [kg/s] 
mdot_NOx = meas.X_NOx * (M_NO2/M_exh) * mdot_exh;          % [kg/s]
mdot_CO2 = meas.X_CO2 * (M_CO2/M_exh) * mdot_exh;          % [kg/s]
isNOx = mdot_NOx ./ P_i;                                   % [kg/J]
isCO2 = mdot_CO2 ./ P_i;                                   % [kg/J]

%% Results and some statistics
toGkWh = 1/g * kWhr;                                       % kg/J -> g/kWh

KPI.Name  = {'W_i'; 'IMEP'; 'P_i'; 'eta_i'; 'isfc'; 'isCO2'; 'isNOx'};
KPI.Unit  = {'J'; 'bar'; 'kW'; '%'; 'g/kWh'; 'g/kWh'; 'g/kWh'};
KPI.Cycle = {W_i; IMEP/bara; P_i/1e3; eta_i*100; isfc*toGkWh; isCO2*toGkWh; isNOx*toGkWh};
KPI.Mean  = cellfun(@mean, KPI.Cycle);
KPI.Std   = cellfun(@std,  KPI.Cycle);

fprintf('\nIndicated KPIs, %d cycles at %d rpm (inputs marked DUMMY in section 2 are not measured):\n', Ncycles, N);
fprintf('  V_d = %.4f L, V_c = %.4f L\n', Vd/dm^3, Vc/dm^3);
for k = 1:numel(KPI.Name)
    fprintf('  %-6s = %9.3f +/- %7.3f %s\n', KPI.Name{k}, KPI.Mean(k), KPI.Std(k), KPI.Unit{k});
end
fprintf('  COV_IMEP = %.2f %%\n', 100*std(IMEP)/mean(IMEP));

%% ===================================================================
%% PLOTS (this one the style is finetuned with Claude with skills)
%% ===================================================================
iSel = 10;                                                 % cycle highlighted in the plots
pMean = mean(p, 2);

cBlue    = [15 77 146]/255;
cBlue2   = [55 117 186]/255;
cRed     = [182 67 66]/255;
cNeutral = [207 206 206]/255;

%% 12.1 Pressure vs crank angle
figure('Name', 'p-CA', 'Position', [100 600 1100 400], 'Color', 'w');
plot(CaVec, p/bara, 'Color', cNeutral, 'LineWidth', 0.5); hold on;
hAll  = plot(nan, nan, 'Color', cNeutral, 'LineWidth', 1.5);         % legend entry for all cycles
hMean = plot(CaVec, pMean/bara, 'Color', cBlue, 'LineWidth', 2.5);
hSel  = plot(CaVec, p(:,iSel)/bara, 'Color', cRed, 'LineWidth', 1.5);
xlabel('Crank angle \theta [deg]'); ylabel('p [bar]');
xlim([-360 360]); set(gca, 'XTick', -360:60:360);
legend([hAll hMean hSel], {'All cycles', 'Mean of all cycles', sprintf('Cycle %d', iSel)}, 'Location', 'northwest');
title('Pressure vs CA');

%% pV diagram (linear and log-log)
figure('Name', 'pV', 'Position', [100 100 1100 450], 'Color', 'w');
subplot(1,2,1);
plot(V/dm^3, pMean/bara, 'Color', cBlue, 'LineWidth', 2.5); hold on;
plot(V/dm^3, p(:,iSel)/bara, 'Color', cRed, 'LineWidth', 1.5);
xlabel('V [L]'); ylabel('p [bar]');
legend('Mean of all cycles', sprintf('Cycle %d', iSel), 'Location', 'northeast');
title('pV diagram (linear)');
subplot(1,2,2);
loglog(V/dm^3, pMean/bara, 'Color', cBlue, 'LineWidth', 2.5); hold on;
loglog(V/dm^3, p(:,iSel)/bara, 'Color', cRed, 'LineWidth', 1.5);
xlabel('V [L]'); ylabel('p [bar]');
xlim([0.03 0.8]); set(gca, 'XTick', [0.03 0.05 0.1 0.2 0.5 0.8], 'YTick', [1 2 5 10 20 50]);
legend('Mean of all cycles', sprintf('Cycle %d', iSel), 'Location', 'northeast');
title('pV diagram (log-log)');

%% Work and IMEP per cycle
figure('Name', 'Work and IMEP', 'Position', [1250 600 600 450], 'Color', 'w');
subplot(2,1,1);
plot(1:Ncycles, W_i, 'o-', 'Color', cBlue2, 'MarkerFaceColor', cBlue2, 'MarkerSize', 4, 'LineWidth', 1.5); hold on;
yline(mean(W_i), '--', 'Mean', 'Color', cRed, 'LineWidth', 1.5);
xlabel('Cycle number'); ylabel('W_i [J]'); title('Indicated work per cycle');
subplot(2,1,2);
plot(1:Ncycles, IMEP/bara, 'o-', 'Color', cBlue2, 'MarkerFaceColor', cBlue2, 'MarkerSize', 4, 'LineWidth', 1.5); hold on;
yline(mean(IMEP)/bara, '--', 'Mean', 'Color', cRed, 'LineWidth', 1.5);
xlabel('Cycle number'); ylabel('IMEP [bar]'); title('IMEP per cycle');

%% Efficiency and specific KPIs per cycle
kpiPlot = {eta_i*100, '\eta_i [%]', 'Indicated efficiency';
           isfc*toGkWh, 'isfc [g/kWh]', 'Indicated specific fuel consumption';
           isCO2*toGkWh, 'isCO_2 [g/kWh]', 'Indicated specific CO_2';
           isNOx*toGkWh, 'isNO_x [g/kWh]', 'Indicated specific NO_x'};
figure('Name', 'KPIs per cycle', 'Position', [1250 100 900 600], 'Color', 'w');
for k = 1:size(kpiPlot,1)
    subplot(2,2,k);
    y = kpiPlot{k,1};
    if all(isnan(y))
        text(0.5, 0.5, '-', 'HorizontalAlignment', 'center'); axis off;
    else
        plot(1:Ncycles, y, 'o-', 'Color', cBlue2, 'MarkerFaceColor', cBlue2, 'MarkerSize', 4, 'LineWidth', 1.5); hold on;
        yline(mean(y), '--', 'Mean', 'Color', cRed, 'LineWidth', 1.5);
        xlabel('Cycle number'); ylabel(kpiPlot{k,2});
    end
    title(kpiPlot{k,3});
end
sgtitle('KPIs per cycle');

set(findall(groot, 'Type', 'axes'), 'FontName', 'Arial', 'FontSize', 13, 'LineWidth', 1.5, ...
    'Box', 'off', 'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off');
set(findall(groot, 'Type', 'legend'), 'Box', 'off', 'FontSize', 12);
