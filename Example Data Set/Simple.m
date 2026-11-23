%% Info
% Toerental: 1500 RPM
% SOA van 4.2º voor TDC
% Resolutie van 0.2º CA
% Data voor 69 cycles (maximale van de Smetec, de OGO gensets kunnen in principe “onbeperkt” aan)
% 
%% init
clear all; clc;close all;
addpath( "Functions","Nasa");
%% Units
mm      = 1e-3;dm=0.1;
bara    = 1e5;
MJ      = 1e6;
kWhr    = 1000*3600;
volperc = 0.01; % Emissions are in volume percentages
ppm     = 1e-6; % Some are in ppm (also a volume- not a mass-fraction)
g       = 1e-3;
s       = 1;
%% Load NASA maybe you need it at some point?
% Global (for the Nasa database in case you wish to use it).
global Runiv
Runiv = 8.314;
[SpS,El]        = myload('Nasa\NasaThermalDatabase.mat',{'Diesel','O2','N2','CO2','H2O'});
%% Engine geom data (check if these are correct)
Cyl.Bore                = 104*mm;
Cyl.Stroke              = 85*mm;
Cyl.CompressionRatio    = 21.5;
Cyl.ConRod              = 136.5*mm;
Cyl.TDCangle            = 0;
% -- Valve closing events can sometimes be seen in fast oscillations in the pressure signal (due
% to the impact when the Valve hits its seat).
CaIVO = -355;
CaIVC = -135;
CaEVO = 149;
CaEVC = -344;
CaSOI = -3.2;
% Write a function [V] = CylinderVolume(Ca,Cyl) that will give you Volume
% for the given Cyl geometry. If you can do that you can create pV-diagrams
%% Load data (if txt file)
FullName        = fullfile('Data','ExampleDataSet.txt');
dataIn          = table2array(readtable(FullName));
[Nrows,Ncols]   = size(dataIn);                    % Determine size of array
NdatapointsperCycle = 720/0.2;                     % Nrows is a multitude of NdatapointsperCycle
Ncycles         = Nrows/NdatapointsperCycle;       % This must be an integer. If not checkwhat is going on
Ca              = reshape(dataIn(:,1),[],Ncycles); % Both p and Ca are now matrices of size (NCa,Ncycles)
p               = reshape(dataIn(:,2),[],Ncycles)*bara; % type 'help reshape' in the command window if you want to know what it does (reshape is a Matlab buit-in command
%% Plotting 
f1=figure(1);
set(f1,'Position',[ 200 800 1200 400]);             % Just a size I like. Your choice
pp = plot(Ca,p/bara,'LineWidth',1);                 % Plots the whole matrix
xlabel('Ca');ylabel('p [bar]');                     % Always add axis labels
xlim([-360 360]);ylim([0 50]);                      % Matter of taste
iselect = 10;                                    % Plot cycle 10 again in the same plot to emphasize it. Just to show how to access individual cycles.
line(Ca(:,iselect),p(:,iselect)/bara,'LineWidth',2,'Color','r');
YLIM = ylim;
% Add some extras to the plot
line([CaIVC CaIVC],YLIM,'LineWidth',1,'Color','b'); % Plot a vertical line at IVC. Just for reference not a particular reason.
line([CaEVO CaEVO],YLIM,'LineWidth',1,'Color','r'); % Plot a vertical line at EVO. Just for reference not a particular reason.
set(gca,'XTick',[-360:60:360],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title('All cycles in one plot.')
%% pV-diagram
V = CylinderVolume(Ca(:,iselect),Cyl);
f2 = figure(2);
set(f2,'Position',[ 200 400 600 800]);              % Just a size I like. Your choice
subplot(2,1,1)
plot(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0 0.8]);ylim([0.5 50]);                      % Matter of taste
set(gca,'XTick',[0:0.1:0.8],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram'})
subplot(2,1,2)
loglog(V/dm^3,p(:,iselect)/bara);
xlabel('V [dm^3]');ylabel('p [bar]');               % Always add axis labels
xlim([0.02 0.8]);ylim([0 50]);                      % Matter of taste
set(gca,'XTick',[0.02 0.05 0.1 0.2 0.5 0.8],...
    'YTick',[0.5 1 2 5 10 20 50],'XGrid','on','YGrid','on');        % I like specific axis labels. Matter of taste
title({'pV-diagram'})


% Heat release analysis

% Parameters

dCa     = 0.2; % Crank angle resolution
gam     = 1.35; % Gamma value (cp/cv)
CaStart = -20; % Lower limit for crank angle integration
CaEnd = 80; % Upper limit for crank angle integration

Ca1 = Ca(:,1); % Crank angle values for one cycle
V1  = CylinderVolume(Ca1, Cyl); % Cylinder volume array
Nc  = size(p, 2);             % Number of recorded cycles
Npts = length(Ca1);           % Total data points per cycle

ps     = movmean(p, 9, 1); % pressure 9 sample moving average

dpdCa = zeros(Npts, Nc); % Initialize change in pressure per change in crank angle array
dVdCa = zeros(Npts, 1); % Initialize change in volume per change in crank angle array

for i = 2:(Npts - 1) % Iterate from second to second last data point
    dVdCa(i) = (V1(i+1) - V1(i-1)) / (2 * dCa); % Using standard formula for dV/dCa
    for j = 1:Nc % Iterate through all cycles
        dpdCa(i, j) = (ps(i+1, j) - ps(i-1, j)) / (2 * dCa); % Using standard formula for dp/dCa
    end
end

AHRR = zeros(Npts, Nc); % Initialize Apparent Heat Release Rate array

for j = 1:Nc % Iterate through all cycles
    for i = 1:Npts % For each cycle iterate through all data points
        term1 = (gam / (gam - 1)) * ps(i, j) * dVdCa(i); % Using part of standard formula for AHRR
        term2 = (1 / (gam - 1)) * V1(i) * dpdCa(i, j); % Using part of standard formula for AHRR
        AHRR(i, j) = term1 + term2; % Adding terms for AHRR formula
    end
end

iWindow = find(Ca1 >= CaStart & Ca1 <= CaEnd); % Index range of interest

% Cropped version of data (cropped to iWindow range)
CaW     = Ca1(iWindow);
AHRRw   = AHRR(iWindow, :);

Nw      = length(CaW); % Number of points in the window

Qnet   = zeros(1, Nc); % Initialize net heat released array
CA10   = zeros(1, Nc); % Initialize CA10 array
CA50   = zeros(1, Nc); % Initialize CA50 array
CA90   = zeros(1, Nc); % Initialize CA90 array
XB     = zeros(Nw, Nc); % Initialize mass fraction burned array

for j = 1:Nc % Iterate through each cycle
    Q = zeros(Nw, 1); % Initialize array for Q
    running_sum = 0; % Initialize variable for sum

    for i = 1:Nw % Iterate through each point in the window
        running_sum = running_sum + AHRRw(i, j) * dCa; % Add AHRR value for the specific point in the specific cycle to the sum
        Q(i) = running_sum; % Add current sum to the Q array
    end
    
    [~, index_peak] = max(AHRRw(:, j)); % Find peak heat release rate point
    

    Q0 = Q(1); % Set minimum Q to first item
    index_start = 1; % Set index for Q0
    for i = 1:index_peak % Iterate from the beggining till peak heat release

        % If current Q is less than Q0 current Q becomes the new Q0
        % And index is saved to index_start
        if Q(i) < Q0
            Q0 = Q(i);
            index_start = i;
        end
    end
    
    Qend = Q(end); % find Q left at the end
    Qnet(j) = Qend - Q0; % Find the Qnet difference between Q at end and at its minimum
    

    for i = 1:Nw % Iterate through all points in the window
        XB(i, j) = (Q(i) - Q0) / (Qend - Q0); % Using mass fraction burned formula
    end
    
    % Step F: Find CA10, CA50, and CA90 crank angles
    CA10(j) = ca_at(XB(:, j), CaW, index_start, 0.10);
    CA50(j) = ca_at(XB(:, j), CaW, index_start, 0.50);
    CA90(j) = ca_at(XB(:, j), CaW, index_start, 0.90);
end


    %% 7. Summary Printout
fprintf('\n================== HEAT RELEASE ANALYSIS ==================\n');
fprintf('Mean Net Heat Release : %.0f J/cycle\n', mean(Qnet));
fprintf('CA10 (10%% Burned)      : %.1f° ATDC\n', mean(CA10));
fprintf('CA50 (50%% Burned)      : %.1f° ATDC\n', mean(CA50));
fprintf('CA90 (90%% Burned)      : %.1f° ATDC\n', mean(CA90));
fprintf('Combustion Duration   : %.1f° CA\n', mean(CA90 - CA10));
fprintf('===========================================================\n');

%% 8. Plots
figure(3);
subplot(1, 2, 1);
plot(CaW, AHRRw, 'Color', [0.8 0.8 0.8]); hold on;
plot(CaW, mean(AHRRw, 2), 'k', 'LineWidth', 2);
xline(CaSOI, 'b--', 'SOI'); yline(0, 'k:');
xlabel('Crank Angle [°]'); ylabel('AHRR [J/deg]');
title('Apparent Heat Release Rate'); grid on; xlim([CaStart CaEnd]);

subplot(1, 2, 2);
plot(CaW, XB, 'Color', [0.8 0.8 0.8]); hold on;
plot(CaW, mean(XB, 2), 'k', 'LineWidth', 2);
yline([0.1 0.5 0.9], 'b:');
plot([mean(CA10) mean(CA50) mean(CA90)], [0.1 0.5 0.9], 'ro', 'MarkerFaceColor', 'r');
xlabel('Crank Angle [°]'); ylabel('Burned Fraction x_b [-]');
title('Mass Fraction Burned'); grid on; xlim([CaStart CaEnd]);


function ca = ca_at(xb, CaW, index_start, target_fraction)
    idx = index_start;
    while idx <= length(xb) && xb(idx) < target_fraction
        idx = idx + 1;
    end
    
    if idx > 1 && idx <= length(xb)
        x1 = xb(idx - 1); x2 = xb(idx);
        y1 = CaW(idx - 1); y2 = CaW(idx);
        ca = y1 + (target_fraction - x1) * (y2 - y1) / (x2 - x1);
    else
        ca = CaW(idx);
    end
end