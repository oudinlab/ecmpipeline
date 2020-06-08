%% Validate PLSR Model
% This script creates a PLS model using the PLSRegress algorithm in Matlab
% Inputs: 
%   - data set

% Outputs:
%   - Cross validation predictions (variable = Ypred)
%   - Q2 values (variable = Q2)

close all
clear all

%% Load Data - Training Set

inputfile = "PLSData.xlsx";
d = uigetdir('Select a folder');
data = readmatrix(d+"/"+inputfile);

% Assign variables

  Y0i = data(1:12, [13:15]);     % cell response variables (n by m matrix)
   
  X0i = data(1:12, [2:12]);      % cell shape variables (n by p matrix)

%% Run Model Cross Validation

[n,~] = size(X0i);
Pcomps = 1:(n-2);
Ypred = CrossValFunc(X0i,Y0i,Pcomps);
    % inputs [cell shape var, cell response vars, # PCs]

%% Q2 Calculation

 pc = Pcomps(end);

 Q2 = q2calc(Ypred,Y0i,pc);
    % inputs [output from CrossValFunc, cell repsonse vars, # PC)
    
%% Make Q2 Plot

figure
plot(1:pc, Q2,'-ro');
xlabel('Number of PLS components');
ylabel('Q2');
ylim([0,1])

 
 