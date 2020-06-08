%% Create PLSR Model
% This script creates a PLS model using the PLSRegress algorithm in Matlab
% Inputs: 
%   - training and test data sets
%   - Labels for conditions and variables
%   - number of PCs (for PLS models and calculating VIP scores)
% Outputs:
%   - Variance explained in X and Y per number of PCs
%   - VIP scores (variable = VIP)
%   - PLS Loadings and Scores Plot
%   - Prediction of new data (variable = Ypred)

close all
clear all

%% Load Data - Training Set

%----- CHOOSE ONE OPTION BELOW TO LOAD DATA -----%
% -- OPTION 1: --

% inputfile = "231-PLSData.xlsx";
% d = uigetdir('Select a folder');
% data = readmatrix(d+"/"+inputfile);

% -- OPTION 2: --

% Manually input file path: (more streamlined)
d = "/Users/jananibaskaran/Box/Oudin Lab/Shared Data/Janani Baskaran/Matlab Resources";
data = readmatrix(d+"/231-PLSData.xlsx", 'Sheet', 'Sheet1');

% Assign variables

  Y0i = data(1:5, [13:15]);     % cell response variables (n by m matrix)
   
  X0i = data(1:5, [2:12]);      % cell shape variables (n by p matrix)

% Name the ECM conditions and X and Y labels (will be useful for labeling
% graphs)

    CondLabs = {'Control'
        'C1'
        'FN'
        'TNC'
        'C4'
            };
    
    Xlabs = {'Area/Cell (sq um)'
        'Aspect Ratio'
        'Compactness'
        'Eccentricity'
        'Extent'
        'Form Factor'
        'Max Feret Diameter'
        'Mean Radius'
        'Min Feret Diameter'
        'Perimeter'
        'Solidity'       
        };
    
    %--- can also input labels from excel sheet:
    %Xlabs = readcell(d+"/Xlabs.xlsx", 'Sheet', 'Sheet1');
    
    Ylabs = {'2D Speed'
        '2D Persistence'
        '3D Invasion'};
    
% Define which variables/parameters to include in model UNCOMMENT if needed
% varinc = 1:11;
% X0i = X0i(:,varinc);
% Xlabs = Xlabs(varinc);


% Scale training variables (sets standard deviation to 1)
 X0 = (X0i) ./(sqrt(var(X0i)));
 Y0 = (Y0i) ./(sqrt(var(Y0i)));

%% Run the PLS
Pcomps = 4;     % define number of PCS (maximum = no. of ECM Conditions - 1)

[Xloadings,Yloadings,Xscores,Yscores,beta,PctVar,mse,stats] = plsregress(X0,Y0,Pcomps);
    % X loading - p-by-pcomp matrix; coefficients define linear comb. of
        % Pcomps that approx. the original predictor variables
    % Y loading - m-by-pcomp matrix; coefficients define linear comb. of
        % Pcomps that approx. the original response variables
    % X scores - n-by-pcomp matrix; coeffs map out each ECM condition on
        % Pcomps based on predictor variables
    % Y scores - n-by-pcomp matrix; map out each ECM condition on Pcomps
        % based on both predictor and response variables
    % beta - regression coefficients, use this to predict
    % PctVar - percent variation explained from the model (row 1 is X and 2
        % is Y)
    % mse - estimated mean square error (row 1 is X and row 2 is Y)
    % stats - structure with different statitics, used to calculate VIP
        % Scores

%% ---- Plot Variance Explained by PC --> this is the R^2:
R2 = cumsum(100*PctVar(2,:));
figure
subplot(1,2,1)
plot(1:Pcomps,R2,'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in Y');

subplot(1,2,2)
plot(1:Pcomps,cumsum(100*PctVar(1,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in X');

%% ----- Plot Loadings and Scores:
% ----- Plot PLS Loadings:
figure
scatter(Xloadings(:,1),Xloadings(:,2),[])
hold on
scatter(Yloadings(:,1), Yloadings(:,2),[])
hold off
title('PLS Loadings Plot')
xlabel('PC 1')
ylabel('PC 2')
hax = gca;
line([0 0],get(hax,'YLim'),'Color','k','LineStyle','--')
hline = refline([0 0]); hline.Color = 'k'; hline.LineStyle = '--';
dx = 0.03; dy = 0.03; % displacement so the text does not overlay the data points
text(Xloadings(:,1)+dx,Xloadings(:,2)+dy, Xlabs, 'Fontsize', 10, 'Interpreter', 'none'); % labeling points
text(Yloadings(:,1)+dx,Yloadings(:,2)-dy, Ylabs, 'Fontsize', 10, 'Interpreter', 'none', 'Color','red'); % labeling points

% ----- Plot PLS X Scores:
figure
scatter(Xscores(:,1), Xscores(:,2))
title('PLS X Scores Plot')
xlabel('PC 1')
ylabel('PC 2')

hax = gca;
line([0 0],get(hax,'YLim'),'Color','k','LineStyle','--')
hline = refline([0 0]); hline.Color = 'k'; hline.LineStyle = '--';
yline = 0;
xline = 0;
dx = 0.02; dy = 0.02; % displacement so the text does not overlay the data points
text(Xscores(:,1)+dx, Xscores(:,2)+dy, CondLabs, 'Fontsize', 10, 'Interpreter', 'none'); % labeling points

% ----- Plot PLS Y Scores:
figure
scatter(Yscores(:,1), Yscores(:,2))
title('PLS Y Scores Plot')
xlabel('PC 1')
ylabel('PC 2')

hax = gca;
line([0 0],get(hax,'YLim'),'Color','k','LineStyle','--')
hline = refline([0 0]); hline.Color = 'k'; hline.LineStyle = '--';
yline = 0;
xline = 0;
dx = 0.02; dy = 0.02; % displacement so the text does not overlay the data points
text(Yscores(:,1)+dx, Yscores(:,2)+dy, CondLabs, 'Fontsize', 10, 'Interpreter', 'none'); % labeling points


%% ---- Calculate VIP Scores: for each Y variable

Pcomps = 2;     % with #PC ideal for model
for k = 1:3         % change so easier to adjust
[Xl,Yl,Xs,Ys,beta,PctVar,mse,stats] = plsregress(X0,Y0(:,k),Pcomps);

W = stats.W;
[n,p] = size(X0);
SSa = zeros(Pcomps,1);
VIPa = zeros(p, 1);


for i = 1: Pcomps
    SSa(i) = Yl(i)^2 * Xs(:,i)' * Xs(:,i);     % Calculate sum of squares for each PC
    W(:,i) = (W(:,i) / norm(W(:,i))).^2;                        % normalize and square PLS weights
end
SStot = sum(SSa);
for j = 1:p 
    for i = 1 : Pcomps
        VIPa(j) = VIPa(j) + p * W(j,i) * SSa(i) / SStot ;
    end
end
VIPa = VIPa.^0.5; % sq rt everything
VIP(:,k) = VIPa;
end

figure
subplot(1,3,1)
bar(VIP(:,1))
title('VIP Score Prediction of 2D Speed');
ylabel('VIP Score')
box on;
set(gca,'xticklabel',Xlabs,'XTickLabelRotation',45);

subplot(1,3,2)
bar(VIP(:,2))
title('VIP Score Prediction of 2D Persistence');
ylabel('VIP Score')
box on;
set(gca,'xticklabel',Xlabs,'XTickLabelRotation',45);

subplot(1,3,3)
bar(VIP(:,3))
title('VIP Score Prediction of 3D Invasion');
ylabel('VIP Score')
box on;
set(gca,'xticklabel',Xlabs,'XTickLabelRotation',45);


%% Predict with PLS

Pcomps = [1 2 3 4];

   X1i = data(6, 2:12);    %shape parameters for new ECM
   X1 = X1i ./ (sqrt(var(X0i)));        % scale
   Ypreds = zeros(4,3);
   Ypred = zeros(4,3);


for i = 1:length(Pcomps)
    [Xl,Yl,Xs,Ys,beta,PctVar,mse,stats] = plsregress(X0,Y0,Pcomps(i));
    [nt,~] = size(X1);
   pred1 = [ones(nt,1) X1];
   Ypreds(i,:) = pred1*beta;
   
   Ypred(i,:) = Ypreds(i,:) .* (sqrt(var(Y0i)));
end
