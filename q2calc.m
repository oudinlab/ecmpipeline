function [Q2] = q2calc(Ypred,Yact,pcomps)

    % This function calculates the Q2 from leave one out cross validation
    % Inputs:
        % Ypred = Matrix of predicted cell responses (from CrossValFunc)
        % Yact = Matrix of actual cell response values
        % pcomps = number of principal components
    % Output:
        % Q2 = vector of Q2 values per each principal component
        
        
% Concatenate Ypred and Yact
Y = [Ypred,Yact];

[YR, YC] = size(Y);
[~,YVAR] = size(Yact);

% Separate out speed, persistence, invasion
for i = 1:YVAR
    sep(i).sep = Y(:,i:YVAR:end);
end
    
% Find Min and Max values for each measurement

minmax = zeros(YVAR,2);

for i = 1:YVAR
    minmax(i,1) = min(sep(i).sep, [], 'all');
    minmax(i,2) = max(sep(i).sep, [], 'all');
end

% Mean Center

MC = zeros(YR, YC);

for i = 1:YVAR
    MC(:,[i:YVAR:end]) = ((sep(i).sep-minmax(i,1))./(minmax(i,2)-minmax(i,1)))+1;
end 
    
% Calculate PRESS
for i = 1:pcomps
    if (YVAR==3)
    press(i).a = (MC(:,[(3*i)-2, (3*i)-1, (3*i)])-MC(:,[end-2,end-1,end])).^2;
    elseif (YVAR==2)
        press(i).a = (MC(:,[(2*i)-1, (2*i)])-MC(:,[end-1,end])).^2;
    elseif (YVAR==1)
        press(i).a = (MC(:,i)-MC(:,end)).^2;
    end
    PRESS(i) = [sum(press(i).a, 'all')];
end
    
% Calculate DEN
    
Ypredsum = sum(MC(:,[1:end-YVAR]).^2, 'all');
Ypredave = (sum(MC(:,[1:end-YVAR]), 'all').^2)./numel(Ypred);

DEN = Ypredsum-Ypredave;

% Calculate Q2

Q2 = 1 - PRESS./DEN;

end