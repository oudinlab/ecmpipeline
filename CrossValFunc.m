function [Ypred] = CrossValFunc(X0i,Y0i,Pcomps)
    % Function performs a leave-one-out cross validation of PLS model
    % Inputs:
        % Unscaled cell shape variable data (X0i)
        % Unscaled cell response data (Y0i)
        % Vector of Principal Components
    % Output:
        % Matrix of predicted cell responses
        
[LLY,VAR]=size(X0i);
[LLY,YVAR]=size(Y0i);
Ypred = zeros(LLY, YVAR*length(Pcomps));

for x = 1 : LLY
    if(x==1)
        sampvec = 2 : LLY;
    else
        sampvec = [1:(x-1) x+1:LLY];
    end
    
    yyi  = Y0i(sampvec,:);
        yy = yyi ./(sqrt(var(yyi)));
    YYTi = Y0i(x,:); % not used in the code
    XXi  = X0i(sampvec,:);
        XX = XXi ./(sqrt(var(XXi)));
    XXTi = X0i(x,:);
        XXT = XXTi ./(sqrt(var(XXi)));
    [n,p] = size(XX);
    [nt,pt] = size(XXT);
    
    for i = 1:length(Pcomps)
    
        [Xloadings,Yloadings,Xscores,Yscores,betaPLS10,PLSPctVar,mse,stats] = plsregress(XX,yy,Pcomps(i));
    
        Ypreds = [ones(nt,1) XXT]*betaPLS10;
        Ypredus = Ypreds .* (sqrt(var(yyi)));
        if (YVAR==3)
            Ypredi(:,[(3*i)-2 (3*i)-1 3*i]) = Ypredus;
        elseif (YVAR==2)
            Ypredi(:,[(2*i)-1 2*i]) = Ypredus;
        elseif (YVAR==1)
            Ypredi(:,i) = Ypredus;
        end
    end
    Ypred(x,:) = Ypredi';
end