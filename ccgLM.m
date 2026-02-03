%% parameter fits to LM - nulling/tracking task
% zeb - a,b,c,sigma, estimated Rmax-X/Y;
% cfx - a,b,c,sigma, estimated Rmax-X/Y -- for both phase & Ydiff
function [LMfit,T] = paramLMfit(funcFit,check)
npars = size(funcFit,1)+3; 
E = [2.5,5,10,15]; SF = ["0.5SF","1SF","2SF"];
% model function table 
varNames = convertCharsToStrings({'Intercept', 'Slope'});
rowNames = convertCharsToStrings({'a','b','c','sigma','estYmax','estX@max'});

% LMfit results structure 
LMfit.est =[];
LMfit.funcexp =[];
LMfit.tStat = [];
LMfit.pval = [];
LMfit.rsq = [];

for sf = 1:3
    for np = 1:npars
        for ecc = 1:4
            if np <4
                rawFit(ecc,np,1) = funcFit(1,np,ecc);
            elseif np==4       
                rawFit(ecc,np,1) = check(1,1,ecc);  
            else
                if size(check,2) ==4
                    rawFit(ecc,np,1) = check(1,np-2,ecc); % GaussFit
                else rawFit(ecc,np,1) = check(1,np,ecc); % Dog1Fit 
                end
            end
        end
        
        % upfate the LMfit structure 
        lmf = fitlm(E,rawFit(:,np));
        LMfit.est(:,np) = lmf.Fitted;
        LMfit.funcexp(:,np) = lmf.Coefficients.Estimate;
        LMfit.tStat(:,np) = lmf.Coefficients.tStat;
        LMfit.pval(:,np) = lmf.Coefficients.pValue;
        LMfit.rsq(:,np) = lmf.Rsquared.Ordinary;
       
        figure(1); % plot the LM with raw params 
        subplot(1,npars,np);
        raw = plot(E,rawFit(:,np),'LineWidth',2); hold on;
        lm = plot(E,LMfit.est(:,np),'r*','LineWidth',2); hold on;
        
        syms f(x)
        f(x) = LMfit.funcexp(1,np) + x .* LMfit.funcexp(2,np);
        func = fplot(f,[0,15],'LineWidth',2); hold on; 
        
        Enew = [0,E];
        [ypred,yci] = predict(lmf, Enew'); pci = plot(Enew, yci, '--r');
        
        xlabel('Eccentricity(deg)'); ylabel('Prameter values');
        title(  strcat(rowNames(np),"  R^2=",sprintf('%.6f',LMfit.rsq(:,np)),...
            "  p=",sprintf('%.6f',LMfit.pval(2,np))));
    end
end
    legend([raw,lm, func],'raw paramter fits', 'model estimates','model function','show','Location','Northeast');
        
        figure(2); % Table of the functions 
        intercept = LMfit.funcexp(1,:); slope = LMfit.funcexp(2,:);
        T = table(intercept',slope','VariableNames',varNames,'RowNames',rowNames);
        uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
        'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
LMfit; T;
end