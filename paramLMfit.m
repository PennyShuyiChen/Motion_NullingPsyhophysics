%% parameter fits to LM - nulling/tracking task
% zeb - a,b,c,sigma, estimated Rmax-X/Y;
% cfx - a,b,c,sigma, estimated Rmax-X/Y -- for both phase & Ydiff
function [LMfit,T,rawdata] = paramLMfit(funcFit,check)
npars = size(funcFit,1)+3; 
E = [2.5,5,10,15]; SF = ["0.5SF","1SF","2SF"];
SUB =["all"];
%SUB =["LKC","PSC"];
% model function table 
varNames = convertCharsToStrings({'Intercept', 'Slope'});
rowNames = convertCharsToStrings({'a','b','c','sigma','estYmax','estX@max'});

% LMfit results structure 
LMfit.est =[];
LMfit.funcexp =[];
LMfit.tStat = [];
LMfit.pval = [];
LMfit.rsq = [];

for sub = 1:1
for sf = 1:3
    figure();
    for np = 1:npars
        for ecc = 1:4
            if np <4
                rawFit(ecc,np,sf,sub) = funcFit(1,np,ecc,sf,sub);
            elseif np==4       
                rawFit(ecc,np,sf,sub) = check(1,1,ecc,sf,sub);  
            else
                if size(check,2) ==4
                    rawFit(ecc,np,sf,sub) = check(1,np-2,ecc,sf,sub); % GaussFit
                else rawFit(ecc,np,sf,sub) = check(1,np,ecc,sf,sub); % Dog1Fit 
                end
            end
        end
        
        % update the LMfit structure 
        lmf{sf} = fitlm(E,rawFit(:,np,sf,sub));
        LMfit.est(:,np,sf,sub) = lmf{sf}.Fitted;
        LMfit.funcexp(:,np,sf,sub) = lmf{sf}.Coefficients.Estimate;
        LMfit.tStat(:,np,sf,sub) = lmf{sf}.Coefficients.tStat;
        LMfit.pval(:,np,sf,sub) = lmf{sf}.Coefficients.pValue;
        LMfit.rsq(:,np,sf,sub) = lmf{sf}.Rsquared.Ordinary;
       
         % plot the LM with raw params 
        subplot(1,npars,np);
        raw(sf) = plot(E,rawFit(:,np,sf,sub),'LineWidth',2); hold on;
        rawdata(:,np,sf,sub) = rawFit(:,np,sf,sub);
        lm(sf) = plot(E,LMfit.est(:,np,sf,sub),'r*','LineWidth',2); hold on;
        
        syms f(x)
        f(x) = LMfit.funcexp(1,np,sf,sub) + x .* LMfit.funcexp(2,np,sf,sub);
        func(sf) = fplot(f,[0,15],'LineWidth',2); hold on; 
        
        Enew = [0,E];
        [ypred,yci] = predict(lmf{sf}, Enew'); 
       
         plot(Enew, yci, '--r');
        
        xlabel('Eccentricity(deg)'); ylabel('Prameter values');
        title(  strcat(SUB(sub)," ",SF(sf),rowNames(np),"  R^2=",sprintf('%.6f',LMfit.rsq(:,np,sf,sub)),...
            "  p=",sprintf('%.6f',LMfit.pval(2,np,sf,sub))));
    end
end
end

    legend([raw,lm, func],'raw paramter fits', 'model estimates','model function','show','Location','Northeast');
    for sub = 1:1   
    for sf = 1:3
        figure(); % Table of the functions
        intercept = LMfit.funcexp(1,:,sf,sub); slope = LMfit.funcexp(2,:,sf,sub);
        T = table(intercept',slope','VariableNames',varNames,'RowNames',rowNames);
        uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,...
            'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
    end
    end
    LMfit; T;
end