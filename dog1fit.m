%% 1st Derivative of Gauss function Fit - mulling/tracking tasks 

function [dog1Fit,dgf,dcheck,damp] = dog1fit(br)

cly =[0.28,0.08,0.02;  0.64,0.08,0.23;  0.8,0.38,0.03; 0.94,0.66,0.24];
clx =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
cl = clx; cond = ["2.5","5","10","15"]; SF = ["0.5SF","1SF","2SF"];
%SUB =["LKC","PSC"];
SUB =["all"];

for sub = 1:1
for sf = 1:3
    figure();
for ecc = 1:4
%     x = br.storeLags(ecc,1:72,sf,sub)';
%     y = br.diffY(ecc,:,sf,sub)';
    x = br.storeLags(ecc,:,sf,sub)';
    y = br.phaseMean(ecc,:,sf,sub)';
    
    ft= fittype('a.*(x-b).*exp(-((x-b)/c).^2)', ...
        'independent', {'x'}, 'dependent', {'y'}, 'coefficient',{'a','b','c'});
    options = fitoptions(ft); % choose better starting points???
    options.StartPoint = [-0.3804    0.4544   -0.1430];
    allci =[];allg = [];
    
    for ite = 1:100
        [myFit G] = fit(x,y,ft,options);
        g = [G.rsquare,G.sse];
        allg = [allg; g];
        allcf(ite,:) = coeffvalues(myFit);
        allci = [allci;confint(myFit)];
    end
    
    for i = 1:100
    if allg(i,1)== max(g(:,1))
        dgf(:,:,ecc,sf,sub) = allg(ite,:);
        cf = allcf(i,:);
        ci = allci(2*ite-1:2*ite,:);
    end
    end
    a = cf(1); b = cf(2); c = cf(3);
    ye = a.*(x-b).*exp(-((x-b)/c).^2);
    yemax = max(ye);
    damp(1,ecc,sf,sub) = yemax;
    for i = 1:length(x)
        if ye(i) == yemax
            xemax = x(i);
        end
    end
    
    sigmafit = abs(c./sqrt(2)); % sigma calculated from c
    intenorm = trapz(x,y);
    sigmanorm = abs((intenorm/(a.*sqrt(2*pi)))^(1/3)); % sigma calculated from a
    FWHM = abs(2*sqrt(log(2)) * c); % FWHM from c
    
    
    dog1Fit(:,:,ecc,sf,sub) = [cf;ci];
    dcheck(:,:,ecc,sf,sub) =[sigmafit,sigmanorm,yemax,xemax];
    
    
    subplot(1,4,ecc); 
    data = plot(x,y,'Color', cl(ecc,:),'LineWidth',2);hold on;
    estimate = plot(x,ye,':','Color', cl(ecc,:),'LineWidth',2);hold on;
    ylim([-0.1, 0.1]);
    legend([data,estimate],'data','estimate');
     title( strcat(SUB(sub)," ",SF(sf)," Eccentricity:  ", cond(ecc),"  s=", sprintf('%.6f',sigmafit),...
         "  r^2=",sprintf('%.6f',dgf(1,1,ecc,sf)))); 

end
end
end 

end
