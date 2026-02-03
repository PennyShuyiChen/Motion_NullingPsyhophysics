%% Gauss function fit

function [gaussFit,gf,check,amp] = GaussFit(br)

cly =[0.28,0.08,0.02;  0.64,0.08,0.23;  0.8,0.38,0.03; 0.94,0.66,0.24];
clx =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
cl = cly; cond = ["2.5","5","10","15"]; SF = ["0.5SF","1SF","2SF"];
%SUB =["LKC","PSC"];
SUB =["all"];

for sub = 1:1
    for sf = 1:3
        figure();
        for ecc = 1:4
            x = br.storeLags(ecc,:,sf,sub)';
            y = br.YMean(ecc,:,sf,sub)';
            %y = br.XMean(ecc,:,sf,sub)';
            %y = br.phaseMean(ecc,:,sf)'; % for phase check
            
            allci =[];allg = [];
            
            for ite = 1:10
                [myFit G] =  fit(x,y,'gauss1');
                g = [G.rsquare,G.sse];
                allg = [allg; g];
                allcf(ite,:) = coeffvalues(myFit);
                allci = [allci;confint(myFit)];
            end
            
            for i = 1:10
                if allg(i,1)== max(g(:,1))
                    gf(:,:,ecc,sf,sub) = allg(ite,:);
                    cf = allcf(i,:);
                    ci = allci(2*ite-1:2*ite,:);
                end
            end
            
            a = cf(1); b = cf(2); c = cf(3);
            ye = a.*exp(-((x-b)/c).^2);
            yemax = max(ye);
            amp(1,ecc,sf,sub) = yemax;
            
            for i = 1:length(x)
                if ye(i) == yemax
                    xemax = x(i);
                end
            end
            
            %% check
            sigmafit = c./sqrt(2); % sigma calculated from c
            intenorm = trapz(x,y);
            sigmanorm = intenorm/(a.*sqrt(2*pi)); % sigma calculated from a
            FWHM = abs(2*sqrt(log(2)) * c); % FWHM from c
            
            hpkx = 1;
            hmax = max(y).*0.5;
            
            for i = 1:length(y)
                if  (abs(hmax-y(i)))<=(abs(hmax-y(hpkx)))
                    hpkx = i;
                end
                if y(i) == max(y)
                    pkx = i;
                end
            end
            fwhm = 2*abs(x(hpkx)-x(pkx));
            check(:,:,ecc,sf,sub) = [sigmafit, sigmanorm, FWHM, fwhm,yemax,xemax];
            gaussFit(:,:,ecc,sf,sub) = [cf;ci];
            
            subplot(1,4,ecc);
            data = plot(x,y,'Color', cl(ecc,:),'LineWidth',2);hold on;
            estimate = plot(x,ye,':','Color', cl(ecc,:),'LineWidth',2);hold on;
            ylim([-0.05, 0.4]);
            %ylim([-0.1, 0.1]);
            legend([data,estimate],'data','estimate');
            title( strcat(SUB(sub)," ",SF(sf)," Ecc:  ", cond(ecc),"  s=",sprintf('%.6f',sigmafit),...
                "  r^2=",sprintf('%.6f',gf(1,1,ecc,sf))));
        end
        
    end
end


end


