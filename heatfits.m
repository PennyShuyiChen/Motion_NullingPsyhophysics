%% the heatmap 
function [dog1Fit,dgf,dcheck,damp] = heatfits(br)
% heatmap(xvalues,yvalues,cdata)
yecc = 1:400; zeroInd = 73;
%% y raw traces 
figure(1);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.yTraces(:,zeroInd:end,1));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(2);h = heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.yTraces(:,zeroInd:end,2));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};


figure(3);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.yTraces(:,zeroInd:end,3));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

%% x raw traces 
figure(1);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.xTraces(:,zeroInd:end,1));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(2);h = heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.xTraces(:,zeroInd:end,2));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};


figure(3);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.xTraces(:,zeroInd:end,3));colormap(hot(256));grid off;
caxis([-0.25,0.45]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

%% phase raw traces - eccentricity dependence 
%caxis([-0.2,0.25]); % used for cfx 
figure(1);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,1));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(2);h = heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,2));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(3);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,3));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

%% phase raw traces - spatial frequency dependence 
%caxis([-0.2,0.25]); % used for cfx 
figure(1);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,1));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(2);h = heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,2));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(3);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,3));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};

figure(4);h=heatmap(rd.storeLags(1,zeroInd:end),yecc,rd.pTraces(:,zeroInd:end,3));colormap(hot(256));grid off;
caxis([-0.18,0.23]);xlabel('Lag (s)');ylabel('Eccentricity (deg)')
idx = ~mod(rd.storeLags(1,zeroInd:end),0.2);h.XDisplayLabels(~idx) = {''};
h.YDisplayLabels(:)= {''};idy = [1,101,201,301]; h.YDisplayLabels(idy)= {2.5,5,10,15};
%% For individual subject plots 
jaw(1:20,:,:) = rd.yTraces(1:20,:,:);
jaw(21:40,:,:) = rd.yTraces(101:120,:,:);
jaw(41:60,:,:) = rd.yTraces(201:220,:,:);
jaw(61:80,:,:) = rd.yTraces(301:320,:,:);
figure(1);imagesc(jaw(:,:,1))
figure(2);imagesc(jaw(:,:,2))
figure(3);imagesc(jaw(:,:,3))

