%CCG  Create a CCG (cross-correlogram) object
%
%
%   XC = CCG(xcdat) creates a tracking data object XC using
%   the data in the xcdat. 
%   The structure xcdat (generally made by 
%   a trkDat object) should contain:
%   a struct of experimental parameters, ep 
%   a struct of cross-correlogram matrices and a time lag vector
%   
%
%
%   The primary thing that CCG does is to compute the tracking cross-correlations
%   and create a CCG object, which can then be used to summarize and plot the ccgs.

classdef ccgSF
   
    %% properties
   properties
      ep        % experimental parameters structure
      xc        % cross correlations
      frameRate
      br        % Bootstrapping Results - mean peaks and lags along with their bootstrapped SEs.
   end
   properties (Constant)
      nBootReps = 100;
   end
    properties (Dependent)
      Lags  % ccg lags in seconds
    end
   
    %% methods
   methods
       % der konstruktor - see also the "create" static method below
      function obj = ccgSF(myfilename)
         if nargin > 0
            load(myfilename, 'ep', 'xc');
            obj.ep = ep;
            obj.xc = xc;
         end
      end
       
      function lags = get.Lags(obj)
          lags = obj.xc.lags./obj.ep.frameRate;
      end
      
      %% money methods 
      
      %% ***** plot raw ccgs with overlaid means *******************************
      function rd = plotraw_sf(obj, type, snum)
          % ccg.plotraw plots individual trial ccgs lightly and overlays
          % the mean ccg - plot per condition
          % optional types:
          %         'c' plots per condition only (default)
          %         's' to plot by subject and condition
          validTypes = 'cs';
          switch nargin  % obj is counted as an input even for the obj.method call format
              case 0
                  error('Must either pass a ccg object, or use ccg.plotraw format');
              case 1
                  type = 'c';
                  sstart = 1;
              case 2
                  if ~any(type == validTypes)
                      error('type argument must be c or s');
                  end
                  sstart = 1;
              case 3
                  sstart = snum;
              otherwise    % nargin > 3
          end
          
          % set up condition sorting - eccentricity 
           condStrVec = [];
          
          for i = 1:size(obj.ep.ecc,1)
              for ii = 1:length(obj.ep.ecc)
                condStrVec(1,(i-1).*length(obj.ep.ecc)+ii)= obj.ep.ecc(i,ii);
              end
          end
          condStrVec = string(condStrVec)
          
          theConds = ["2.5","5","10","15"];
          nConds = length(theConds);
        
          cond = ["2.5","5","10","15"];

          % set up the SF conditions 
          theSFs = unique(obj.ep.freq);
          for i = 1:size(obj.ep.freq,1)
              for ii = 1:length(obj.ep.freq)
                sfVec(1,(i-1).*length(obj.ep.freq)+ii)= obj.ep.freq(i,ii);
              end
          end
             
          for ii = 1:length(sfVec)
              
              if sfVec(ii) == theSFs(1) || sfVec(ii) == theSFs(2) ||sfVec(ii) == theSFs(5)
                  condSFStrVec(ii) = string(0.5);
              elseif sfVec(ii) == theSFs(3) || sfVec(ii) == theSFs(4) ||sfVec(ii) == theSFs(8)
                  condSFStrVec(ii) = string(1);
              else
                  condSFStrVec(ii) = string(2);
              end
          end
           
          theSFcond = unique(condSFStrVec);
          nSFs = length(theSFcond);
          
          % set up the eccentricity condition 
          %theConds = ["39","79","157","236"];
          nConds = length(theConds);
          cond = ["2.5","5","10","15"]
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStims = length(theStims);
          
          % vGab
          gabStrVec = string(obj.ep.vGabs);
          theGabs = unique(gabStrVec, 'rows');
          nGabs = length(theGabs);
          
          % and now for something completely bootstrap
          nBtReps = obj.nBootReps;
          PkMaxs = zeros(nBtReps, 1);
          PkLags = zeros(size(PkMaxs));
          
          if type == 'c' % plot by condition
              nSubs = 1;
         
          end
          
          theLags = obj.Lags;   % compute dependent property Lags (which is in secs)
          zeroInd = ceil(length(theLags)./2);
          
          
          % set up condition sorting
         
          
          if type == 'c'  % plot by condition
              sstop = 1;
          elseif type == 's' && exist('snum', 'var')
              sstop = sstart;
              nSubs = 1;
          else
              sstop = nSubs;
          end
          
          theLags = obj.Lags;
          
          % store the raw data for all five subjects:
          rd.cond = string;
          rd.subj = string;
          rd.storeLags = theLags;%(zeroInd:end); % the v's are universal
          rd.xTraces = zeros(400,length(theLags),nSFs);
          rd.pTraces = zeros(400,length(theLags),nSFs);
          rd.yTraces = zeros(400,length(theLags),nSFs);
  
                
          % current task - skip plotting in a particular subplot if there are no data to plot there
        for ii = 1:nSFs  
          figure;
          count = 1;
          for i=1:nConds
              for j = sstart:sstop
                  
                  if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
                  sf = strcmp(condSFStrVec,theSFcond(ii));
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs& sf';
                  
                  
                  if any(thisCombo)  % ***** try to plot only if there are data to plot
                      therows = find(thisCombo);
                      
                      subplot(nConds, nSubs, i)
                      
                      plot(theLags, obj.xc.xCorY(therows,:)','Color', [1 0.8 1], 'LineWidth', 0.5);
                      
                      hold on
                      if (strcmp(theConds(i), "2.5") ||strcmp(theConds(i), "5")||strcmp(theConds(i), "10")||...
                              strcmp(theConds(i), "15"))
%                           plot(theLags, obj.xc.xCorX(therows,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
%                           ylim([-0.2,0.4]);
%            
%                           rd.xTraces((i-1)*100+1:100*i,:,ii) = obj.xc.xCorX(therows,:);
%                           mx = mean(obj.xc.xCorX(therows,:));
%                           phx = plot(theLags, mx, 'Color', [.2 .2 1], 'LineWidth', 2);
%                           xx = theLags(find(mx==max(mx)));
%                           plot([xx,xx],[-0.2,max(mx)],'--','Color', [.2 .2 1], 'LineWidth', 1);
                          
%                           rd.pTraces((i-1)*100+1:100*i,:,ii) = -obj.xc.xCorDP(therows,:);
%                           plot(theLags, -obj.xc.xCorDP(therows,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
%                           ylim([-0.15,0.15]);
%                           phdp = plot(theLags, -mean(obj.xc.xCorDP(therows,:)), 'Color', [.2 .2 1], 'LineWidth', 2);
%                           hold on;
                      else
%                           plot(theLags, -obj.xc.xCorDP(therows,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
%                           phdp = plot(theLags, -mean(obj.xc.xCorDP(therows,:)), 'Color', [.2 .2 1], 'LineWidth', 2);
                      end
%                       my = mean(obj.xc.xCorY(therows,:));
%                       phy = plot(theLags, my, 'Color', [1 0.2 1], 'LineWidth', 2);
%                       yy = theLags(find(my==max(my)));
%                       rd.yTraces((i-1)*100+1:100*i,:,ii) = obj.xc.xCorY(therows,:);
%                       plot([yy,yy],[-0.2,max(my)],'--','Color', [1 .2 1], 'LineWidth', 1);
%                       hold on;
%                       
%                       hold off
%                       
%                       xlabel('seconds');
%                       ylabel('ccg value');
                      
                      
                      if type == 'c' % plot by condition
                          substr = "all";
                      else
                          substr = theSubs(j);
                      end
                      titstr = strcat("Condition = ", string(theConds(i)), " Subject = ", substr);
                      title(titstr);
                      
                      if strcmp(theConds(i), "2.5") || strcmp(theConds(i), "5") || strcmp(theConds(i), "10")||strcmp(theConds(i), "15")
                          %legend([phx, phy], {'H', 'V'});
                      else
                          %legend([phdp, phy], {'Phase', 'V'});
                      end
                  end % ***** end plotting if block
              end 
              end % end subj loop
          end % end condition loop
      end
      %***** end plotraw method *****************************************************************************
      %% ***** plot the mean ccgs per cond and (opt) subj *******************************
      function plotmeans_sf(obj, type)
          % ccg.plotmeans plots the mean ccgs
          % and also plots a schmear around zero based on the standard error of the +lag noise
          % optional types:
          %         'c' plots per condition only (default)
          %         's' to plot by subject and condition
          validTypes = 'cs';
          switch nargin  % obj is counted as an input even for the obj.method call format
              case 0
                  error('Must either pass a ccg object, or use ccg.plotraw format');
              case 1
                  type = 'c';
              case 2
                  if ~any(type == validTypes)
                      error('type argument must be c or s');
                  end
              otherwise    % nargin > 3
          end
          
          % set up condition sorting
          condStrVec = [];
          for i = 1:size(obj.ep.stimRad,1)
              condStrVec= [condStrVec,obj.ep.stimRad(i,:)];
              condStrVec = string(condStrVec);
          end % string vector for easier comparison
          theConds = ["39","79","157","236"];%unique(condStrVec);
          nConds = length(theConds);
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          if type == 'c' % plot by condition
              nSubs = 1;
          end
          
          theLags = obj.Lags;   % compute dependent property Lags (which is in secs)
          zeroInd = ceil(length(theLags)./2);
          
          figure;
          for i=1:nConds
              for j = 1:nSubs
                  if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs;
                  if any(thisCombo)  % ***** try to plot only if there are data to plot
                      therows = find(thisCombo);
                      subplot(nConds,nSubs,(i-1).*nSubs + j)
                      leMean = nanmean(obj.xc.xCorY(therows,:));
                      leNoise = rms(leMean(1:zeroInd));
                      yy = theLags(find(leMean==max(leMean)));
                      
                      plot([yy,yy],[-0.1,max(leMean)],'--','Color', [1 .2 1], 'LineWidth', 1);hold on;

                      phy = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [1 0.2 1], 'LineWidth', 2);
                      ylim([-0.1,0.3]);
                      % add schmear
                      nScl = 1.5; % how many rms errors to schmear
                      myxlim = xlim;
                      xe = [myxlim, fliplr(myxlim)];
                      ye = [leNoise, leNoise, -leNoise, -leNoise];
                      hold on;
                      sHandle = fill(xe,nScl.*ye,'r');
                      hold off;
                      set(sHandle, 'EdgeColor', 'none');
                      set(sHandle, 'FaceAlpha', 0.1);
                      
                      hold on
                      if strcmp(theConds(i), "39") || strcmp(theConds(i), "79") || strcmp(theConds(i), "157")|| strcmp(theConds(i), "236")
                          leMean = nanmean(obj.xc.xCorX(therows,:));
                          xx = theLags(find(leMean==max(leMean)));
                          plot([xx,xx],[-0.2,max(leMean)],'--','Color', [.2 .2 1], 'LineWidth', 1);
                          phx = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [.2 .2 1], 'LineWidth', 2);
                          
                          
                      else
                          leMean = nanmean(obj.xc.xCorDP(therows,:));
                          phdp = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', [.2 .2 1], 'LineWidth', 2);
                      end
                      hold off
                      
                      xlabel('seconds');
                      ylabel('ccg value');
                      
                      if type == 'c' % plot by condition
                          substr = "all";
                      else
                          substr = theSubs(j);
                      end
                      titstr = strcat("Condition = ", string(theConds(i)), " Subject = ", substr);
                      title(titstr);
                      
                      if strcmp(theConds(i), "39") || strcmp(theConds(i), "79") || strcmp(theConds(i), "157")|| strcmp(theConds(i), "236")
                          legend([phx, phy], {'H', 'V'});
                      else
                          legend([phdp, phy], {'Phase', 'V'});
                      end
                  end % ***** end plotting if block
                  
              end % subjects
          end % conditions
      end
      %***** end plotmeans method *****************************************************************************
      
      %% ***** bootstrap the peak and lag errors per cond and (opt) subj *******************************
      function br = bootanal_sf(obj, type)  %**** oh what oh what should yon function return?
          % ccg.bootanal bootstraps distributions on the peaks and lags
          % perhaps it optionally makes a bar plot or something
          % optional types:
          %         'c' plots per condition only (default)
          %         's' to plot by subject and condition
          validTypes = 'cs';
          switch nargin  % obj is counted as an input even for the obj.method call format
              case 0
                  error('Must either pass a ccg object, or use ccg.plotraw format');
              case 1
                  type = 'c';
              case 2
                  if ~any(type == validTypes)
                      error('type argument must be c or s');
                  end
              otherwise    % nargin > 3
          end
                    
          
          % set up condition sorting - eccentricity 
           condStrVec = [];
          
          for i = 1:size(obj.ep.ecc,1)
              for ii = 1:length(obj.ep.ecc)
                condStrVec(1,(i-1).*length(obj.ep.ecc)+ii)= obj.ep.ecc(i,ii);
              end
          end
          condStrVec = string(condStrVec)
          
          theConds = ["2.5","5","10","15"];
          nConds = length(theConds);
        
          cond = ["2.5","5","10","15"];

          % set up the SF conditions 
          theSFs = unique(obj.ep.freq);
          for i = 1:size(obj.ep.freq,1)
              for ii = 1:length(obj.ep.freq)
                sfVec(1,(i-1).*length(obj.ep.freq)+ii)= obj.ep.freq(i,ii);
              end
          end
             
          for ii = 1:length(sfVec)
              
              if sfVec(ii) == theSFs(1) || sfVec(ii) == theSFs(2) ||sfVec(ii) == theSFs(5)
                  condSFStrVec(ii) = string(0.5);
              elseif sfVec(ii) == theSFs(3) || sfVec(ii) == theSFs(4) ||sfVec(ii) == theSFs(8)
                  condSFStrVec(ii) = string(1);
              else
                  condSFStrVec(ii) = string(2);
              end
          end
           
          theSFcond = unique(condSFStrVec);
          nSFs = length(theSFcond);
          
          % set up the eccentricity condition 
          %theConds = ["39","79","157","236"];
          nConds = length(theConds);
          cond = ["2.5","5","10","15"]
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStims = length(theStims);
          
          % vGab
          gabStrVec = string(obj.ep.vGabs);
          theGabs = unique(gabStrVec, 'rows');
          nGabs = length(theGabs);
          
          % and now for something completely bootstrap
          nBtReps = obj.nBootReps;
          PkMaxs = zeros(nBtReps, 1);
          PkLags = zeros(size(PkMaxs));
          
          if type == 'c' % plot by condition
              nSubs = 1;
         
          end
          
          theLags = obj.Lags;   % compute dependent property Lags (which is in secs)
          zeroInd = ceil(length(theLags)./2);
          pLags = theLags(zeroInd:end);
          
          % set up struct elements for saving boot results
          br.cond = string;
          br.subj = string;
          br.vLags = zeros(nConds * nSubs, 3); % the v's are universal
          br.vPks = br.vLags;
          br.vLagErr = br.vLags;
          br.vPksErr = br.vLags;
          br.hLags = br.vLags;                  % the h's are h for zebra,and dPhi for the others
          br.hPks = br.vLags;
          br.hLagErr = br.vLags;
          br.hPksErr = br.vLags;
          
          loopCount = 0;
          figHans = zeros(nConds.*nSubs, 1);
          %figure 
          nConds 
          cond = ["2.5","5","10","15"];
          
          for i=1:nConds
              for j = 1:nSubs
                  figure(j);
                  for e = 1:nStims
                   %figure(e);
                   if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
              for ii = 1:nSFs
                  
                  loopCount = loopCount + 1;
                 
                 %thisCombo = (strcmp(condStrVec, theConds(i))) &
                 %theseSubs; % zeb conditions etc 
                  sf = strcmp(condSFStrVec,theSFcond(ii));
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs & sf' &(strcmp(stimStrVec, theStims(e)));
                 subplot(1,3,ii);
                  if any(thisCombo)  % ***** is there anything to bootstrap?
                      
                      deezRows = find(thisCombo); % if there is, den deez r da rows ta uze
                      
                      %thisIter = (i-1).*nSubs + j;
                      thisIter = i;
                      if  (strcmp(theConds(i), "2.5") || strcmp(theConds(i), "5") || strcmp(theConds(i), "10")||...
                                  strcmp(theConds(i), "15"))
                      %                       subplot(nConds,nSubs,thisIter)
                      
                      %figHans(loopCount) = figure;  % why is this here?  Possibly an ouput arg at some point
                      %figHans(loopCount) = nConds;
                      
                      
                      if strcmp(theStims(e),"cfx") && strcmp(theGabs,"0")
                            leMean = nanmean(-obj.xc.xCorDP(thisCombo,:));
                            vrmtrials = -obj.xc.xCorDP(thisCombo,:);
                      else
                         leMean = nanmean(obj.xc.xCorY(thisCombo,:));
                         vrmtrials = obj.xc.xCorY(thisCombo,:);
                         leMeanD = nanmean(obj.xc.xCorYD(thisCombo,:)); 
                      end 
                      leNoise = rms(leMean(1:zeroInd));
                      leSNR = max(leMean)./leNoise;
                      
                      if leSNR >= 1.5
                          cly =[0.28,0.08,0.02;  0.64,0.08,0.23;  0.8,0.38,0.03; 0.94,0.66,0.24];
                          %*****    dew da bewtstrappin  *********
                          %*****    first for the y - direction, which is univesal ********
                          if strcmp(theStims(e),"zeb")||strcmp(theStims(e),"apt") || (strcmp(theStims(e),"cfx") && strcmp(theGabs,"1"))||(strcmp(theStims(e),"cfc") && strcmp(theGabs,"1"))
                              YMean = leMean((zeroInd:end));
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorY(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.01, [1 0.5 1]);
                              end % end lil' boot loop
%                               hold on;
%                               
%                               phy(i) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10,'Color',cly(i,:));
%                               legend(phy,cond);
%                               ylim([-0.05,0.25]);
%                               yticks(-0.05:0.05:0.25)
%                               hold off

                               for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorYD(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.01, [1 0.5 1]);
                              end % end lil' boot loop
%                               hold on;
%                               phyd(i) = plot(theLags(zeroInd:end-1), leMeanD(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10,'Color',cly(i,:));
%                               legend(phyd,cond);
%                               ylim([-0.05,0.05]);
%                               hold off
                               
                          elseif (strcmp(theStims(e),"cfx") && strcmp(theGabs,"0"))
                              leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
            
                              %*****   for the phase, if in a cuttlefish condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                  % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              
%                               phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',cly(i,:));
%                               legend(phdp,cond);
%                               ylim([-0.1,0.1]);
%                               hold off
                          end
                          
                          % and store the lags, the peaks, and their error estimates
                          br.vLags(thisIter,ii,j) = mean(PkLags);
                          br.vPks(thisIter,ii,j) = mean(PkMaxs);
                          br.vLagErr(thisIter,ii,j) = std(PkLags);
                          br.vPksErr(thisIter,ii,j) = std(PkMaxs);
                          
                          % add schmear for ESP noise along the bottom
                          nScl = 1.5; % how many rms errors to schmear
                          myxlim = xlim;
                          xe = [myxlim, fliplr(myxlim)];
                          ye = [leNoise, leNoise, -leNoise, -leNoise];
                          hold on;
                          sHandle = fill(xe,nScl.*ye,'r');
                          hold off;
                          set(sHandle, 'EdgeColor', 'none');
                          set(sHandle, 'FaceAlpha', 0.1);
                          
%                          yticks(-0.1:0.1:0.4);
%                          legend(phy,cond);
                          hold on
                          
                          if strcmp(theStims(e),"zeb") ||strcmp(theStims(e),"apt")|| (strcmp(theStims(e),"cfx") && strcmp(theGabs,"0"))
                              
                              leMean = nanmean(obj.xc.xCorX(thisCombo,:));
                              XMean = leMean((zeroInd:end));
                              %*****    dew moar bewtstrappin  *********
                              %*****    now for the x - direction, if in the zeb condition ********
                              clx =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorX(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.02, clx(i,:));
                              end % end lil' boot loop
                              hold on
                              
                              phx(i) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', clx(i,:), 'LineWidth', 2);
                              errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                  'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                  'LineWidth', 1.5, 'CapSize', 10,'Color',clx(i,:));   
                              legend(phx,cond);
                             ylim([-0.05,0.25]);
                             yticks(-0.05:0.05:0.25)
                              hold off
                              
                          elseif  (strcmp(theStims(e),"cfx") && strcmp(theGabs,"1"))||(strcmp(theStims(e),"cfc") && strcmp(theGabs,"1"))
                              
                              clp =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
                              leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
                              leMeanS = nanmean(obj.xc.xCorS(thisCombo,:));
                              leMeanTVSA = nanmean(obj.xc.xCorTVSA(thisCombo,:));
                              phaseMean = -leMean((zeroInd:end));
                              %*****    dew eevin moar bewtstrappin  *********
                              %*****    or for the phase, if in a cuttlefish condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                 % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              % **** plot phdp ****
%                               phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', clp(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',clp(i,:));
%                               legend(phdp,cond);
%                               ylim([-0.03,0.065]);
%                               yticks(-0.03:0.03:0.065);
%                               hold off
%                               for k = 1:nBtReps
%                                   deezBootRows = randsample(deezRows, length(deezRows), true);
%                                   disBootMean = nanmean(obj.xc.xCorS(deezBootRows,zeroInd:end));
%                                   [PkMaxs(k), li] = max(-disBootMean);
%                                   PkLags(k) = pLags(li);
%                                  % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
%                               end % end lil' boot loop
%                               hold on
                              
%                               phsp(i) = plot(theLags(zeroInd:end), -leMeanS(zeroInd:end), 'Color', clp(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',clp(i,:));
%                               legend(phsp,cond);
                              
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorTVSA(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                 % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
%                               hold on
%                               phtvsa(i) = plot(theLags(zeroInd:end), -leMeanTVSA(zeroInd:end), 'Color', clp(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',clp(i,:));
%                               legend(phtvsa,cond);
%                               ylim([-0.02,0.02]);%yticks(-0.05:0.01:0.05);
%                               hold off                               
                              
                              
                          end % if zeb or cuttlefish or indie1D
                          hold off
                          
                          storeLags = theLags(zeroInd:end);
                          xLabStr = 'Lag (seconds)';
                          yLabStr = 'CCG value';
                          
                          
                          if type == 'c'  && (numel(theseSubs) > 1) % more that one subj in data file
                              substr = "all";
                          else
                              substr = theSubs(j);      % plotting by condition, or single S in file
                          end
                          
                          %titsf = string(theSFcond(ii));
                          titsf = ["Low","Mid","High"];
                          titstr = strcat(titsf(ii), " SF ");
                          %titstr = strcat("Stimulus Acceleration ",titsf, " SF ", theStims(e)," ",substr);

                          %formatFigure(xLabStr, yLabStr, titstr, 0, 0, 32, 32);
                          formatFigure(xLabStr, yLabStr, titstr);
                          
                          
                          xlim([0,1.2]);
                          set(gcf, 'Name', strcat(string(theConds(i)), substr));
                          
                          % and store the horizontal or phase results
                          br.hLags(thisIter,ii,j) = mean(PkLags);
                          br.hPks(thisIter,ii,j) = mean(PkMaxs);
                          br.hLagErr(thisIter,ii,j) = std(PkLags);
                          br.hPksErr(thisIter,ii,j) = std(PkMaxs);
                          
                          br.cond(thisIter, :,j) = theConds(i);
                          br.subj(thisIter, :,j) = substr;
                          
                          br.figHans = figHans;
                          
                          %br.phx(thisIter, :) = phx(i);
                          %br.phy(thisIter, :) = phy(i);
                          
                          diffY = diff(YMean,1);
                          %br.YMean(thisIter, :,ii,j) = YMean;
                          br.XMean(thisIter, :,ii,j) = XMean;
%                          br.phaseMean(thisIter, :,ii,j) = phaseMean;
                          br.storeLags(thisIter, :,ii,j) = storeLags;
                          br.diffY(thisIter, :,ii,j) = diffY;
                          
                      else % SNR critereon block
                          
                          br.hLags(thisIter) = NaN;
                          br.hPks(thisIter) = NaN;
                          br.hLagErr(thisIter) = NaN;
                          br.hPksErr(thisIter) = NaN;
                          
                          br.cond(thisIter, :) = theConds(i);
                          br.subj(thisIter, :) = substr;
                          
                          br.figHans = figHans;
                          
                         
                      end % end of SNR critereon block
                      end  
                  end % ***** end if(data in this condition) block
                  
              end 
              end 
          end% subjects
          end % conditions
          hold on;
          %leg = legend([phy,phdp],[cond,cond]);
          %leg = legend([phy,phx],[cond,cond]);
          legend boxoff;
%          title(leg,'Eccentricity (deg)');
      
          hold on;
      end
      %***** end bootanal method *****************************************************************************
    
      

      
      %% ***** bootstrap the peak and lag errors per cond and (opt) subj *******************************
      function br = bootanal_ecc(obj, type)  %**** oh what oh what should yon function return?
          % ccg.bootanal bootstraps distributions on the peaks and lags
          % perhaps it optionally makes a bar plot or something
          % optional types:
          %         'c' plots per condition only (default)
          %         's' to plot by subject and condition
          validTypes = 'cs';
          switch nargin  % obj is counted as an input even for the obj.method call format
              case 0
                  error('Must either pass a ccg object, or use ccg.plotraw format');
              case 1
                  type = 'c';
              case 2
                  if ~any(type == validTypes)
                      error('type argument must be c or s');
                  end
              otherwise    % nargin > 3
          end
          
          % set up condition sorting - eccentricity 
           condStrVec = [];
          
          for i = 1:size(obj.ep.ecc,1)
              for ii = 1:length(obj.ep.ecc)
                condStrVec(1,(i-1).*length(obj.ep.ecc)+ii)= obj.ep.ecc(i,ii);
              end
          end
          condStrVec = string(condStrVec)
          
          theConds = ["2.5","5","10","15"];
          nConds = length(theConds);
        
          cond = ["2.5","5","10","15"];

          % set up the SF conditions 
          theSFs = unique(obj.ep.freq);
          for i = 1:size(obj.ep.freq,1)
              for ii = 1:length(obj.ep.freq)
                sfVec(1,(i-1).*length(obj.ep.freq)+ii)= obj.ep.freq(i,ii);
              end
          end
             
          for ii = 1:length(sfVec)
              
              if sfVec(ii) == theSFs(1) || sfVec(ii) == theSFs(2) ||sfVec(ii) == theSFs(5)
                  condSFStrVec(ii) = string(0.5);
              elseif sfVec(ii) == theSFs(3) || sfVec(ii) == theSFs(4) ||sfVec(ii) == theSFs(8)
                  condSFStrVec(ii) = string(1);
              else
                  condSFStrVec(ii) = string(2);
              end
          end
           
          theSFcond = unique(condSFStrVec);
          nSFs = length(theSFcond);
          
          % set up the eccentricity condition 
          %theConds = ["39","79","157","236"];
          nConds = length(theConds);
          cond = ["2.5","5","10","15"]
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStims = length(theStims);
          
          % vGab
          gabStrVec = string(obj.ep.vGabs);
          theGabs = unique(gabStrVec, 'rows');
          nGabs = length(theGabs);
          
          % and now for something completely bootstrap
          nBtReps = obj.nBootReps;
          PkMaxs = zeros(nBtReps, 1);
          PkLags = zeros(size(PkMaxs));
          
          if type == 'c' % plot by condition
              nSubs = 1;
         
          end
          
          theLags = obj.Lags;   % compute dependent property Lags (which is in secs)
          zeroInd = ceil(length(theLags)./2);
          pLags = theLags(zeroInd:end);
          
          % set up struct elements for saving boot results
          br.cond = string;
          br.subj = string;
          br.vLags = zeros(nConds * nSubs, 3); % the v's are universal
          br.vPks = br.vLags;
          br.vLagErr = br.vLags;
          br.vPksErr = br.vLags;
          br.hLags = br.vLags;                  % the h's are h for zebra,and dPhi for the others
          br.hPks = br.vLags;
          br.hLagErr = br.vLags;
          br.hPksErr = br.vLags;
          
          loopCount = 0;
          figHans = zeros(nConds.*nSubs, 1);
          %figure 
          nConds 
          cond = ["2.5","5","10","15"];
          
          for i=1:nConds
              for j = 1:nSubs
                  figure(j);
                  for e = 1:nStims
                   %figure(e);
                   if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
              for ii = 1:nSFs
                  
                  loopCount = loopCount + 1;
                 
                 %thisCombo = (strcmp(condStrVec, theConds(i))) &
                 %theseSubs; % zeb conditions etc 
                  sf = strcmp(condSFStrVec,theSFcond(ii));
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs & sf' &(strcmp(stimStrVec, theStims(e)));
                 subplot(1,4,i);
                  if any(thisCombo)  % ***** is there anything to bootstrap?
                      
                      deezRows = find(thisCombo); % if there is, den deez r da rows ta uze
                      
                      %thisIter = (i-1).*nSubs + j;
                      thisIter = i;
                      if  (strcmp(theConds(i), "2.5") || strcmp(theConds(i), "5") || strcmp(theConds(i), "10")||...
                                  strcmp(theConds(i), "15"))
                      %                       subplot(nConds,nSubs,thisIter)
                      
                      %figHans(loopCount) = figure;  % why is this here?  Possibly an ouput arg at some point
                      %figHans(loopCount) = nConds;
                      
                      
                      if strcmp(theStims(e),"cfx") && strcmp(theGabs,"0")
                            leMean = nanmean(-obj.xc.xCorDP(thisCombo,:));
                      else
                         leMean = nanmean(obj.xc.xCorY(thisCombo,:)); 
                      end 
                      leNoise = rms(leMean(1:zeroInd));
                      leSNR = max(leMean)./leNoise;
                      YMean = leMean((zeroInd:end));
                      if leSNR >= 1.5
                          cly(:,:,ii) =[0.28,0.08,0.02,1-(ii-1)*0.3;  0.64,0.08,0.23,1-(ii-1)*0.3;  0.8,0.38,0.03,1-(ii-1)*0.3; 0.94,0.66,0.24,1-(ii-1)*0.4];
                          
                          %*****    dew da bewtstrappin  *********
                          %*****    first for the y - direction, which is univesal ********
                          if strcmp(theStims(e),"zeb") ||strcmp(theStims(e),"apt")|| (strcmp(theStims(e),"cfx") && strcmp(theGabs,"1"))||(strcmp(theStims(e),"cfc") && strcmp(theGabs,"1"))
                              YMean = leMean((zeroInd:end));
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorY(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.01, [1 0.5 1]);
                              end % end lil' boot loop
%                               hold on;
%                               
%                               phy(i,:,ii) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', cly(i,:,ii), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10,'Color',cly(i,:,ii));
%                               
%                               ylim([-0.1,0.4]);
%                               hold off
                              
                          elseif (strcmp(theStims(e),"cfx") && strcmp(theGabs,"0"))
                              leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
            
                              %*****   for the phase, if in a cuttlefish condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                  % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              
%                               phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',cly(i,:));
%                               legend(phdp,cond);
%                               ylim([-0.1,0.1]);
%                               hold off
                          end
                          
                          % and store the lags, the peaks, and their error estimates
                          br.vLags(thisIter,ii,j) = mean(PkLags);
                          br.vPks(thisIter,ii,j) = mean(PkMaxs);
                          br.vLagErr(thisIter,ii,j) = std(PkLags);
                          br.vPksErr(thisIter,ii,j) = std(PkMaxs);
                          
                          % add schmear for ESP noise along the bottom
                          nScl = 1.5; % how many rms errors to schmear
                          myxlim = xlim;
                          xe = [myxlim, fliplr(myxlim)];
                          ye = [leNoise, leNoise, -leNoise, -leNoise];
                          hold on;
                          sHandle = fill(xe,nScl.*ye,'r');
                          hold off;
                          set(sHandle, 'EdgeColor', 'none');
                          set(sHandle, 'FaceAlpha', 0.1);
                          
                         %yticks(-0.1:0.1:0.4);
                         %legend(phy,cond);
                          hold on
                          
                          if strcmp(theStims(e),"zeb") ||strcmp(theStims(e),"apt")|| (strcmp(theStims(e),"cfx") && strcmp(theGabs,"0"))
                              
                              leMean = nanmean(obj.xc.xCorX(thisCombo,:));
               
                              XMean = leMean((zeroInd:end));
                              %*****    dew moar bewtstrappin  *********
                              %*****    now for the x - direction, if in the zeb condition ********
                              clx(:,:,ii) =[0.5,0.05,0.93,1-(ii-1)*0.3;  0.4,0.33,0.93,1-(ii-1)*0.3;  0.3,0.52,0.93,1-(ii-1)*0.3; 0.2,0.95,0.93,1-(ii-1)*0.3];
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorX(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.02, clx(i,:));
                              end % end lil' boot loop
                              hold on
                              
%                               phx(i) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', clx(i,:,ii), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10,'Color',clx(i,:,ii));   
%                               legend(phx,cond);
%                              ylim([-0.1,0.4]);
%                               hold off
                              
                          elseif  (strcmp(theStims(e),"cfx") && strcmp(theGabs,"1"))||(strcmp(theStims(e),"cfc") && strcmp(theGabs,"1"))
                              
                              clp(:,:,ii) =[0.5,0.05,0.93,1-(ii-1)*0.3;  0.4,0.33,0.93,1-(ii-1)*0.3;  0.3,0.52,0.93,1-(ii-1)*0.3; 0.2,0.95,0.93,1-(ii-1)*0.3];
                              leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
                              leMeanS = nanmean(obj.xc.xCorS(thisCombo,:));
                              phaseMean = -leMean((zeroInd:end));
                              %*****    dew eevin moar bewtstrappin  *********
                              %*****    or for the phase, if in a cuttlefish condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                 % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              
                              phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', clp(i,:,ii), 'LineWidth', 2);
                              errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                  'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                  'LineWidth', 1.5, 'CapSize', 10, 'Color',clp(i,:,ii));
                              legend(phdp,cond);
                              ylim([-0.03,0.07]);
                              yticks(-0.03:0.03:0.07);
                              hold off
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorS(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                 % ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
%                               phsp(i) = plot(theLags(zeroInd:end), -leMeanS(zeroInd:end), 'Color', clp(i,:,ii), 'LineWidth', 2);
%                               errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
%                                   'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
%                                   'LineWidth', 1.5, 'CapSize', 10, 'Color',clp(i,:,ii));
%                               legend(phsp,cond);
%                               ylim([-0.5,0.5]);yticks(-0.5:0.1:0.5);
%                               hold off                               
                              
                              
                          end % if zeb or cuttlefish or indie1D
                          hold off
                          
                           %ylim([-0.1,0.4]);
                           %yticks(-0.1:0.1:0.4);
                          storeLags = theLags(zeroInd:end);
                          xLabStr = 'Lag (seconds)';
                          yLabStr = 'CCG value';
                          
                          
                          if type == 'c'  && (numel(theseSubs) > 1) % more that one subj in data file
                              substr = "all";
                          else
                              substr = theSubs(j);      % plotting by condition, or single S in file
                          end
                          
                          titecc = string(cond(i));
                          titstr = strcat("Y direction ",titecc, " deg ", theStims(e)," ",substr);
                          %titstr = strcat("Stimulus Acceleration ",titsf, " SF ", theStims(e)," ",substr);

                          %formatFigure(xLabStr, yLabStr, titstr, 0, 0, 32, 32);
                          formatFigure(xLabStr, yLabStr, titstr);
                          
                          
                          xlim([0,1.2]);
                          set(gcf, 'Name', strcat(string(theConds(i)), substr));
                          
                          % and store the horizontal or phase results
                          br.hLags(thisIter,ii,j) = mean(PkLags);
                          br.hPks(thisIter,ii,j) = mean(PkMaxs);
                          br.hLagErr(thisIter,ii,j) = std(PkLags);
                          br.hPksErr(thisIter,ii,j) = std(PkMaxs);
                          
                          br.cond(thisIter, :,j) = theConds(i);
                          br.subj(thisIter, :,j) = substr;
                          
                          br.figHans = figHans;
                          
                          %br.phx(thisIter, :) = phx(i);
                          %br.phy(thisIter, :) = phy(i);
                          
                          diffY = diff(YMean,1);
                          br.YMean(thisIter, :,ii,j) = YMean;
                          br.XMean(thisIter, :,ii,j) = XMean;
                          %br.phaseMean(thisIter, :,ii,j) = phaseMean;
                          %br.storeLags(thisIter, :,ii,j) = storeLags;
                          br.diffY(thisIter, :,ii,j) = diffY;
                          
                      else % SNR critereon block
                          
                          br.hLags(thisIter) = NaN;
                          br.hPks(thisIter) = NaN;
                          br.hLagErr(thisIter) = NaN;
                          br.hPksErr(thisIter) = NaN;
                          
                          br.cond(thisIter, :) = theConds(i);
                          br.subj(thisIter, :) = substr;
                          
                          br.figHans = figHans;
                          
                         
                      end % end of SNR critereon block
                      end  
                  end % ***** end if(data in this condition) block
                  
              end 
              end 
          end% subjects
          end % conditions
          hold on;
          %legend([phy,phdp],[cond,cond]);
          %legend([phy,phx],[cond,cond]);
          
      
          hold on;
      end
      %***** end bootanal method *****************************************************************************
    
      
       %% Acceleration
        function Acc = Acc_sf(obj)
           
            E = [2.5, 5, 10, 15];
            SF = ["0.5SF","1SF","2SF"];
            cond = ["2.5","5","10","15"];
            cly =[0.28,0.08,0.02;  0.64,0.08,0.23;  0.8,0.38,0.03; 0.94,0.66,0.24];
            clp =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
            nConds = size(obj.br.storeLags,1);
            xplot = size(obj.br.storeLags,2);
            
            for sf = 1:3
                figure();
            for i = 1:nConds
                
                subplot(1,4,i)
                pphaseMean(i,sf) = plot(obj.br.storeLags(i,2:xplot,sf),obj.br.phaseMean(i,1:xplot-1,sf),...
                    'Color', clp(i,:), 'LineWidth', 2); hold on;
                phaseMax(i,sf) = max(obj.br.phaseMean(i,1:xplot-1,sf));
                pdiffY(i,sf) = plot(obj.br.storeLags(i,2:xplot,sf),obj.br.diffY(i,:,sf),...
                    'Color', cly(i,:), 'LineWidth', 2); hold on;
                diffYMax(i,sf) = max(obj.br.diffY(i,:,sf));
                scalar(i,sf) = phaseMax(i,sf) ./ diffYMax(i,sf);
                sdiffY(i,:,sf) = obj.br.diffY(i,:,sf) .* scalar(i,sf);
                psdiffY(i,sf) = plot(obj.br.storeLags(i,2:xplot,sf),sdiffY(i,:,sf),...
                    '--','Color', cly(i,:), 'LineWidth', 1.5); hold on;
                
                [rho(i,sf),pval(i,sf)] = corr(obj.br.phaseMean(i,1:xplot-1,sf)',sdiffY(i,:,sf)','Tail','right','Rows','pairwise');
                
                %ylim([-0.1, 0.1]);
                 legend([pphaseMean(i,sf),pdiffY(i,sf),psdiffY(i,sf)],'phase','dY','scaled-dY');
                title( strcat(SF(sf),"- "," Eccentricity:  ", cond(i), "  rho=",sprintf('%.6f',rho(i,sf)),...
                    "  p=",sprintf('%.6f',pval(i,sf)), "  scalar=",sprintf('%.6f',scalar(i,sf))));
                
                ylabel('ccg value');
                xlabel('Time (s)');
            end
            
            end 
            
            for sf = 1:3
                lm{sf} = fitlm(E,scalar(:,sf)');
                    
            end
            
            Acc.lm = lm;
             Acc.phaseMax = phaseMax;
             Acc.diffYMax = diffYMax;
             Acc.scalar   = scalar;
             Acc.sdiffY   = sdiffY;
          
             Acc.rho = rho;
             Acc.pval = pval;
             
             Enew = [0,E]'; 
             figure();
             for sf = 1:3
                subplot(1,3,sf);
                ps(sf)  = plot(E,scalar(:,sf),'LineWidth',2); hold on;
                plm(sf) = plot(E,Acc.lm{1,sf}.Fitted,'r*','LineWidth',2); hold on;
                syms f(x)
              f(x) = Acc.lm{1,sf}.Coefficients.Estimate(1) + x .* Acc.lm{1,sf}.Coefficients.Estimate(2);
              pfunc(sf) = fplot(f,[0,15],'LineWidth',2); hold on;
              
              [ypred,yci] = predict(Acc.lm{1,sf}, Enew); 
              Acc.ci(:,:,sf) = yci
               plot(Enew, Acc.ci(:,:,sf), '--r');
               %ylim([0,10]);
               xlabel('Eccentricity(deg)'); ylabel('Parameter values');
               title(  strcat(SF(sf),"Scalar ","  R^2=",sprintf('%.6f',Acc.lm{1,sf}.Rsquared.Ordinary),...
                  "  p=",sprintf('%.6f',Acc.lm{1,sf}.Coefficients.pValue(2))));
              legend([ps(sf),plm(sf),pfunc(sf)],'scalar','fitted scalar','function','Location','Northwest');
             end 
             
              
%             
        end
      %% display
      % display object
      function disp(obj)
         disp('experimental parameters (ep): ');
         disp(obj.ep);
         disp('ccgs (xc): ');
         disp(obj.xc);
      end
   end % end regular methods

    %% static constructor method
   methods (Static)
      function obj = create % this is a function you can call from command line, as in myTD = CCG.createObj
            [myfilename, pathname] = uigetfile('*.mat', 'Pick the xc mat file made by trkDat'); 
            myfullfile = fullfile(pathname, myfilename);
            obj = ccg(myfullfile);
      end
   end % end static methods
   
end % end classdef