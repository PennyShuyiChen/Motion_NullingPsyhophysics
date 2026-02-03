%CCG  Create a CCG (cross-correlogram) object
%    PSC 
%    Last modified: 12202020 =) 
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

classdef ccgLOOP
   
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
      function obj = ccgLOOP(myfilename)
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
      function plotraw_lp(obj, type, snum)
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
          
          % set up condition sorting
          condStrVec = [];
          for i = 1:size(obj.ep.stimRad,1)
              condStrVec= [condStrVec,obj.ep.stimRad(i,:)];
              condStrVec = string(condStrVec);
          end
          %condStrVec = string(obj.ep.stimRad); % string vector for easier comparison
          %theConds = unique(condStrVec);
          theConds = ["39","79","157","236"];
          nConds = length(theConds);
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStim = length(theStims);
          
          % vGab
          theGabs = string(obj.ep.vGabs); 
          
          if type == 'c'  % plot by condition
              sstop = 1;
          elseif type == 's' && exist('snum', 'var')
              sstop = sstart;
              nSubs = 1;
          else
              sstop = nSubs;
          end
          
          theLags = obj.Lags;
          
          % current task - skip plotting in a particular subplot if there are no data to plot there
          figure;
          count = 0;
          for i=1:nConds
              for j = sstart:sstop
                  count = count + 1;
                  if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs;
                  
                  if any(thisCombo)  % ***** try to plot only if there are data to plot                   
                      therows = find(thisCombo);
                      
                      %subplot(nConds, nSubs, (i-1).*nSubs + count)
                      subplot(nConds, nSubs, i)
                      
                      if (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79"))
                      
                             % *** Y-direction ***
                             if  strcmp(theStims,"zeb") || (strcmp(theStims,"cfx") && strcmp(theGabs,"1"))
                                 plot(theLags, obj.xc.xCorY(therows,:)','Color', [1 0.8 1], 'LineWidth', 0.5);
                                 hold on;
                                 ylim([-0.2,0.4]);
                                 my = mean(obj.xc.xCorY(therows,:));
                                 phy = plot(theLags, my, 'Color', [1 0.2 1], 'LineWidth', 2);
                                 yy = theLags(find(my==max(my)));
                                 plot([yy,yy],[-0.2,max(my)],'--','Color', [1 .2 1], 'LineWidth', 1);
                                 hold on;
                             elseif (strcmp(theStims,"cfx") && strcmp(theGabs,"0"))
                                 plot(theLags, -obj.xc.xCorDP(therows,:)','Color', [1 0.8 1], 'LineWidth', 0.5);
                                 hold on;
                                 ylim([-0.1,0.1]);
                                 mdp = -mean(obj.xc.xCorDP(therows,:));
                                 phdp = plot(theLags, mdp, 'Color', [1 0.2 1], 'LineWidth', 2);
                                 dp = theLags(find(mdp==max(mdp)));
                                 plot([dp,dp],[-0.2,max(mdp)],'--','Color', [1 0.2 1], 'LineWidth', 1);
                                 hold on;
                             end
                             
                             % *** Now for the X-direction *** 
                  
                      if  strcmp(theStims,"zeb") || (strcmp(theStims,"cfx") && strcmp(theGabs,"0"))
                          plot(theLags, obj.xc.xCorX(therows,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
                          hold on;
                          ylim([-0.2,0.4]);
                          mx = mean(obj.xc.xCorX(therows,:));
                          phx = plot(theLags, mx, 'Color', [.2 .2 1], 'LineWidth', 2);
                          xx = theLags(find(mx==max(mx)));
                          plot([xx,xx],[-0.2,max(mx)],'--','Color', [.2 .2 1], 'LineWidth', 1);
                          hold on;
                      elseif  strcmp(theStims,"cfx")&& strcmp(theGabs,"1")
                          plot(theLags, -obj.xc.xCorDP(therows,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
                          hold on;
                          ylim([-0.1,0.1]);
                          mdp = -mean(obj.xc.xCorDP(therows,:));
                          phdp = plot(theLags, mdp, 'Color', [.2 .2 1], 'LineWidth', 2);
                          dp = theLags(find(mdp==max(mdp)));
                          plot([dp,dp],[-0.2,max(mdp)],'--','Color', [.2 .2 1], 'LineWidth', 1);
                          hold on;
                      end
                  end
                     
                      hold off
                      
                      xlabel('seconds');
                      ylabel('ccg value');
                      
                      
                      if type == 'c' % plot by condition
                          substr = "all";
                      else
                          substr = theSubs(j);
                      end
                      titstr = strcat("Condition = ", string(theConds(i)), " Subject = ", substr, " Task = ", theStims);
                      title(titstr);
                      
                      if strcmp(theStims,"zeb")
                         legend([phx, phy], {'H', 'V'});
                      elseif  strcmp(theStims,"cfx")
                          if strcmp(theGabs,"0")
                              %legend([phx, phdp], {'X', 'Y-Phase'});
                              %legend(phx,'X');
                              legend(phdp,'Y-Phase');
                          elseif strcmp(theGabs,"1")
                              legend([phdp, phy], {'X-Phase', 'Y'});
                          end 
                      end
                  end % ***** end plotting if block
              
              end % end subj loop
          end % end condition loop
      end
      %***** end plotraw method *****************************************************************************
      %% ***** plot the mean ccgs per cond and (opt) subj *******************************
      function plotmeans_lp(obj, type)
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
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStim = length(theStims);
          
          % vGab
          theGabs = string(obj.ep.vGabs);
          
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
                      if (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79"))
                      therows = find(thisCombo);
                      subplot(nConds,nSubs,(i-1).*nSubs + j)
                      
                      if  strcmp(theStims,"zeb") || (strcmp(theStims,"cfx") && strcmp(theGabs,"1"))
                          
                          leMean = nanmean(obj.xc.xCorY(therows,:));
                          leNoise = rms(leMean(1:zeroInd));
                          yy = theLags(find(leMean==max(leMean)));
                          
                          plot([yy,yy],[-0.1,max(leMean)],'--','Color', [1 .2 1], 'LineWidth', 1);hold on;
                          
                          phy = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [1 0.2 1], 'LineWidth', 2);
                          ylim([-0.1,0.3]);
                      elseif (strcmp(theStims,"cfx") && strcmp(theGabs,"0"))
                           leMean = -nanmean(obj.xc.xCorDP(therows,:));
                           leNoise = rms(leMean(1:zeroInd));
                          dp = theLags(find(leMean==max(leMean)));
                          plot([dp,dp],[-0.2,max(leMean)],'--','Color', [1 .2 1], 'LineWidth', 1);hold on;
                          phdp = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [1 0.2 1], 'LineWidth', 2);
                          hold on;
                          ylim([-0.1,0.1]);
                      end
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
                      if  strcmp(theStims,"zeb") || (strcmp(theStims,"cfx") && strcmp(theGabs,"0"))
%                           leMean = nanmean(obj.xc.xCorX(therows,:));
%                           xx = theLags(find(leMean==max(leMean)));
%                           plot([xx,xx],[-0.2,max(leMean)],'--','Color', [.2 .2 1], 'LineWidth', 1); hold on;
%                           phx = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [.2 .2 1], 'LineWidth', 2);
%                           ylim([-0.1,0.3]);
%                           hold on;
                          
                      elseif  (strcmp(theStims,"cfx") && strcmp(theGabs,"1"))
                          leMean = -nanmean(obj.xc.xCorDP(therows,:));
                          dp = theLags(find(leMean==max(leMean)));
                          plot([dp,dp],[-0.2,max(leMean)],'--','Color', [.2 .2 1], 'LineWidth', 1); hold on;
                          phdp = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [.2 .2 1], 'LineWidth', 2);
                          ylim([-0.1,0.3]);
                          hold on;
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
                      
                      if strcmp(theStims,"zeb")
                         legend([phx, phy], {'H', 'V'});
                      elseif  strcmp(theStims,"cfx")
                          if strcmp(theGabs,"0")
                              %legend([phx, phdp], {'X', 'Y-Phase'});
                              %legend(phx,'X');
                              legend(phdp,'Y-Phase');
                          elseif strcmp(theGabs,"1")
                              legend([phdp, phy], {'X-Phase', 'Y'});
                          end 
                      end
                      
                  end % ***** end plotting if block
                  end 
              end % subjects
          end % conditions
      end
      %***** end plotmeans method *****************************************************************************
      
      %% ***** bootstrap the peak and lag errors per cond and (opt) subj *******************************
      function br = bootanal_lp(obj, type)  %**** oh what oh what should yon function return?
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
          
          % set up condition sorting
          condStrVec = [];
          for i = 1:size(obj.ep.stimRad,1)
              condStrVec= [condStrVec,obj.ep.stimRad(i,:)];
              condStrVec = string(condStrVec);
          end % string vector for easier comparison
          %theConds = unique(condStrVec);
          theConds = ["39","79","157","236"];
          nConds = length(theConds);
          cond = ["2.5","5","10","15"]
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
          % stimulus type
          stimStrVec = string(obj.ep.stimType);
          theStims = unique(stimStrVec, 'rows');
          nStim = length(theStims);
          
          % vGab
          theGabs = string(obj.ep.vGabs); 
          
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
          br.vLags = zeros(nConds * nSubs, 1); % the v's are universal
          br.vPks = br.vLags;
          br.vLagErr = br.vLags;
          br.vPksErr = br.vLags;
          br.hLags = br.vLags;                  % the h's are h for zebra,and dPhi for the others
          br.hPks = br.vLags;
          br.hLagErr = br.vLags;
          br.hPksErr = br.vLags;
          
          loopCount = 0;
          figHans = zeros(nConds.*nSubs, 1);
          figure 
          nConds 
          for i=1:nConds
              for j = 1:nSubs
                  
                  loopCount = loopCount + 1;
                  
                  if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
                 %thisCombo = (strcmp(condStrVec, theConds(i))) &
                 %theseSubs; % zeb conditions etc 
                  thisCombo = (strcmp(condStrVec, theConds(i))') & theseSubs;
                  if any(thisCombo)  % ***** is there anything to bootstrap?
                      
                      deezRows = find(thisCombo); % if there is, den deez r da rows ta uze
                      
                      thisIter = (i-1).*nSubs + j;
                      %                       subplot(nConds,nSubs,thisIter)
                      
                      %figHans(loopCount) = figure;  % why is this here?  Possibly an ouput arg at some point
                      %figHans(loopCount) = nConds;
                      
                      if strcmp(theStims,"cfx") && strcmp(theGabs,"0")
                            leMean = nanmean(-obj.xc.xCorDP(thisCombo,:));
                      else
                         leMean = nanmean(obj.xc.xCorY(thisCombo,:)); 
                      end 
                      leNoise = rms(leMean(1:zeroInd));
                      leSNR = max(leMean)./leNoise;
                      
                      if leSNR >= 1.5
                          
                          %*****    dew da bewtstrappin  *********
                          %*****    first for the y - direction, which is univesal ********
                          if  (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79"))
                              cly =[0.28,0.08,0.02;  0.64,0.08,0.23;  0.8,0.38,0.03; 0.94,0.66,0.24];
                              
                              if strcmp(theStims,"cfx") && strcmp(theGabs,"0")
                                  for k = 1:nBtReps
                                      deezBootRows = randsample(deezRows, length(deezRows), true);
                                      disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                      [PkMaxs(k), li] = max(-disBootMean);
                                      PkLags(k) = pLags(li);
                                      %ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                                  end % end lil' boot loop
                                  hold on
                                  
                                  phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
                                  errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                      'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                      'LineWidth', 1.5, 'CapSize', 10,'Color',cly(i,:));
                                  legend(phdp,cond);
                                  hold off
                                  
                              else
                                  
                                  for k = 1:nBtReps
                                      deezBootRows = randsample(deezRows, length(deezRows), true);
                                      disBootMean = nanmean(obj.xc.xCorY(deezBootRows,zeroInd:end));
                                      [PkMaxs(k), li] = max(disBootMean);
                                      PkLags(k) = pLags(li);
                                      %ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.01, [1 0.5 1]);
                                  end % end lil' boot loop
                                  hold on
                                  
                                  phy(i) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', cly(i,:), 'LineWidth', 2);
                                  errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                      'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
                                      'LineWidth', 1.5, 'CapSize', 10,'Color',cly(i,:));
                                  legend(phy,cond);
                                  hold off
                              end
                          end
                          
                          % and store the lags, the peaks, and their error estimates
                          br.vLags(thisIter) = mean(PkLags);
                          br.vPks(thisIter) = mean(PkMaxs);
                          br.vLagErr(thisIter) = std(PkLags);
                          br.vPksErr(thisIter) = std(PkMaxs);
                          
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
                          ylim([-0.1,0.3]);
                          yticks(-0.1:0.1:0.3);
                          
                          hold on
                          
                          if (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79")) 
                              if  strcmp(theStims,"zeb") || ((strcmp(theStims, "cfx"))&& strcmp(theGabs,"0"))
                                  leMean = nanmean(obj.xc.xCorX(thisCombo,:));
                                  
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
                                  
                                  
                                  %subplot(1,2,1);
                                  phx(i) = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', clx(i,:), 'LineWidth', 2);
                                  errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                      'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                      'LineWidth', 1.5, 'CapSize', 10,'Color',clx(i,:));
                                  
                                  legend(phx,cond);
                                  hold off
                                  
                              elseif  (strcmp(theStims, "cfx"))&& strcmp(theGabs,"1")
                                  leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
                                  %*****    dew eevin moar bewtstrappin  *********
                                  %*****    or for the phase, if in a cuttlefish condition ********
                                  for k = 1:nBtReps
                                      deezBootRows = randsample(deezRows, length(deezRows), true);
                                      disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                      [PkMaxs(k), li] = max(-disBootMean);
                                      PkLags(k) = pLags(li);
                                      %ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                                  end % end lil' boot loop
                                  hold on
                                  
                                  clp =[0.5,0.05,0.93;  0.4,0.33,0.93;  0.3,0.52,0.93; 0.2,0.95,0.93];
                                  phdp(i) = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', clp(i,:), 'LineWidth', 2);
                                  errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                      'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                      'LineWidth', 1.5, 'CapSize', 10,'Color',clp(i,:));
                                  legend(phdp,cond);
                                  hold off
                              end
                              
                          end % if the condition of the eccentricity 
                          hold off
                          
                          xLabStr = 'Lag (seconds)';
                          %xLabStr = 'Lag'
                          yLabStr = 'CCG value';
                          
                          if type == 'c'  && (numel(theseSubs) > 1) % more that one subj in data file
                              substr = "all";
                          else
                              substr = theSubs(j);      % plotting by condition, or single S in file
                          end
                          %titstr = strcat("Condition = ", string(theConds(i)), " Subject = ", substr); 
                          % PSC
                          titstr = strcat("X-direction (hGab)");
             
                          
                          if (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79")) &&(strcmp(theStims, "zeb"))
                              %legend([phx, phy], {'H', 'V'});
                              %legend(phx(i),["157","236","39","79"])
                              hold on;
                          elseif (strcmp(theConds(i), "157") || strcmp(theConds(i), "236") || strcmp(theConds(i), "39")||...
                                  strcmp(theConds(i), "79")) &&(strcmp(theStims, "cfx"))
                              %legend([phdp, phy], {'Phase', 'V'});
                          end
                          ylim([-0.1,0.3]); %% YLIM control
                          %ylim([-0.1, 0.1]);
                          set(gcf, 'Name', strcat(string(theConds(i)), substr));
                          
                          % and store the horizontal or phase results
                          br.hLags(thisIter) = mean(PkLags);
                          br.hPks(thisIter) = mean(PkMaxs);
                          br.hLagErr(thisIter) = std(PkLags);
                          br.hPksErr(thisIter) = std(PkMaxs);
                          
                          br.cond(thisIter, :) = theConds(i);
                          br.subj(thisIter, :) = substr;
                          
                          br.figHans = figHans;
                          
                          %br.phx(thisIter, :) = phx(i);
                          %br.phy(thisIter, :) = phy(i);
                          
                      else % SNR critereon block
                          
                          br.hLags(thisIter) = NaN;
                          br.hPks(thisIter) = NaN;
                          br.hLagErr(thisIter) = NaN;
                          br.hPksErr(thisIter) = NaN;
                          
                          br.cond(thisIter, :) = theConds(i);
                          br.subj(thisIter, :) = substr;
                          
                          br.figHans = figHans;
                          
                          
                      end % end of SNR critereon block
                      
                  end % ***** end if(data in this condition) block
                  
              end % subjects
              
          end % conditions
          hold on;
%           cond = ["2.5","5","10","15","2.5","5","10","15"];
%           if (strcmp(theStims, "cfx"))&& strcmp(theGabs,"0")
%               legend([phdp,phx],cond, 'Location','northeast');
%           else
%           legend([phy,phdp],cond, 'Location','northeast');
%           end

            cond = ["2.5","5","10","15"];
            legend(phx,cond);
          formatFigure(xLabStr, yLabStr, titstr);
          hold on;
      end
      %***** end bootanal method *****************************************************************************
      
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