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

classdef ccg
   
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
      function obj = ccg(myfilename)
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
      function plotraw(obj, type, snum)
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
          condStrVec = string(obj.ep.stimType); % string vector for easier comparison
          theConds = unique(condStrVec, 'rows');
          nConds = length(theConds);
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
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
                  
                  thisCombo = (strcmp(condStrVec, theConds(i))) & theseSubs;
                  if any(thisCombo)  % ***** try to plot only if there are data to plot
                      subplot(nConds, nSubs, (i-1).*nSubs + count)
                      plot(theLags, obj.xc.xCorY(thisCombo,:)','Color', [1 0.8 1], 'LineWidth', 0.5);
                      hold on
                      if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || strcmp(theConds(i), "pld")
                          plot(theLags, obj.xc.xCorX(thisCombo,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
                          phx = plot(theLags, mean(obj.xc.xCorX(thisCombo,:)), 'Color', [.2 .2 1], 'LineWidth', 2);
                      else
                          plot(theLags, -obj.xc.xCorDP(thisCombo,:)','Color', [.8 .8 1], 'LineWidth', 0.5);
                          phdp = plot(theLags, -mean(obj.xc.xCorDP(thisCombo,:)), 'Color', [.2 .2 1], 'LineWidth', 2);
                      end
                      phy = plot(theLags, mean(obj.xc.xCorY(thisCombo,:)), 'Color', [1 0.2 1], 'LineWidth', 2);
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
                      
                      if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || strcmp(theConds(i), "pld")
                          legend([phx, phy], {'H', 'V'});
                      else
                          legend([phdp, phy], {'Phase', 'V'});
                      end
                  end % ***** end plotting if block
              
              end % end subj loop
          end % end condition loop
      end
      %***** end plotraw method *****************************************************************************
      
      %% ***** plot the mean ccgs per cond and (opt) subj *******************************
      function plotmeans(obj, type)
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
          condStrVec = string(obj.ep.stimType); % string vector for easier comparison
          theConds = unique(condStrVec, 'rows');
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
                  thisCombo = (strcmp(condStrVec, theConds(i))) & theseSubs;
                  if any(thisCombo)  % ***** try to plot only if there are data to plot
                      
                      subplot(nConds,nSubs,(i-1).*nSubs + j)
                      leMean = nanmean(obj.xc.xCorY(thisCombo,:));
                      leNoise = rms(leMean(1:zeroInd));

                      phy = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [1 0.2 1], 'LineWidth', 2);
                      
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
                      if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || strcmp(theConds(i), "pld")
                          leMean = nanmean(obj.xc.xCorX(thisCombo,:));
                          phx = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [.2 .2 1], 'LineWidth', 2);
                      else
                          leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
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
                      
                      if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || strcmp(theConds(i), "pld")
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
      function br = bootanal(obj, type)  %**** oh what oh what should yon function return?
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
          condStrVec = string(obj.ep.stimType); % string vector for easier comparison
          theConds = unique(condStrVec, 'rows');
          nConds = length(theConds);
          
          % set up subject sorting
          subStrVec = string(obj.ep.subj);
          theSubs = unique(subStrVec, 'rows');
          nSubs = length(theSubs);
          
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
          
          for i=1:nConds
              for j = 1:nSubs
                  
                  loopCount = loopCount + 1;
                  
                  if type == 'c'
                      theseSubs = ones(size(subStrVec));
                  else
                      theseSubs = strcmp(subStrVec, theSubs(j));
                  end
                  
                  thisCombo = (strcmp(condStrVec, theConds(i))) & theseSubs;
                  if any(thisCombo)  % ***** is there anything to bootstrap?
                      
                      deezRows = find(thisCombo); % if there is, den deez r da rows ta uze
                      
                      thisIter = (i-1).*nSubs + j;
                      %                       subplot(nConds,nSubs,thisIter)
                      
                      figHans(loopCount) = figure;  % why is this here?  Possibly an ouput arg at some point
                      
                      leMean = nanmean(obj.xc.xCorY(thisCombo,:));
                      leNoise = rms(leMean(1:zeroInd));
                      leSNR = max(leMean)./leNoise;
                      
                      if leSNR >= 1.5
                          
                          %*****    dew da bewtstrappin  *********
                          %*****    first for the y - direction, which is univesal ********
                          for k = 1:nBtReps
                              deezBootRows = randsample(deezRows, length(deezRows), true);
                              disBootMean = nanmean(obj.xc.xCorY(deezBootRows,zeroInd:end));
                              [PkMaxs(k), li] = max(disBootMean);
                              PkLags(k) = pLags(li);
                              ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.01, [1 0.5 1]);
                          end % end lil' boot loop
                          hold on
                          phy = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [1 0.5 1], 'LineWidth', 2);
                          errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                              'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [1, .2, 1], ...
                              'LineWidth', 1.5, 'CapSize', 10);%16% upper and lower quantile
                          hold off
                          
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
                          
                          hold on
                          if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || strcmp(theConds(i), "pld")
                              
                              leMean = nanmean(obj.xc.xCorX(thisCombo,:));
                              
                              %*****    dew moar bewtstrappin  *********
                              %*****    now for the x - direction, if in the zeb condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorX(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(disBootMean);
                                  PkLags(k) = pLags(li);
                                  ghostline(theLags(zeroInd:end), disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              phx = plot(theLags(zeroInd:end), leMean(zeroInd:end), 'Color', [.8 .8 1], 'LineWidth', 2);
                              errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                  'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                  'LineWidth', 1.5, 'CapSize', 10);
                              hold off
                              
                          else
                              leMean = nanmean(obj.xc.xCorDP(thisCombo,:));
                              %*****    dew eevin moar bewtstrappin  *********
                              %*****    or for the phase, if in a cuttlefish condition ********
                              for k = 1:nBtReps
                                  deezBootRows = randsample(deezRows, length(deezRows), true);
                                  disBootMean = nanmean(obj.xc.xCorDP(deezBootRows,zeroInd:end));
                                  [PkMaxs(k), li] = max(-disBootMean);
                                  PkLags(k) = pLags(li);
                                  ghostline(theLags(zeroInd:end), -disBootMean, 0.01, 0.02, [.8 .8 1]);
                              end % end lil' boot loop
                              hold on
                              phdp = plot(theLags(zeroInd:end), -leMean(zeroInd:end), 'Color', [.8 .8 1], 'LineWidth', 2);
                              errorbar(mean(PkLags), mean(PkMaxs), std(PkMaxs), std(PkMaxs), std(PkLags), std(PkLags), ...
                                  'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', [.2 .2 1], ...
                                  'LineWidth', 1.5, 'CapSize', 10);
                              hold off
                          end % if zeb or cuttlefish
                          hold off
                          
                          xLabStr = 'Lag (seconds)';
                          yLabStr = 'CCG value';
                          
                          if type == 'c'  && (numel(theseSubs) > 1) % more that one subj in data file
                              substr = "all";
                          else
                              substr = theSubs(j);      % plotting by condition, or single S in file
                          end
                          titstr = strcat("Condition = ", string(theConds(i)), " Subject = ", substr); 
                          
                          formatFigure(xLabStr, yLabStr, titstr, 0, 0, 32, 32);
                          
                          if strcmp(theConds(i), "zeb") || strcmp(theConds(i), "apt") || ...
                                  strcmp(theConds(i), "pld") || strcmp(theConds(i), "cfx")
                              legend([phx, phy], {'H', 'V'});
                          else
                              legend([phdp, phy], {'Phase', 'V'});
                          end
                          
                          set(gcf, 'Name', strcat(string(theConds(i)), substr));
                          
                          % and store the horizontal or phase results
                          br.hLags(thisIter) = mean(PkLags);
                          br.hPks(thisIter) = mean(PkMaxs);
                          br.hLagErr(thisIter) = std(PkLags);
                          br.hPksErr(thisIter) = std(PkMaxs);
                          
                          br.cond(thisIter, :) = theConds(i);
                          br.subj(thisIter, :) = substr;
                          
                          br.figHans = figHans;
                          
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