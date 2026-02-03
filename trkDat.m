%TRKDAT  Create a tracking data object
%
%
%   TD = TRKDAT(FILENAME) creates a tracking data object TD using
%   THE data in the file FILENAME.
%   The file FILENAME should contain:
%   a struct of experimental parameters, ep, and
%   a struct of tracking data, d
%
%
%   The primary thing that TRKDAT does is to compute the tracking cross-correlations
%   and create a CCG object, which can then be used to summarize and plot the ccgs.

classdef trkDat
    
    %% trkDat properties
    properties
        ep        % experimental parameters structure
        d         % data structure
        frameRate
        xc    % cross correlations
    end
    properties (Constant)
        maxLagSec =  1.2;    % specify max lag (secs) over which to compute cross-correlations
        frontPorchSec = 1;     % initial time (secs) to drop from the analysis
    end
    properties (Dependent)
        maxLag
        frontPorch = round(frontPorchSec*ep.frameRate);
    end
    
    %% trkDat methods 
    methods
        %%  der konstruktor - see also the "create" static method below
        function obj = trkDat(myfilename)
            if nargin > 0
                load(myfilename, 'ep', 'd');
                obj.ep = ep;
                obj.d = d;
                obj.frameRate = ep.frameRate;
            end
        end
        
        %% methods for dependent properties
        % get the maximum lag for the cross-correlation, and how much initial data in a trial to ignore
        function maxLag = get.maxLag(obj)
            maxLag = round(obj.maxLagSec*obj.frameRate);
        end
        function frontPorch = get.frontPorch(obj)
            frontPorch = round(obj.frontPorchSec*obj.frameRate);
        end
        
        %% methods for plot and display
        % plot selected stimulus and response traces
        function plotthis(obj)
            % dummey for now
            plot(obj.d.normTargCoords(1,:, 1));
            hold on
            plot(obj.d.normRespCoords(1,:,1));
            hold off
        end
        
        % display object
        function disp(obj)
            disp('experimental parameters (ep): ');
            disp(obj.ep);
            disp('data (d): ');
            disp(obj.d);
        end
        
        %% money methods
        %*****************************************************
        %***** Kompute the ccgs and save the ccg object. *****
        %***** This is the big-kahuna method for trkDat  *****
        function xc = makeccgs(obj)
            
            %***** compute velocities (from raw data) for the envelope motion
            respVel = -diff(obj.d.normRespCoords, 1, 2);  % first difference, taken along the rows
            windVel = diff(obj.d.normTargCoords, 1, 2);   % could compute from temp and ept.stepSize, but this is easier...
            respAcc = diff(respVel,1 ,2); 
            PhaseAcc = diff(obj.d.dPhases,1,2); % stimulus phase acceleration 
            dPhases = obj.d.dPhases;
%             targVel = diff(dt.targCoordsActual, 1, 2);   % the actual target position, i.e. the visual stimulus

            fp = obj.frontPorch;
            mxl = obj.maxLag;
            [nTrials, ~, ~] = size(obj.d.normTargCoords);
            
            % allocate the xc arrays
            xCorX = zeros(nTrials, 2.*mxl+1);
            xCorY = xCorX;
            xCorYD = xCorX;
            xCorDP = xCorX;
            xCorS = xCorX;
            xCorTVSA = xCorX;
            
            % make an index vector to get the right response orientations for correlating
            % with the carrier motion
            theDirec = ones(nTrials, 1);
            theDirec(obj.ep.stAngle == 90) = 2; % 90 is horizontal bar 
            
            % okay, here we go
            for j = 1:nTrials
                
                %***** analysis for x and y envelope motion and response *****
                % x (horizontal) responses are page 1 (:,:,1), y = 2 (:,:,2)
                [xCorX(j,:), lags] = xcorr(respVel(j,fp:end,1), windVel(j,fp:end,1), mxl, 'coeff');
                [xCorY(j,:), ~] = xcorr(respVel(j,fp:end,2), windVel(j,fp:end,2), mxl, 'coeff');
                
                [xCorS(j,:), ~] = xcorr(respVel(j,fp:718,1), respAcc(j,fp:718,1), mxl, 'coeff'); 
                
                %**** and the carrier speed  *******
                [xCorDP(j,:), ~] = xcorr(respVel(j,fp:end,theDirec(j)), dPhases(j, fp:end-1), mxl, 'coeff');
                [xCorTVSA(j,:), ~] = xcorr(respVel(j,fp:718,theDirec(j)), PhaseAcc(j,fp:718,1), mxl, 'coeff');            
            end
            
            xCorYD = diff(xCorY,1,2);
            
            xc.xCorX = xCorX;
            xc.xCorY = xCorY;
            xc.xCorYD = xCorYD;
            xc.xCorDP = xCorDP;
            xc.xCorS = xCorS;
            xc.xCorTVSA = xCorTVSA; 
            xc.lags = lags;
            
            ep = obj.ep;    % naming this local variable the same as the property, so it has that name in saved file
            
            save('xc.mat', 'xc', 'ep');
            
            % for debugging
            save('vels.mat', 'respVel', 'windVel', 'dPhases');
            
        end % end of makeccgs function
        %*****************************************************
        
        
    end % end regular methods
    
    %% static constructor method
    methods (Static)
        function obj = create % this is a function you can call from command line, as in myTD = TRKDAT.createObj
            [myfilename, pathname] = uigetfile('*.mat', 'Pick the tracking mat file');
            myfullfile = fullfile(pathname, myfilename);
            obj = trkDat(myfullfile);
        end
    end % end static methods
    
end % end classdef