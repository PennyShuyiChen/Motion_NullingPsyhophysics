% WindyNullingDemo.m 
% ___________________________________________________________________
% 
% Try to stabilize a target with the mouse as it's trying to move randomly.
% 
% After the trial, some plots of your performance will be shown, including
% cross-correlation functions (CCGs), which are quite useful.
% See:
% Mulligan, Stevenson, Cormack (2013) SPIE: Human Vision and Electronic Imaging
% Bonnen et al. (2015) Journal of Vision
% Bonnen, Huk, and Cormack (2017) Journal of Neurophysiology
% _______________________________________________________________________
%
%% simulate the target move le to right - the subject try to follow the motion using a mouse
% 1. the position of the target as a function generate white noise as the
% velocityie. Integrate the velocities, the vector of cumulative sum of the
% briownian noise 
% 2. perfect subjecy copy the subject response  velocities of the response
% vector - derivative if the 1st difference 
% stimulus response vector and subject response vectir
% 3. cross-correlation: for loop - 20s * 60Hz = 1200 data points; maxLag
% 2s, 2*60 - 120 samples; v
% 4. plot the vector for the perfect subject 

% HISTORY
% 18-Apr-2018   lkc Wrote it based on BasicTrackingDemo by the same author

% NOTES
% You must have PsychToolbox installed (psychtoolbox.org)
% WARNING: This is a script, and it's going to start by clearing your workspace!!!
% 18-Apr-2018    lkc cleaning it up for the Githubs 

clear

%% set up initial parameters

cont = .08; % .32, .16, .08, .04, .02, .01  see a pattern?
gabHoriz = 0;  % 1 for horizontal stripes, 0 for vertical
spFreq = 1;     % in cyc/deg
bw = 1;         % in octaves
gabsize = 100;
halfg = round(gabsize/2);
mygray = 120;

trialTime = 10;      % trial duration in s
stepSize = 4;      % ave stepsize in pixels - determines target speed

smooveSD = 1;   % set to zero for no smoothing

% some analysis params
maxLag =  1.2;    % specify max lag (secs) over which to compute cross-correlations
frontPorch = 1;     % initial time (secs) to drop from the analysis

if smooveSD
    kerSize = ceil(6.*smooveSD);
	x=-kerSize/2:kerSize/2-1;
    SmooveKer=(1/(smooveSD*sqrt(2*pi))).*exp(-x.^2./(2*smooveSD^2));
else
    SmooveKer = [];  % in case there's an old one hanging around
end

% define some colornames
% bkCol = 128;    % background color
bkCol = mygray;    % background color % PSC: 120 

% dot color & size
dotSizes = 10; % first one is target
dotColors = [20; 200; 20];

% make a Gabor matlab matrix
cinc = round(cont.*mygray);

targGab = ggab(gabsize, spFreq, bw);
targGab = round(normalize(targGab, mygray+cinc, mygray-cinc)); % change contrast here - ha!
if gabHoriz
    targGab = targGab';
end

%% stimulus presentation
try
    Screen('Preference', 'SkipSyncTests', 1);
    % Open up a window on the screen and clear it.
    whichScreen = max(Screen('Screens'));
    [theWindow,theRect] = Screen(whichScreen,'OpenWindow',[bkCol, bkCol, bkCol],[],[],2);

    mfi = Screen('GetFlipInterval', theWindow);  % should work always, unlike Screen('FrameRate')
    frameRate = round(1./mfi);     % get nominal frame rate
    totFrames = trialTime*frameRate;    % yup
    
    % now make the gabor texture and rects
    gabTex = Screen('MakeTexture', theWindow, targGab);


    %***** Generate the target sequence *****
    % generate random velocities, smooove them, and then integrate them to get positions.
    temp = randn(2,totFrames);
    if smooveSD
        temp(1,:) = conv(temp(1,:), SmooveKer, 'same');
        temp(2,:) = conv(temp(2,:), SmooveKer, 'same');
        temp(:,1) = [0; 0];
    end

    targCoordsNoResp = stepSize*cumsum(temp, 2);    % what target positions would be if wind only (no response from S)
    respCoords = zeros(size(targCoordsNoResp));     % response coordinate array
    targCoordsActual = respCoords;                  % actual target positions given wind and response

    % Move the cursor to the center of the screen
    centerX = theRect(RectRight)/2;
    centerY = theRect(RectBottom)/2;
    targCoordsNoResp(1,:) = targCoordsNoResp(1,:) + centerX;
    targCoordsNoResp(2,:) = targCoordsNoResp(2,:) + centerY;
    xTargOld = centerX;
    yTargOld = centerY;
    
%     refRectSz = round(0.8.*dotSizes.*stepSize);
    refRectSz = 200;
    refRect = [centerX-refRectSz, centerY-refRectSz, centerX+refRectSz, centerY+refRectSz];
    
    HideCursor();

    % Wait for a click
    Screen(theWindow,'FillRect',[bkCol, bkCol, bkCol]);
    Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
    Screen(theWindow,'DrawText','Klick to start trial',50,50,255);
    Screen('Flip', theWindow);
    while(1)
        [x,y,buttons] = GetMouse(whichScreen);
        if buttons(1)
            break
        end
    end
    
    Screen('DrawTexture', theWindow, gabTex, [], [centerX-halfg, centerY-halfg, centerX+halfg, centerY+halfg]);
%     Screen('DrawDots', theWindow, [centerX, centerX; centerY, centerY], dotSizes, dotColors);
    Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
    vbl = Screen('Flip', theWindow);
    
    pause(0.4); % ready... set...
    
    SetMouse(centerX,centerY);

    tic % PSC start the timer here 
    for i = 1:totFrames
        
        % save last mouse position and get new position, compute change
        xMouseOld = x;
        yMouseOld = y;
        [x,y,buttons] = GetMouse(whichScreen);
        
        deltaMouse = [x-xMouseOld; y-yMouseOld]; 
        if i ==1; deltaMouse = [0; 0]; end
        
        % new position is old position + random jump + mouse movement
        xTargNew = xTargOld + stepSize.*temp(1, i) + deltaMouse(1);
        yTargNew = yTargOld + stepSize.*temp(2, i) + deltaMouse(2);
                
        % change color of target when close to make more entertaining
        if abs(xTargNew-centerX)<refRectSz && abs(yTargNew-centerY)<refRectSz
            dotColors = [20; 200; 20];
        else
            dotColors = [200; 20; 20];
        end
        
        Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
%         Screen('DrawDots', theWindow, [xTargNew; yTargNew], dotSizes, dotColors);
        Screen('DrawTexture', theWindow, gabTex, [], [xTargNew-halfg, yTargNew-halfg, xTargNew+halfg, yTargNew+halfg]);

        vbl = Screen('Flip', theWindow, vbl + 0.5.*mfi);  
        
        respCoords(:,i) = [x; y];               % where the mouse or finger actually is
        targCoordsActual(:,i) = [xTargNew ; yTargNew];  % where the target ended up, hopefully near the center
        
        xTargOld = xTargNew;
        yTargOld = yTargNew; % PSC: x,y updated for every frame ?

        
    end % end trial for loop
    actualtime = toc;
    % Close up
    ShowCursor(0);
    Screen(theWindow,'Close');

catch
    Screen('CloseAll');
    ShowCursor;
    psychrethrow(psychlasterror);
end %try..catch..

%% data analysis!
% let's look at x and y correlation peaks and latencies as our index of
% adaptation

% convert to samples for lag and front bumper
maxLag = round(maxLag*frameRate); % PSC: 72 ?
frontPorch = round(frontPorch*frameRate); % PSC: 60 ?

% center coords - since this is a demo deployed in the wild, leave everything in pixels
% respCoords = where the mouse is on the desk (or the finger on the pad)
normRespCoords(1,:) = (respCoords(1,:) - centerX);
normRespCoords(2,:) = (respCoords(2,:) - centerY);

normTargCoords(1,:) = (targCoordsNoResp(1,:) - centerX);
normTargCoords(2,:) = (targCoordsNoResp(2,:) - centerY);

% compute velocities (from raw data) % PSC: velocities or position
% (coordinates) ?
respVel = -diff(respCoords, 1, 2); % PSC: dimension 2 - columns 
windVel = diff(targCoordsNoResp, 1, 2);   % could compute from temp and stepSize, but this is easier...
targVel = diff(targCoordsActual, 1, 2);   % the actual target position, i.e. the visual stimulus

%% do the cross-correlation analysis
% compute velocity cross correlations
[xCorX, xLagsVel] = xcorr(respVel(1,frontPorch:end), windVel(1,frontPorch:end), maxLag, 'coeff');
[xCorY, yLagsVel] = xcorr(respVel(2,frontPorch:end), windVel(2,frontPorch:end), maxLag, 'coeff');

timeZeroIndex = ceil(length(xCorX)/2); % PSC ?

% find peak lags and heights for velocity
[xMax, xLagAtMax] = max(xCorX(timeZeroIndex:end));
[yMax, yLagAtMax] = max(xCorY(timeZeroIndex:end));

xRelLag = xLagAtMax-1; %-maxLag-1; % lag relative to zero
yRelLag = yLagAtMax-1; %-maxLag-1; % lag relative to zero

xMaxLagSecs = xRelLag./frameRate;
yMaxLagSecs = yRelLag./frameRate;

if xRelLag < 0 
    xRelLag = 0; 
    fprintf('Warning: x lag less than zero; S has ESP.');
end
if yRelLag < 0 
    yRelLag = 0; 
    fprintf('Warning: y lag less than zero; S has ESP.');
end

corrNoise = rms([xCorX(1:timeZeroIndex), xCorY(1:timeZeroIndex)]);

%% running spatial errors 
% compute absolute (pythagorian) error vs. time
cursErr = targCoordsActual - [centerX; centerY];
absCursErr = sqrt(cursErr(1,:).^2 +cursErr(2,:).^2);
aveAbsErr = mean(absCursErr);

% using hist instead of histcounts for all the pre-2014a people out there
[aveCnts, aveBs] = hist(absCursErr, 15);


%% save coords adn  for offline playing
save demodat.mat respCoords targCoordsNoResp xCorX xCorY xLagsVel
    

%%
%****************************************************************************
%***** That is the end of the actual data collection and analysis       *****
%***** Everything below is plotting plotting plotting                   *****
%****************************************************************************

%% plotting!  Yay!
    
% raw data traces, i.e. the walks in x, y for both stimulus and response
rawFig = figure;
plot(normTargCoords(1,:),-normTargCoords(2,:), 'go:');
hold on;%    
plot(-normRespCoords(1,:),normRespCoords(2,:), 'ko-');
% mark the beginning and end points
plot(normTargCoords(1,1),-normTargCoords(2,1), 'bd', 'MarkerSize', 20, ...
    'MarkerFaceColor', 'b');
plot(normTargCoords(1,end),-normTargCoords(2,end), 'rh', 'MarkerSize', 20, ...
    'MarkerFaceColor', 'r');
hold off;
legend('the wind', 'your response', 'start', 'finish');
xlabel('Horizontal', 'FontSize', 18);
ylabel('Vertical', 'FontSize', 18);
set(gca,'FontSize',14);
axis equal;

% ***** plots of the positions and velocities as a function of time *****
% The instances of 1:end below can be replace by frontPorch:end to omit
% plotting the first n=frontPorch samples

% make a time axis for our plots
realTime = 0:1/frameRate:trialTime;
realTime(end) = [];
if length(realTime) ~= totFrames    % length reality check
    fprintf('warning: number of frames and time axis not lining up');
end

tSeriesFig = figure;

subplot(2,2,1); % H pos
plot(realTime, -normRespCoords(1,1:end), 'b')
hold on
plot(realTime, normTargCoords(1,1:end), 'g')
hold off;
legend('you', 'the wind');
title('Horizontal');
xlabel('Time', 'FontSize', 14);
ylabel('H Mag.', 'FontSize', 14);
set(gca,'FontSize',12);

subplot(2,2,2); % V pos
plot(realTime, -normRespCoords(2,1:end), 'k')
hold on
plot(realTime, normTargCoords(2,1:end), 'g')
hold off;
legend('you', 'the wind');
title('Vertical');
xlabel('Time', 'FontSize', 14);
ylabel('V Mag.', 'FontSize', 14);
set(gca,'FontSize',12);

subplot(2,2,3); % H vel
plot(realTime(2:end), windVel(1,1:end), 'g')
hold on
plot(realTime(2:end), respVel(1,1:end), 'b')
hold off;
legend('you', 'the wind');
title('Horizontal Velocity');
xlabel('Time', 'FontSize', 14);
ylabel('H Velocity', 'FontSize', 14);
set(gca,'FontSize',12);

subplot(2,2,4); % V vel
plot(realTime(2:end), windVel(2,1:end), 'g')
hold on
plot(realTime(2:end), respVel(2,1:end), 'k')
hold off;
legend('you', 'the wind');
title('Vertical Velocity');
xlabel('Time', 'FontSize', 14);
ylabel('V Velocity', 'FontSize', 14);
set(gca,'FontSize',12);

%% error vs. time figure
errVsTimeFig = figure;

shiftedTime = realTime(xRelLag+1:end); % ooohh weeeee ooohh

posvec = [.1 .1 .6 .8];
subplot('Position', posvec);
plot(realTime, absCursErr, 'r')
myxlim = xlim;
myylim = ylim;
text(mean(myxlim), 1.2*aveAbsErr, ['ave. error: ', num2str(aveAbsErr, '%1.3f')]);
line([myxlim(1), myxlim(2)], [aveAbsErr, aveAbsErr], ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
line([myxlim(1), myxlim(2)], [refRectSz, refRectSz], ...
    'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);
legend('absolute cursor error', 'ave. error', 'target zone boundary');
title('Absolute (Pythagorean) Error');
xlabel('Time', 'FontSize', 14);
ylabel('Absolute Position Error', 'FontSize', 14);
set(gca,'FontSize',12);

% Add marginal distribution of errors, which should be approximately Rayleigh.  It is.
posvec = [.75 .1 .2 .8];
subplot('Position', posvec);
barh(aveBs, aveCnts, 'r');
ylim(myylim);
myxlim = xlim;
line([myxlim(1), myxlim(2)], [aveAbsErr, aveAbsErr], ...
    'Color', 'k', 'LineStyle', ':', 'LineWidth', 2);
line([myxlim(1), myxlim(2)], [refRectSz, refRectSz], ...
    'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);
title('Error Distribution');
xlabel('Count', 'FontSize', 14);


%% ******  plots of the cross-correlations ******
xcorFig = figure;
% horiz velocity
plot(xLagsVel./frameRate, xCorX);
hold on;   % vert velocity
plot(yLagsVel./frameRate, xCorY, 'k');
hold off;
myylim = ylim;
line([xMaxLagSecs, xMaxLagSecs], [myylim(1), xMax], ...
    'Color', 'b', ...
    'LineStyle', ':');
text(xMaxLagSecs, xMax, [num2str(xMaxLagSecs, '%1.3f'), '  ', num2str(xMax, '%1.2f')]);
line([yMaxLagSecs, yMaxLagSecs], [myylim(1), yMax], ...
    'Color', 'k', ...
    'LineStyle', ':');
text(yMaxLagSecs, yMax, [num2str(yMaxLagSecs, '%1.3f'), '  ', num2str(yMax, '%1.2f')]);
xlabel('Lag (s)', 'FontSize', 18);
ylabel('Correlation', 'FontSize', 18);
set(gca,'FontSize',14);

myxlim = xlim;
xe = [myxlim, fliplr(myxlim)];
ye = [corrNoise, corrNoise, -corrNoise, -corrNoise];
hold on;
sHandle = fill(xe,1.5.*ye,'r');
hold off;
set(sHandle, 'EdgeColor', 'none');
set(sHandle, 'FaceAlpha', 0.1);

legend('Horizontal', 'Vertical', ...
    'H Peak', 'V Peak', '1.5xRMS error', ...
    'Location', 'NorthWest');


%***** now a cross-correlation plot showing only postive lags *****
xLagsVel = xLagsVel(timeZeroIndex:end);
xCorX = xCorX(timeZeroIndex:end);

yLagsVel = yLagsVel(timeZeroIndex:end);
xCorY = xCorY(timeZeroIndex:end);

anotherxcorFig = figure;
% horiz velocity
plot(xLagsVel./frameRate, xCorX);
hold on;   % vert velocity
plot(yLagsVel./frameRate, xCorY, 'k');
hold off;
myylim = ylim;
line([xMaxLagSecs, xMaxLagSecs], [myylim(1), xMax], ...
    'Color', 'b', ...
    'LineStyle', ':');
line([yMaxLagSecs, yMaxLagSecs], [myylim(1), yMax], ...
    'Color', 'k', ...
    'LineStyle', ':');

text(xMaxLagSecs, xMax, [num2str(xMaxLagSecs, '%1.3f'), '  ', num2str(xMax, '%1.2f')]);
text(yMaxLagSecs, yMax, [num2str(yMaxLagSecs, '%1.3f'), '  ', num2str(yMax, '%1.2f')]);
xlabel('Lag (s)', 'FontSize', 18);
ylabel('Correlation', 'FontSize', 18);
set(gca,'FontSize',14);

myxlim = xlim;
xe = [myxlim, fliplr(myxlim)];
hold on;
sHandle = fill(xe,1.5.*ye,'r');
hold off;
set(sHandle, 'EdgeColor', 'none');
set(sHandle, 'FaceAlpha', 0.1);

legend('Horizontal', 'Vertical', ...
    'H Peak', 'V Peak', '1.5xRMS error', ...
    'Location', 'NorthEast');


%% put the figs in a sensible order for user to click through
figure(anotherxcorFig);
figure(xcorFig)
figure(errVsTimeFig);
figure(tSeriesFig)
figure(rawFig)

