function myout = WindyMoIndPosShiftNulling_v1_psc(subj, stimType, vGabs, fname)
% WindyMoIndPosShiftNulling.m 
% PSC - 10.27.2020
% Not in a proper loop form, just for testing gabors 
% TO RUN TESTS: WindyMoIndPosShiftNulling_v1_psc('psc', 'zeb', 1, 'ptest10date.mat')


% Windy Motion Induced Position Shift Nulling
% valid stimulus types:  'zeb', 'cfo', 'cfx', 'cfh', 'itd', 'iod', 'apt', 'pld'
% 'zebra', 'cuttlefishOrtho', 'cuttlefish Ortho with x-mouse on', 'cuttlefish', 'indieTwoD', 'indieOneD', 'aperature', 'plaid'

% $$$ conditions are zeb, cfx, and iod.  pld is the plaid control condition with both H and V gratings superimposed

% vGabs = 1 for vertical Gabors, 0 for horizontal, 2 for plaids
% ___________________________________________________________________
% 
% Try to stabilize a global target with the mouse as it's trying to move randomly.
% 
% This is the experimental version - saves data, plotting optional
% 
% See:
% Mulligan, Stevenson, Cormack (2013) SPIE: Human Vision and Electronic Imaging
% Bonnen et al. (2015) Journal of Vision
% Bonnen, Huk, and Cormack (2017) Journal of Neurophysiology
% _______________________________________________________________________
%
% Stimulus type descriptions:
% zebra (zeb) - A herd of zebras. The apertures (zebras) jitter in both dimensions, and 
%               the Gabor carrier ("stripes") are fixed with respect to the aperature,
%               like the stripes on a zebra.
% cuttlefishOrtho (cfo) - The aperature jitters only in 1D parallel to the stripes, and the 
%               carrier does a  1D jitter re. the body (ortho to the stripes/envelope).
%               Like a cuttlefish moving vertically while it makes stripes move
%               along its body horizontally.  The horizontal motion of the mouse is disabled, so the 
%               responses to the MIPS are "open loop".
% cuttlefish (cfx) - A hypnosis of cuttlefish (school). The aperatures (cuttlefish bodies) jitter 
%               in both dimensions, while the stripes jitter re. the body in 1D.  The
%               mouse drives the envelop in both directions, but does not effect the carrier motion
%               re. the aperature.
% cuttlefish (cfh) - A hypnosis of cuttlefish (school). The aperatures (cuttlefish bodies) jitter 
%               in both dimensions, while the stripes jitter re. the body in 1D.
% indie2D and indie1D (itd, iod) - Like the cuttlefish, except the carrier brownian walk is
%               independent of the envelop motion instead of added to it.  Like aperatures
%               jittering over an independetly jittering carrier.
% aperature (apt) - Aperature motion only.  Apertures jittering over a carrier that is stationary
%               in world coordinates (whereas for zebras, the stripes are stationary in zebra
%               coordinates.
%
% Zebra and CuttlefishOrtho are kinda the money conditions.

% HISTORY
% 18-Apr-2018   lkc Wrote it based on BasicTrackingDemo by lkc

% NOTES
% You must have PsychToolbox installed (psychtoolbox.org)
% 
% 18-June-2018    Windy MIPS Nulling!
% Changing data saving strategy somewhat - I'm going to save out just the data, 
% and then have an object handle the analysis via methods.  

%% Trials & duration 
nTrials = 1; %5
ept.trialTime = 12 %12;      % trial duration in s 5x12 gives one minute of exp. time
                             % keep it 12 
makeMovie = 0;
showCursorAndCentroid = 0;  % only for movies and demos


%% argument handling
ni = nargin;
myout = 0;
% need to add checking for valid input arguments
ept.noHorzMouse = 0;

switch ni
    case 0
        error('must enter subj initials');
    case 1
        ept.subj = subj;
        ept.stimType = 'cfo';
        ept.vGabs = 1;
        fname = 'windy301.mat';
    case 2
        ept.subj = subj;
        ept.stimType = stimType;
        ept.vGabs = 1;
        fname = 'windy301.mat';
    case 3
        ept.subj = subj;
        ept.stimType = stimType;
        ept.vGabs = vGabs;
        fname = 'windy301.mat';
    case 4
        ept.subj = subj;
        ept.stimType = stimType;
        ept.vGabs = vGabs;
    otherwise    % nargin > 4
        error('args are subj, stimType, vertical flag, filename');
end


%% set up initial parameters
% The key flags that get set depending upon stimulus type are
        % ***           ept.TwoDimJitter [0, 1] does the carrier move in both dirs [1], or only parallel to stripes [0]
        % ***           ept.carrierMotion [0, 1] does the carrier phase change [1] or remain fixed relative to the envelope [0]
        % ***           ept.correctCarrierMotionForEnvlope [0, 1] is the phase walk independent of the evelope walk [1]
        % ***           ept.subjectChangesPhase [0, 1] does the subject's mouse control carrier phase [1]
        % ***           ept.supjectChangesEnvelop [0, 1] does the subject's mouse move the envelop also too?

switch ept.stimType
    case 'zeb'
        ept.TwoDimJitter = 1;
        ept.carrierMotion = 0;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 0;
        ept.subjectChangesPhase = 1;
        ept.supjectChangesEnvelop = 1;
    case 'cfo'
        ept.TwoDimJitter = 0;
        ept.carrierMotion = 1;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 0;
        ept.subjectChangesPhase = 0;
        ept.supjectChangesEnvelop = 0;
    case 'cfx'
        ept.TwoDimJitter = 0;
        ept.carrierMotion = 1;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 0;
        ept.subjectChangesPhase = 0;
        ept.supjectChangesEnvelop = 1;
    case 'cfh'
        ept.TwoDimJitter = 1;
        ept.carrierMotion = 1;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 0;
        ept.subjectChangesPhase = 0;
        ept.supjectChangesEnvelop = 1;
    case 'itd'
        ept.TwoDimJitter = 1;
        ept.carrierMotion = 1;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 1;
        ept.subjectChangesPhase = 1;
        ept.supjectChangesEnvelop = 0;
    case 'iod'
        ept.TwoDimJitter = 0;
        ept.carrierMotion = 1;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 1;
        ept.subjectChangesPhase = 1;
        ept.supjectChangesEnvelop = 0;  % need to let the envelope be adjusted in ortho direction!
    case 'apt'
        ept.TwoDimJitter = 2;
        ept.carrierMotion = 0;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 1;
        ept.subjectChangesPhase = 0;
        ept.supjectChangesEnvelop = 1;
    case 'pld'  % same as zeb, just force vGabs to be 2 for plaid control stimulus
        ept.TwoDimJitter = 1;
        ept.carrierMotion = 0;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.correctCarrierMotionForEnvlope = 0;
        ept.vGabs = 2;
        ept.subjectChangesPhase = 1;
        ept.supjectChangesEnvelop = 1;
    otherwise
        error('please specify a valid stimulus type');
end



%% Rig params - PSC 
vdist = 50.0; % in cm (closer vdist help to present the stimuli (?)) % 57.3/56.7
H_res = 1920; V_res = 1080; % monitor resolution #info from Rect 
H_im = 54.34; V_im = 30.56; % effective image size (?)# info online

VA = rad2deg(atan((H_im/2) / vdist)); % Visual angle for 1/2 screen (H)
VA_pix = VA/(H_res/2); % how much visual degree a pixel extends 

%% V1 & MT RFs - PSC
% change the parameter to ept.para for proper data storing af

% Eccentricities 
E_min = 2.5;
E_max = 15;
E_ncond = 4; % # eccentricity conditions

%E = linspace(E_min,E_max,E_ncond); % Eccentricity used
E = [2.5, 5, 10, 15];
% ring radius in pix(change to ept.stimRadius af testing) - PSC
r_R = round( E/VA_pix); % ring radius in pix

%% set Gabor stimulus parameters - PSC 
% 0 deg is vertical bars and increasing angle is CCW.  Postive phase change = rightward mo at zero deg
% fun Params

% Gabor orientation
if ept.vGabs
    ept.stAngle = 0;
else
    ept.stAngle = 90;
end

% Gabor carrier (MIPS) motions
ept.carrierSpeed = 30;     
drftSpeedDegPerFrame = ept.carrierMotion .* ept.carrierSpeed;    

% global form parameters
test = 4;
ept.ngabors = 16;    % 16 for experiment
ept.stimRad = r_R(test); %300;  % radius of the circle of Gabors in pixels, 300 for experiment
gabOffsets = ept.stimRad.*unitcirc(ept.ngabors);

refRectSz = 10;     % reference rectangle inside circle of Gabors (fixation center)


% Gabor initial parameters:
% Phase of underlying sine grating in degrees:
phase = 0;
% Frequency of sine grating:
pSF = exp(-0.49 * log(E(test)) + log(1.99*3.1)); % cpd
ept.freq = (pSF) * VA_pix; % dpp for the rig -> cpp (?)
% Spatial constant of the exponential "hull"
bw = (0.6808 + 0.0930 * log(E(test))); % bw in octave 
sc = bw2sig(ept.freq, bw);
% Size of support in pixels
tw = round(6 * sc);
th = round(6 * sc);

% Contrast of grating:
ept.contrast = 10;  %10
% Aspect ratio width vs. height:
aspectratio = 1.0;


% Initialize matrix with spec for all 'ept.ngabors' patches to start off
% identically:
mypars = repmat([phase+180, ept.freq, sc, ept.contrast, aspectratio, 0, 0, 0]', 1, ept.ngabors);

%%
% set up mouse and circle center dots if necessary
if showCursorAndCentroid
    % make two dots to draw, one corresponding to the walk, and the other to the response
    % and perhaps a dot for unwrapped phase (i.e. carrier walk)
    % dot colorz & sizes
    switch ept.stimType
        case 'zeb'
            % draw 2 dots, one for cursor, one for centroid
            dotSizes = [10, 1]; % first two are targets
            walkCol = [.1; 1]; % green(ish) is the new black
            phasCol = [1; .1];  % red is the new black
            cursCol = [0.1; 0.1];   % blue is the new black
            dotColors = [walkCol, cursCol];
        case 'cfx'
            % draw 3 dots; cursor, centoid, and centroid + unwrapped phase
            dotSizes = [10, 10, 10]; % first two are targets
            walkCol = [.1; 1; .1]; % green(ish) is the new black
            phasCol = [1; .1; .1];  % red is the new black
            cursCol = [0.1; 0.1; 1];   % blue is the new black
            dotColors = [walkCol, phasCol, cursCol];
        otherwise
    end
end


%%

ept.stepSize = 2;      % ave ept.stepSize in pixels - determines target speed - set to 0 for no envelope motion

ept.smooveSD = 1;   % set to zero for no smoothing

if ept.smooveSD
    kerSize = ceil(6.*ept.smooveSD);
	x=-kerSize/2:kerSize/2-1;
    SmooveKer=(1/(ept.smooveSD*sqrt(2*pi))).*exp(-x.^2./(2*ept.smooveSD^2));
else
    SmooveKer = [];  % in case there's an old one hanging around
end

% define some colornames

bkCol = 0.5;    % background color


%% set up PTB stuff
HideCursor();

Screen('Preference', 'SkipSyncTests', 1);

% Setup defaults and unit color range:
PsychDefaultSetup(2);

% Open up a window on the screen and clear it.
whichScreen = max(Screen('Screens'));
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange');
[theWindow,theRect] = PsychImaging('OpenWindow', whichScreen, 0.5);
gamma = 1/2.2;
PsychColorCorrection('SetEncodingGamma', theWindow, gamma); % gamma correction


centerX = theRect(RectRight)/2;
centerY = theRect(RectBottom)/2;

refRect = [centerX-refRectSz, centerY-refRectSz, centerX+refRectSz, centerY+refRectSz];
    
mfi = Screen('GetFlipInterval', theWindow);  % should work always, unlike Screen('FrameRate')
ept.frameRate = round(1./mfi);     % get nominal frame rate
totFrames = ept.trialTime*60;    %ept.frameRate;    

% Enable alpha-blending, set it to a blend equation useable for linear
% superposition with alpha-weighted source. This allows to linearly
% superimpose gabor patches in the mathematically correct manner, should
% they overlap. Alpha-weighted source means: The 'globalAlpha' parameter in
% the 'DrawTextures' can be used to modulate the intensity of each pixel of
% the drawn patch before it is superimposed to the framebuffer image, ie.,
% it allows to specify a global per-patch ept.contrast value:
Screen('BlendFunction', theWindow, GL_ONE, GL_ONE);

% Build a procedural gabor texture for a gabor with a support of tw x th
% pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
% change of aspect-ratio:
gabortex = CreateProceduralGabor(theWindow, tw, th, 1); % theWindow, tw, th, 1

% will prob need to add additional arg vals to get ept.contrast to play nice as per the help
     
texrect = Screen('Rect', gabortex);
inrect = repmat(texrect', 1, ept.ngabors);

dstRects = zeros(4, ept.ngabors);
for i=1:ept.ngabors
%     scale(i) = 1;   % *(0.1 + 0.9 * randn);  % use to mult. by texrect below to scale sizes
    dstRects(:, i) = CenterRectOnPoint(texrect, centerX+gabOffsets(i,1), centerY+gabOffsets(i,2))';
end
startRects = dstRects;

%% ***** Generate the target sequences *****

% Preallocate array with rotation angles:
% zero is rightward, increasing angles go anticlockwise
% tempRotAng = randn(totFrames, nTrials);
%     for j = 1:nTrials
%         tempRotAng(:,j) = conv(tempRotAng(:,j), SmooveKer, 'same');
%         tempRotAng(1,j) = 0;
%     end
% rotAngles = angStepSize*cumsum(tempRotAng, 1);

% %***** for debugging angular calculations   *****
% rotAngles = linspace(0, 720, totFrames);    %*****     
% rotAngles = repmat(rotAngles', 1, nTrials); %***** 
% %***** for debugging angular calculations   *****

rotAngles = ept.stAngle.*ones(totFrames, nTrials);    %*****     
dPhases = drftSpeedDegPerFrame .* randn(totFrames, nTrials); % 30 * normal, so most jumps < 30 deg, very rarely >= 90
carDrift = cumsum(dPhases.*(1./360)./(ept.freq));   % deg * (1/(deg/cyc)) / (cyc/pix) = pix 
% carDrift is the unwrapped carrier drift over the whole trial, expressed in pixels.  If you imagine the carrier as
% stripes painted on a tablecloth, it is the drift of the tablecloth itself.

% generate random velocities, smooove them, and then integrate them to get positions.
% doing this after PTB setup because we need the framerate
temp = randn(2,totFrames, nTrials);
if ept.smooveSD
    for j = 1:nTrials
        temp(1,:,j) = conv(temp(1,:,j), SmooveKer, 'same');
        temp(2,:,j) = conv(temp(2,:,j), SmooveKer, 'same');
        temp(:,1,j) = [0; 0];
    end
end

targCoordsNoResp = zeros(2,totFrames,nTrials);
targCoordsNoResp(:,:,:) = ept.stepSize*cumsum(temp, 2);    % what target positions would be if wind only (no response from S)

if ~ept.TwoDimJitter
    if ept.vGabs
        temp(1,:,:) = 0;
    else
        temp(2,:,:) = 0;
    end
end


respCoords = zeros(size(targCoordsNoResp));     % response coordinate array
dt.targCoordsActual = respCoords;                  % actual target positions given wind and response

%% start trial outer loop  
for j = 1:nTrials
    
    % make an nTrials x ept.ngabors array of rotation angles for this trial

    theseRotAngs = repmat(rotAngles(:,j), 1, ept.ngabors);
    
    % Move the cursor to the center of the screen
    targCoordsNoResp(1,:,j) = targCoordsNoResp(1,:,j) + centerX;
    targCoordsNoResp(2,:,j) = targCoordsNoResp(2,:,j) + centerY;
    xTargOld = centerX;
    yTargOld = centerY;
    
    dstRects = startRects;
    
    % Wait for a click
    Screen(theWindow,'FillRect',[bkCol, bkCol, bkCol]);
    Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
    tstrg = char(cellstr(strcat('Klick to start trial ', {' '}, num2str(j), {' of '}, num2str(nTrials))));
    Screen(theWindow,'DrawText',tstrg,50,50,255);
    Screen('Flip', theWindow);
    
    while(1)
        [x,y,buttons] = GetMouse(whichScreen);
        if buttons(1)
            break
        end
    end
    
    Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:), [], [], [], [], kPsychDontDoRotation, mypars);
    if ept.vGabs ==2
       Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:)+90, [], [], [], [], kPsychDontDoRotation, mypars);
    end
%     Screen('DrawDots', theWindow, [centerX, centerX; centerY, centerY], dotSizes, dotColors);
    Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
    vbl = Screen('Flip', theWindow);
    
    pause(0.4); % ready... set...
    
    SetMouse(centerX,centerY);

    tic
    for i = 1:totFrames
        
        % save last mouse position and get new position, compute change
        xMouseOld = x;
        yMouseOld = y;
        [x,y,~] = GetMouse(whichScreen);
        
        deltaMouse = [x-xMouseOld; y-yMouseOld]; 
        if i ==1; deltaMouse = [0; 0]; end
        
        % ********************************************************************************************************************
        % ********* okay, here is where we do all the adjustments of aperature position and carrier phase based upon *********
        % *********  we have/need flags for 1) what moves [Gabor, envelope, carrier] and what the subject controls   *********
        % ********* with the mouse [Gabor, envelope, carrier, nothing].
        % ********* The current relevant flags are:
        % ***           ept.TwoDimJitter [0, 1] does the carrier move in both dirs [1], or only parallel to stripes [0]
        % ***           ept.carrierMotion [0, 1] does the carrier phase change [1] or remain fixed relative to the envelope [0]
        % ***           ept.correctCarrierMotionForEnvlope [0, 1] is the phase walk independent of the evelope walk [1]
        % ***           ept.subjectChangesPhase [0, 1] does the subjects mouse control carrier phase [1]
        % ***           ept.supjectChangesEnvelop [0, 1] does mouse control envelope position
        % ***           note: when both previous flags are 1, the mouse just moves the Gabor
        % ***           
        % *** Steps:
        % ***    Do the wind motion
        % ***       if TwoDimJitter, then move the Gabor perpendicular to the stripes (in addition to parallel)
        % ***       if correctCarrierMotionForEvelope, then "undo" the (absolute) phase change caused by the aperature
        % ***               this results in "aperture only" motion
        % ***       if "open loop" then we're done (but this open loop idea isn't really working out anyway)
        % ***    Now do the subject mouse motion
        % ***       if subjectChangesPhase AND supjectChangesEnvelop, the add deltaMouse to dstRects and we're done.
        % ***       else if subjectChangesPhase (only) then recompute phase, and only add deltaMouse to dstRects in one D,
        % ***               depeding on the orientation of the grating
        % ***       else if subjectChangesEvelope (only) then change dstRects AND undo the absolute phase change caused by the
        % ***       aperature motion

        
        % new position is old position + random jump due do wind 
        xTargNew = xTargOld + ept.stepSize.*temp(1, i, j); % + deltaMouse(1);
        yTargNew = yTargOld + ept.stepSize.*temp(2, i, j); % + deltaMouse(2);
                
        % Increment phase-shift of each gabor by dPhases optionally correcting for envelope motion:
        % mypars is 8 params (last 3 dummy) x ept.ngabors
        % non-dummys are [phase+180, ept.freq, sc, ept.contrast, aspectratio]
        phaseCorrection = 0;
        if ept.correctCarrierMotionForEnvlope
         % phase correction is jump (pix) x ept.freq (cyc/pix) x 360 (deg/cyc)
         % this is correcting the phase for any envelop motion, be it from the trajectory or the subject's response.
            if ept.vGabs
                phaseCorrection = -(ept.stepSize.*temp(1, i, j)).* ept.freq .* 360; %  + deltaMouse(1)) 
            else
                phaseCorrection = (ept.stepSize.*temp(2, i, j)) .* ept.freq .* 360; %  + deltaMouse(2)
            end
        end

        % *** done with the wind, now incoporate the subject's mouse movements as necessary ***
        % first, handle case in which S moves the whole Gabor (as in Zebra condition)
        sPhaseChange = 0;
        sEnvPhaseCorrection = 0;
        if ept.supjectChangesEnvelop && ept.subjectChangesPhase
            xTargNew = xTargNew + deltaMouse(1);
            yTargNew = yTargNew + deltaMouse(2);
        % or, now, if S only changes phase
        elseif ept.subjectChangesPhase
            if ept.vGabs
                sPhaseChange = -deltaMouse(1) .* ept.freq .* 360;
                yTargNew = yTargNew + deltaMouse(2);
            else
                sPhaseChange = deltaMouse(2) .* ept.freq .* 360;
                xTargNew = xTargNew + deltaMouse(1);
            end
        % or if they are changing the envelope only
        elseif ept.supjectChangesEnvelop % need to move Gabors and then undo absolute phase change due to motion
            % so we move the envelope...
            xTargNew = xTargNew + deltaMouse(1);
            yTargNew = yTargNew + deltaMouse(2);
            % and then undo the resulting absolute phase change...
            if ept.vGabs
                sEnvPhaseCorrection = deltaMouse(1).* ept.freq .* 360; %  + deltaMouse(1)) 
            else
                sEnvPhaseCorrection = deltaMouse(2) .* ept.freq .* 360; %  + deltaMouse(2)
            end
        else 
            % mouse has no effect, so do nothing
        end
        
        % and, finally, update dstRects and phases
        dstRects = CenterRectOnPointd(inrect, xTargNew+gabOffsets(:,1)', yTargNew+gabOffsets(:,2)');
        mypars(1,:) = mypars(1,:) + dPhases(i, j) - phaseCorrection + sPhaseChange + sEnvPhaseCorrection;
        
        % ********* Okay, done with all the postion and phase adjusting, and ready get our drawing on!

        Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
        Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(i,:), [], [], [], [], kPsychDontDoRotation, mypars);
        if ept.vGabs ==2
            Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:)+90, [], [], [], [], kPsychDontDoRotation, mypars);
        end
        
        
        % For demos or really for making figures, draw the would-be walk and the subjects mouse input
        if showCursorAndCentroid
            switch ept.stimType
                case 'zeb'
                    % draw 2 dots, one for cursor, one for centroid
                    % dotCoords = [[walk x; walk y], [mouse x; mouse y]]; 
                    dotCoords = [[targCoordsNoResp(1,i,j); targCoordsNoResp(2,i,j)], [x; y]];
                case 'cfx'
                    % draw 3 dots; cursor, centoid, and centroid + unwrapped phase
                    % dotCoords = [[walk x; walk y], [phase x; phase y], [mouse x; mouse y]]; 
                    dotCoords = [[centerX; targCoordsNoResp(2,i,j)], ...
                        [centerX - carDrift(i,j); targCoordsNoResp(2,i,j)], [x; y]];
                otherwise
            end
            Screen('DrawDots', theWindow, dotCoords, dotSizes, dotColors);
        end

        vbl = Screen('Flip', theWindow, vbl + 0.5.*mfi);  
        
        respCoords(:,i,j) = [x; y];               % where the mouse or finger actually is
        dt.targCoordsActual(:,i,j) = [xTargNew ; yTargNew];  % where the target ended up, hopefully near the center
        
        xTargOld = xTargNew;
        yTargOld = yTargNew;

        if makeMovie && j == nTrials   % only make movie on last trial if you make movie accedentally
            % here is the bit to make the movie demo
            stimFrameCapture = Screen('GetImage', theWindow, [], [], 1, 3);
            frameName = strcat('cfxWithMouse', num2str(i), '.jpg');
            imwrite(stimFrameCapture, frameName, 'jpg');
        end

    end % end single trial stimulus for-loop

end % end 1:nTrials loops

% Close up
Screen(theWindow,'Close');
ShowCursor('Arrow');


% center coords - since this is a demo deployed in the wild, leave everything in pixels
% respCoords = where the mouse is on the desk (or the finger on the pad)
dt.normRespCoords(1,:,:) = (respCoords(1,:,:) - centerX);
dt.normRespCoords(2,:,:) = (respCoords(2,:,:) - centerY);

dt.normTargCoords(1,:,:) = (targCoordsNoResp(1,:,:) - centerX);
dt.normTargCoords(2,:,:) = (targCoordsNoResp(2,:,:) - centerY);

dt.dPhases = dPhases;

saveTheStuff(ept, dt, nTrials, fname);

myout = 1;

end % main function end

%% saving to file
function saveTheStuff(ept, dt, nTrials, fname)
    
    % for analysis, we want one row of everything per trial
    ept.subj = repmat(ept.subj, nTrials, 1);
    ept.stimType = repmat(ept.stimType, nTrials, 1);
    ept.stAngle = repmat(ept.stAngle, nTrials, 1);
    ept.carrierSpeed = repmat(ept.carrierSpeed, nTrials, 1);
    ept.freq = repmat(ept.freq, nTrials, 1);
    ept.contrast = repmat(ept.contrast, nTrials, 1);
    ept.ngabors = repmat(ept.ngabors, nTrials, 1);
    ept.stimRad = repmat(ept.stimRad, nTrials, 1);
    
    % for the data, rather than redo the above code, I'm going to permute
    % the arrays so that trials are the rows, and H and V are pages 1 and 2, repectively
    dt.normTargCoords = permute(dt.normTargCoords, [3 2 1]);
    dt.targCoordsActual = permute(dt.targCoordsActual, [3 2 1]);
    dt.normRespCoords = permute(dt.normRespCoords, [3 2 1]);
    dt.dPhases = dt.dPhases';

    if exist(fname, 'file')
        
        load(fname, 'ep', 'd');   % load the old data
        
        ep.subj = [ep.subj; ept.subj];  % append the new data...
        ep.stimType = [ep.stimType; ept.stimType];
        ep.stAngle = [ep.stAngle; ept.stAngle];
        ep.carrierSpeed = [ep.carrierSpeed; ept.carrierSpeed];
        ep.freq = [ep.freq; ept.freq];
        ep.contrast = [ep.contrast; ept.contrast];
        ep.ngabors = [ep.ngabors; ept.ngabors];
        ep.stimRad = [ep.stimRad; ept.stimRad];
        
        disp(['old: ', num2str(size(d.normTargCoords))]);   % fucking array size bug
        disp(['new: ', num2str(size(dt.normTargCoords))]);
        
        d.normTargCoords = [d.normTargCoords; dt.normTargCoords];
        d.targCoordsActual = [d.targCoordsActual; dt.targCoordsActual];
        d.normRespCoords = [d.normRespCoords; dt.normRespCoords];
        d.dPhases = [d.dPhases; dt.dPhases];
        
        save(fname, 'ep', 'd');     % and save
        
    else
        
        ep = ept;
        d = dt;
        save(fname, 'ep', 'd');
        
    end
    
end % end local save function
