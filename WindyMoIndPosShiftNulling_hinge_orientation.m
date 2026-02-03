function myout = WindyMoIndPosShiftNulling_hinge_orientation(subj, stimType, vGabs, fname)
% WindyMoIndPosShiftNulling.m
% PSC - 18-Jul-2022 orientation sensitivity test 

%% EXP - RUN!
 %WindyMoIndPosShiftNulling_hinge_orientation('oritest', 'ori', 1, 'oritest.mat');
%%

% Windy Motion Induced Position Shift Nulling
% valid stimulus types:  'zeb', 'cfo', 'cfx', 'cfh', 'itd', 'iod', 'apt', 'pld'
% 'zebra', 'cuttlefishOrtho', 'cuttlefish Ortho with x-mouse on', 'cuttlefish', 'indieTwoD', 'indieOneD', 'aperature', 'plaid'

% $$$ conditions are zeb, cfx, and iod.  pld is the plaid control condition with both H and V gratings superimposed

% vGabs = 1 for vertical Gabors, 0 for horizontal, 2 for plaids
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
%               in 1D, while the stripes jitter re. the body in 1D.  The
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
% 22-Jan-2021   psc modified it for eccentricity dependence 

% NOTES
% You must have PsychToolbox installed (psychtoolbox.org)
%
% 18-June-2018    Windy MIPS Nulling!
% Changing data saving strategy somewhat - I'm going to save out just the data,
% and then have an object handle the analysis via methods.


%% argument handling
deviceIds = PsychPowerMate('List');
handle = PsychPowerMate('Open',deviceIds,10000);

ni = nargin;
myout = 0;
% need to add checking for valid input arguments
ept.noHorzMouse = 0;

switch ni
    case 0
        error('must enter subj initials'); 
    case 4
        ept.subj = subj;
        ept.stimType = stimType;
        ept.vGabs = vGabs;
end


%% set up initial parameters
% The key flags that get set depending upon stimulus type are
% ***           ept.TwoDimJitter [0, 1] does the carrier move in both dirs [1], or only parallel to stripes [0]
% ***           ept.carrierMotion [0, 1] does the carrier phase change [1] or remain fixed relative to the envelope [0]
% ***           ept.correctCarrierMotionForEnvlope [0, 1] is the phase walk independent of the evelope walk [1]
% ***           ept.subjectChangesPhase [0, 1] does the subject's mouse control carrier phase [1]
% ***           ept.supjectChangesEnvelop [0, 1] does the subject's mouse move the envelop also too?

switch ept.stimType
    case 'ori'
        ept.TwoDimJitter = 0;
        ept.carrierMotion = 0;  % 1 = carrier Motion for MIPS, 0 = non-drifting Gabors control
        ept.carrierOriMotion = 1;
        ept.correctCarrierMotionForEnvlope = 0;
        ept.subjectChangesPhase = 0;
        ept.supjectChangesEnvelop = 0;
        ept.subjectChangesOrientation = 1;
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
        ept.TwoDimJitter = 1;
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
vdist = 25.0; % 15 original in cm (closer vdist help to present the stimuli (?)) % 57.3/56.7
H_res = 1920; V_res = 1080; % monitor resolution #info from Rect
H_im = 54.34; V_im = 30.56; % effective image size in cm(?)# info online

VA = rad2deg(atan((H_im/2) / vdist)); % Visual angle for 1/2 screen (H)
VA_pix = VA/(H_res/2); % how much visual degree a pixel extends

%% Eccentricity - PSC
% Eccentricities
E_min = 2.5;
E_max = 15;
E_ncond = 4; % # eccentricity conditions

%E = linspace(E_min,E_max,E_ncond); % Eccentricity used
E = [2.5, 5, 10, 15];
% ring radius in pix(change to ept.stimRadius af testing) - PSC
r_R = round( E/VA_pix); % ring radius in pix

%% Trials & duration

ept.trialTime = 12; % trial duration 12s * keep it 12 
nrepeat = 1;%20; % 20 repeat / condition 
ntrialsblock = E_ncond * nrepeat;
nTrials = nrepeat * E_ncond;

makeMovie = 0;
showCursorAndCentroid = 0;  % only for movies and demos

%% RF calculation
% RF size in visual angle - diameter 
for i = 1:E_ncond
    if E(i) <= 5
        RF_V(i) = 1;
        RF_MT(i) = 5;
    else
        RF_V(i) = 0.26 * E(i); %psychophysics-0.26, physiology-0.21
        RF_MT(i) = 1* E(i);
    end
end

% RF diamter in pix
d_V = RF_V/VA_pix;

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
ept.ngabors = 1;    % 1 for orientation nulling experiment
refRectSz = repmat([5],1,4);     % reference rectangle inside circle of Gabors (fixation center)


%% Gabor initial parameters:
% Phase of underlying sine grating in degrees:
phase = 0;
% assme RF = 6 * sc
SC = d_V./6;
sc = RF_V./6;
% fixed bandwidth in octave
BW = 1.4;
% Frequency of sine grating in cpp:
pSF = bwsig2freq(BW,sc);
pSF = (pSF) * VA_pix;
% Contrast of grating:
ept.contrast = 0.75; % 10% contrast% original = 10; 
% Aspect ratio width vs. height:
aspectratio = 1.0;

% Size of support in pixels
TW = round(d_V);
TH = round(d_V);

%% Experiment trial sequence and conditions - emat
        % emat columns: 1)trial sequence  2)peak Spatial Frequency
        % 3)Spatial Constant  4)ring radius (eccentricity)  5&6) Tw & Th 
        % eseq: refer to the row number in emat column 1 
emat = transpose([pSF;SC;r_R;TW;TH;E;refRectSz]);
[eseq, emat] = randeseq(emat, nrepeat, ntrialsblock); % eseq - run which row in emat
% Initialize matrix with spec for all 'ept.ngabors' patches to start off
% identically:
%mypars = repmat([phase+180, ept.freq, sc, ept.contrast, aspectratio, 0, 0, 0]', 1, ept.ngabors);

%%

ept.stepSize = 15;      % ave ept.stepSize in pixels - determines target speed - set to 0 for no envelope motion

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

mfi = Screen('GetFlipInterval', theWindow);  % should work always, unlike Screen('FrameRate')
ept.frameRate = round(1./mfi);     % get nominal frame rate
totFrames = ept.trialTime*60;    %ept.frameRate;

Screen('BlendFunction', theWindow, GL_ONE, GL_ONE);


%% ***** Generate the target sequences *****


rotAngles = ept.stAngle.*ones(totFrames, nTrials);    %*****
dPhases = drftSpeedDegPerFrame .* randn(totFrames, nTrials); % 30 * normal, so most jumps < 30 deg, very rarely >= 90

% generate random velocities, smooove them, and then integrate them to get positions.
% doing this after PTB setup because we need the framerate
temp = randn(1,totFrames, nTrials);
if ept.smooveSD
    for j = 1:nTrials
        temp(1,:,j) = conv(temp(1,:,j), SmooveKer, 'same');
        
        temp(:,1,j) = [0];
    end
end

targCoordsNoResp = zeros(1,totFrames,nTrials);
targCoordsNoResp(:,:,:) = ept.stepSize*temp;    % what target orientation will be 

respOri = zeros(size(targCoordsNoResp));     % response coordinate array
dt.targOriActual = respOri;                  % actual target positions given wind and response

%% start trial outer loop
% Build a procedural gabor texture for a gabor with a support of tw x th
% pixels and the 'nonsymetric' flag set to 1 == Gabor shall allow runtime
% change of aspect-ratio:

for j = 1:nTrials
    seq = eseq(j);
    ept.trial(j) = j;
    ept.freq(j) = emat(seq,2);
    sc = emat(seq,3);
    ept.stimRad(j) = emat(seq,4);
    tw = emat(seq,5);
    th = emat(seq,6);
    ept.ecc(j) = emat(seq,7);
    refRect = [centerX-emat(seq,8), centerY-emat(seq,8), centerX+emat(seq,8), centerY+emat(seq,8)];

    gabOffsets = ept.stimRad(j).*unitcirc(ept.ngabors);
    backgroundOffset = [0,0,0,0];
disableNorm = 1;
preContrastMultiplier = 0.5;
    gabortex = CreateProceduralGabor(theWindow, tw, th,0,backgroundOffset,disableNorm,preContrastMultiplier);
    
    mypars = repmat([phase+180, ept.freq(j), sc, ept.contrast, aspectratio, 0, 0, 0]', 1, ept.ngabors);
    
    texrect = Screen('Rect', gabortex);
    inrect = repmat(texrect', 1, ept.ngabors);
    
    dstRects = zeros(4, ept.ngabors);
    for i=1:ept.ngabors
        %     scale(i) = 1;   % *(0.1 + 0.9 * randn);  % use to mult. by texrect below to scale sizes
        dstRects(:, i) = CenterRectOnPoint(texrect, centerX+gabOffsets(i,1), centerY+gabOffsets(i,2))';
    end
    startRects = dstRects;
    % make an nTrials x ept.ngabors array of rotation angles for this trial
    
    theseRotAngs = repmat(rotAngles(:,j), 1, ept.ngabors);
    
    % Move the cursor to the center of the screen
    targCoordsNoResp(1,:,j) = targCoordsNoResp(1,:,j);

    oriTargOld = ept.stAngle;  
    dstRects = startRects;
    
    Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:), [], [], [1 1 1 0], [], kPsychDontDoRotation, mypars);
    if ept.vGabs ==2
        Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:)+90, [], [], [], [], kPsychDontDoRotation, mypars);
    end
    %     Screen('DrawDots', theWindow, [centerX, centerX; centerY, centerY], dotSizes, dotColors);
    Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
    vbl = Screen('Flip', theWindow);
    
    pause(0.4); % ready... set...
    
    tic
    for i = 1:totFrames
        
        % save last mouse position and get new position, compute change
        oriOld = oriTargOld;
        [button, dialPos] = PsychPowerMate('Get', handle)
       
        deltaDial = -[dialPos-oriOld];
        if i ==1; deltaDial = [0]; end
        
        % new position is old position + random jump due do wind
        
        % Increment phase-shift of each gabor by dPhases optionally correcting for envelope motion:
        % mypars is 8 params (last 3 dummy) x ept.ngabors
        % non-dummys are [phase+180, ept.freq, sc, ept.contrast, aspectratio]
        phaseCorrection = 0;
        if ept.correctCarrierMotionForEnvlope
            % phase correction is jump (pix) x ept.freq (cyc/pix) x 360 (deg/cyc)
            % this is correcting the phase for any envelop motion, be it from the trajectory or the subject's response.
            if ept.vGabs
                phaseCorrection = -(ept.stepSize.*temp(1, i, j)).* ept.freq(j) .* 360; %  + deltaMouse(1))
            else
                phaseCorrection = (ept.stepSize.*temp(2, i, j)) .* ept.freq(j) .* 360; %  + deltaMouse(2)
            end
        end
        
        % *** done with the wind, now incoporate the subject's mouse movements as necessary ***
        % first, handle case in which S moves the whole Gabor (as in Zebra condition)
        
        oriTargNew = oriTargOld + ept.stepSize.*temp(1, i, j);
        oriTargNew = oriTargNew + deltaDial;
       
        
        % and, finally, update dstRects and phases
        
        % ********* Okay, done with all the postion and phase adjusting, and ready get our drawing on!
        
        Screen('FrameRect', theWindow, [200, 200, 200], refRect, 1);
        Screen('DrawTextures', theWindow, gabortex, [], dstRects, oriTargNew, [], [], [], [], kPsychDontDoRotation, mypars);
        if ept.vGabs ==2
            Screen('DrawTextures', theWindow, gabortex, [], dstRects, oriTargNew+90, [], [], [], [], kPsychDontDoRotation, mypars);
        end
        
        vbl = Screen('Flip', theWindow, vbl + 0.5.*mfi);
        respOri(:,i,j) = dialPos;
        dt.targOriActual(:,i,j) = oriTargNew;
        
        oriTargOld = oriTargNew;
    
        
    end % end single trial stimulus for-loop
    
end % end 1:nTrials loops

% Close up
Screen(theWindow,'Close');
ShowCursor('Arrow');

% center coords - since this is a demo deployed in the wild, leave everything in pixels
% respCoords = where the mouse is on the desk (or the finger on the pad)
dt.RespOris(1,:,:) = respOri(1,:,:);
dt.TargOris(1,:,:) = targCoordsNoResp(1,:,:);

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
ept.freq = ept.freq;%repmat(ept.freq, nTrials, 1);
ept.contrast = repmat(ept.contrast, nTrials, 1);
ept.ngabors = repmat(ept.ngabors, nTrials, 1);
ept.stimRad = ept.stimRad; %repmat(ept.stimRad, nTrials, 1);
ept.ecc = ept.ecc;
ept.vGabs = repmat(ept.vGabs,nTrials,1);

% for the data, rather than redo the above code, I'm going to permute
% the arrays so that trials are the rows, and H and V are pages 1 and 2, repectively
dt.TargOris = permute(dt.TargOris, [3 2 1]);
dt.targOriActual = permute(dt.targOriActual, [3 2 1]);%given wind and response
dt.RespOris = permute(dt.RespOris, [3 2 1]);
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
    ep.ecc = [ep.ecc;ept.ecc];
    ep.vGabs = [ep.vGabs;ept.vGabs];
    
    disp(['old: ', num2str(size(d.TargOris))]);   % fucking array size bug
    disp(['new: ', num2str(size(dt.TargOris))]);
    
    d.TargOris = [d.TargOris; dt.TargOris];
    d.targOriActual = [d.targOriActual; dt.targOriActual];
    d.RespOris = [d.RespOris; dt.RespOris];
    d.dPhases = [d.dPhases; dt.dPhases];
    
    save(fname, 'ep', 'd');     % and save
    
else
    
    ep = ept;
    d = dt;
    save(fname, 'ep', 'd');
    
end



end % end local save function
