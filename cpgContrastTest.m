% CreateProcedualGabor contrast test

% Initial parameters of gabors:

phase = 0;
SC = 400;
pSF = 0.01 * 0.0494;
contrast = 0.2; 
aspectratio = 1.0;
TW = 6*SC;
TH = 6*SC;

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
totFrames = 12*60;    %ept.frameRate;
Screen('BlendFunction', theWindow, GL_ONE, GL_ONE);

%% present the Gabor 

%backgroundOffset = [0.5 0.5 0.5 0.0];
backgroundOffset = [0,0,0,0];
disableNorm = 1;
preContrastMultiplier = 0.5;


gabortex = CreateProceduralGabor(theWindow, TW, TH,0,backgroundOffset,disableNorm,preContrastMultiplier);
mypars = repmat([phase+180, pSF, SC, contrast, aspectratio, 0, 0, 0]', 1, 1);

texrect = Screen('Rect', gabortex);
dstRects(1:4, 1) = CenterRectOnPoint(texrect, centerX, centerY)';
theseRotAngs=0.*ones(12*60, 1)

Screen('DrawTextures', theWindow, gabortex, [], dstRects, theseRotAngs(1,:), [], [], [1 1 1 0], [], kPsychDontDoRotation, mypars);
vbl = Screen('Flip', theWindow);
pause(50); 
SetMouse(centerX,centerY);
    
Screen(theWindow,'Close');
ShowCursor('Arrow');

