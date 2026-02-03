function JustBkMouse(wptr, secs,text)
    if exist('text', 'var') == 0 
		text = 'Need Coffee? Break Time! % 2d:%02d';
    end

	rect = Screen('Rect', wptr);
	[cx, cy] = RectCenter(rect);
	[smin, ssec] = sec2min(secs);
	stext = sprintf(text, smin, ssec);
	tbound = Screen('TextBounds', wptr, stext);
	[tcx, tcy] = RectCenter(tbound);
	tx = cx - tcx;
	ty = cy - tcy;

	Screen('DrawText', wptr, stext, tx, ty);
	[timestamp, breakOnset] = Screen('Flip', wptr);
	breakEnd = breakOnset + secs + 1;

    while GetSecs < breakEnd
        [min, sec] = sec2min(secs - floor(GetSecs - breakOnset));
%         [clicks, x,y, button] = GetClicks;
%         if button == 1
%             break;
%         end
        Screen('DrawText', wptr, sprintf(text, min, sec), tx, ty);
        Screen('Flip', wptr);
    end
% 	WaitSecs(0.5);
end


function [min, sec] = sec2min(sec)
	min = floor(sec / 60);
	sec = sec - min * 60;
end