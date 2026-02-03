function w = bwsig2freq(b,s)
% bwsig2freq
% returns the carrier frequency of a gabor function given an octave bandwidth
% and the corresponding spatial constant (sigma)
% usage: w = bw2sig(b, s) where
%				w is the frequency in cyc/unit
%				b is the full width at half-height in octaves
%				s is returned in the unit of the frequecy denominator

a = sqrt(log(2)/2);
w = (a./(pi.*s)) .* ((2.^b +1)./(2.^b-1));

