

phase = 0;
% Frequency of sine grating:
pSF = exp(-0.49 * log(E(test)) + log(1.99*3.1)); % cpd
ept.freq = (pSF) * VA_pix; % dpp for the rig -> cpp (?)
% Spatial constant of the exponential "hull"
bw = 3.7268 + 0.6458 * log(E(test)); % bw in octave 
sc = 3.1 * bw2sig(ept.freq, bw);
% Size of support in pixels
tw = round(6 * sc);
th = round(6 * sc)
pSF



phase = 0;
% Frequency of sine grating:
pSF = exp(-0.49 * log(E(test)) + log(1.99*2.5)); % cpd
ept.freq = (pSF) * VA_pix; % dpp for the rig -> cpp (?)
% Spatial constant of the exponential "hull"
bw = 1.35* (0.6808 + 0.0930 * log(E(test))); % bw in octave 
sc = bw2sig(ept.freq, bw);
% Size of support in pixels
tw = round(6 * sc);
th = round(6 * sc)
pSF


