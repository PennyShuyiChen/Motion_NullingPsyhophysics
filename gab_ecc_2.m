function [gb, pSF, s, bw] = gab_ecc_2(ecc)
% input eccentricity in visual deg 
% From the paper
    % pSF: linear relationship on log-log scale
        % b0 = log(1.99), b1 = -0.49
    % bw: eyeballed fits

pSF = exp(-0.49 * log(ecc) + log(1.99)); % cpd
%pSF = (pSF) * 0.0297; % dpp for the rig -> cpp (?)

bw = 3.7268 + 0.6458 * log(ecc); % bw in octave 
s = bw2sig(pSF, bw);

gb = gab(200,pSF,s);
%plot(gb);

