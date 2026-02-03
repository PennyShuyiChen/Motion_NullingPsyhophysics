function gb = gab_ecc(ecc)
% input eccentricity in visual deg 
% From the paper
    % pSF: linear relationship on log-log scale
        % b0 = log(1.99), b1 = -0.49
    % bw: eyeballed fits

pSF = exp(-0.49 * log(ecc) + log(1.99)); % cpd
%pSF = (pSF) * 0.0297; % dpp for the rig -> cpp (?)

%bw = 3.7268 + 0.6458 * log(ecc); % bw in octave 
bw =  0.5977 * log(ecc)+ 3.3380;
%bw = 2.5 + 0.5*ecc
s = bw2sig(pSF, bw);

gb = gab(20,pSF,s,90);
hold on;
plot(gb);
hold off;

