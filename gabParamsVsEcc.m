% gab params vs. eccentricity

ecc = 1:15;
nEcc = numel(ecc);

sfs = zeros(nEcc, 1);
ss = sfs;
bws = sfs;

for i = 1:nEcc
    [gb, sfs(i), ss(i), bws(i)] = gab_ecc_2(ecc(i));
end
    
subplot(3, 1, 1)
plot(ecc, sfs)
formatFigure(" ", "spatial frequency");
subplot(3, 1, 2)
plot(ecc, ss)
formatFigure(" ", "space constant");
subplot(3, 1, 3)
plot(ecc, bws)
formatFigure(" ", "bandwidth");