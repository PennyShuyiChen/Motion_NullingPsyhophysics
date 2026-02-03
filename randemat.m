function [eseq,emat] = randemat(emat, nrepeat,ntrialsblock)


emat = repmat(emat,nrepeat,1);
[eseq, emat] = randeseq(emat, ntrialsblock);

