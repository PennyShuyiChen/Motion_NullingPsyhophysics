function [rseq, emat] = randeseq(emat, nrepeat, ntrialsblock)
% exp mat with randomized trials 	
emat = repmat(emat,nrepeat,1);

n = size(emat, 1);
rseq = [];
    
	while length(rseq) < n
		rseq = [rseq, length(rseq) + randperm(ntrialsblock)];
        
	end
	rseq = rseq(1:n); % the row in the emat

	idseq(rseq) = 1:length(rseq); 
	emat = [idseq', emat]; % col 1 = trial number 
end