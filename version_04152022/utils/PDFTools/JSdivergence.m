function [D_js] = JSdivergence(P,Q)

%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1

if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

PDFSUM = (Q+P)/2;

KL1 = KLdivergence(P, PDFSUM);
KL2 = KLdivergence(Q, PDFSUM);

D_js = (1/2) * (KL1 + KL2);

end
