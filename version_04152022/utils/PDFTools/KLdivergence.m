function [D_kl] = KLdivergence(P,Q)

%  dist = KLDiv(P,Q) Kullback-Leibler divergence of two discrete probability
%  distributions
%  P and Q  are automatically normalised to have the sum of one on rows
% have the length of one at each 
% P =  n x nbins
% Q =  1 x nbins or n x nbins(one to one)
% dist = n x 1
% figure;
% subplot(2,1,1)
% plot(P,'r-'); hold on; plot(Q,'k-'); title('input P/Q')

if size(P,2)~=size(Q,2)
    error('the number of columns in P and Q should be the same');
end

if sum(~isfinite(P(:))) + sum(~isfinite(Q(:)))
   error('the inputs contain non-finite values!') 
end

% normalizing the P and Q
if size(Q,1)==1
    Q = Q ./sum(Q);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    
    % DO KL
    temp =  P.*log(P./repmat(Q,[size(P,1) 1]));
    temp(isnan(temp))=0;% resolving the case when P(i)==0
    dist = sum(temp,2);
    
elseif size(Q,1)==size(P,1)
    
    Q = Q ./repmat(sum(Q,2),[1 size(Q,2)]);
    P = P ./repmat(sum(P,2),[1 size(P,2)]);
    
%     temp =  P.*log(P./Q);
%     temp(isnan(temp))=0; % resolving the case when P(i)==0
%     dist = sum(temp,2);
    
    M=log(P./Q); 
    M(isnan(M))=0; 
    dist = sum(P.*M,2);
end

D_kl = dist;
% subplot(2,1,2)
% plot(P,'r-'); hold on; plot(Q,'k-'); title('normalized P/Q')
end

