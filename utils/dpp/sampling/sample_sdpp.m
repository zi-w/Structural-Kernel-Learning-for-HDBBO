function Y = sample_sdpp(model,C,k)
% sample from an SDPP. 
% model contains the details of the SDPP (see helpers/bp.m for spec)
% C is the decomposed covariance matrix, computed using:
%   C = decompose_kernel(bp(model,'covariance'));
% k is (optionally) the size of the set to return.

if ~exist('k','var')  
  % choose eigenvectors randomly
  D = C.D ./ (1+C.D);
  v = find(rand(length(D),1) <= D);
else
  % k-DPP
  v = sample_k(C.D,k);
end
k = length(v);
V = C.V(:,v);

% rescale eigenvectors so they normalize in the projected space
V = bsxfun(@times,V,1./sqrt(C.D(v)'));

% iterate
Y = zeros(model.T,k);
for i = k:-1:1

  % choose a labeling
  Y(:,i) = bp(model,'sample2',V);

  % choose a vector to eliminate
  S = sum(model.G(Y(:,i),:),1) * V;
  j = find(S,1);
  Vj = V(:,j);
  Sj = S(j);
  V = V(:,[1:j-1 j+1:end]);
  S = S(:,[1:j-1 j+1:end]);

  % update V
  V = V - bsxfun(@times,Vj,S/Sj);

  % orthogonalize in the projected space
  for a = 1:i-1
    for b = 1:a-1
      V(:,a) = V(:,a) - (V(:,a)'*C.M *V(:,b))*V(:,b);
    end
    V(:,a) = V(:,a) / sqrt(V(:,a)'*C.M*V(:,a));
  end

end
