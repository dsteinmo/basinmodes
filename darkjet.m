function c = darkjet(N)
% Through a jet, darkly -- darkened version of colormap jet for enchanced
% contrast near the middle

if nargin < 1
   N = size(get(gcf,'colormap'),1);
end
c = jet(N);
c = c.*repmat(1-0.5./cosh(linspace(-10,10,N)'),[1 3]);