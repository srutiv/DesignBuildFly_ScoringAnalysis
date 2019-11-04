% Make the scattered interpolant.
F = scatteredInterpolant(xi, yi, dataValues);
  
% Get a grid of points at everypixel location in the RGB image.
[xGrid, yGrid] = meshgrid(1:columns, 1:rows);
xq = xGrid(:);
yq = yGrid(:);
  
% Evaluate the interpolant at query locations (xq,yq).
vq = F(xq, yq);
fittedImage = reshape(vq, rows, columns);
imshow(fittedImage, []);