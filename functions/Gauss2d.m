function F = Gauss2d(a, data)
% a(1) - A
% a(2) - Xpos
% a(3) - Ypos
% a(4) - sigmaX
% a(5) - sigmaY
% a(6) - Rotation Angle
% a(7) - B
% a(6) = 0;
x0rot = a(2)*cos(a(6)) - a(3)*sin(a(6));
y0rot = a(2)*sin(a(6)) + a(3)*cos(a(6));

[~,sizex] = size(data);
sizex= sizex/2;
X = data(:,1:sizex);
Y = data(:,sizex+1:end);

Xrot = X*cos(a(6)) - Y*sin(a(6));
Yrot = X*sin(a(6)) + Y*cos(a(6));

expPart = exp(-(Xrot-x0rot).^2./(2.*(a(4).^2))-(Yrot-y0rot).^2./(2.*(a(5).^2)));
F =  a(1)*expPart + a(7);
end