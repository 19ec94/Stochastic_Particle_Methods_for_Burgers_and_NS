function [r,theta] = cart2plr(x,y)
	%   cart2plr  Convert Cartesian coordinates to polar coordinates
	%
	%   [r,theta] = cart2plr(x,y) computes r and theta with
	%
	%       r = sqrt(x^2 + y^2);
	%       theta = atan2(y,x);
	
	r = sqrt(x^2 + y^2);
	theta = atan2(y,x);