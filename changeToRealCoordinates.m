function new = changeToRealCoordinates(old)
pkg load symbolic;
syms theta phi;
RyTheta = [[cos(theta), 0, sin(theta)];[0,1,0];[-sin(theta), 0, cos(theta)]];
RzPhi = [[cos(phi), sin(phi), 0]; [-sin(phi), cos(phi), 0]; [0,0,1]];
new = RyTheta'*old;
new = RzPhi'*new;
