function new = changeToPendulumCoordinates(old)
pkg load symbolic;
syms alpha beta psi;
RzPsi = [[cos(psi), sin(psi), 0]; [-sin(psi), cos(psi), 0]; [0,0,1]] ;
RyBeta = [[cos(beta), 0, -sin(beta)];[0,1,0];[sin(beta), 0, cos(beta)]];
RxAlpha = [[1,0,0];[0, cos(alpha), sin(alpha)]; [0, -sin(alpha), cos(alpha)]];
new = RyBeta'*old;
new = RxAlpha'*new;
new = RzPsi'*new;
