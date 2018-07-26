function reflectedRay = getReflectionX(incidentRay,mirrorToSuspensionPoint)
pkg load symbolic;
syms theta alpha beta psi phi;
mirrorNormal = [-1;0;0];
mirrorNormal = changeToRealCoordinates(changeToPendulumCoordinates(mirrorNormal));
reflectedRay = incidentRay - 2*(dot(mirrorNormal,incidentRay)*mirrorNormal);
syms r;
posMirror = r*mirrorNormal; 
syms l;
posOrigin = [l*sin(theta)*cos(phi);l*sin(theta)*sin(phi);l*cos(theta)]
posMirror = posOrigin + posMirror;
screenPosition = [-mirrorToSuspensionPoint;0;0];
magnitude = posMirror - screenPosition;
magnitude = magnitude(1);
cosGamma = dot([-1;0;0],reflectedRay);
magnitude = magnitude/cosGamma;
reflectedRay = magnitude*reflectedRay;

