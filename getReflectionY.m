function reflectedRay = getReflectionY(incidentRay,mirrorToSuspensionPoint)
pkg load symbolic;
syms theta alpha psi phi;
mirrorNormal = [1;0;0];
mirrorNormal = change(mirrorNormal);
reflectedRay = incidentRay - 2*(dot(mirrorNormal,incidentRay)*mirrorNormal);
syms r;
posMirror = r*mirrorNormal; 
syms l;
posOrigin = [l*sin(theta)*cos(phi)+l*tan(alpha)*cos(theta)*cos(phi);l*sin(theta)*sin(phi)+l*tan(alpha)*cos(theta)*sin(phi);-l*cos(theta)+ l*tan(alpha)*sin(theta)];
posMirror = posOrigin + posMirror;
screenPosition = [mirrorToSuspensionPoint;0;0];
magnitude = posMirror - screenPosition;
magnitude = magnitude(1);
cosGamma = dot([1;0;0],reflectedRay);
magnitude = magnitude/cosGamma;
reflectedRay = magnitude*reflectedRay;

