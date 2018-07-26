function new = change(old)
pkg load symbolic;
syms phi psi theta alpha;
R_psi = [[cos(psi),sin(psi),0];[-sin(psi),cos(psi),0];[0,0,1]];
R_beta = [[cos(phi),sin(phi),0];[-sin(phi),cos(phi),0];[0,0,1]];
R_theta = [[cos(theta+2*alpha),0,sin(theta+2*alpha)];[0,1,0];[-sin(theta+2*alpha),0,cos(theta+2*alpha)]];
new = R_psi'*old;
new = R_theta'*new;
new = R_beta'*new;
end