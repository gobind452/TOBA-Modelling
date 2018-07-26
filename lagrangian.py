import sympy as sp
import math

class Tensor: ## Tensor Class
    def __init__(self,rows,columns):
        self.rows = rows
        self.columns = columns
        self.contents = []
    def add(self,tensor):
        assert self.rows == tensor.rows
        assert self.columns == tensor.columns
        answer = Tensor(self.rows,self.columns)
        for i in range(self.rows):
            currentRow = []
            for j in range(self.columns):
                currentRow.append(self.contents[i][j]+tensor.contents[i][j])
            answer.contents.append(currentRow[:])
        return answer
    def multiply(self,tensor):
        assert self.columns == tensor.rows
        answer = Tensor(self.rows,tensor.columns)
        for i in range(self.rows):
            currentRow = []
            for j in range(tensor.columns):
                sumd = 0
                for m in range(self.columns):
                    sumd = sumd + self.contents[i][m]*tensor.contents[m][j]
                currentRow.append(sumd)
            answer.contents.append(currentRow[:])
        return answer
    def transpose(self):
        answer = Tensor(self.columns,self.rows)
        for i in range(self.columns):
            currentRow = []
            for j in range(self.rows):
                currentRow.append(self.contents[j][i])
            answer.contents.append(currentRow[:])
        return answer
    
    def getScalar(self):
        assert self.rows == 1
        assert self.columns == 1
        return self.contents[0][0]

    def changeToRealCoordinates(self):
        assert self.rows == 3
        assert self.columns == 1
        RyTheta = Tensor(3,3)
        RyTheta.contents = [[1-theta(t)*theta(t)/2, 0, theta(t)],[0,1,0],[-theta(t), 0, 1-theta(t)*theta(t)/2]]
        RzPhi = Tensor(3,3)
        RzPhi.contents = [[sp.cos(phi(t)), sp.sin(phi(t)), 0],[-sp.sin(phi(t)), sp.cos(phi(t)), 0], [0,0,1]];
        new = (RyTheta.transpose()).multiply(self)
        new = (RzPhi.transpose()).multiply(new)
        return new
    
    def changeToPendulumCoordinates(self):
        assert self.rows == 3
        assert self.columns == 1
        RzPsi = Tensor(3,3)
        RzPsi.contents = [[1- psi(t)*psi(t)/2, psi(t), 0],[-psi(t), 1-psi(t)*psi(t)/2, 0],[0,0,1]]
        RyBeta = Tensor(3,3)
        RyBeta.contents = [[1-beta(t)*beta(t)/2, 0, -beta(t)],[0,1,0],[beta(t), 0, 1-beta(t)*beta(t)/2]]
        RxAlpha = Tensor(3,3)
        RxAlpha.contents = [[1,0,0],[0, 1- alpha(t)*alpha(t)/2, alpha(t)],[0, -alpha(t),1-alpha(t)*alpha(t)/2]]
        new = (RyBeta.transpose()).multiply(self)
        new = (RxAlpha.transpose()).multiply(new)
        new = (RzPsi.transpose()).multiply(new)
        return new

    def changeToRealCoordinates2(self):
        assert self.rows == 3
        assert self.columns == 1
        RzPsi = Tensor(3,3)
        RzPhi = Tensor(3,3)
        RzPsi.contents = [[sp.cos(psi(t)), sp.sin(psi(t)), 0],[-sp.sin(psi(t)), sp.cos(psi(t)), 0],[0,0,1]]
        RxAlphaTheta = Tensor(3,3)
        RzPhi.contents = [[sp.cos(phi(t)),sp.sin(phi(t)),0],[-sp.sin(phi(t)),sp.cos(phi(t)),0],[0,0,1]]
        RxAlphaTheta.contents = [[sp.cos(2*alpha(t)+theta(t)),0, sp.sin(2*alpha(t)+theta(t))],[0,1,0],[-sp.sin(2*alpha(t)+theta(t)),0,sp.cos(2*alpha(t)+theta(t))]]
        new = (RzPsi.transpose()).multiply(self)
        new = (RxAlphaTheta.transpose()).multiply(new)
        new = (RzPhi.transpose()).multiply(new)
        return new

class Vector(Tensor): ##Vector Class
    def __init__(self,x,y,z):
        self.rows = 3
        self.columns = 1
        self.contents = [[x],[y],[z]]
        
    def dot(self,vector):
        answer = Tensor(1,1)
        sumd = 0
        for i in range(3):
            sumd = sumd + self.contents[i][0]*vector.contents[i][0]
        answer.contents = [[sumd]]
        return answer
    
    def cross(self,vector):
        answer = Vector(0,0,0)
        for i in range(3):
            j = i+1
            k = i+2
            answer.contents[i][0] = (self.contents[j%3][0]*vector.contents[k%3][0] - vector.contents[j%3][0]*self.contents[k%3][0])
        return answer

def changeToVector(tensor):
     assert tensor.rows == 3
     assert tensor.columns == 1
     return Vector(tensor.contents[0][0],tensor.contents[1][0],tensor.contents[2][0])
    
def getLagrangian1(): ## Model without bending of wire. Assumes free rotation of the block about suspension point
    InertiaTensor = Tensor(3,3)
    InertiaTensor.contents = [[I_xx,I_xy,I_xz],[I_xy,I_yy,I_yz],[I_xz,I_yz,I_zz]]
    AngularVelocity = Vector(alphaD(t),betaD(t),psiD(t))
    Rotational = (((AngularVelocity.transpose()).multiply(InertiaTensor)).multiply(AngularVelocity)).getScalar()
    Pos = Vector(0,-offSet,-h)
    VelocitySuspensionPoint = Vector(l*thetaD(t),l*theta(t)*phiD(t),0)
    VelocityCOMFromSuspension = AngularVelocity.cross(Pos)
    VelocityCOMFromSuspension = changeToVector(VelocityCOMFromSuspension.changeToPendulumCoordinates())
    CrossTerm = (VelocityCOMFromSuspension.dot(VelocitySuspensionPoint)).getScalar()
    Pos = changeToVector((Pos.changeToPendulumCoordinates()).changeToRealCoordinates())
    Pos = Pos.add(Vector(l*theta(t)*sp.cos(phi(t)),l*theta(t)*sp.sin(phi(t)),-l*(1-theta(t)*theta(t)/2)))
    BendingEnergy = (Modulus)*(InertiaWire)*(alpha(t)*alpha(t) + beta(t)*beta(t))/(l-lamda)
    lagrangian = M*((l*thetaD(t))*(l*thetaD(t)) + (l*theta(t))*phiD(t)*(l*theta(t))*phiD(t))+ Rotational+ 2*M*CrossTerm - k*psi(t)*psi(t) - 2*M*g*Pos.contents[2][0] - BendingEnergy
    return lagrangian


def getLagrangian2(Hamiltonian): ## Model with bending of wire ( More accurate)
    Pos = generatePosition()
    Velocity = differentiate(Pos)
    BendingEnergy = BendingK*(alpha(t)*alpha(t))/2
    PosCOM = Vector(offSet,0,-h)
    AngularVelocity = Vector(0,0,psiD(t))
    VelocityCOM = AngularVelocity.cross(PosCOM)
    VelocityCOM = changeToVector(VelocityCOM.changeToRealCoordinates2())
    CrossTerm = (VelocityCOM.dot(Velocity)).getScalar()
    PosCOM = changeToVector(PosCOM.changeToRealCoordinates2())
    Pos = Pos.add(PosCOM)
    KineticEnergy = 0.5*M*Velocity.dot(Velocity).getScalar() + 0.5*I_zz*psiD(t)*psiD(t)
    PotentialEnergy = M*g*Pos.contents[2][0] + BendingEnergy + 0.5*TorsionK*psi(t)*psi(t)
    ##print PotentialEnergy
    ##print KineticEnergy + M*CrossTerm
    Hamiltonian = PotentialEnergy + KineticEnergy + M*CrossTerm
    lagrangian= KineticEnergy + M*CrossTerm - PotentialEnergy
    return lagrangian
 
def generatePosition(): ## Helper function
    Pos = Vector(0,0,0)
    Pos.contents[0][0] = l*sp.sin(theta(t))*sp.cos(phi(t)) + l*sp.sin(alpha(t))*sp.cos(theta(t))*sp.cos(phi(t))
    Pos.contents[1][0] = l*sp.sin(theta(t))*sp.sin(phi(t)) + l*sp.sin(alpha(t))*sp.cos(theta(t))*sp.sin(phi(t))
    Pos.contents[2][0] = l*sp.sin(alpha(t))*sp.sin(theta(t)) - l*sp.cos(theta(t))
    return Pos

def KE(): ## Helper function
    Pos1 = Vector(sp.sin(theta(t)),0,sp.cos(theta(t))) ##-l*sin*alphadot
    Velocity1 = differentiate(Pos1) ##l*cosalpha 
    Velocity2 = Vector(sp.cos(theta(t)),0,sp.sin(theta(t)))##l cos (alpha) alphadot
    Velocity3 = differentiate(Vector(sp.cos(theta(t)),0,sp.sin(theta(t)))) ## lsin alpha
    Term = l*l*sp.sin(alpha(t))*sp.sin(alpha(t))*alphaD(t)*alphaD(t)+SwingingLength*SwingingLength*(thetaD(t)*thetaD(t)) + sp.cos(alpha(t))*sp.cos(alpha(t))*alphaD(t)*alphaD(t)*l*l + ProtudingLength*ProtudingLength*(thetaD(t)*thetaD(t))
    CrossTerm = SwingingLength*alphaD(t)*l*(Velocity1.dot(Velocity2).getScalar()) + SwingingLength*ProtudingLength*(Velocity1.dot(Velocity3).getScalar()) - l*sp.sin(alpha(t))*alphaD(t)*ProtudingLength*(Pos1.dot(Velocity3).getScalar())
    return Term + 2*CrossTerm

def differentiate(Pos): ##Helper
    Velocity = Vector(0,0,0)
    Velocity.contents[0][0] = Pos.contents[0][0].diff(t)
    Velocity.contents[0][0] = Velocity.contents[0][0].subs([(alpha(t).diff(t),alphaD(t)),(theta(t).diff(t),thetaD(t)),(psi(t).diff(t),psiD(t)),(phi(t).diff(t),phiD(t))])
    Velocity.contents[1][0] = Pos.contents[1][0].diff(t)
    Velocity.contents[1][0] = Pos.contents[1][0].diff(t).subs([(alpha(t).diff(t),alphaD(t)),(theta(t).diff(t),thetaD(t)),(psi(t).diff(t),psiD(t)),(phi(t).diff(t),phiD(t))])
    Velocity.contents[2][0] = Pos.contents[2][0].diff(t)
    Velocity.contents[2][0] = Velocity.contents[2][0].subs([(alpha(t).diff(t),alphaD(t)),(theta(t).diff(t),thetaD(t)),(psi(t).diff(t),psiD(t)),(phi(t).diff(t),phiD(t))])
    return Velocity

def diffT(term):
    return (term.diff(t)).subs([(phi(t).diff(t),phiD(t)),(phiD(t).diff(t),phiDD(t)),(alpha(t).diff(t),alphaD(t)),(psi(t).diff(t),psiD(t)),(theta(t).diff(t),thetaD(t)),(alphaD(t).diff(t),alphaDD(t)),(psiD(t).diff(t),psiDD(t)),(thetaD(t).diff(t),thetaDD(t))])
   
def smallAngle(term): ## Applies the small angle approximation
    return term.subs([(sp.sin(phi(t)),phi(t)),(sp.cos(phi(t)),1),(sp.sin(2*alpha(t)+theta(t)),(2*alpha(t)+theta(t))),(sp.cos(2*alpha(t)+theta(t)),1),(sp.sin(theta(t)),theta(t)),(sp.cos(theta(t)),1),(sp.sin(alpha(t)),alpha(t)),(sp.cos(alpha(t)),1),(sp.sin(psi(t)),psi(t)),(sp.cos(psi(t)),1),(sp.sin(2*alpha(t)),2*alpha(t)),(sp.cos(2*alpha(t)),1),(sp.sin(2*theta(t)),2*theta(t)),(sp.cos(2*theta(t)),1)])

def getEquations(lagrangian): ## Extracts the equations
    eqn1 = sp.Eq(lagrangian.diff(psi(t)),diffT(lagrangian.diff(psiD(t)))+TorsionY*psiD(t))
    eqn2 = sp.Eq(lagrangian.diff(alpha(t)),diffT(lagrangian.diff(alphaD(t))))
    #eqn3 = sp.Eq(lagrangian.diff(beta(t)),diffT(lagrangian.diff(betaD(t))))
    eqn4 = sp.Eq(lagrangian.diff(theta(t)),diffT(lagrangian.diff(thetaD(t)))+ PendulumY*thetaD(t))
    eqn5 = sp.Eq(lagrangian.diff(phi(t)),diffT(lagrangian.diff(phiD(t))))
    eqns = [eqn1,eqn2,eqn4,eqn5]
    for i in range(4):
        arg1 = smallAngle(eqns[i].args[0])
        arg2 = smallAngle(eqns[i].args[1])
        eqns[i] = sp.Eq(arg1,arg2)
    return eqns

def partial(expr,func,value): ##Partial derivative the term to return only linear terms
    if value is 0:
        return ((expr.diff(func)).subs([(alpha(t),value),(beta(t),value),(psi(t),value),(theta(t),value),(alphaD(t),value),(betaD(t),value),(psiD(t),value),(thetaD(t),value)]))*func
    else:
        return ((expr.diff(func)).subs([(phi(t),0),(phiD(t),0),(phiDD(t),0),(alpha(t),0),(psi(t),0),(theta(t),0),(alphaD(t),0),(psiD(t),0),(thetaD(t),0),(psiDD(t),0),(thetaDD(t),0),(alphaDD(t),0)]))*func

def lineariseEq(eqns,value): ## lineariseEq
     answer = []
    for eqn in eqns:
        arg1 = eqn.args[0]
        arg2 = eqn.args[1]
        arg1 = partial(arg1,alpha(t),value) + partial(arg1,psi(t),value) + partial(arg1,theta(t),value) + partial(arg1,alphaD(t),value) + partial(arg1,thetaD(t),value) + partial(arg1,psiD(t),value) + arg1.diff(psiDD(t))*psiDD(t) + arg1.diff(thetaDD(t))*thetaDD(t) + arg1.diff(alphaDD(t))*alphaDD(t) + partial(arg1,phiD(t),value) + partial(arg1,phiDD(t),value)
        arg2 = partial(arg2,alpha(t),value) + partial(arg2,psi(t),value) + partial(arg2,theta(t),value) + partial(arg2,alphaD(t),value) + partial(arg2,thetaD(t),value) + partial(arg2,psiD(t),value) + arg2.diff(psiDD(t))*psiDD(t) + arg2.diff(thetaDD(t))*thetaDD(t) + arg2.diff(alphaDD(t))*alphaDD(t) + partial(arg2,phiD(t),value) + partial(arg2,phiDD(t),value)
        answer.append(sp.Eq(arg1,arg2))
    return answer
    
TorsionK,BendingK,TorsionY = sp.symbols('TorsionK,BendingK,TorsionY')
t = sp.symbols('t')
alpha,beta,psi,phi,theta = sp.symbols('alpha,beta,psi,phi,theta',cls = sp.Function)
alphaD,betaD,psiD,thetaD,phiD = sp.symbols('alphaD,betaD,psiD,thetaD,phiD',cls = sp.Function)
alphaDD,betaDD,psiDD,thetaDD,phiDD = sp.symbols('alphaDD,betaDD,psiDD,thetaDD,phiDD',cls = sp.Function)

def lineariseTerm(term,value):
    q = [sp.Eq(0,term)]
    q = lineariseEq(q,value)
    return q[0].args[1]

def secondOrderEq(eqn,value): ## Second Order Eq (For the azimuthal angle)
    func1 = [alpha(t),beta(t),psi(t),theta(t),alphaD(t),betaD(t),psiD(t),thetaD(t),alphaD(t).diff(t),betaD(t).diff(t),psiD(t).diff(t),thetaD(t).diff(t),phiD(t),phiD(t).diff(t)]
    arg = eqn.args[1]
    finalArg = 0
    if value is 0:
        for i in range(len(func1)):
            for j in range(len(func1)):
                func = func1[i]*func1[j]
                finalArg = finalArg+ (arg.diff(func1[i])).diff(func1[j]).subs([(alpha(t),value),(beta(t),value),(psi(t),value),(theta(t),value),(alphaD(t),value),(betaD(t),value),(psiD(t),value),(thetaD(t),value)])*func
    else:
        for i in range(len(func1)):
            for j in range(len(func1)):
                func = func1[i]*func1[j]
                finalArg = finalArg + (arg.diff(func1[i])).diff(func1[j]).subs([(alpha(t),0),(psi(t),0),(theta(t),0),(alphaD(t),0),(psiD(t),0),(thetaD(t),0)])*func
    return sp.Eq(0,finalArg)

def secondOrderTerm(term,value):
    q = sp.Eq(0,term)
    q = secondOrderEq(q,value)
    if len(q.args) <= 1:
        return 0
    return q.args[1]

'''m = 0.048
R = 0.03
a = 0.006
b = 0.065
c = 0.025
x = 0.065*0.9063077870366499
y = 0.065*0.42261826174069944
I_yy = M*R*R/4 + 4*(m*(a*a*0.42261826174069944*0.42261826174069944+b*b*0.9063077870366499*0.9063077870366499+c*c)/12 + m*x*x)
I_xx = M*R*R/4+4*(m*(a*a*0.9063077870366499*0.9063077870366499+b*b*0.42261826174069944*0.42261826174069944+c*c)/12 + m*y*y)
I_xy = 4*(m*(a*a-b*b)*0.42261826174069944*0.9063077870366499/12 - x*y*m)
I_yz = 4*(m*y*c/2)
I_xz = 4*(m*x*c/2)
M = M+4*m'''
l=0.125
I_zz = 0.0016
M = 1
g = 9.8
Q = 100000
D = 0.001
offSet = 0.001
h = 0.05
InertiaWire = (math.pi)*D*D*D*D/32
PendulumQ = 4*math.pow(10,4)
PendulumY = 2*M*l*l/PendulumQ
SwingingLength = l*sp.cos(alpha(t))
ProtudingLength = l*sp.sin(alpha(t))
Hamiltonian = 0
lagrangian = getLagrangian2(Hamiltonian)
##lagrangian = lagrangian.subs([(phi(t),0),(beta(t),0),(phiD(t),0),(betaD(t),0)])
eqns = getEquations(lagrangian)
new = lineariseEq(eqns,1)



'''P,F,psi,theta,vb = sp.symbols('P,F,psi,theta,vb' , cls = sp.Function)
eqn1 = sp.Eq(theta(t).diff(t),P(t))
eqn2 = sp.Eq(psi(t).diff(t),F(t))
eqn3 = sp.Eq(P(t).diff(t)*I_zz + g*I_zz*theta(t)/l,Y*b*F(t)/l + k*psi(t)*b/l+I_zz*psi(t)*b*g/(l*l))
eqn4 = sp.Eq(I_zz*F(t).diff(t)+Y*F(t)+(k-M*g*b*b/l)*psi(t),M*g*b*theta(t))
'''

