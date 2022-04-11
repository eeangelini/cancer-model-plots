from pylab import *

# pharmacodynamic parameters for dS and dR
# logistic growth rate (drug "sensitivity")
rS = 0.2 
rR = 0.2
# half-maximum point (drug "potency")
PS = 20
PR = 20 
# amplitude (drug "efficacy")
ES = 0.3
ER = 0.04
# shift constant (~ basal rate)
delS = 0.02 
delR = 0.02 

# pharmacodynamic parameters for kSR
rSR = 0.2 # sensitivity
PSR = 50 # potency
ESR = 0.5 # efficacy
kappa = 0.01 # ~ basal rate

# parameters that do not depend on drug dose
# cell division rates
bS = 0.2 
bR = 0.1 
# phenotype switching rate sentitive to resistant
kRS = 0.01 
# initial conditions
N0 = 1e9 # tumor size at detection
xi = 0.1 # fraction of resistant cells at detection
x0 = array([[(1-xi)*N0], [xi*N0]])

# set drug dose range
mvec = array(range(101))

# set range of times (for plotting trajectories)
t = linspace(0,200,num = 10000)

# for plotting total population at 4 fixed drug doses
mvals = [23,35,45,85]
stylevec = ['solid','dashdot',(0,(3,2,1,2,1,2)),(0,(5,1))]

# define logistic curves
dS = lambda m: delS + ES/(1+exp(-rS*(m-PS))) # sens death rate
dR = lambda m: delR + ER/(1+exp(-rR*(m-PR))) # res death rate
kSR = lambda m: kappa + ESR/(1+exp(-rSR*(m-PSR))) # transition S->R

# define net growth rates
gS = lambda m: bS - dS(m)
gR = lambda m: bR - dR(m)

# miscellaneous constants & eigenvalues/vectors
alpha = lambda m: gR(m)-gS(m)+kSR(m)-kRS
mybeta = lambda m: sqrt(alpha(m)**2 + 4*kSR(m)*kRS)
lam1 = lambda m: 0.5*(gS(m)+gR(m)-kSR(m)-kRS+mybeta(m))
lam2 = lambda m: 0.5*(gS(m)+gR(m)-kSR(m)-kRS-mybeta(m))
vec1 = lambda m: (-alpha(m)+mybeta(m))/(2*kSR(m))
vec2 = lambda m: (-alpha(m)-mybeta(m))/(2*kSR(m))
