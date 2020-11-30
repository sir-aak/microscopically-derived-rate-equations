# version: 30.05.2019

import numpy as np
#~ import sympy as sy
from scipy.optimize import fsolve

import mdre_parameters
#import mdre_processes as pr
import mdre_equations as eq
import mdre_jacobian as jac
#import mdre_fixed_point as fp

p = mdre_parameters.Parameters()


# fixed point
#~ A              = 2.48381020498
#~ phi            = -1.64037139145
#~ rho_GS_e_act   = 0.878195637535
#~ rho_GS_h_act   = 0.338997190197
#~ rho_GS_e_inact = 0.942771088427
#~ rho_GS_h_inact = 0.375689959511
#~ rho_ES_e       = 0.702003359083
#~ rho_ES_h       = 0.215928876075
#~ w_e            = 4.48559364045
#~ w_h            = 7.53617114873

#A              = 2.48
#phi            = -1.64
#rho_GS_e_act   = 0.88
#rho_GS_h_act   = 0.34
#rho_GS_e_inact = 0.94
#rho_GS_h_inact = 0.38
#rho_ES_e       = 0.70
#rho_ES_h       = 0.22
#w_e            = 4.49

# fixed point state vector
#x = np.array([A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e])






# calculation of fixed point (does not work)

#~ FP = fsolve(eq.mdre_RHS, x0)
#~ FP = fp.fixedPointNewton(eq.mdre_RHS, jac.JmdreMat, x0)
#~ FP = fp.fixedPointSecant(eq.mdre_RHS, x0)
#~ FP = fp.fixedPointBroyden(eq.mdre_RHS, jac.JmdreMat, x0)

#~ print("A              = " + str(FP[0]))
#~ print("phi            = " + str(FP[1]))
#~ print("rho_GS_e_act   = " + str(FP[2]))
#~ print("rho_GS_h_act   = " + str(FP[3]))
#~ print("rho_GS_e_inact = " + str(FP[4]))
#~ print("rho_GS_h_inact = " + str(FP[5]))
#~ print("rho_ES_e       = " + str(FP[6]))
#~ print("rho_ES_h       = " + str(FP[7]))
#~ print("w_e            = " + str(FP[8]))




# write Jacobian matrix to file
#~ def makeJacobian ():
    
    #~ J = jac.JmdreSym()
    
    #~ filename = "mdre_jacobian.txt"
    #~ file     = open(filename, "w")
    #~ file.close()
    #~ file = open(filename, "a")
    #~ i = 0
    
    #~ for entry in J:
        
        #~ if (i + 1) % 9 != 0:
            #~ file.write(str(J[i]) + ", ")
        
        #~ else:
            #~ file.write(str(J[i]) + "\n")
        
        #~ i += 1
    
    #~ file.close()


#~ makeJacobian()

A              = 2.48
phi            = -1.64
rho_GS_e_act   = 0.88
rho_GS_h_act   = 0.34
rho_GS_e_inact = 0.94
rho_GS_h_inact = 0.38
rho_ES_e       = 0.70
rho_ES_h       = 0.22
w_e            = 4.49

x = np.array([A, phi, rho_GS_e_act, rho_GS_h_act, rho_GS_e_inact, rho_GS_h_inact, rho_ES_e, rho_ES_h, w_e])


print(jac.Jmdre(x))






