# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 10:25:11 2020

@author: Saumya Dholakia
"""

# Importing the requisite files and libraries:
import numpy as np
import matplotlib.pyplot as plt
from Riemann_mod import Riemann
import math

# Extracting the initial values of the conservative variales from the shock tube problem:
if __name__ == '__main__':
    R = Riemann()
    rhoL,uL,pL,rhoR,uR,pR,max_time = R.get_cases()
    case = 1
    gam = 1.4

# Geometrical parameters:   
# When the number of nodes are given:
#L=1
#nx = 25
#ny = 25
#dx = L/(nx-1)
#dy = L/(ny-1)
#x = np.linspace(0,L,nx)
#y = np.linspace(0,L,ny)
# When the interval length is given:
L=1
xc = L/2
yc = L/2
Nx = 50
Ny = 50
dx = L/Nx
dy = L/Ny
dA = dx*dy
x = np.arange(0,L+dx,dx)
y = np.arange(0,L+dy,dy)
nx = np.size(x)
ny = np.size(y)
Rc = L/10
cmax = 0.8

# Generic constants:
Rgas = 287.
p0 = 101300
T0 = 298
yint = int(ny/2)

# Project 6:
# Farfield values:
uinf = uR[case]
vinf = 0
pinf = pR[case]
Tinf = pR[case]/(Rgas*rhoR[case])
rhoinf = rhoR[case]
# Thermal parameters:
aL=((gam*pL[case])/rhoL[case])**(1/2)
aR=((gam*pR[case])/rhoR[case])**(1/2)
SR = R.find_star_state(gam,rhoL[case],pL[case],aL,uL[case],rhoR[case],pR[case],aR,uR[case])[9]
Mac = 0.3
#Max = uR[case]/aR
temp1 = 0.5*(gam - 1)*(Mac**2)
temp2 = temp1/(1 + temp1)
Tc = Tinf/(1 + temp1)
ac = (gam*Rgas*Tc)**0.5
# Spatial parameters:
xth = L/4 # It controls the initialization in space
xs = L/4 # It controls the end time
ns = xs/dx
# Temporal paramters:
time_end = xs/SR
#time_end = max_time[case]
# Sod tube:
#rhoexact,pexact,uexact,eexact = R.plot(xs,gam,time_end,case,"rusanov"+str(case),nx)

# Project 5:
# Baseline case:
# Thermal parameters:
#Mac = 0.3
#Max = 0
#temp1 = 0.5*(gam - 1)*(Mac**2)
#temp2 = temp1/(1 + temp1)
#Tc = T0/(1 + temp1)
#ac = (gam*Rgas*Tc)**0.5
#a0 = (gam*Rgas*T0)
## Farfield values:
#uinf = 0
#vinf = 0
## Temporal parameters:
#time_end = Rc/ac

# X convection:
# Thermal parameters:
#Mac = 0.3
#Max = 0.3
#May = 0
#temp1 = 0.5*(gam - 1)*(Mac**2)
#temp2 = temp1/(1 + temp1)
#Tc = T0/(1 + temp1)
#ac = (gam*Rgas*Tc)**0.5
#a0 = (gam*Rgas*T0)
## Farfield values:
#uinf = Max*a0
#vinf = 0
## Temporal parameters:
#time_end = Rc/(Max*a0)

# Y convection:
# Thermal parameters:
#Mac = 0.3
#Max = 0
#May = 0.3
#temp1 = 0.5*(gam - 1)*(Mac**2)
#temp2 = temp1/(1 + temp1)
#Tc = T0/(1 + temp1)
#ac = (gam*Rgas*Tc)**0.5
#a0 = (gam*Rgas*T0)
## Farfield values:
#uinf = 0
#vinf = May*a0
## Temporal parameters:
#time_end = Rc/(May*a0)

# Diagonal convection:
# Thermal parameters:
#Mac = 0.3
#Max = 0.3/math.sqrt(2) 
#May = 0.3/math.sqrt(2) 
#temp1 = 0.5*(gam - 1)*(Mac**2)
#temp2 = temp1/(1 + temp1)
#Tc = T0/(1 + temp1)
#ac = (gam*Rgas*Tc)**0.5
#a0 = (gam*Rgas*T0)
## Farfield values:
#uinf = Max*a0
#vinf = May*a0
## Temporal parameters:
#time_end = (math.sqrt(2)*L)/(math.sqrt(uinf**2 + vinf**2))
## Normalized time:
##time_end = (2.*(math.pi)*Rc)/ac

# Problem initialization:
def initial(rhoL,uL,pL,rhoR,uR,pR,nx,ny,xth,x,y,xc,yc,Rc,gam,Rgas,Mac,ac,uinf,vinf,pinf,Tinf,temp1,temp2,T0,p0):
    u = np.zeros((nx,ny))
    v = np.zeros((nx,ny))
    p = np.zeros((nx,ny))
    T = np.zeros((nx,ny))
    rho = np.ones((nx,ny))
    xstar = np.zeros(nx)
    ystar = np.zeros(ny)
    rstar = np.zeros((nx,ny))
    for jj in range(ny):
        for ii in range(nx):
            if x[ii] <= xth: # Condition only required for the project 6:
                # LHS initial conditions for the Sod tube and vortex addition:
                rho[jj][ii] = rhoL
                u[jj][ii] = uL
                p[jj][ii] = pL 
                # These geometrical parameters shall remain common to all:
                xstar[ii] = (x[ii]-xc)/Rc
                ystar[jj] = (y[jj]-yc)/Rc
                rstar[jj][ii] = (((x[ii] - xc)**2 + (y[jj] - yc)**2)**0.5)/Rc
                # Project 5:
#                u[jj][ii] = uinf - (Mac*ac*ystar[jj]*math.exp((1-(rstar[jj][ii]**2))/2))
#                v[jj][ii] = vinf + (Mac*ac*xstar[ii]*math.exp((1-(rstar[jj][ii]**2))/2))
#                T[jj][ii] = T0*(1 - (temp2*(rstar[jj][ii]**2)*(math.exp(1-(rstar[jj][ii]**2)))))
#                p[jj][ii] = p0*((((T[jj][ii])/T0))**(gam/(gam - 1)))
#                rho[jj][ii] = p[jj][ii]/(Rgas*T[jj][ii])
            else: # Condition only required for the project 6:
                # RHS initial conditions for the Sod tube problem:
                rho[jj][ii] = rhoR
                u[jj][ii] = uR
                p[jj][ii] = pR
                # These geometrical parameters shall remain common to all:
                xstar[ii] = (x[ii]-xc)/Rc
                ystar[jj] = (y[jj]-yc)/Rc
                rstar[jj][ii] = (((x[ii] - xc)**2 + (y[jj] - yc)**2)**0.5)/Rc
                # RHS initial conditions for vortex addition:
                u[jj][ii] = uinf - (Mac*ac*ystar[jj]*math.exp((1-(rstar[jj][ii]**2))/2))
                v[jj][ii] = vinf + (Mac*ac*xstar[ii]*math.exp((1-(rstar[jj][ii]**2))/2))
                T[jj][ii] = Tinf*(1 - (temp2*(rstar[jj][ii]**2)*(math.exp(1-(rstar[jj][ii]**2)))))
                p[jj][ii] = pinf*((((T[jj][ii])/Tinf))**(gam/(gam - 1)))
                rho[jj][ii] = p[jj][ii]/(Rgas*T[jj][ii])
    return(rho,u,v,p,xstar,ystar,rstar)
rhoi,ui,vi,pi,xstar,ystar,rstar = initial(rhoL[case],uL[case],pL[case],rhoR[case],uR[case],pR[case],nx,ny,xth,x,y,xc,yc,Rc,gam,Rgas,Mac,ac,uinf,vinf,pinf,Tinf,temp1,temp2,T0,p0)

# Reconstructing the Q,F,G and a matrices:
def cons(rhoi,ui,vi,pi,gam):
    p = np.copy(pi)
    u = np.copy(ui)
    v = np.copy(vi)
    rho = np.copy(rhoi)
    rho_et = p/(gam-1) + 0.5*rho*(u**2 + v**2)
    Q = [rho,rho*u,rho*v,rho_et]
    F = [rho*u,rho*u*u + p ,rho*u*v, (rho_et + p)*u]
    G = [rho*v,rho*u*v ,rho*v*v + p, (rho_et + p)*v]
    a = np.sqrt(np.abs((gam*p)/rho))
    return(Q,F,G,a)
Q,F,G,a = cons(rhoi,ui,vi,pi,gam)

# Reconstructing the Q,F,G and a matrices:(Alternative method)
#def cons(rhoi,ui,vi,pi,gam,nx,ny):
#    p = np.copy(pi)
#    u = np.copy(ui)
#    v = np.copy(vi)
#    rho = np.copy(abs(rhoi))
#    rho_et = p/(gam-1) + 0.5*rho*(u**2 + v**2)
#    a = np.sqrt(abs((gam*p)/rho))
#    Q = np.zeros((4,nx,ny))
#    F = np.zeros((4,nx,ny))
#    G = np.zeros((4,nx,ny))
#    for kk in range(4):
#        if kk == 0: # Mass
#           for jj in range(ny):
#               for ii in range(nx):
#                   Q[kk][jj][ii] = rho[jj][ii]
#                   F[kk][jj][ii] = rho[jj][ii]*u[jj][ii]
#                   G[kk][jj][ii] = rho[jj][ii]*v[jj][ii]
#        elif kk == 1: # X momentum
#             for jj in range(ny):
#               for ii in range(nx):
#                     Q[kk][jj][ii] = rho[jj][ii]*u[jj][ii]
#                     F[kk][jj][ii] = rho[jj][ii]*u[jj][ii]**2 + p[jj][ii]
#                     G[kk][jj][ii] = rho[jj][ii]*u[jj][ii]*v[jj][ii]
#        elif kk == 2: # Y momentum
#             for jj in range(ny):
#               for ii in range(nx):
#                     Q[kk][jj][ii] = rho[jj][ii]*v[jj][ii]
#                     F[kk][jj][ii] = rho[jj][ii]*u[jj][ii]*v[jj][ii]
#                     G[kk][jj][ii] = rho[jj][ii]*v[jj][ii]**2 + p[jj][ii]
#        elif kk == 3: # Energy
#             for jj in range(ny):
#               for ii in range(nx):
#                     Q[kk][jj][ii] = rho_et[jj][ii]
#                     F[kk][jj][ii] = (rho_et[jj][ii] + p[jj][ii])*u[jj][ii]
#                     G[kk][jj][ii] = (rho_et[jj][ii] + p[jj][ii])*v[jj][ii]
#    return(Q,F,G,a)
#Q,F,G,a = cons(rhoi,ui,vi,pi,gam,nx,ny)

# Decoding the conservative variables:
def decons(Q,nx,ny,gam):
    e = np.zeros((nx,ny))
    rho_e = np.zeros((nx,ny))
    rho_et = np.zeros((nx,ny))
    p = np.zeros((nx,ny))
    u = np.zeros((nx,ny))
    v = np.zeros((nx,ny))
    rho = np.zeros((nx,ny))
    for jj in range(ny):
        for ii in range(nx):
            rho[jj][ii] = Q[0][jj][ii]
            u[jj][ii] = Q[1][jj][ii]/rho[jj][ii]
            v[jj][ii] = Q[2][jj][ii]/rho[jj][ii]
            rho_et[jj][ii] = Q[3][jj][ii]
            rho_e[jj][ii] = rho_et[jj][ii] - rho[jj][ii]*(u[jj][ii]**2 + v[jj][ii]**2)*0.5
            p[jj][ii] = (gam - 1)*rho_e[jj][ii]
            e[jj][ii] = rho_e[jj][ii]/rho[jj][ii]
    return(rho,u,v,p,e)
rho,u,v,p,e = decons(Q,nx,ny,gam)

# Decoding the conservative variables:(Alternative method)
#def decons(Q,nx,ny,gam):
#    e = np.zeros((nx,ny))
#    rho_et = np.zeros((nx,ny))
#    p = np.zeros((nx,ny))
#    u = np.zeros((nx,ny))
#    rho = np.zeros((nx,ny))
#    rho = (Q[0])
#    u = Q[1]/rho
#    v = Q[2]/rho
#    rho_et = Q[3]
#    e = rho_et - rho*(u*u + v*v)*0.5
#    p = (gam - 1)*e
#    return(rho,u,v,p,e)
#rho,u,v,p,e = decons(Q,nx,ny,gam)

# Function for dynamic time stepping:
#def dtmin(rho,p,u,v,gam,cmax,dx,dy):
#    a = np.sqrt(np.abs((gam*p)/rho))
#    velmag = np.sqrt(u**2 + v**2)
#    dtx = (cmax*dx)/(np.max(a + velmag))
#    dty = (cmax*dy)/(np.max(a + velmag))
#    dt = min(dtx,dty)
#    return(dt)
#dt = dtmin(rho,p,u,v,gam,cmax,dx,dy)

# Function for dynamic time stepping:(Alternative method)
#def dtmin(rho,p,u,v,gam,cmax,dx,dy):
#    dt = 1.e+30
#    a = np.sqrt(np.abs((gam*p)/rho))
#    dtx = (cmax*dx)/(np.max(a + np.abs(u)))
#    dty = (cmax*dy)/(np.max(a + np.abs(v)))
#    dt = min(dtx,dt)
#    dt = min(dty,dt)
#    return(dt)
#dt = dtmin(rho,p,u,v,gam,cmax,dx,dy)

# The Maccormack method:
def mac(rho,u,v,p,time_end,nx,ny,gam,cmax,dx,dy,Q):
    # Initializating parameters:
    timemac = []
    rhot_mac = []
    ut_mac = []
    vt_mac = []
    pt_mac = []
    et_mac = []
    rhonew,unew,vnew,pnew,enew = decons(Q,nx,ny,gam)
    Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam) # Initial values
    Qpred,Fpred,Gpred,apred = cons(rhonew,unew,vnew,pnew,gam) # Predicted values
    Qcorr,Fcorr,Gcorr,acorr = cons(rhonew,unew,vnew,pnew,gam) # Corrected values
    # Initializing time:
    dtx = (cmax*dx)/(np.max(np.sqrt(u**2+v**2)+a))
    dty = (cmax*dy)/(np.max(np.sqrt(u**2+v**2)+a))
    dt = min(dtx,dty)
    # Lists for storing values in time:
    timemac.append(dt)
    rhot_mac.append(rhonew)
    ut_mac.append(unew)
    vt_mac.append(vnew)
    pt_mac.append(pnew)
    et_mac.append(enew)
    t0 = 0
    it = 0
    # Initializing time loop:
    for t in range(3000):
        if t0 + dt > time_end:
            break
        # This contains the most updated values of Q,F,G and a at time level 'n+1'
        # They are equal to the ones defined below at the end of the time loop
        Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam)
        # We can also use the following configuration for better accuracy:
        # Predictor: F: Forward, G: Forward ; Corrector: F: Backward, G: Backward
        # Predictor: F: Forward, G: Backward ; Corrector: F: Backward, G: Forward
        # Predictor: F: Backward, G: Forward ; Corrector: F: Forward, G: Backward
        # Predictor: F: Backward, G: Backward ; Corrector: F: Backward, G: Backward
        # Predictor: Forward step
        for kk in range(4):
            for jj in range(0,ny-1):
                for ii in range(0,nx-1): 
                    Qpred[kk][jj][ii] = Q[kk][jj][ii] - (dt/dx)*(F[kk][jj][ii+1] - F[kk][jj][ii]) - (dt/dy)*(G[kk][jj+1][ii] - G[kk][jj][ii])
            # Boundary conditions:    
            for jj in range(0,ny-1): # Neumann BC - Right boundary
                Qpred[kk][jj][nx-1] = Qpred[kk][jj][nx-2]
###            for jj in range(0,ny-1): # Periodic BC - Right boundary
###                Qpred[kk][jj][nx-1] = Qpred[kk][jj][0]
            for ii in range(0,nx): # Top boundary
                Qpred[kk][ny-1][ii] = Qpred[kk][0][ii]
        # Updating the intermediate values:
        rhostar,ustar,vstar,pstar,estar = decons(Qpred,nx,ny,gam)
        Qpred,Fpred,Gpred,apred = cons(rhostar,ustar,vstar,pstar,gam)
        # Corrector: Backward step
        for kk in range(4):
            for jj in range(1,ny):
                for ii in range(1,nx):
                    Qcorr[kk][jj][ii] = 0.5*(Q[kk][jj][ii] + Qpred[kk][jj][ii] - (dt/dx)*(Fpred[kk][jj][ii] - Fpred[kk][jj][ii-1]) - (dt/dy)*(Gpred[kk][jj][ii] - Gpred[kk][jj-1][ii]))
            # Boundary conditions:    
            for jj in range(1,ny): # Neumann BC - Left boundary
                Qcorr[kk][jj][0] = Qcorr[kk][jj][1]
###            for jj in range(1,ny): # Periodic BC - Left boundary
###                Qcorr[kk][jj][0] = Qcorr[kk][jj][nx-1]
            for ii in range(0,nx): # Bottom boundary
                Qcorr[kk][0][ii] = Qcorr[kk][ny-1][ii]
        # Updating the final corrected values:
        rhonew,unew,vnew,pnew,enew = decons(Qcorr,nx,ny,gam)
        # The values of Q,F,G and a defined here are equal to the ones defined above after the time loop
        Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam)
        # Recompute dt dynamically to maintain stability
        # Remember to use new values of u,v,p and rho
        dtx = (cmax*dx)/(np.max(np.sqrt(unew**2+vnew**2)+a))
        dty = (cmax*dy)/(np.max(np.sqrt(unew**2+vnew**2)+a))
        dt = min(dtx,dty)
        it = it + 1
        print('for iteration =', it, 'for time step =',dt,end = '\n')
        # Update time:
        t0 = t0 + dt
        timemac.append(t0)
        rhot_mac.append(rhonew)
        ut_mac.append(unew)
        vt_mac.append(vnew)
        pt_mac.append(pnew)
        et_mac.append(enew)
    # Storing the final values at time level 'n+1'
    rhomac,umac,vmac,pmac,emac = decons(Qcorr,nx,ny,gam)
    return(rhomac,umac,vmac,pmac,emac,timemac,rhot_mac,ut_mac,vt_mac,pt_mac,et_mac)
rhomac,umac,vmac,pmac,emac,timemac,rhot_mac,ut_mac,vt_mac,pt_mac,et_mac = mac(rho,u,v,p,time_end,nx,ny,gam,cmax,dx,dy,Q)
print('Maccormack done')

def Rusanov(rho,u,v,p,time_end,nx,ny,gam,cmax,dx,dy,Q):
    # Initializating parameters:
    timerus = []
    rhot_rus = []
    ut_rus = []
    vt_rus = []
    pt_rus = []
    et_rus = []
    rhonew,unew,vnew,pnew,enew = decons(Q,nx,ny,gam)
    Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam)
    Qnew,Fnew,Gnew,anew = cons(rhonew,unew,vnew,pnew,gam)
    dtx = (cmax*dx)/(np.max(np.sqrt(u**2+v**2)+a))
    dty = (cmax*dy)/(np.max(np.sqrt(u**2+v**2)+a))
    dt = min(dtx,dty)
    # Lists for storing values in time:
    timerus.append(dt)
    rhot_rus.append(rhonew)
    ut_rus.append(unew)
    vt_rus.append(vnew)
    pt_rus.append(pnew)
    et_rus.append(enew)
    t0 = 0
    it = 0
    shx = 0
    shy = 0
    # Initializing time loop:
    for t in range(3000):
        if t0 + dt > time_end:
            break
        # This contains the most updated values of Q,F,G and a at time level 'n+1'
        Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam)
        # Calculating the F matrix:
        for kk in range(4):
            for jj in range(ny):
                for ii in range(0,nx-1):
                    # Remember to use new values of u for shx
                    shx = max(abs(unew[jj][ii] + a[jj][ii]),abs(unew[jj][ii+1] + a[jj][ii+1]))
                    Fnew[kk][jj][ii] = 0.5*((F[kk][jj][ii] + F[kk][jj][ii+1]) - shx*(Q[kk][jj][ii+1] - Q[kk][jj][ii]))
            # Boundary conditions:
            for jj in range(ny): # Neumann BC : Right BC
                Fnew[kk][jj][nx-1] = Fnew[kk][jj][nx-2] 
#            for jj in range(ny): # Periodic BC : Right BC
#                Fnew[kk][jj][nx-1] = Fnew[kk][jj][0]
        # Calculating the G matrix:
        for kk in range(4):
            for ii in range(nx):
                for jj in range(0,ny-1):
                    # Remember to use new values of v for shy
                    shy = max(abs(vnew[jj][ii] + a[jj][ii]),abs(vnew[jj+1][ii] + a[jj+1][ii]))
                    Gnew[kk][jj][ii] = 0.5*((G[kk][jj][ii] + G[kk][jj+1][ii]) - shy*(Q[kk][jj+1][ii] - Q[kk][jj][ii]))
            # Boundary conditions:
            for ii in range(nx): # Periodic in y: Top BC
                Gnew[kk][ny-1][ii] = Gnew[kk][0][ii] 
        # Calculating the Q matrix:
        for kk in range(4):   
            for jj in range(1,ny):
                for ii in range(1,nx):
                    Qnew[kk][jj][ii] = Q[kk][jj][ii] - (dt/dx)*(Fnew[kk][jj][ii] - Fnew[kk][jj][ii-1]) - (dt/dy)*(Gnew[kk][jj][ii] - Gnew[kk][jj-1][ii])
            # Boundary conditions:
            for jj in range(1,ny): # Neumann BC : Left BC 
                Qnew[kk][jj][0] = Qnew[kk][jj][1]
##            for jj in range(1,ny): # Periodic BC : Left BC 
##                Qnew[kk][jj][0] = Qnew[kk][jj][nx-1]
            for ii in range(0,nx): # Periodic in y: Bottom BC
                Qnew[kk][0][ii] = Qnew[kk][ny-1][ii]
        # For the Sod tube:
#        for kk in range(4):   
#            for jj in range(1,ny-1):
#                for ii in range(1,nx-1):
#                    Qnew[kk][jj][ii] = Q[kk][jj][ii] - (dt/dx)*(Fnew[kk][jj][ii] - Fnew[kk][jj][ii-1]) - (dt/dy)*(Gnew[kk][jj][ii] - Gnew[kk][jj-1][ii])
        # Updating the final corrected values:
        rhonew,unew,vnew,pnew,enew = decons(Qnew,nx,ny,gam)
        Q,F,G,a = cons(rhonew,unew,vnew,pnew,gam)
        # Recompute dt dynamically to maintain stability
        # Remember to use new values of u,v,p and rho
        dtx = (cmax*dx)/(np.max(np.sqrt(unew**2+vnew**2)+a))
        dty = (cmax*dy)/(np.max(np.sqrt(unew**2+vnew**2)+a))
        dt = min(dtx,dty)
        it = it + 1
        print('for iteration =', it, 'for time step =',dt,end = '\n')
        # Update time:
        t0 = t0 + dt
        timerus.append(t0)
        rhot_rus.append(rhonew)
        ut_rus.append(unew)
        vt_rus.append(vnew)
        pt_rus.append(pnew)
        et_rus.append(enew)
    # Storing the final values at time level 'n+1'
    rhorus,urus,vrus,prus,erus = decons(Qnew,nx,ny,gam)
    print('Rusanov done')
    return(rhorus,urus,vrus,prus,erus,timerus,rhot_rus,ut_rus,vt_rus,pt_rus,et_rus)
rhorus,urus,vrus,prus,erus,timerus,rhot_rus,ut_rus,vt_rus,pt_rus,et_rus = Rusanov(rho,u,v,p,time_end,nx,ny,gam,cmax,dx,dy,Q)

# Calculation of vorticity:
def vort(u,v,nx,ny,dx,dy):
    vortz = np.zeros((nx,ny))
    dvdx, dudy = np.zeros((nx,ny)), np.zeros((nx,ny))
    for jj in range(0,ny):
        for ii in range(1,nx-1):
            dvdx[jj][ii] = (v[jj][ii+1] - v[jj][ii-1])/(2*dx)
        dvdx[jj][0] = (v[jj][1] -v[jj][nx-2])/(2*dx)
        dvdx[jj][nx-1] = dvdx[jj][nx-2] # Neumann boundaries
#        dvdx[jj][nx-1] = dvdx[jj][0] # Periodic boundaries
    for ii in range(0,nx):
        for jj in range(1,ny-1):
            dudy[jj][ii] = (u[jj+1][ii] - u[jj-1][ii])/(2*dy)
        dudy[0][ii] = (u[1][ii] - u[ny-2][ii])/(2*dy)
        dudy[ny-1][ii] = dudy[0][ii] # Periodic boundaries
    vortz = dvdx - dudy
    return(vortz)
vortnum = vort(u,v,nx,ny,dx,dy)

# Calculation of analytical vorticity:
def vort_exact(nx,ny,Mac,Rc,rstar,ac):
    vortexact = np.zeros((nx,ny))
    for jj in range(ny):
        for ii in range(nx):
            vortexact[jj][ii] = ((Mac*ac)/Rc)*math.exp((1-(rstar[jj][ii]**2))*0.5)*(2-rstar[jj][ii]**2)
    return(vortexact)
vortexact = vort_exact(nx,ny,Mac,Rc,rstar,ac)

# Calculation of the dilatation term:
def dilatation(u,v,dx,dy,nx,ny):
    dil = np.zeros((nx,ny))
    dudx, dvdy = np.zeros((nx,ny)), np.zeros((nx,ny))
    for jj in range(0,ny):
        for ii in range(1,nx-1):
            dudx[jj][ii] = (u[jj][ii+1] - u[jj][ii-1])/(2*dx)
        # Boundary conditions:
        dudx[jj][0] = (u[jj][1] - u[jj][nx-2])/(2*dx)
        dudx[jj][nx-1] = dudx[jj][nx-2] # Neumann boundaries
#        dudx[jj][nx-1] = dudx[jj][0] # Periodic boundaries
    for ii in range(0,nx):
        for jj in range(1,ny-1):
            dvdy[jj][ii] = (v[jj+1][ii] - v[jj-1][ii])/(2*dy)
        # Boundary conditions:
        dvdy[0][ii] = (v[1][ii] - v[ny-2][ii])/(2*dy)
        dvdy[ny-1][ii] = dvdy[0][ii] # Periodic boundaries
    dil = dudx + dvdy
    return(dil)
dil = dilatation(u,v,dx,dy,nx,ny)

# Calculation of Baroclinic torque:
def btorque(nx,ny,dx,dy,rho,p):
    barotorque = np.zeros((nx,ny))
    drhodx = np.zeros((nx,ny))
    drhody = np.zeros((nx,ny))
    dpdx = np.zeros((nx,ny))
    dpdy = np.zeros((nx,ny))
    for jj in range(0,ny):
        for ii in range(1,nx-1):
            drhodx[jj][ii] = (rho[jj][ii+1] - rho[jj][ii-1])/(2*dx)
            dpdx[jj][ii] = (p[jj][ii+1] - p[jj][ii-1])/(2*dx)
        # Boundary conditions:
        drhodx[jj][0] = (rho[jj][1] - rho[jj][nx-2])/(2*dx)
        drhodx[jj][nx-1] = drhodx[jj][nx-2] # Neumann boundaries
#        drhodx[jj][nx-1] = drhodx[jj][0] # Periodic boundaries
        dpdx[jj][0] = (p[jj][1] - p[jj][nx-2])/(2*dx)
        dpdx[jj][nx-1] = dpdx[jj][nx-2] # Neumann boundaries
    for ii in range(0,nx):
        for jj in range(1,ny-1):
            drhody[jj][ii] = (rho[jj+1][ii] - rho[jj-1][ii])/(2*dy)
            dpdy[jj][ii] = (p[jj+1][ii] - p[jj-1][ii])/(2*dy)
        # Boundary conditions:
        drhody[0][ii] = (rho[1][ii] - rho[ny-2][ii])/(2*dy)
        drhody[ny-1][ii] = drhody[0][ii] # Periodic boundaries
        dpdy[0][ii] = (p[1][ii] - p[ny-2][ii])/(2*dy)
        dpdy[ny-1][ii] = dpdy[0][ii] # Periodic boundaries
    barotorque = (drhodx*dpdy - drhody*dpdx)
    return(barotorque)
barotorque = btorque(nx,ny,dx,dy,rho,p)

# Calculation of enstrophy:
def enstrophy(nx,ny,dx,dy,time,dA,ut,vt):
    ens = 0
    eps = []
    ud = []
    vd = []
    for tt in range(len(time)):
        ud = ut[tt]
        vd = vt[tt]
        vortnum = vort(ud,vd,nx,ny,dx,dy)
        ens = np.sum(((vortnum)**2)*dA)
        eps.append(ens)
    return(eps)
#eps = enstrophy(nx,ny,dx,dy,time,dA,ut,vt)

# Calculation of source terms:
def source_dil(nx,ny,dx,dy,time,ut,vt,dA):
    Sdiltemp = 0
    Sdil = []
    ud = []
    vd = []
    for tt in range(len(time)):
        ud = ut[tt]
        vd = vt[tt]
        vortnum = vort(ud,vd,nx,ny,dx,dy)
        dil = dilatation(ud,vd,dx,dy,nx,ny)
        Sdiltemp = np.sum(-2*(vortnum**2)*dil*dA)
        Sdil.append(Sdiltemp)
    return(Sdil)
#Sdil = source_dil(nx,ny,dx,dy,time,ut,vt)

def source_tor(nx,ny,dx,dy,time,ut,vt,rhot,pt,dA):
    Stortemp = 0
    Stor = []
    ud = []
    vd = []
    rhod = []
    pd = []
    for tt in range(len(time)):
        ud = ut[tt]
        vd = vt[tt]
        rhod = rhot[tt]
        pd = pt[tt]
        vortnum = vort(ud,vd,nx,ny,dx,dy)
        barotorque = btorque(nx,ny,dx,dy,rhod,pd)
        Stortemp = np.sum(((2*vortnum)/(rho**2))*barotorque*dA)
        Stor.append(Stortemp)
    return(Stor)
#Stor = source_tor(nx,ny,dx,dy,time,ut,vt,rhot,pt)
    
# Calculation of shadowgraphs: 
#def shadow(nx,ny,dx,dy,rho):
#    shade1 = np.zeros((nx,ny))
#    shade2 = np.zeros((nx,ny))
#    dx2 = dx*dx
#    dy2 = dy*dy
#    for jj in range(0,ny):
#        for ii in range(1,nx-1):
#            shade1[jj][ii] = (rho[jj][ii+1] - 2*rho[jj][ii] + rho[jj][ii-1])/dx2
#        shade1[jj][nx-1] = (rho[jj][1] - 2*rho[jj][nx-1] + rho[jj][nx-2])/dx2
##        shade1[jj][0] = shade1[jj][1] # Neumann boundaries
#        shade1[jj][0] = shade1[jj][nx-1] # Periodic boundaries
#    for ii in range(0,nx):
#        for jj in range(1,ny-1):
#            shade2[jj][ii] = (rho[jj+1][ii] - 2*rho[jj][ii] + rho[jj-1][ii])/dy2
#        shade2[ny-1][ii] = (rho[1][ii] - 2*rho[ny-1][ii] + rho[ny-2][ii])/dy2
#        shade2[0][ii] = shade2[ny-1][ii] # Periodic boundaries
#    shade = shade1 + shade2
#    return(shade)
#shade = shadow(nx,ny,rho)  

# Calculation of circulation and error:
#def error(nx,ny,eps,vortnum):
#    circ = 0
#    err = 0
#    for jj in range(ny):
#        for ii in range(nx):
#            circ = circ + vortnum[jj][ii]*dx*dy
#            circ = circ + vortnum[jj][ii]
#    err = circ/math.sqrt(eps)
#    return(err)
#err = error(nx,ny,eps,vortnum)

# Post processing:
# Project 6:
# Baseline cases: (Sod tube)
# 1) Maccormack for all the grid resolutions and cmax = 1.0 and xth = L/4,time_end = L/4/SR
#R.plot_compare(x,xth,L,rhomac[yint],pmac[yint],umac[yint],emac[yint],gam,time_end,case,"rusanov"+str(case),nx)
# 2) Rusanov for all the grid resolutions and cmax = 0.7 and xth = L/4,time_end = L/4/SR
#R.plot_compare(x,xth,L,rhorus[yint],prus[yint],urus[yint],erus[yint],gam,time_end,case,"rusanov"+str(case),nx)

# Vortex addition:
# This is a generic plotting function for plotting contours of different cases:
def ccompare(x,y,var,varexact, casename, titlename):
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x, y, var)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    cp.collections[0].set_label(casename)
    cp = ax.contour(x, y, varexact,linestyles = 'dashed')
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    cp.collections[0].set_label('exact')
    plt.colorbar(cp)
    ax.set_title(titlename)
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()

# Definiing variables for the grid:
if Nx == 25 and Ny == 25:
    x25 = x
    y25 = y
    yint25 = int(ny/2)
    vortmac25 = vort(umac,vmac,nx,ny,dx,dy)
    vortrus25 = vort(urus,vrus,nx,ny,dx,dy)
    vortexact25 = vort_exact(nx,ny,Mac,Rc,rstar,ac)
    timemac25 = timemac
    timerus25 = timerus
    epsmac25 = enstrophy(nx,ny,dx,dy,timemac,dA,ut_mac,vt_mac)
    epsrus25 = enstrophy(nx,ny,dx,dy,timerus,dA,ut_rus,vt_rus)
    Sdilmac25 = source_dil(nx,ny,dx,dy,timemac,ut_mac,vt_mac,dA)
    Sdilrus25 = source_dil(nx,ny,dx,dy,timerus,ut_rus,vt_rus,dA)
    dilmac25 = dilatation(umac,vmac,dx,dy,nx,ny)
    dilrus25 = dilatation(urus,vrus,dx,dy,nx,ny)
    Stormac25 = source_tor(nx,ny,dx,dy,timemac,ut_mac,vt_mac,rhot_mac,pt_mac,dA)
    Storrus25 = source_tor(nx,ny,dx,dy,timerus,ut_rus,vt_rus,rhot_rus,pt_rus,dA)
    btormac25 = btorque(nx,ny,dx,dy,rhomac,pmac)
    btorrus25 = btorque(nx,ny,dx,dy,rhorus,prus)
elif Nx == 50 and Ny == 50:
    x50 = x
    y50 = y
    yint50 = int(ny/2)
    vortmac50 = vort(umac,vmac,nx,ny,dx,dy)
    vortrus50 = vort(urus,vrus,nx,ny,dx,dy)
    vortexact50 = vort_exact(nx,ny,Mac,Rc,rstar,ac)
    timemac50 = timemac
    timerus50 = timerus
    epsmac50 = enstrophy(nx,ny,dx,dy,timemac,dA,ut_mac,vt_mac)
    epsrus50 = enstrophy(nx,ny,dx,dy,timerus,dA,ut_rus,vt_rus)
    Sdilmac50 = source_dil(nx,ny,dx,dy,timemac,ut_mac,vt_mac,dA)
    Sdilrus50 = source_dil(nx,ny,dx,dy,timerus,ut_rus,vt_rus,dA)
    dilmac50 = dilatation(umac,vmac,dx,dy,nx,ny)
    dilrus50 = dilatation(urus,vrus,dx,dy,nx,ny)
    Stormac50 = source_tor(nx,ny,dx,dy,timemac,ut_mac,vt_mac,rhot_mac,pt_mac,dA)
    Storrus50 = source_tor(nx,ny,dx,dy,timerus,ut_rus,vt_rus,rhot_rus,pt_rus,dA)
    btormac50 = btorque(nx,ny,dx,dy,rhomac,pmac)
    btorrus50 = btorque(nx,ny,dx,dy,rhorus,prus)
elif Nx == 100 and Ny == 100:
    x100 = x
    y100 = y
    yint100 = int(ny/2)
    vortmac100 = vort(umac,vmac,nx,ny,dx,dy)
    vortrus100 = vort(urus,vrus,nx,ny,dx,dy)
    vortexact100 = vort_exact(nx,ny,Mac,Rc,rstar,ac)
    timemac100 = timemac
    timerus100 = timerus
    epsmac100 = enstrophy(nx,ny,dx,dy,timemac,dA,ut_mac,vt_mac)
    epsrus100 = enstrophy(nx,ny,dx,dy,timerus,dA,ut_rus,vt_rus)
    Sdilmac100 = source_dil(nx,ny,dx,dy,timemac,ut_mac,vt_mac,dA)
    Sdilrus100 = source_dil(nx,ny,dx,dy,timerus,ut_rus,vt_rus,dA)
    dilmac100 = dilatation(umac,vmac,dx,dy,nx,ny)
    dilrus100 = dilatation(urus,vrus,dx,dy,nx,ny)
    Stormac100 = source_tor(nx,ny,dx,dy,timemac,ut_mac,vt_mac,rhot_mac,pt_mac,dA)
    Storrus100 = source_tor(nx,ny,dx,dy,timerus,ut_rus,vt_rus,rhot_rus,pt_rus,dA)
    btormac100 = btorque(nx,ny,dx,dy,rhomac,pmac)
    btorrus100 = btorque(nx,ny,dx,dy,rhorus,prus)
    
# Profile and Contour plots:
if Nx == 25 and Ny == 25:  
    # Maccormack:
    # Vorticity:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(x25,vortmac25[yint25],label='MacCormack Nx = 25',color='red')
    plt.plot(x50,vortmac50[yint50],label='MacCormack Nx = 50',color='blue')
    plt.plot(x100,vortmac100[yint100],label='MacCormack Nx = 100',color='black')
    plt.plot(x25,vortexact25[yint25],label='Exact Nx = 25',color='red',linestyle='--')
    plt.plot(x50,vortexact50[yint50],label='Exact Nx = 50',color='blue',linestyle='--')
    plt.plot(x100,vortexact100[yint100],label='Exact Nx = 100',color='black',linestyle='--')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('Vorticity',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Enstrophy:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timemac25,epsmac25,label='Maccormack Nx = 25',color='red')
    plt.plot(timemac50,epsmac50,label='Maccormack Nx = 50',color='blue')
    plt.plot(timemac100,epsmac100,label='Maccormack Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Enstrophy',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Dilatation:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timemac25,Sdilmac25,label='Maccormack Nx = 25',color='red')
    plt.plot(timemac50,Sdilmac50,label='Maccormack Nx = 50',color='blue')
    plt.plot(timemac100,Sdilmac100,label='Maccormack Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Dilatation',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Baro torque:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timemac25,Stormac25,label='Maccormack Nx = 25',color='red')
    plt.plot(timemac50,Stormac50,label='Maccormack Nx = 50',color='blue')
    plt.plot(timemac100,Stormac100,label='Maccormack Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Baroclinic torque',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Rusanov:
    # Vorticity:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(x25,vortrus25[yint25],label='Rusanov Nx = 25',color='red')
    plt.plot(x50,vortrus50[yint50],label='Rusanov Nx = 50',color='blue')
    plt.plot(x100,vortrus100[yint100],label='Rusanov Nx = 100',color='black')
    plt.plot(x25,vortexact25[yint25],label='Exact Nx = 25',color='red',linestyle='--')
    plt.plot(x50,vortexact50[yint50],label='Exact Nx = 50',color='blue',linestyle='--')
    plt.plot(x100,vortexact100[yint100],label='Exact Nx = 100',color='black',linestyle='--')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('Vorticity',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Enstrophy:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timerus25,epsrus25,label='Rusanov Nx = 25',color='red')
    plt.plot(timerus50,epsrus50,label='Rusanov Nx = 50',color='blue')
#    plt.plot(timerus100,epsrus100,label='Rusanov Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Enstrophy',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show() 

    # Dilatation:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timerus25,Sdilrus25,label='Rusanov Nx = 25',color='red')
    plt.plot(timerus50,Sdilrus50,label='Rusanov Nx = 50',color='blue')
#    plt.plot(timerus100,Sdilrus100,label='Rusanov Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Dilatation',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
    
    # Baro torque:
    f, ax = plt.subplots(1,1,figsize=(8,5))
    plt.plot(timerus25,Storrus25,label='Rusanov Nx = 25',color='red')
    plt.plot(timerus50,Storrus50,label='Rusanov Nx = 50',color='blue')
#    plt.plot(timerus100,Storrus100,label='Rusanov Nx = 100',color='black')
    ax.set_xlabel('$Time history$',size=20)
    ax.set_ylabel('Baroclinic torque',size=20)
    ax.grid()
    ax.legend(fontsize=16)
    ax.plot()
    plt.show()
        
if Nx == 25 and Ny == 25:    
# Contours for vorticity, dilatation and bar torque:
    # Macormack:
    ccompare(x25,y25,vortmac25,vortexact25,'Numerical','Maccormack Nx = 25')
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x25, y25, dilmac25)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 25')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x25, y25, btormac25)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 25')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    # Rusanov:
    ccompare(x25,y25,vortrus25,vortexact25,'Numerical','Rusanov Nx = 25')
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x25, y25, dilrus25)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 25')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x25, y25, btorrus25)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 25')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
elif Nx == 50 and Ny == 50:
    # Macormack:
    ccompare(x50,y50,vortmac50,vortexact50,'Numerical','Maccormack Nx = 50')
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x50, y50, dilmac50)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 50')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x50, y50, btormac50)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 50')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    # Rusanov:
    ccompare(x50,y50,vortrus50,vortexact50,'Numerical','Rusanov Nx = 50') 
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x50, y50, dilrus50)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 50')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x50, y50, btorrus50)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 50')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
elif Nx == 100 and Ny == 100:
    # Macormack:
    ccompare(x100,y100,vortmac100,vortexact100,'Numerical','Maccormack Nx = 100')
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x100, y100, dilmac100)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 100')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x100, y100, btormac100)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Maccormack Nx = 100')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    # Rusanov:
    ccompare(x100,y100,vortrus100,vortexact100,'Numerical','Rusanov Nx = 100') 
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x100, y100, dilrus100)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 100')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    
    fig, ax = plt.subplots(1,1,figsize=(8,8))
    cp = ax.contour(x100, y100, btorrus100)
    ax.clabel(cp,inline=True, fontsize=10, colors = 'black')
    plt.colorbar(cp)
    cp.collections[0].set_label('Rusanov Nx = 100')
    ax.set_xlabel('$x$',size=20)
    ax.set_ylabel('$y$',size=20)
    ax.legend(fontsize=12)
    ax.grid()
    plt.show()
    
    
    
    
    
    
    
    