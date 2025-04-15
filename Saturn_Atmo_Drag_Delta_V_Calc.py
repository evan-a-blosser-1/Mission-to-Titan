################################################
# Author: Evan A. Blosser                      #
# Email: evan.a.blosser-1@ou.edu               #
#                                              #
# Program: Saturn_Atmo_Drag_Delta_V_Calc       #
#           -Calculating amount of time needed #
#            to Achieve a User defined Delta-V #
#             to slow down on Saturn arrival   # 
################################################
#%% Packages & Constants
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
import math  as math                #
# System Time                       #
import time                         #
# System Error/Exit                 #
from sys import exit                #
# ODE solver                        #
from scipy.integrate import odeint  #
# Plotting & Animation              #
import matplotlib.pyplot as plt     #
#####################################
#########################
# Constants Declaration #
##############################################################
# Gravitational Parameter for Saturn and Titan in km^3/s^2   #
# [If working in Canonical Units replace with 1 (DU^3/TU^2)] #
mu_saturn = 3.986004418e5                                    #
mu_titan  = 8978.13710                                       #
# Equitorial Radius of the Earth and Titan in km             #
# [If working in Canonical Units replace with 1 (DU)]        #
r_saturn = 60270                                             #
r_titan  = 2575.5                                            #
##############################################################
############
### Main ###
############
#%% User Menu
welcome =f"""
{'-'*42}
            Welcome to Saturn Delta-V
             Atmo. Drag Calculator
{'-'*42}
| Author: Evan A. Blosser
{'-'*42}
| Given the initial position and velocity vecotrs of the 
|   space craft upon entry near saturn/ Orbital Elements 
|    around: can calculate the Delta-V breaking recieved 
|     from Saturn  
{'-'*42}
|  Note: Please enter simulation length 
|        as in reasonable units.
|        
{'-'*42}
"""                              
print(welcome)         

###############
# User Inputs # Hardwired
####################################################################
r_p_a    = 0
e        = 0
i_in     = 0
nu_in    = 0
Omega_in = 0
w_in     = 0                                    
Graph_input = 'h'             
t_inp      = 24
####################################################################
#
##################
# User Selection #
##############################################
# Setting for Hours                          #
if Graph_input == 'h' or Graph_input == 'H': #
# In seconds (60*60) from hours              #
    Graph_time = 3600                        #                                    
else:                                          ###############
    exit('ERROR: input for time unit is incorrect, exit...') #
##############################################################

#%% Initial Calculations of State Vectors From User Input
#######################
# Initial Calculation #
# Cpu Clock           #################
Initial_Calc_Start_Time = time.time() #
#######################################
#
##########################################
# Time step, t_inp times the chosen unit #
n     = t_inp*Graph_time                 #
# Graph points, set to the chosen unit   # 
points = Graph_time                      #
##########################################
#
######################## 
# Angle Inputs deg2Rad #
###############################
i_0  = np.deg2rad(i_in)        #
nu_0 = np.deg2rad(nu_in)        #
Omega_0 = np.deg2rad(Omega_in)   #
w_0  = np.deg2rad(w_in)         #
################################
#
##################
# Rotaton matrix #
###############################################################################
R_11 =  np.cos(Omega_0)*np.cos(w_0) - np.sin(Omega_0)*np.sin(w_0)*np.cos(i_0) #
R_12 = -np.cos(Omega_0)*np.sin(w_0) - np.sin(Omega_0)*np.cos(w_0)*np.cos(i_0) #
R_13 =  np.sin(Omega_0)*np.sin(i_0)                                           #
R_21 =  np.sin(Omega_0)*np.cos(w_0) + np.cos(Omega_0)*np.sin(w_0)*np.cos(i_0) #
R_22 = -np.sin(Omega_0)*np.sin(w_0) + np.cos(Omega_0)*np.cos(w_0)*np.cos(i_0) #
R_23 = -np.cos(Omega_0)*np.sin(i_0)                                          #
R_31 =  np.sin(w_0)*np.sin(i_0)                                             #
R_32 =  np.cos(w_0)*np.sin(i_0)                                           #
R_33 =  np.cos(i_0)                                                     #
# Assign Matrix ########################################################
R_matrix = np.array([[R_11, R_12, R_13], [R_21, R_22, R_23],   # 
                     [R_31, R_32, R_33]])                      #
################################################################
#
##########
# Saturn #
############################################
# Set altitude to radius relative to Earth #
r_p = r_p_a + r_saturn                     #
# Calculate the Semi-Latus Rectum          #
p = r_p*(1 + e)                            #
# Position magnitue for elements           #
r = p / (1 + e*np.cos(nu_0))               #
############################################
# Calculate the perifocal position vector #
r_pf_p = r*np.cos(nu_0)                   #
r_pf_q = r*np.sin(nu_0)                   #
# Set perifocal position vector in array   #
r_pf = np.array([r_pf_p, r_pf_q, 0])        #
# Calculate perifocal velocity vector         #
v_pf_p = np.sqrt(mu_saturn/p)*(-np.sin(nu_0))  #
v_pf_q = np.sqrt(mu_saturn/p)*(e+np.cos(nu_0)) #
# Set perifocal velocity vector arry          #
v_pf = np.array([v_pf_p, v_pf_q, 0])         #
############################################
######################################
# Geocentric Position Vector         #
r_geo = np.matmul(R_matrix,r_pf)     #
# Geocentric Velocity Vector         #
v_geo = np.matmul(R_matrix,v_pf)     #
######################################
#
##########
# Titan  #
############################################
# Set altitude to radius relative to Earth #
r_p_t = r_p_a + r_titan                    #
# Calculate the Semi-Latus Rectum          #
p_t = r_p_t*(1 + e)                        #
# Position magnitue for elements           #
r_t = p_t / (1 + e*np.cos(nu_0))           #
###################################################
# Calculate the perifocal position vector         #
r_pf_p_t = r_t*np.cos(nu_0)                       #
r_pf_q_t = r_t*np.sin(nu_0)                       #
# Set perifocal position vector in array          #
r_pf_t = np.array([r_pf_p_t, r_pf_q_t, 0])        #
# Calculate perifocal velocity vector             #
v_pf_p_t = np.sqrt(mu_titan/p_t)*(-np.sin(nu_0))  #
v_pf_q_t = np.sqrt(mu_titan/p_t)*(e+np.cos(nu_0)) #
# Set perifocal velocity vector arry              #
v_pf_t = np.array([v_pf_p_t, v_pf_q_t, 0])        #
###################################################
######################################
# Geocentric Position Vector         #
r_geo_t = np.matmul(R_matrix,r_pf_t) #
# Geocentric Velocity Vector         #
v_geo_t = np.matmul(R_matrix,v_pf_t) #
######################################
#############
# Cpu Clock #########################
Initial_Calc_End_Time = time.time() #
#######################################
# Initial Calculations Execution Time ####################################
Initial_Calc_Execution = Initial_Calc_End_Time - Initial_Calc_Start_Time #
##########################################################################
#%% Orbit ODE
###################
# ODE Calculation #
# Cpu Clock       #################
ODE_Calc_Start_Time = time.time() #
####################################
#
# Settings #
##############################################################
# Time and Points                                            #
t_span  = np.linspace(0,n,points)                            #
##############################################################
#%% Saturn
###############
## ODE Setup ##
###############
v_rel_sat = [((2*np.pi)/10.656*3600),((2*np.pi)/10.656*3600),0]  
# kg/km^3
rho_saturn = 0.19e-9
C_drag = 1.5
Craft_Area_mass = 0.02
#####################
# Define ODE System #
#####################
# Define ODE System #
def orbit(a,t):
      #
      # Main Fix!
      ###########################################
      r_geo_vec = [0]*3                         # Calcualte new position 
      v_geo_vec = [0]*3                         #
      for j in range(0,3):                      #  norm/magnitude with
          r_geo_vec[j] = a[j]                   #  the current position!
          r_mag = np.linalg.norm(r_geo_vec)     #
          # Delta-V Calculations by change in   #
          #  velocity at each new position      #
          v_geo_vec[j] = a[j + 3]               #
          v_rel_vec = v_geo_vec[j] - v_rel_sat  #
          V_rel = np.linalg.norm(v_rel_vec)     #
      ###########################################
      rx  = a[0] # Define Initial Position
      ry  = a[1] #   x, y, z component
      rz  = a[2] #
      vx  = a[3] # Define Initial Velocity
      vy  = a[4] #   x, y, z component
      vz  = a[5] #
      ##################################################
      # Drag Force Equations       #####################
      Dragx = - 1/2*rho_saturn*V_rel*(C_drag*Craft_Area_mass)*v_rel_vec[0]
      Dragy = - 1/2*rho_saturn*V_rel*(C_drag*Craft_Area_mass)*v_rel_vec[1]
      Dragz = - 1/2*rho_saturn*V_rel*(C_drag*Craft_Area_mass)*v_rel_vec[2]
      ##################################################
      # Derivatives/ 6 Equations of Motion #############
      drxdt = vx 
      drydt = vy 
      drzdt = vz 
      dvxdt = -mu_saturn*(rx/r_mag**3) + Dragx
      dvydt = -mu_saturn*(ry/r_mag**3) + Dragy
      dvzdt = -mu_saturn*(rz/r_mag**3) + Dragz
      dadt  = [drxdt,drydt,drzdt,dvxdt,dvydt,dvzdt]
      return dadt
# End ODE Definition ###################################                                                    
######################

# Initial Conditions                                         #
a0 = [r_geo[0],r_geo[1],r_geo[2],v_geo[0],v_geo[1],v_geo[2]] #
##############################################################
################
# ODE Solution ##############
a = odeint(orbit,a0,t_span) #
#############################
#%%Titan
###############
## ODE Setup ##
###############
spin_titan = 15.945421*60*60*24
v_rel_titan = [((2*np.pi)/spin_titan),((2*np.pi)/spin_titan),0]  
# kg/km^3
rho_titan = 0.19e-9
C_drag = 1.5
Craft_Area_mass = 0.02
#####################
# Define ODE System #
#####################
# Define ODE System #
def orbit(a_t,t):
      ################################################
      r_geo_vec_t = [0]*3                            # Calcualte new position 
      v_geo_vec_t = [0]*3                            #
      for j in range(0,3):                           #  norm/magnitude with
          r_geo_vec_t[j] = a_t[j]                    #  the current position!
          r_mag = np.linalg.norm(r_geo_vec_t)        #
          # Delta-V Calculations by change in Vel    #
          #  at each new position                    #
          v_geo_vec_t[j] = a_t[j + 3]                #
          v_rel_vec_t = v_geo_vec_t[j] - v_rel_titan #
          V_rel_t = np.linalg.norm(v_rel_vec_t)      #
      ################################################
      rx  = a_t[0] # Define Initial Position
      ry  = a_t[1] #   x, y, z component
      rz  = a_t[2] #
      vx  = a_t[3] # Define Initial Velocity
      vy  = a_t[4] #   x, y, z component
      vz  = a_t[5] #
      ##################################################
      # Drag Force Equations       #####################
      Dragx = - 1/2*rho_titan*V_rel_t*(C_drag*Craft_Area_mass)*v_rel_vec_t[0]
      Dragy = - 1/2*rho_titan*V_rel_t*(C_drag*Craft_Area_mass)*v_rel_vec_t[1]
      Dragz = - 1/2*rho_titan*V_rel_t*(C_drag*Craft_Area_mass)*v_rel_vec_t[2]
      ##################################################
      # Derivatives/ 6 Equations of Motion #############
      drxdt = vx 
      drydt = vy 
      drzdt = vz 
      dvxdt = -mu_titan*(rx/r_mag**3) + Dragx
      dvydt = -mu_titan*(ry/r_mag**3) + Dragy
      dvzdt = -mu_titan*(rz/r_mag**3) + Dragz
      da_tdt  = [drxdt,drydt,drzdt,dvxdt,dvydt,dvzdt]
      return da_tdt
# End ODE Definition ###################################                                                    
######################

# Initial Conditions                                         #
a_t0 = [r_geo_t[0],r_geo_t[1],r_geo_t[2],v_geo_t[0],v_geo_t[1],v_geo_t[2]] #
##############################################################
################
# ODE Solution ##############
a_t = odeint(orbit,a_t0,t_span) #
#############################
#%% Using ODE Solution to Solve for Orbital Elements
###################################
# Solve for elements at each time #
###################################
# Setting empty matrices ##############################################
rcl, vcl = [],[]                                                      #
r_m,v_ml = [],[]                                                      #
rc_tl, vc_tl = [],[]                                                      #
r_m_t,v_m_tl = [],[]                                                      #
delta_v_total = 0
delta_v_l     = []
#######################################################################
#
#####################
# Start Calculating #
#####################
for i in range(0,points):
    rc = np.array([a[i,0],a[i,1],a[i,2]])
    rcl.append(rc)
    vc = np.array([a[i,3],a[i,4],a[i,5]])
    vcl.append(vc)
    r_m = np.linalg.norm(rc)
    v_m = np.linalg.norm(vc)
    v_ml.append(v_m)
    rc_t = np.array([a_t[i,0],a_t[i,1],a_t[i,2]])
    rc_tl.append(rc_t)
    vc_t = np.array([a_t[i,3],a_t[i,4],a_t[i,5]])
    vc_tl.append(vc_t)
    r_m_t = np.linalg.norm(rc_t)
    v_m_t = np.linalg.norm(vc_t)
    v_m_tl.append(v_m_t)
###################################
#
########################
# ODE Calculations End #
#  Cpu Clock           ##########
ODE_Calc_End_Time = time.time() #
#######################################
# ODE Calculations Execution Time ############################
ODE_Calc_Execution = ODE_Calc_End_Time - ODE_Calc_Start_Time #
##############################################################
ODE_message = f"""
{'-'*42}
| ODE Solved in:
|    {ODE_Calc_Execution} Seconds
{'-'*42}
| Plotting data, Please wait...
{'-'*42}
"""
print(ODE_message)
#%% Plots & Animation
##########
# Colors #
##########
#
Space      = "#000000"
Background = "#d3d4d3"
earth      = "#0000ff"
orbit_line = "#ff06b5"
#############################
# Orbital Elements Plotting #
################################################
# Set plot                                     #
fig = plt.figure('Saturn Atmosphere Drag')     #
# Set axis                                     #
ax1 = plt.axes()                               #
# Set Window Size                              #
fig.set_size_inches(5,5)                       #
ax1.plot(t_span/Graph_time,v_ml, color='red')  #
ax1.set_title("Saturn Atmosphere Drag")        #
ax1.set_xlabel('Time (hours)')                 #
ax1.set_ylabel('Magnitude of velocity (km/s)') # 
plt.show()                                     #
################################################
#
####################
# 3D Plot of Orbit #
#########################################
# Set plot                              #
fig = plt.figure('Orbit')               #
# Set axis                              #
ax1 = plt.axes(projection='3d')         #  
# Set Window Size                       #
fig.set_size_inches(12,10)              #
# Set x data to position, i. given by a #
xline = a[:,0]                          #
# Set y data to position, j. given by a #
yline = a[:,1]                          #
# Set z data to position, k. given by a #
zline = a[:,2]                          #
# Axis settings                         #
ax1.plot3D(xline, yline, zline, 'red')  #
# Axis Labels                           #
ax1.set_xlabel('x')                     #
ax1.set_ylabel('y')                     #
ax1.set_zlabel('z')                     #
# Axis Colors                           #
ax1.tick_params(axis='x', colors='red') #
ax1.tick_params(axis='y', colors='red') #
ax1.tick_params(axis='z', colors='red') #
ax1.yaxis.label.set_color('red')        #  
ax1.xaxis.label.set_color('red')        #  
ax1.zaxis.label.set_color('red')        #
# Background Color                      # 
fig.set_facecolor(Space)                #
ax1.set_facecolor(Space)                #
# Grid Pane Color/set to clear          #
ax1.xaxis.set_pane_color((0.0, 0.0,     #
                            0.0, 0.0))  #
ax1.yaxis.set_pane_color((0.0, 0.0,     #
                            0.0, 0.0))  #
ax1.zaxis.set_pane_color((0.0, 0.0,     #
                            0.0, 0.0))  #  
# Grid Color                            #
plt.rcParams['grid.color'] = 'red'      #                       
#########################################

# Show The Plots & Animation ###################################
plt.show()                                                     #
################################################################
