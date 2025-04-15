#########################################
# Author: Evan A. Blosser               #
# Email: evan.a.blosser-1@ou.edu        #
#                                       #
# Program: Mission to Titan             #
#            - X-fer calculations and   #
#              charting of a patched    #
#               conic trajectory        # 
#########################################
############
# Packages #
#####################################
# Mathmatical!!                     #
import numpy as np                  #
# ODE solver                        #
from scipy.integrate import odeint  #
# Plotting & Animation              #
import matplotlib.pyplot as plt     #
from matplotlib import animation    #
#####################################
#
#############
# Constants #
##############################
# Grav. Parameter for Earth  #  
mu_e = 3.986004418e5         #
# Grav. Parameter for Saturn #  
mu_sat = 37931187.0          #
# Equi.Radius of Earth in km #
r_e = 6378.136               #
# Equi. Radius of Saturn km  #
# Table A.1 Curtis           #
r_sat = 60270.0              #
##############################
# Earth Smi-major axis       #
R1 = 149.6e6                 #
# Saturn Smi-major axis      #
R2 = 1.433e9                 #
# Grav. P of Sun in km^3/s^2 #
mu_sun = 132712440018        #
# Accel. of Gravity (km/s^2) #
g_0 = 9.81e-3                #
##############################
#########
# Titan #
##################################
# Titan's Eccentricity less then
# 0.05 so it is approximated as zero 
e = 0.0288
# Grav. Parameter of Titan   #
#  Horizens System
mu_Titan = 8978.13710
##################################
############
### Main ###
############
#####################################
Welcome_Message = f"""
{'-'*42}
| Welcome to The Mission to Titan: KOMODO
{'-'*42}
"""                                 #
print(Welcome_Message)              #
#####################################
##################
## Calculations ##
##################
#%% Time
###################
# User Time Input #
###################
#   Calculation Errors:
#       Keeps second digit is always off by 1
#       ie. 2464819.0 days, was calculated as 
#           2464809.0 Days???
##############################################################################
# Time Inputs                                                                #
# Mission Time May 5th 2036 12:00:00                                         #
year  = 2036                                                                 #
month = 10                                                                   #
day   = 1                                                                    #
UT_H = 12                                                                    #
UT_M = 0                                                                     #
UT_S = 0                                                                     #
# Universal Time Calculations                                                #
UT = UT_H + (UT_M/60) + (UT_S/3600)                                          #
##############################################################################
# Julian/Sidereal Time #
##############################################################################
# Simplest Formula fo obtaining J0, Boulet(1991)  ############################
Julian_0 = 367*year - int((7*(year + int((month + 9)/12)))/4) +\
    int((257*month)/9) + day + 1721013.5
# Julian Date given Year, Month, and Day #####################################
Julian_D = Julian_0 + UT/24                                                  #
# T0 in Julian Cenbturies                                                    #
G_T0 = (Julian_0 - 2451545)/36525                                            #
# Greenwich Sidereal Time at 0h UT                                           #
Green_T0 = 100.4606184 + 36000.77004*G_T0 + 0.000387933*G_T0**2 -\
    2.583*10**(-8)*G_T0**3
# Greenwich Sidreal Time #####################################################
Green_Sidreal = Green_T0 + 360.98564724*(UT/24)                              #
##############################################################################
#
#################
# Synodic Period #######################################
T_Earth  = 365.256                                     #
T_Saturn = 29.46*T_Earth                               #
T_Synodic = (T_Earth*T_Saturn)/abs(T_Earth - T_Saturn) #
Mission_Window = f"""                               
{'-'*42}                                               
| The synodic Period between Earth and Saturn is:      
|    {T_Synodic} Days                                  
{'-'*42}                                               
"""                                                    #
print(Mission_Window)                                  #
########################################################
#%% Transfer to Saturn
######################################
# Altitude, Isp, and Spacecraft Mass #
R1_park_alt = 600
R2_park_alt = 591600.0
I_sp    = 311
Craft_mass = 1600
###################################
#################################
# Transfer From Earth To Jupiter #
#################################
RJ = 778.6e6
####################################################################
# Heliocentric velocity of craft on departure from earth           #
v_helio_craft_j = np.sqrt(mu_sun*((2/R1)-(2/(R1+RJ))))             #
# Earth's Heliocentric Velocity                                    #
v_helio_planet_j = np.sqrt(mu_sun/R1)                              #
# The velocity of the departure hyperbola                          #
v_inf_dep_j =  v_helio_craft_j - v_helio_planet_j                  #
# Geocentric Spacecraft velocity                                   #
v_geo_j = np.sqrt(mu_e/(R1_park_alt+r_e))                          #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.          #
v_geo_dep_j = np.sqrt(v_inf_dep_j**2 + (2*mu_e)/(R1_park_alt+r_e)) #
# Departure Delta-V                                                #
Delta_V_dep_j = v_geo_dep_j - v_geo_j                              #
####################################################################
#
#################################
# Transfer From Earth To Saturn #
#################################
####################################################################
# Heliocentric velocity of craft on departure from earth           #
v_helio_craft = np.sqrt(mu_sun*((2/R1)-(2/(R1+R2))))               #
# Earth's Heliocentric Velocity                                    #
v_helio_planet = np.sqrt(mu_sun/R1)                                #
# The velocity of the departure hyperbola                          #
v_inf_dep =  v_helio_craft - v_helio_planet                        #
# Geocentric Spacecraft velocity                                   #
v_geo = np.sqrt(mu_e/(R1_park_alt+r_e))                            #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.          #
v_geo_dep = np.sqrt(v_inf_dep**2 + (2*mu_e)/(R1_park_alt+r_e))     #
# Departure Delta-V                                                #
Delta_V_dep = v_geo_dep - v_geo                                    #
####################################################################
#####################
# Arrival at Saturn #
#####################
#######################################################################
# Heliocentric velocity of craft on departure from earth              #
v_helio_craft_arrive = np.sqrt(mu_sun*((2/R2)-(2/(R1+R2))))           #
# Earth's Heliocentric Velocity                                       #
v_helio_arrive_planet = np.sqrt(mu_sun/R2)                            #
# The velocity of the departure hyperbola                             #
v_inf_ari = v_helio_craft_arrive - v_helio_arrive_planet              #
# Geocentric Spacecraft velocity                                      #
v_geo_Rende = np.sqrt(mu_sat/(R2_park_alt+r_sat))                     #
# Geocentric Space craft Velocity at perigee of Dep. Hyp.             #
v_geo_ari = np.sqrt(v_inf_ari**2 + (2*mu_sat)/(R2_park_alt+r_sat))    #
# Departure Delta-V                                                   #
Delta_V_ari = v_geo_ari - v_geo_Rende                                 #
#######################################################################
############################################################
# Time of Filght from Earth to Saturn for Hohman           #
TOF_Sec = (np.pi)/(np.sqrt(mu_sun))*((R1 + R2)/2)**(3/2)   #
# Time of Flight in years                                  #
TOF_days = TOF_Sec/(60*60*24)                              #
TOF_years = TOF_Sec/(60*60*24*365.256)                     # 
############################################################
####################
# TOTAL DELTA-V    #
# Before Maneuvers #
######################
# Total Delta-V Req. ######################################
Delta_V_To_Saturn_0 = abs(Delta_V_dep) + abs(Delta_V_ari) #
###########################################################
# Propellent Req. for the departure ##################################
Fuel_Mass_Ratio_0 = 1 - np.exp(-Delta_V_To_Saturn_0/(I_sp*g_0))      #
# Spacecraft Initial mass Equation                                   #
Mass_Prop_0 = (Fuel_Mass_Ratio_0*Craft_mass)/(1 - Fuel_Mass_Ratio_0) #
######################################################################
#
#############
# Maneuvers #
######################################################
# Delta-V imparted to craft by Jupiter               #
Delta_V_from_Jupiter = Delta_V_dep - Delta_V_dep_j   #
# Delta-V from Atmo drag with Saturn and Titan       #
Delta_V_ari_atmo = 0                                 #
######################################################
# TOTAL DELTA-V   #
# After Maneuvers #
######################
# Total Delta-V Req. ###########################################
Delta_V_To_Saturn = abs(Delta_V_dep_j) + abs(Delta_V_ari_atmo) #
################################################################
# Propellent Req. for the departure ##############################
Fuel_Mass_Ratio = 1 - np.exp(-Delta_V_To_Saturn/(I_sp*g_0))      #
# Spacecraft Initial mass Equation                               #
Mass_Prop = (Fuel_Mass_Ratio*Craft_mass)/(1 - Fuel_Mass_Ratio)   #
##################################################################
##########################
# Mission Time in Julian ########################
# From Horizon's systems as min is not accurate #
Mission_Start_Time_Horizons = 2464968.0         ##########
# Add initial to final time                              #
Mission_Time = Mission_Start_Time_Horizons + TOF_days    # 
##########################################################
Mission_Time_Data_Message = f"""
{'-'*42}
| Mission TOF: {TOF_days} Days
|            : {TOF_years} Years
{'-'*42}
| Mission Start Time: {Mission_Start_Time_Horizons} Days
{'-'*42}
| Mission End Time: {Mission_Time} Days
{'-'*42}
"""                                    #
print(Mission_Time_Data_Message)       #
########################################
#%%Saturn Parking Radius Optimization
#
# Calcualte all Delta-V's varying them with respect to parking orbit
#  within the equations that the Pakring orbit is dependant upon.
###################################
# Set Optimization list for Earth #
Delta_V_To_Earth_opt = []
# Parking orbit as linespace
R1_park_opt = np.linspace(100, 1000, 100)
# Calcualte Delta-V with respect to Earth's Parking Orbit
for i in range(0,100):
    #################################
    # Transfer From Earth To Saturn #
    #################################
    ####################################################################
    # Geocentric Spacecraft velocity                                   #
    v_geo_opt = np.sqrt(mu_e/(R1_park_opt+r_e))                        #
    # Geocentric Space craft Velocity at perigee of Dep. Hyp.          #
    v_geo_dep_opt = np.sqrt(v_inf_dep**2 + (2*mu_e)/(R1_park_opt+r_e)) #
    # Departure Delta-V                                                #
    Delta_V_dep_opt = v_geo_dep_opt - v_geo_opt                        #
    # Total Delta-V Req. ##################################################
    Delta_V_To_Earth_optimize = abs(Delta_V_dep_opt) + abs(Delta_V_ari)   #
    #######################################################################
    Delta_V_To_Earth_opt.append(Delta_V_To_Earth_optimize)                #
###########################################################################
# Set plot                                                    #
fig_opt_e = plt.figure('Optimized Earth Parking Radius')      #
# Set axis                                                    #
axopt_e = plt.axes()                                          #
# Set Window Size and Titles                                  #
fig_opt_e.set_size_inches(5,5)                                #
axopt_e.set_title("Earth Parking Radius vs. Arrival Delta-V") #
axopt_e.set_xlabel('Parking Radius (km)')                     #
axopt_e.set_ylabel(' Delta-V from Earth to Saturn (km/s)')    # 
# Plot Data                                        ############
Del_line_Earth = Delta_V_To_Earth_opt[0]                    #
axopt_e.plot(R1_park_opt,Del_line_Earth)                    #
plt.legend()                                                #
#############################################################
# Set Optimization list for Saturn
Delta_V_To_Sat_opt = []
# Parking orbit as linespace
R2_park_opt = np.linspace(122187, 1221870, 100)
# Calcualte Delta-V with respect to Saturn's Parking Orbit
for i in range(0,100):
    #####################
    # Arrival at Saturn #
    #####################
    #######################################################################
    # Geocentric Spacecraft velocity                                      #
    v_geo_Rende_opt = np.sqrt(mu_sat/(R2_park_opt+r_sat))                 #
    # Geocentric Space craft Velocity at perigee of Dep. Hyp.             #
    v_geo_ari_opt = np.sqrt(v_inf_ari**2 + (2*mu_sat)/(R2_park_opt+r_sat))#
    # Departure Delta-V                                                   #
    Delta_V_ari_opt = v_geo_ari_opt - v_geo_Rende_opt                     #
    # Total Delta-V Req. ##################################################
    Delta_V_To_Sat_optimize = abs(Delta_V_dep) + abs(Delta_V_ari_opt)     #
    #######################################################################
    Delta_V_To_Sat_opt.append(Delta_V_To_Sat_optimize)                    #
###########################################################################
# Set plot                                                    #
fig_opt_s = plt.figure('Optimized Saturn Parking Radius')     #
# Set axis                                                    #
axopt_s = plt.axes()                                          #
# Set Window Size and Titles                                  #
fig_opt_s.set_size_inches(5,5)                                #
axopt_s.set_title("Saturn Parking Radius vs. Arrival Delta-V")#
axopt_s.set_xlabel('Parking Radius (km)')                     #
axopt_s.set_ylabel(' Delta-V from Earth to Saturn (km/s)')    # 
# Plot Data                                        ############
Del_line_Earth = Delta_V_To_Earth_opt[0]                    #
Del_line_sat = Delta_V_To_Sat_opt[0]                        #
axopt_s.plot(R2_park_opt,Del_line_sat)                      #
axopt_s.plot(R1_park_opt,Del_line_Earth)                    #
plt.legend()                                                #
plt.show() ##################################################
#%% Ephemeris
#############
# Ephemeris #
########################################################################
# October 1st
# Earth's Velocity vector at Transfer Start Time   #
v_E_ephem = np.array([-4.822006259693160e0, 
                      2.934625595766267e1, 
                      -1.780127279500832e-3])
# Saturn's Velocity vector at TRansfer Start Time
v_Sat_ephem = np.array([5.529729354646601e0, 
                        -7.549336343850271e0, 
                        -9.000393743228008e-2])
###############################################
#  Lambert's Output
# Craft Departure Velocity Vector
v_c_dep_ephem = np.array([-15.536719157854575,  
                          36.86395166588286,  
                          -3.164186373826042])
# Craft Fianl Velocity Vector
v_c_final_ephem = np.array([1.9635319022677116,
                            -3.4718844295142666,
                            0.30368160207068684])
################################################
# Delta-V Calculations                                                 #
V_ephem_depart_vec =  v_c_dep_ephem - v_E_ephem                        #
V_ephem_depart     = np.linalg.norm(V_ephem_depart_vec)                # 
V_ephem_arrive_vec = v_c_final_ephem - v_Sat_ephem                     #
V_ephem_arriv      = np.linalg.norm(V_ephem_arrive_vec)                #
Delta_V_To_Saturn_Ephemeris = abs(V_ephem_arriv) + abs(V_ephem_depart) #
########################################################################
#########################################################################
# May 5th
# Earth's Velocity vector at Transfer Start Time   #
v_E_ephem_may = np.array([2.062296080826162e1,    
                      -2.114282897288339e1, 
                      2.673577113423420e4])
# Saturn's Velocity vector at TRansfer Start Time
v_Sat_ephem_may  = np.array([4.891747698073565, 
                        -8.021447995695407, 
                        -5.565149276508308e2])
###############################################
#  Lambert's Output
# Craft Departure Velocity Vector
v_c_dep_ephem_may  = np.array([75.66103209904955, 
                          75.27936917353489,  
                          -0.19546888914148836])
# Craft Fianl Velocity Vector
v_c_final_ephem_may  = np.array([-81.9988374737535,  
                            -55.52613670496863, 
                            4.21041003238867])         
################################################                             
# Delta-V Calculations                                                       #
V_ephem_depart_vec_may =  v_c_dep_ephem_may - v_E_ephem_may                  #
V_ephem_depart_may     = np.linalg.norm(V_ephem_depart_vec_may)              #
V_ephem_arrive_vec_may = v_c_final_ephem_may - v_Sat_ephem_may               #
V_ephem_arriv_may      = np.linalg.norm(V_ephem_arrive_vec_may)              #
Delta_V_To_Saturn_Eph_may = abs(V_ephem_arriv_may) + abs(V_ephem_depart_may) #
##############################################################################
#%% Transfer to Titan
####################################
## Transfer/Rendezvous with Titan ##
####################################
# This transfer is from  Saturn to Titan
#  
#
#  For the finally delta-v we can use Saturn's Atmospheric
#   drag to slow the spacecraft down
#
#  Flyby Mercury instead of venus?
#   Maybe this will get us there faster and more efficiantly... 
#
#####################################
# Hohmann to Titan (Saturn Central) #
#####################################
# Parking orbit of spacecraft
R_Craft = R2_park_alt+r_sat
#Titan's Semi-major Axis
A_Titan = 1221870.0           
# 18 Days? Titan's Orbital Period is 16 days around Saturn 
Period_Titan = 16*24*60*60
# Orbital period at the R2_parking_alt
Period_SpaceCraft = 2*np.pi*np.sqrt(R_Craft**3/mu_sat)   
# Hohman Transfer Window for Satrun orbit to Titan #
Window_Titan = (Period_SpaceCraft*Period_Titan)/\
    abs(Period_SpaceCraft - Period_Titan)
Window_Titan_Days = Window_Titan /(24*60*60)
####################################################                                
##################################################################
# Delta-V equations to leave parking orbit around Saturn         #
vel_Tit_craft = np.sqrt(mu_sat/R_Craft)                          #
E_Craft = -mu_sat/(2*R_Craft)                                    #
vel_Tit_craft_dep = np.sqrt(2*(mu_sat/R_Craft + E_Craft))        #
Delta_V_Tit_a = vel_Tit_craft_dep - vel_Tit_craft                #
# Delta-V equations to arrive at Titan                           #
vel_Tit_craft_f = np.sqrt(mu_sat/A_Titan)                        #
E_Craft_f = -mu_sat/(2*A_Titan)                                  #
vel_Tit_craft_arri = np.sqrt(2*(mu_sat/A_Titan + E_Craft))       #
Delta_V_Tit_b = vel_Tit_craft_f - vel_Tit_craft_arri             #
# Delta-V for Hohmann transfer to Titan                          #
Delta_V_Titan_Transfer = abs(Delta_V_Tit_a) + abs(Delta_V_Tit_b) #
# Total Delta-V for Missin                                       #
Delta_V_Total = Delta_V_To_Saturn + Delta_V_Titan_Transfer       #
#################################################################
############################
# Semi major axis of Titan ##################################
a_Titan = ((Period_Titan*np.sqrt(mu_sat))/(2*np.pi))**(2/3) #
#############################################################
#%%Jupiter Flyby
#########################
# Question 5 Homework 5 #
#########################
##################################
# Fly by Altitude around Jupiter #
r_flyby_5 = 200000               #
##################################
# Earth's Semi_Major Axis     #
R_E_5 = 149.6e6               #
# Jupiter's Semi_Major Axis   #
R_J_5 = 778.6e6               #
# Jupiter's Grav. Per.        #
mu_j = 126686534              #
# Jupiter's Radius            #
r_j_5  = 71490                #
###########################################
# Assumed Hohmann transfer from Earth     #
#   - in order to define angular momentum #
###########################################
# Semi-Major Axis of Transfer       #
semi_major_5_hohm = (R_E_5+R_J_5)/2 # 
# Jupiter's Speed                   #
V_Jup_5 = np.sqrt(mu_sun/R_J_5)     #
#  Heliocentric Spacecraft Velocity ###########################
V_Craft_5 = np.sqrt(mu_sun*((2/R_J_5)-(1/semi_major_5_hohm))) #
# Excess Speed upon Jupiter Arrival                           #
V_inf_5 = V_Jup_5 - V_Craft_5                                 #
# Eccentricity of Fly By Entry                                #
e_5 = 1 + ((r_j_5 + r_flyby_5)*(V_inf_5)**2)/mu_j             #
# Turning Angle                                               #
turn_ang_rad = 2*np.arcsin(1/e_5)                             #
turn_ang = np.rad2deg(turn_ang_rad)                           #
# Delta-V Transfered to the spacecraft (!Use Rad!)            #
Delta_V_5 = 2*V_inf_5*np.sin(turn_ang_rad/2)                  #
###############################################################
###################################
# Calculate Exit Orbital Elements #
############################################
# Turning angle for the outbound direction ################## 
Outbound_ang = 180 + turn_ang              
Outbound_ang_rad = np.deg2rad(Outbound_ang)                 
# Outbound Velocity Vector
V_outbound_Uv = V_Jup_5 + V_inf_5*np.cos(Outbound_ang_rad)
V_outbound_Us = V_inf_5*np.sin(Outbound_ang_rad)
# Magnitude/Normal of the outbound Velocity
V_5_outbound = np.sqrt(V_outbound_Uv**2 + V_outbound_Us**2)
#########
# Angular Momentum at exit
#    -assuming we have made a new apogee/perigee
h_5   = R_J_5*V_outbound_Uv
######################
# Eccentricity is Calculated as shown in Lec. 20230410
#  were we use relations between the eccentricity and
#   true anomaly at exit
#####################################################
# e*cos relation                                    #
Ecos = h_5**2/(mu_sun*R_J_5) - 1                    #
# e*sin relation                                    #
Esin = (V_outbound_Us*h_5)/mu_sun                   #
# divide both relations to cancel out eccentricity  #
Tan_theta_5 = Esin/Ecos                             #
# Solve for True anomaly at exit                    #
Theta_5_2   =np.arctan(Tan_theta_5)                 #
# New Eccenmtricity calculated with Theta_2         #
Ecc_5 =  Ecos/np.cos(Theta_5_2)                     #
# New Perihelion of Spacecraft                      #
R_perihelion_5 = (h_5**2/mu_sun)*(1/(1 + Ecc_5))    #
# New Semi-Major Axis Calculations using new        #
#   eccentricity and angular momentumn              #
a_5 =  R_perihelion_5/(1 - Ecc_5)                   #
#####################################################
#%% Output and Ephemris Plot
###############################################################################
Data_Message = f"""
{'-'*42}
| Window_Titan : {Window_Titan_Days} Days
| assuming the period of Titan is 16 days
{'-'*42}
{'#'*42}
{'-'*42}
| Hohmann Transfer
{'-'*42}
| Initial Calculations Without Maneuvers
{'-'*42} 
| Delta_V_To_Saturn: {Delta_V_To_Saturn_0} (km/s)
| Fuel_Mass_Ratio: 
|    {Fuel_Mass_Ratio_0}
| Propellant Mass:
|   {Mass_Prop_0} (kg)
|{'-'*42}
| Earth to Jupiter: {Delta_V_dep_j} (km/s)
|{'-'*42}
| --To Saturn--
| Earth Departure: {Delta_V_dep} (km/s) 
| 
| Saturn Arrival: {Delta_V_ari} (km/s) 
{'-'*42}
| Delta-V imparted to craft by jupiter:
|    {Delta_V_from_Jupiter} (km/s)
|
{'-'*42}
| Final Calculations With Maneuvers 
|  for a craft with Isp = {I_sp} (seconds)
|  with a mass of:  {Craft_mass} (kg) 
{'-'*42}
| Delta_V_To_Saturn: {Delta_V_To_Saturn} (km/s)
|
| Delta-V to Titan : {Delta_V_Titan_Transfer} (km/s)
|
| Total Mission Delta-V to Titan:
|    {Delta_V_Total} (km/s)
|
| Fuel_Mass_Ratio: 
|    {Fuel_Mass_Ratio}
| Propellant Mass:
|   {Mass_Prop} (kg)
|
{'-'*42}
|  Jupiter Flyby
{'-'*42}
|    The Semi-Major axis of the post
|        flyby trajectory is:
|            {a_5} (km) 
{'-'*42}
|    The Eccentricity of the post 
|       flyby trajectory is:
|           {Ecc_5}   
{'-'*42}
|    The Delta-V imparted on the
|       spacecraft is:
|         {Delta_V_5} (km/s)
{'-'*42}
{'#'*42}
{'-'*42}
| Ephemris Calculations
|{'-'*42}
| October 1st, 2036
| Earth Departure: {V_ephem_depart} (km/s)
|
| November 1st,2036
| Saturn Arrival: {V_ephem_arriv} (km/s)
|
| Total Delta-V: {Delta_V_To_Saturn_Ephemeris} (km/s)
{'-'*42}
| May 5th, 2036
|
| Earth Departure: {V_ephem_depart_may} (km/s)
|
| Saturn Arrival: {V_ephem_arriv_may} (km/s)
|
| Total Delta-V: {Delta_V_To_Saturn_Eph_may} (km/s)
{'-'*42}
"""
print(Data_Message)
#########################################
#####################
# Ephemris Plotting ########################################
Earth_Loc_Start = [148101352.8509693, 
                   20706303.64949361,
                   1837.641222303733]
Satrun_Loc_End  = [-1145097770.344292, 
                   -919599022.8293214,
                   61545527.29974103]
Earth_Loc_mag   = np.linalg.norm(Earth_Loc_Start)
Satrun_Loc_mag  = np.linalg.norm(Satrun_Loc_End)
############
# Settings #
##########################################
# Time and Points                        #
t_inp      =  6.08                       #
Graph_time = 60*60*24*365                #
n     = t_inp*Graph_time                 #
# Graph points, set to the chosen unit   # 
points = Graph_time                      #
t_span  = np.linspace(0,n,points)        #
v_c_dep_ephem                            #
##########################################
###############
## ODE Setup ##
#########################
r_geo = Earth_Loc_Start #
v_geo = v_c_dep_ephem   #
#########################
#####################
# Define ODE System #
def orbit(a,t):
      #
      ##########################################
      r_geo_vec = [0]*3                        # Calcualte new position 
      for j in range(0,3):                     #  norm/magnitude with
          r_geo_vec[j] = a[j]                  #  the current position!
          r_mag = np.linalg.norm(r_geo_vec)    #
      ##########################################
      rx  = a[0] # Define Initial Position
      ry  = a[1] #   x, y, z component
      rz  = a[2] #
      vx  = a[3] # Define Initial Velocity
      vy  = a[4] #   x, y, z component
      vz  = a[5] #
      ##################################################
      # Derivatives/ 6 Equations of Motion #############
      drxdt = vx 
      drydt = vy 
      drzdt = vz 
      dvxdt = -mu_sun*(rx/r_mag**3) 
      dvydt = -mu_sun*(ry/r_mag**3) 
      dvzdt = -mu_sun*(rz/r_mag**3) 
      dadt  = [drxdt,drydt,drzdt,dvxdt,dvydt,dvzdt]
      return dadt
# End ODE Definition ###################################                                                    
######################
######################
# Initial Conditions #########################################
a0 = [r_geo[0],r_geo[1],r_geo[2],v_geo[0],v_geo[1],v_geo[2]] #
##############################################################
################
# ODE Solution ##############
a = odeint(orbit,a0,t_span) #
#############################
#
####################
# 2D Plot of Orbit #
#########################################
# Set plot                              #
fig = plt.figure('Earth to Saturn')     #
# Set axis                              #
ax1 = plt.axes()                        #
# Set Window Size                       #
fig.set_size_inches(5,5)                 #
# Set Theta for drawing circles           #
theta = np.linspace( 0 , 2 * np.pi , 1000) #
# Earth's Assumed circular orbit           #
Earth_x = R1*np.cos(theta)                 #
Earth_y = R1*np.sin(theta)                 #
# Saturn's Assumed circular orbit          #
Sat_x = R2*np.cos(theta)                   #
Sat_y = R2*np.sin(theta)                   #
# Plot Data ###################################################
ax1.plot( Earth_y, Earth_x,linestyle='dashed',label='Earth')  #
ax1.plot(Sat_y, Sat_x,linestyle='dashed',label='Saturn')      #
ax1.plot(0, 0, 'yo',label='Sun')                              #
ax1.plot(a[:,0],a[:,1],label='Komodo')                        #
plt.legend()                                                  #
plt.show()                                                    #
###############################################################