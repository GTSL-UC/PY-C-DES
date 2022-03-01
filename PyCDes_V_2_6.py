# PYCDES Script             --- Version 2.6         #
# University of Cincinnati  --- GTSL                #
# By Kaden Wells            --- wellskc@mail.uc.edu #
#                                                   #
#                                                   #
# sCO2 Properties Incorporated                      #
# Actual Gas                                        #
# Last Date Edited          --- 02/28/2022          #
#####################################################
#####################################################
''' *************** TABLE OF CONTENTS *************** 
PART A) IMPORT PYTHON PACKAGES
PART B) MAIN FUNCTION
        i)    FILE CREATION
        ii)   READ INPUT FILES
            1)      TCI
            2)      TCS
            3)      OPTIONAL IGV
            4)      IMPORT TABLES BASED ON FLUID TYPE
        iii)  CALL INLET FUNCTION
        iv)   GENERATE ARRAYS TO HOLD OUTPUT VALUES
        v)    FOR LOOP GOES THROUGH EACH COMPRESSOR STAGE
            1)      Call Rotor Outlet/Stator Inlet Function
            2)      Rotor Geometry 
            3)      Call Stator Outlet/Rotor Inlet Function
            4)      Stator Geometry
            5)      Cross Section Geometry
            6)      Exporting Blade Data to .TINF File
        vi)   EXPORT GEOMETRIC VECTORS 
        vii)  EXPORT FLUID + COMPRESSOR PROPERTIES
        viii) INSENTROPIC CALCULATIONS
        ix)   FINAL COMPRESSOR DATA OUTPUT
        x)    REFERENCE VALUES FOR .TSTK and .TWAL files
        xi)   CALL WALL FILE FUNCTION
        xii)  CALL STACK FILE FUNCTION
        xiii) CALL PLOTS FUNCTION
PART C) INLET FUNCTION
        i)    WITH AN IGV
            1)      IGV Calculations
            2)      Rotor Inlet Calculations
            3)      IGV Geometry
            4)      IGV Cross Section Geometry
        ii)   WITHOUT AN IGV
            1)      Rotor Inlet Calculations
            2)      (no) IGV Geometry
        iii)  DISPLAY DATA
PART D) ROTOR OUTLET FUNCTION
        i)    ROTOR OUTLET CALCULATIONS (REAL GAS)
        ii)   ROTOR OUTLET CALCULATIONS (PERFECT GAS)
        iii)  MISCELLANEOUS CALCULATIONS
        iv)   DISPLAY DATA
PART E) STATOR OUTLET FUNCTION
        i)    CHECK IF LAST STAGE
        ii)   STATOR OUTLET CALCULATIONS (REAL GAS)
        iii)  STATOR OUTLET CALCULATIONS (PERFECT GAS)
        iv)   MISCELLANEOUS CALCULATIONS
        v)    DISPLAY DATA
PART F) RADIUS TYPE FUNCTION
PART G) VELOCITY TRIANGLES FUNCTION
PART H) STATIC PROPERTIES CONVERSION FUNCTION
PART I) TOTAL PROPERTIES CONVERSION FUNCTION
PART J) WALLS FUNCTION
PART K) STACK FUNCTION
PART L) EXTRAPOLATION FUNCTION
PART M) sCO2 PROPERTIES FUNCTION
PART N) AIR PROPERTIES FUNCTION
PART O) CH4 PROPERTIES FUNCTION
PART P) 'OTHER' PROPERTIES FUNCTION
PART Q) PLOTTING FUNCTION
        i)      CREATE AND WRITE TO PDF FILE
END
                                                '''
#####################################################
#####################################################
# PART A - IMPORT PYTHON PACKAGES 
#####################################################
#####################################################
import sys
import math
import numpy
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
import os
from datetime import date
today = date.today()
# mm/dd/y
day = today.strftime("%m_%d_%y")
print("date ran =", day)

#####################################################
#####################################################
# PART B - MAIN FUNCTION
#####################################################
#####################################################

def main(a):
    
#####################################################
#####################################################
# PART B.i - FILE CREATION 
#####################################################
#####################################################
    
#####################################################
#####################################################
# PART B.ii - READ IMPORT FILES
#####################################################
#####################################################

    #Name of INPUT
    INPUT = str(a)

#####################################################
#####################################################
# PART B.ii.1 - READ .TCI FILE
#####################################################
#####################################################
    
    tci    = open(f'{INPUT}' + '.tci', 'r')
    lines  = tci.readlines()
    result = [] 
    
    for x in lines:
    # 
        #      
        result.append(x.split(' ')[0])
        #
    #
    tci.close()

    typer    = int(result[1])            #Rad Type, 1- Hub, 2- Pitch, 3- Tip, 4- Arbitrary
    n        = float(result[2])          #Num of Stages
    MFR      = float(result[3])          #Mass Flow Rate {kg/s}
    Rev      = float(result[4])          #Angular Velocity {rpm}
    Poi      = float(result[5])          #Total Pressure {Pa}
    Toi      = float(result[6])          #Total Temperature {K}
    alphaend = float(result[7])          #Flow Angle {deg}
    Mend     = float(result[8])          #Mach at Outlet {-}
    GammaCUS = float(result[9])          #{-}     For custom fluid
    RgasCUS  = float(result[10]) * 1000  #{J/kgK} For custom fluid
    TCLC     = float(result[11])         #Ratio of Clearance {-}
    Radius   = float(result[12])         #{m} Radius used if Constant Hub, Pitch, or Tip selected
    Fluid     = float(result[13])        #Fluid --- Real == 1, Perfect Gas Eqns == 0
    Fluidtype = float(result[14])        #Fluidtype --- CO2 == 1, Air == 0, CH4 == 2, H2 == 3, O2 == 4, Other3 == 5

    KKI      =  0                        #Tracks if the radius is 'arbitrary'
    if typer == 1 : 
        typert = 'Hub'
    elif typer == 2 : 
        typert = 'Pitch'
    elif typer == 3 : 
        typert = 'Tip'
    elif typer == 4 : 
        typert = 'Arbitrary'
        typer  = 3                       #Arbitrary uses 'constant tip'
        KKI    = 1
    print(f'Radius Type Constant <{typert}> selected \n')
    print(' ')
    if Fluid == 0: 
        print('Perfect Gas Assumption Selected \n')
    elif Fluid == 1: 
        print('Real Gas Properties Selected \n')
    print(f'Number of Stages: {n}')
        
#####################################################
#####################################################
# PART B.ii.2 - READ .TCS FILE
#####################################################
#####################################################

    tcs          = open(f'{INPUT}' + '.tcs', 'r')
    lines2       = tcs.readlines()
    L            = len(lines2)
    result2      = []
    Stage_Number = []
    alpha1       = []
    MACH1        = []
    DH           = []
    R_Phi        = []
    S_Phi        = []
    R_Sol        = []
    S_Sol        = []
    Lambda       = []
    Bleed        = []
    R_AR         = []
    S_AR         = []
    VZR          = []
    R_Coef_Space = []
    S_Coef_Space = []
    Rad_tip      = []
    b            = 0
    
    while b < L : 
    #
        #
        Le  = []
        Le2 = []
        Le  = lines2[b]
        Le2 = Le.split()
        result2.append(Le2) 
        b += 1
        #
    #
    tcs.close()
    t = 0
    while t < (n):
    #
        #
        #If there's going to be 10+ stages, make sure their is data for n+1 stages in the compressor. 
        Stage_Number.append(float(result2[1][t]))    #Num of Stages
        alpha1.append      (float(result2[2][t]))    #Flow Angle {deg}
        MACH1.append       (float(result2[3][t]))    #Mach Number
        DH.append          (float(result2[4][t]))    #{kJ/kg} or {K} Depending on Setting
        R_Phi.append       (float(result2[5][t]))    #Rotor Phi Coef
        S_Phi.append       (float(result2[6][t]))    #Stator Phi Coef
        R_Sol.append       (float(result2[7][t]))    #Rotor Solidity
        S_Sol.append       (float(result2[8][t]))    #Stator Solidity
        Lambda.append      (float(result2[9][t]))    #Stage Lambda
        Bleed.append       (float(result2[10][t]))   #Stage Bleed {%}
        R_AR.append        (float(result2[11][t]))   #Rotor Aspect Ratio
        S_AR.append        (float(result2[12][t]))   #Stator Aspect Ratio
        VZR.append         (float(result2[13][t]))   #Rotor Axial Velo Ratio
        R_Coef_Space.append(float(result2[14][t]))   #Rotor Row Spacing Coef
        S_Coef_Space.append(float(result2[15][t]))   #Stator Row Spacing Coef
        Rad_tip.append     (float(result2[16][t]))   #Tip Radius [Arbitrary] {m}
        
        if Fluid == 1:
            Lambda[t] = Lambda[t]# - 0.03 * Lambda[t]
            R_Phi[t]  = R_Phi[t] #* 3.5
            S_Phi[t]  = S_Phi[t] #* 3.5
            DTT = DH
        t += 1
        
        if Fluid == 0: 
            DTT = DH
        
        #Changing it to work with current datasets below
        #DTT          = [32.625, 32.125, 31.125 ,30.125, 29.625, 29.125, 28.625, 28.625, 28.625, 48.255] #New temp data
        #DTT          = [37, 36, 35.5, 34.5, 35.0, 34.0, 33.0, 47.269] #Old temp data
        #
    #
    
#####################################################
#####################################################
# PART B.ii.3 - READ .IGV FILE
#####################################################
#####################################################
    
    Q        = 0
    while Q == 0 :
    #
        #
        #Check if IGV file exists in Folder
        ANS = os.path.isfile(f'{INPUT}' + '.igv')
        
        if ANS == True :
        #
            # 
            Q       = 1
            IGV_Q   = 1
            igv     = open(f'{INPUT}' + '.igv', 'r')
            lines3  = igv.readlines()
            result3 = []
            for x in lines3:
            # 
                #      
                result3.append(x.split(' ')[0])
                #
            #
            igv.close()
            igv_S       = float(result3[0])    #Solidity
            igv_AR      = float(result3[1])    #Aspect Ratio
            igv_Phi     = float(result3[2])    #Phi Loss Coef
            igv_a       = float(result3[3])    #Inlet Angle {deg}
            igv_M       = float(result3[4])    #Inlet Mach
            igv_L       = float(result3[5])    #Lambda 
            igv_R_Coef  = float(result3[6])    #Row Space Coef
            igv_rad_tip = float(result3[7])    #Tip Radius {m}
            if Fluid == 1:
                igv_L = igv_L - 0.03 * igv_L 
                #
            #
        #
        elif ANS == False :
        #
            #
            Q           = 2
            IGV_Q       = 0
            igv         = '**** No IGV Selected ****'
            igv_S       = None   #Solidity
            igv_AR      = None   #Aspect Ratio
            igv_Phi     = None   #Phi Loss Coef
            igv_a       = None   #Inlet Angle {deg}
            igv_M       = None   #Inlet Mach
            igv_L       = Lambda[0]   #Lambda 
            igv_R_Coef  = None   #Row Space Coef
            igv_rad_tip = None   #Tip Radius {m}
            print(igv)
            print(' ')
            #
        #
    #

#####################################################
#####################################################
# PART B.ii.4 - IMPORT TABLES BASED ON FLUID TYPE
#####################################################
#####################################################
    
    if Fluidtype == 1:
        
        '''IMPORT sCO2 TABLES'''
        HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT = sCO2_Properties()
        Fluidtypename = 'sCO2'

        Gamma    = 1.289  #{-}
        Rgas     = 188.9  #{J/kgK}
        
    elif Fluidtype == 0 or Fluidtype == 5:
        '''IMPORT AIR TABLES'''
        HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT = AIR_Properties()
        Fluidtypename = 'Air'
        
        Gamma    = 1.4    #{-}
        Rgas     = 287    #{J/kgK}
                
        if Fluidtype == 5: 
            Gamma    = GammaCUS  #{-}
            Rgas     = RgasCUS   #{J/kgK}
        
    elif Fluidtype == 2: 
        
        '''IMPORT METHANE TABLES'''
        HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT = CH4_Properties()
        Fluidtypename = 'CH4'
        
        Gamma    = 1.299  #{-}
        Rgas     = 518    #{J/kgK}
        
    else: 
        if Fluidtype == 3: 
            Fluidtypename = 'H2'
            
            Gamma    = 1.405  #{-}
            Rgas     = 4124   #{J/kgK}
        
        elif Fluidtype == 4: 
            Fluidtypename = 'O2'
            
            Gamma    = 1.395  #{-}
            Rgas     = 259.8  #{J/kgK}
        
        # elif Fluidtype == 5: 
        #     Fluidtypename = 'H2'
        #     Gamma    = GammaCUS  #{-}
        #     Rgas     = RgasCUS   #{J/kgK}
        
        '''IMPORT OTHER TABLES'''
        HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT = OTHER_Properties(Fluidtypename)  
    
    Output = f'{Fluidtypename}_' + day
    print(f'Fluid Type Name = {Fluidtypename} \n')
    
    #First Stage == 0 in Python Language
#####################################################
#####################################################
# PART B.iii - CALL INLET FUNCTION
#####################################################
#####################################################
    if n <= 0: 
        print('Please put in a real amount of stages')
    STG = int(Stage_Number[0] - 1)
    
    '''INLET CALCULATIONS - Call Inlet_Calculations Function'''
    Poi_initial      = Poi / 1000   #{kPa}
    Toi_initial      = Toi          #{K}
    CP, Tsi, Psi, Rsi, Ai, Vi, ALPHA, Vzi, A_comp, Radt, Height, Radp, RAV, U, Vti, Wti, Wzi, BETA, Wi, TTR, PTR, MR, Q, Poi, H, S, Po_h, S_constant, ALPHAd, BETAd, Toi_h, Poi_h, Radh, Psi_h, M_h, Vti_h, Vzi_h, ALPHA_h, BETA_h, Toi_t, Poi_t, Tsi_t, Psi_t, Rsi_t, M_t, Vi_t, Vti_t, Vzi_t, ALPHA_t, BETA_t, Tsi_p, Toi_p, Psi_p, Poi_p, Rsi_p, Mi_p, Vi_p, Vzi_p, Vti_p, ALPHAd, BETAd, Rsi_p, UU, Wi_p, Wti_p, Wzi_p, h, GammaValue, Viscosity, Wi_t, Wti_t, Wzi_t, CCPP, Height_1, XRIF, X_igv, Y_igv, X_igv_h, Y_igv_h, X_igv_t, Y_igv_t, X_igv_blades, Y_igv_blades, SNB, H_igv, S_igv = Inlet_Calculations(STG, MFR, Rev, Poi/1000, Toi, Gamma, Rgas, TCLC, Radius, Stage_Number, MACH1, alpha1, Q, igv_S, igv_AR, igv_Phi, igv_a, igv_M, igv_L, igv_R_Coef, igv_rad_tip, typer, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, Fluid, VPT, Output, Rad_tip, KKI)
    Current_Pressure = Poi
    H_initial        = H            #{J/kg}
    Height_1         = Height       #{m}
    print(f'Toi Initial    {Toi}')
    print(f'Toi_h          {Toi_h}')

#####################################################
#####################################################
# PART B.iv - GENERATE ARRAYS FOR OUTPUTS
#####################################################
#####################################################
    
    '''FOR PLOTS'''
    Pplot            = []
    Tplot            = []
    Rplot            = []
    Splot            = []
    Hplot            = []
    LOADING_C        = []
    FLOW_C           = []
    RELATIVE_M       = []
    Pplot.append(Poi)
    Tplot.append(Toi)
    Rplot.append(Rsi)
    Hplot.append(H)
    Splot.append(S) 
    LOADING_C.append(DH[0]*1000 / UU**2)
    FLOW_C.append(Vzi_p / UU)
    RELATIVE_M.append(MR)
    
    '''FOR EXPORT - PITCH LINE'''
    RAD              = []
    PT               = []
    PS               = []
    TT               = []
    TS               = []
    HH               = []
    HS               = []
    DHH              = []
    SS               = []
    RS               = []
    VV               = []
    VZV              = []
    VTV              = []
    WW               = []
    WZW              = []
    WTW              = []
    ALP              = []
    BET              = []
    MACHSP           = []
    CAMBER           = []
    DIFFUSION        = []
    Value_Reaction   = []
    GAMMAVALUE       = []
    VISCOSITY        = []
    SPECIFICHEAT     = []
    RAD.append(round(Radp,4))
    PT.append(round(Poi_p,2))
    PS.append(round(Psi_p,2))
    TT.append(round(Toi_p,2))
    TS.append(round(Tsi_p,2))
    HH.append(round(H,2))
    HS.append(round(h,2))
    SS.append(round(S,2))
    RS.append(round(Rsi_p,2))
    VV.append(round(Vi_p,2))
    VZV.append(round(Vzi_p,2))
    VTV.append(round(Vti_p,2))
    WW.append(round(Wi_p,2))
    WZW.append(round(Wzi_p,2))
    WTW.append(round(Wti_p,2))
    ALP.append(round(ALPHAd,2))
    BET.append(round(BETAd,2))
    MACHSP.append(round(Mi_p, 3))
    CAMBER.append(0)
    DIFFUSION.append(0)
    GAMMAVALUE.append(round(GammaValue, 5))
    VISCOSITY.append(round(Viscosity*10**6, 5))
    SPECIFICHEAT.append(round(CCPP, 5))
    
    '''FOR EXPORT - EFFICIENCY'''
    S_isentropic = S
    H_initial_2  = H
    AD_EFF_perstage = []
    
    '''FOR EXPORT - HUB'''
    RADh             = []
    PTh              = []
    PSh              = []
    TTh              = []
    TSh              = []
    RSh              = []
    VVh              = []
    VZVh             = []
    VTVh             = []
    WWh              = []
    WZWh             = []
    WTWh             = []
    ALPh             = []
    BETh             = []
    MACHSH           = []
    RADh.append(round(Radh,4))
    PTh.append(round(Poi_h,2))
    PSh.append(round(Psi_h,2))
    TTh.append(round(Toi_h,2))
    TSh.append(round(Tsi,2))
    RSh.append(round(Rsi,2))
    VVh.append(round(Vi,2))
    VZVh.append(round(Vzi_h,2))
    VTVh.append(round(Vti_h,2))
    WWh.append(round(Wi,2))
    WZWh.append(round(Wzi,2))
    WTWh.append(round(Wti,2))
    ALPh.append(round(ALPHA_h,2))
    BETh.append(round(BETA_h,2))
    MACHSH.append(round(M_h, 3))
    
    '''FOR EXPORT - TIP'''
    RADt             = []
    PTt              = []
    PSt              = []
    TTt              = []
    TSt              = []
    RSt              = []
    VVt              = []
    VZVt             = []
    VTVt             = []
    WWt              = []
    WZWt             = []
    WTWt             = []
    ALPt             = []
    BETt             = []
    MACHST           = []
    RADt.append(round(Radt,4))
    PTt.append(round(Poi_t,2))
    PSt.append(round(Psi_t,2))
    TTt.append(round(Toi_t,2))
    TSt.append(round(Tsi_t,2))
    RSt.append(round(Rsi_t,2))
    VVt.append(round(Vi_t,2))
    VZVt.append(round(Vzi_t,2))
    VTVt.append(round(Vti_t,2))
    WWt.append(round(Wi_t,2))
    WZWt.append(round(Wzi_t,2))
    WTWt.append(round(Wti_t,2))
    ALPt.append(round(ALPHA_t,2))
    BETt.append(round(BETA_t,2))
    MACHST.append(round(M_t, 3))
    
    '''FOR EXPORT - ROTOR GEOMETRY VECTORS'''
    STAGGERANGLE1_h = []      #GAMMA ANGLE AT ROTOR INLET (HUB) RAD
    STAGGERANGLE1_p = []      #GAMMA ANGLE AT ROTOR INLET (PITCH) RAD
    STAGGERANGLE1_t = []      #GAMMA ANGLE AT ROTOR INLET (TIP) RAD 
    STAGGERANGLE2_h = []      #GAMMA ANGLE AT ROTOR OUTLET (HUB) RAD
    STAGGERANGLE2_p = []      #GAMMA ANGLE AT ROTOR OUTLET (PITCH) RAD
    STAGGERANGLE2_t = []      #GAMMA ANGLE AT ROTOR OUTLET (TIP) RAD
    ZETAANGLE_h     = []      #ZETA ANGLE AT ROTOR (HUB) RAD
    ZETAANGLE_p     = []      #ZETA ANGLE AT ROTOR (PITCH) RAD
    ZETAANGLE_t     = []      #ZETA ANGLE AT ROTOR (PITCH) RAD
    ZETAANGLE_hD    = []      #ZETA ANGLE AT ROTOR (HUB) DEGREES
    ZETAANGLE_pD    = []      #ZETA ANGLE AT ROTOR (PITCH) DEGREES
    ZETAANGLE_tD    = []      #ZETA ANGLE AT ROTOR (PITCH) DEGREES
    AXIAL_WIDTH_R   = []      #millimeters
    CHORD_LENGTH_R  = []      #millimeters
    RAR_R           = []      #Real Aspect Ratio
    ROTORWIDTH_h    = []      #Width of Rotor (HUB) (mm)
    ROTORWIDTH_p    = []      #Width of Rotor (PITCH) (mm)
    ROTORWIDTH_t    = []      #Width of Rotor (TIP) (mm)
    NBR             = []      #NUMBER OF BLADES ON THE ROTOR
    
    '''FOR EXPORT - STATOR GEOMETRY VECTORS'''
    STAGGERANGLE1_Sh = []      #GAMMA ANGLE AT ROTOR INLET (HUB) RAD
    STAGGERANGLE1_Sp = []      #GAMMA ANGLE AT ROTOR INLET (PITCH) RAD
    STAGGERANGLE1_St = []      #GAMMA ANGLE AT ROTOR INLET (TIP) RAD 
    STAGGERANGLE2_Sh = []      #GAMMA ANGLE AT ROTOR OUTLET (HUB) RAD
    STAGGERANGLE2_Sp = []      #GAMMA ANGLE AT ROTOR OUTLET (PITCH) RAD
    STAGGERANGLE2_St = []      #GAMMA ANGLE AT ROTOR OUTLET (TIP) RAD
    ZETAANGLE_Sh     = []      #ZETA ANGLE AT ROTOR (HUB) RAD
    ZETAANGLE_Sp     = []      #ZETA ANGLE AT ROTOR (PITCH) RAD
    ZETAANGLE_St     = []      #ZETA ANGLE AT ROTOR (PITCH) RAD
    ZETAANGLE_ShD    = []      #ZETA ANGLE AT ROTOR (HUB) DEGREES
    ZETAANGLE_SpD    = []      #ZETA ANGLE AT ROTOR (PITCH) DEGREES
    ZETAANGLE_StD    = []      #ZETA ANGLE AT ROTOR (PITCH) DEGREES
    AXIAL_WIDTH_S    = []      #millimeters
    CHORD_LENGTH_S   = []      #millimeters
    RAR_S            = []      #Real Aspect Ratio
    STATORWIDTH_h    = []      #Width of Rotor (HUB) (mm)
    STATORWIDTH_p    = []      #Width of Rotor (PITCH) (mm)
    STATORWIDTH_t    = []      #Width of Rotor (TIP) (mm)
    NBS              = []      #NUMBER OF BLADES ON THE STATOR
    
    '''FOR EXPORT - ROTOR & STATOR CROSS SECTION GEOMETRY VECTORS'''
    X_CROSS_SECTION  = []
    Y_CROSS_SECTION  = []
    X_BLADES         = []      #Getting Arrays for the Blades for Graphing
    Y_BLADES         = []      #Getting Arrays for the Blades for Graphing
    X_ROTOR          = []      #Getting Arrays for the Rotors for Graphing
    Y_ROTOR          = []      #Getting Arrays for the Rotors for Graphing
    X_STATOR         = []      #Getting Arrays for the Stators for Graphing
    Y_STATOR         = []      #Getting Arrays for the Stators for Graphing
    X_HUB            = []      #Getting Arrays for the Hub and Tip for Graphing
    Y_HUB            = []      #Getting Arrays for the Hub and Tip for Graphing
    X_TIP            = []      #Getting Arrays for the Hub and Tip for Graphing
    Y_TIP            = []      #Getting Arrays for the Hub and Tip for Graphing
    X_CROSS_SECTION.extend(X_igv)
    Y_CROSS_SECTION.extend(Y_igv)
    X_BLADES.extend(X_igv_blades)
    Y_BLADES.extend(Y_igv_blades)
    X_STATOR.extend(X_igv_blades)
    Y_STATOR.extend(Y_igv_blades)
    X_HUB.extend(X_igv_h)
    Y_HUB.extend(Y_igv_h)
    X_TIP.extend(X_igv_t)
    Y_TIP.extend(Y_igv_t)
    NBS.append(SNB)
    
    '''FOR EXPORT - RV THETA'''
    RVTHI = []
    BLADEROWS = []
    BLADEROWS.append(SNB)
    
#####################################################
#####################################################
# PART B.v - FOR LOOP TO GO THROUGH EACH STAGE
#####################################################
#####################################################
    
    ''' GOING INTO THE LOOP TO GET PROPERTIES AND GEOMETRY FOR EVERY STAGE '''
    Pi = 3.1415926
    for i in Stage_Number : 
    #
        #
        STG = int(i) - 1
        print(' ')
        print(f'**** STAGE NUMBER {STG + 1} ****')
        
#####################################################
#####################################################
# PART B.v.1 - CALL ROTOR OUTLET FUNCTION
#####################################################
#####################################################
    
        '''Rotor/Stator Interface - Call Rotor_Outlet_Calculations Function'''
        To_2, Ts_2, Ps_2, Po_2, A_2, M_2, Radt_2, Radp_2, Rad_2, V_2, H_2, S_2, V_2p, Vz_2p, Vt_2p, ALPHA_2pd, BETAR_d, Ps_2h, Vz_2h, Vt_2h, ALPHA_2, BETA_2, Ts_2t, To_2t, Ps_2t, Po_2t, M_2t, V_2t, Vz_2t, Vt_2t, Rs_2t, ALPHA_2t, BETA_2t, Ts_2p, M_2p, Rs_2, Rs_2p, DF, h_2p, GammaValue2, Viscosity2, W_2p, Wz_2p, Wt_2p, W_2, Wz_2, Wt_2, W_2t, Wz_2t, Wt_2t, CCPP, Height_2, ROTS, BETA_2pd  = Rotor_Outlet_Calculations(n, STG, Stage_Number, DH, DTT, VZR, R_Phi, Lambda, Radius, Radp, RAV, Toi, Vti, CP, Vzi, Gamma, Rgas, PTR, Psi, TTR, MFR, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, typer, H, S, Po_h, Fluid, R_Sol, Wi, Wti, Wzi, VPT, Wi_p, Wti_p, Wzi_p, Vti_p, Vti_t, Vzi_p, Vzi_t, Radh, Radt, Toi_p, Toi_t, KKI, Rad_tip)
        
        #Add values to compressor export
        #Pitch Line
        RAD.append(round(Radp_2,4))
        PT.append(round(Po_2,2))
        PS.append(round(Ps_2,2))
        TT.append(round(To_2,2))
        TS.append(round(Ts_2p,2))
        HH.append(round(H_2,2))
        HS.append(round(h_2p,2))
        SS.append(round(S_2,2))
        RS.append(round(Rs_2p,2))
        VV.append(round(V_2p,2))
        VZV.append(round(Vz_2p,2))
        VTV.append(round(Vt_2p,2))
        WW.append(round(W_2p,2))
        WZW.append(round(Wz_2p,2))
        WTW.append(round(Wt_2p,2))
        ALP.append(round(ALPHA_2pd,2))
        BET.append(round(BETAR_d,2))
        MACHSP.append(round(M_2p, 3))
        DIFFUSION.append(round(DF, 4))
        GAMMAVALUE.append(round(GammaValue2, 5))
        VISCOSITY.append(round(Viscosity2*10**6, 5))
        SPECIFICHEAT.append(round(CCPP, 5))
        
        if STG == 0: 
            CAMBER.append(round(BETAR_d, 2) - round(BETAd, 2))
        elif STG > 0:
            CAMBER.append(round(BETAR_d, 2) - round(BETA_3d, 2))
        
        #Hub
        RADh.append(round(Rad_2,4))
        PTh.append(round(Po_2,2))
        PSh.append(round(Ps_2h,2))
        TTh.append(round(To_2,2))
        TSh.append(round(Ts_2,2))
        RSh.append(round(Rs_2,2))
        VVh.append(round(V_2,2))
        VZVh.append(round(Vz_2h,2))
        VTVh.append(round(Vt_2h,2))
        WWh.append(round(W_2,2))
        WZWh.append(round(Wz_2,2))
        WTWh.append(round(Wt_2,2))
        ALPh.append(round(ALPHA_2,2))
        BETh.append(round(BETA_2,2))
        MACHSH.append(round(M_2, 3))
        
        #Tip
        RADt.append(round(Radt_2,4))
        PTt.append(round(Po_2t,2))
        PSt.append(round(Ps_2t,2))
        TTt.append(round(To_2t,2))
        TSt.append(round(Ts_2t,2))
        RSt.append(round(Rs_2t,2))
        VVt.append(round(V_2t,2))
        VZVt.append(round(Vz_2t,2))
        VTVt.append(round(Vt_2t,2))
        WWt.append(round(W_2t,2))
        WZWt.append(round(Wz_2t,2))
        WTWt.append(round(Wt_2t,2))
        ALPt.append(round(ALPHA_2t,2))
        BETt.append(round(BETA_2t,2))
        MACHST.append(round(M_2t, 3))
        
#####################################################
#####################################################
# PART B.v.2 - ROTOR GEOMETRY
#####################################################
#####################################################
        
        '''ROTOR GEOMETRY [BASED ON CURRENT STAGE ROTOR, PREVIOUS STAGE STATOR]'''
        #Converting to Radians in order to use equations
        #Then converting back to degrees for Zetad
        
        '''CASE 1 -- STAGE 1, BASE GEOMETRY OFF INLET'''
        if STG == 0: 
            #HUB
            #if IGV_Q == 1: 
                #BETA_h = 0   #TC-DES HAS IT WRONG, BUT TO GET THE SAME VALUES AND COMPARE WE DO THIS SWAP
                
            GAMMA1_h_ROTOR = BETA_h*(Pi/180)                           #STAGGER ANGLE
            GAMMA2_h_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2*(Pi/180) - BETA_h*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_h_ROTOR   = (GAMMA1_h_ROTOR + GAMMA2_h_ROTOR) / 2     #ZETA
            ZETAd_h_ROTOR  = ZETA_h_ROTOR*(180/Pi)
                    
            #PITCH
            GAMMA1_p_ROTOR = BETAd*(Pi/180)                            #STAGGER ANGLE                             
            GAMMA2_p_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2pd*(Pi/180) - BETAd*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_p_ROTOR   = (GAMMA1_p_ROTOR + GAMMA2_p_ROTOR) / 2     #ZETA
            ZETAd_p_ROTOR  = ZETA_p_ROTOR*(180/Pi)
        
            #TIP
            GAMMA1_t_ROTOR = BETA_t*(Pi/180)                           #STAGGER ANGLE
            GAMMA2_t_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2t*(Pi/180) - BETA_t*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_t_ROTOR   = (GAMMA1_t_ROTOR + GAMMA2_t_ROTOR) / 2     #ZETA
            ZETAd_t_ROTOR  = ZETA_t_ROTOR*(180/Pi)
        
            #OVERALL
            WIDTH_ROTOR    = (Height_1 + Height_2) / 2 / R_AR[STG]     #AXIAL WIDTH
            CHORD_ROTOR    = WIDTH_ROTOR / math.cos(ZETA_p_ROTOR)      #CHORD LENGTH
        
            #REAL ASPECT RATIO
            REAL_AR_R      = (Height_1 + Height_2) / 2 / CHORD_ROTOR   #BLADE HEIGHT / CHORD
            WIDTH_h_ROTOR  = CHORD_ROTOR*math.cos(ZETA_h_ROTOR)        #WIDTH AT HUB
            WIDTH_p_ROTOR  = CHORD_ROTOR*math.cos(ZETA_p_ROTOR)        #WIDTH AT PITCH
            WIDTH_t_ROTOR  = CHORD_ROTOR*math.cos(ZETA_t_ROTOR)        #WIDTH AT TIP
            NUM_BLADES_R   = math.floor(Pi*(Radp + Radp_2) / (CHORD_ROTOR/R_Sol[STG]) +1) #Number of Blades in ROTOR
        
            #INDEX VALUES INTO GEOMETRY VECTORS
            STAGGERANGLE1_h.append(round(GAMMA1_h_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (HUB) RAD
            STAGGERANGLE1_p.append(round(GAMMA1_p_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (PITCH) RAD
            STAGGERANGLE1_t.append(round(GAMMA1_t_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (TIP) RAD 
            STAGGERANGLE2_h.append(round(GAMMA2_h_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (HUB) RAD
            STAGGERANGLE2_p.append(round(GAMMA2_p_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (PITCH) RAD
            STAGGERANGLE2_t.append(round(GAMMA2_t_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (TIP) RAD
            ZETAANGLE_h.append(round(ZETA_h_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (HUB) RAD
            ZETAANGLE_p.append(round(ZETA_p_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
            ZETAANGLE_t.append(round(ZETA_t_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
            ZETAANGLE_hD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (HUB) DEGREES
            ZETAANGLE_pD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
            ZETAANGLE_tD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
            AXIAL_WIDTH_R.append(round(WIDTH_ROTOR, 6)*1000)           #millimeters
            CHORD_LENGTH_R.append(round(CHORD_ROTOR, 6)*1000)          #millimeters
            RAR_R.append(round(REAL_AR_R, 6))                          #Real Aspect Ratio
            ROTORWIDTH_h.append(round(WIDTH_h_ROTOR, 6)*1000)          #Width of Rotor (HUB) (mm)
            ROTORWIDTH_p.append(round(WIDTH_p_ROTOR, 6)*1000)          #Width of Rotor (PITCH) (mm)
            ROTORWIDTH_t.append(round(WIDTH_t_ROTOR, 6)*1000)          #Width of Rotor (TIP) (mm)
            NBR.append(NUM_BLADES_R)                                   #NUMBER OF BLADES ON THE ROTOR
            BLADEROWS.append(NUM_BLADES_R)
            
            ''' Adds Geometry Info on Stage By Stage Basis '''
            FG = open(f'{Output}_Geometry.txt', "a") #Appends the Current File Created in Inlet_Calculations Function
            FG.write("**** STAGE NUMBER  : " + f'{STG + 1}' " ************************ \n")
            FG.write("**** ROTOR GEOMETRY **** \n")
            FG.write("Chord Length            : " + f'{round(CHORD_ROTOR, 6)*1000} mm' + "\n")
            FG.write("Axial Width             : " + f'{round(WIDTH_ROTOR, 6)*1000} mm' + "\n")
            FG.write("Number of Blades        : " + f'{NUM_BLADES_R}' + "\n")
            FG.write("Real Aspect Ratio       : " + f'{REAL_AR_R}' + "\n")
            FG.write("                                            \n") 
            FG.write("Hub Stagger Angle 1     : " + f'{round(GAMMA1_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Stagger Angle 2     : " + f'{round(GAMMA2_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Zeta Angle          : " + f'{round(ZETA_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Zeta Angle          : " + f'{round(ZETAd_h_ROTOR, 3)} deg' + "\n")
            FG.write("Hub Width               : " + f'{round(WIDTH_h_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n")  
            FG.write("Pitch Stagger Angle 1   : " + f'{round(GAMMA1_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Stagger Angle 2   : " + f'{round(GAMMA2_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Zeta Angle        : " + f'{round(ZETA_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Zeta Angle        : " + f'{round(ZETAd_p_ROTOR, 3)} deg' + "\n")
            FG.write("Pitch Width             : " + f'{round(WIDTH_p_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n") 
            FG.write("Tip Stagger Angle 1     : " + f'{round(GAMMA1_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Stagger Angle 2     : " + f'{round(GAMMA2_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Zeta Angle          : " + f'{round(ZETA_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Zeta Angle          : " + f'{round(ZETAd_t_ROTOR, 3)} deg' + "\n")
            FG.write("Tip Width               : " + f'{round(WIDTH_t_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n") 
            FG.close
        
        '''CASE 2 -- STAGE > 1, BASE GEOMETRY OFF LAST STAGES STATOR'''
        if STG > 0: 
            #HUB
            GAMMA1_h_ROTOR = BETA_h3*(Pi/180)                           #STAGGER ANGLE
            GAMMA2_h_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2*(Pi/180) - BETA_h3*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_h_ROTOR   = (GAMMA1_h_ROTOR + GAMMA2_h_ROTOR) / 2     #ZETA
            ZETAd_h_ROTOR  = ZETA_h_ROTOR*(180/Pi)
        
            #PITCH
            GAMMA1_p_ROTOR = BETA_3d*(Pi/180)                         #STAGGER ANGLE                             
            GAMMA2_p_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2pd*(Pi/180) - BETA_3d*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_p_ROTOR   = (GAMMA1_p_ROTOR + GAMMA2_p_ROTOR) / 2     #ZETA
            ZETAd_p_ROTOR  = ZETA_p_ROTOR*(180/Pi)
        
            #TIP
            GAMMA1_t_ROTOR = BETA_t3*(Pi/180)                          #STAGGER ANGLE
            GAMMA2_t_ROTOR = (4*(R_Sol[STG]**0.5)*BETA_2t*(Pi/180) - BETA_t3*(Pi/180)) / (4*(R_Sol[STG]**0.5) - 1)
            ZETA_t_ROTOR   = (GAMMA1_t_ROTOR + GAMMA2_t_ROTOR) / 2     #ZETA
            ZETAd_t_ROTOR  = ZETA_t_ROTOR*(180/Pi)
        
            #OVERALL
            WIDTH_ROTOR    = (Height_1 + Height_2) / 2 / R_AR[STG]     #AXIAL WIDTH
            CHORD_ROTOR    = WIDTH_ROTOR / math.cos(ZETA_p_ROTOR)      #CHORD LENGTH
        
            #REAL ASPECT RATIO
            REAL_AR_R      = (Height_1 + Height_2) / 2 / CHORD_ROTOR   #BLADE HEIGHT / CHORD
            WIDTH_h_ROTOR  = CHORD_ROTOR*math.cos(ZETA_h_ROTOR)        #WIDTH AT HUB
            WIDTH_p_ROTOR  = CHORD_ROTOR*math.cos(ZETA_p_ROTOR)        #WIDTH AT PITCH
            WIDTH_t_ROTOR  = CHORD_ROTOR*math.cos(ZETA_t_ROTOR)        #WIDTH AT TIP
            NUM_BLADES_R   = math.floor(Pi*(Radp + Radp_2) / (CHORD_ROTOR/R_Sol[STG]) +1) #Number of Blades in ROTOR
            
            #INDEX VALUES INTO GEOMETRY VECTORS
            STAGGERANGLE1_h.append(round(GAMMA1_h_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (HUB) RAD
            STAGGERANGLE1_p.append(round(GAMMA1_p_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (PITCH) RAD
            STAGGERANGLE1_t.append(round(GAMMA1_t_ROTOR, 3))           #GAMMA ANGLE AT ROTOR INLET (TIP) RAD 
            STAGGERANGLE2_h.append(round(GAMMA2_h_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (HUB) RAD
            STAGGERANGLE2_p.append(round(GAMMA2_p_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (PITCH) RAD
            STAGGERANGLE2_t.append(round(GAMMA2_t_ROTOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (TIP) RAD
            ZETAANGLE_h.append(round(ZETA_h_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (HUB) RAD
            ZETAANGLE_p.append(round(ZETA_p_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
            ZETAANGLE_t.append(round(ZETA_t_ROTOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
            ZETAANGLE_hD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (HUB) DEGREES
            ZETAANGLE_pD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
            ZETAANGLE_tD.append(round(ZETAd_h_ROTOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
            AXIAL_WIDTH_R.append(round(WIDTH_ROTOR, 6)*1000)           #millimeters
            CHORD_LENGTH_R.append(round(CHORD_ROTOR, 6)*1000)          #millimeters
            RAR_R.append(round(REAL_AR_R, 6))                          #Real Aspect Ratio
            ROTORWIDTH_h.append(round(WIDTH_h_ROTOR, 6)*1000)          #Width of Rotor (HUB) (mm)
            ROTORWIDTH_p.append(round(WIDTH_p_ROTOR, 6)*1000)          #Width of Rotor (PITCH) (mm)
            ROTORWIDTH_t.append(round(WIDTH_t_ROTOR, 6)*1000)          #Width of Rotor (TIP) (mm)
            NBR.append(NUM_BLADES_R)                                   #NUMBER OF BLADES ON THE ROTOR
            BLADEROWS.append(NUM_BLADES_R)
        
            ''' Adds Geometry Info on Stage By Stage Basis '''
            FG = open(f'{Output}_Geometry.txt', "a") #Appends the Current File Created in Inlet_Calculations Function
            FG.write("**** STAGE NUMBER  : " + f'{STG + 1}' " ************************ \n")
            FG.write("**** ROTOR GEOMETRY **** \n")
            FG.write("Chord Length            : " + f'{round(CHORD_ROTOR, 6)*1000} mm' + "\n")
            FG.write("Axial Width             : " + f'{round(WIDTH_ROTOR, 6)*1000} mm' + "\n")
            FG.write("Number of Blades        : " + f'{NUM_BLADES_R}' + "\n")
            FG.write("Real Aspect Ratio       : " + f'{REAL_AR_R}' + "\n")
            FG.write("                                            \n") 
            FG.write("Hub Stagger Angle 1     : " + f'{round(GAMMA1_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Stagger Angle 2     : " + f'{round(GAMMA2_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Zeta Angle          : " + f'{round(ZETA_h_ROTOR, 3)} rad' + "\n")
            FG.write("Hub Zeta Angle          : " + f'{round(ZETAd_h_ROTOR, 3)} deg' + "\n")
            FG.write("Hub Width               : " + f'{round(WIDTH_h_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n")  
            FG.write("Pitch Stagger Angle 1   : " + f'{round(GAMMA1_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Stagger Angle 2   : " + f'{round(GAMMA2_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Zeta Angle        : " + f'{round(ZETA_p_ROTOR, 3)} rad' + "\n")
            FG.write("Pitch Zeta Angle        : " + f'{round(ZETAd_p_ROTOR, 3)} deg' + "\n")
            FG.write("Pitch Width             : " + f'{round(WIDTH_p_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n") 
            FG.write("Tip Stagger Angle 1     : " + f'{round(GAMMA1_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Stagger Angle 2     : " + f'{round(GAMMA2_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Zeta Angle          : " + f'{round(ZETA_t_ROTOR, 3)} rad' + "\n")
            FG.write("Tip Zeta Angle          : " + f'{round(ZETAd_t_ROTOR, 3)} deg' + "\n")
            FG.write("Tip Width               : " + f'{round(WIDTH_t_ROTOR, 6)*1000} mm' + "\n")
            FG.write("                                            \n") 
            FG.close
    
#####################################################
#####################################################
# PART B.v.3 - CALL STATOR OUTLET FUNCTION
#####################################################
#####################################################            
    
        '''Stator Outlet = Next Stage's Rotor Inlet'''
    
        '''Stator Outlet - Call Stator_Outlet_Calculations Function'''
        Radius, Radp, Toi, Poi, Tsi, Psi, Vti, Vzi, PTR, TTR, MFR, H, S, PR, Po_h, Rsi, V_3, Vz_3, Vt_3, ALPHA_3d, BETA_3d, Rad_h3, Ts_h3, Ps_h3, Rs_h3, To_h3, Po_h3, M_h3, V_h3, Vz_h3, Vt_h3, ALPHA_h3, BETA_h3, Radt_3, Ts_t3, Ps_t3, Rs_t3, To_t3, Po_t3, M_t3, V_t3, Vz_t3, Vt_t3, ALPHA_t3, BETA_t3, Rs_3p, Ts_3p, M_3p, UU, Wi, Wti, Wzi, DF, h_3, GammaValue, Viscosity, Wi_p, Wti_p, Wzi_p, W_t3, Wt_t3, Wz_t3, CCPP, Height_3 = Stator_Outlet_Calculations(n, alphaend, Mend, STG, alpha1, MACH1, Radius, To_2, H_2, Gamma, S_Phi, Po_2, Ps_2, Rgas, MFR, Lambda, RAV, CP, Bleed, Poi, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, Fluid, V_2p, Vt_2p, S_Sol, VPT, typer, KKI, Rad_tip)
        
        #Add values to compressor export
        #Pitch Line
        RAD.append(round(Radp,4))
        PT.append(round(Poi,2))
        PS.append(round(Psi,2))
        TT.append(round(Toi,2))
        TS.append(round(Ts_3p,2))
        HH.append(round(H,2))
        HS.append(round(h_3,2))
        SS.append(round(S,2))
        RS.append(round(Rs_3p,2))
        VV.append(round(V_3,2))
        VZV.append(round(Vz_3,2))
        VTV.append(round(Vt_3,2))
        WW.append(round(Wi_p,2))
        WZW.append(round(Wzi_p,2))
        WTW.append(round(Wti_p,2))
        ALP.append(round(ALPHA_3d,2))
        BET.append(round(BETA_3d,2))
        MACHSP.append(round(M_3p, 3))
        CAMBER.append(round(ALPHA_3d, 2) - round(ALPHA_2pd, 2))
        DIFFUSION.append(round(DF, 4))
        GAMMAVALUE.append(round(GammaValue, 5))
        VISCOSITY.append(round(Viscosity*10**6, 5))
        SPECIFICHEAT.append(round(CCPP, 5))
        
        '''Efficiency Per Stage'''
        P_isentropic = Poi
        H_isentropic = HSP(S_isentropic, P_isentropic)
        AD_EFF_perstage.append((H_isentropic - H_initial_2) / (H - H_initial_2)*100)
        '''Update for next loop'''
        H_initial_2 = H
        S_isentropic = S
        
        #Hub
        RADh.append(round(Rad_h3,4))
        PTh.append(round(Po_h3,2))
        PSh.append(round(Ps_h3,2))
        TTh.append(round(To_h3,2))
        TSh.append(round(Ts_h3,2))
        RSh.append(round(Rs_h3,2))
        VVh.append(round(V_h3,2))
        VZVh.append(round(Vz_h3,2))
        VTVh.append(round(Vt_h3,2))
        WWh.append(round(Wi,2))
        WZWh.append(round(Wzi,2))
        WTWh.append(round(Wti,2))
        ALPh.append(round(ALPHA_h3,2))
        BETh.append(round(BETA_h3,2))
        MACHSH.append(round(M_h3, 3))
        
        #Tip
        RADt.append(round(Radt_3,4))
        PTt.append(round(Po_t3,2))
        PSt.append(round(Ps_t3,2))
        TTt.append(round(To_t3,2))
        TSt.append(round(Ts_t3,2))
        RSt.append(round(Rs_t3,2))
        VVt.append(round(V_t3,2))
        VZVt.append(round(Vz_t3,2))
        VTVt.append(round(Vt_t3,2))
        WWt.append(round(W_t3,2))
        WZWt.append(round(Wz_t3,2))
        WTWt.append(round(Wt_t3,2))
        ALPt.append(round(ALPHA_t3,2))
        BETt.append(round(BETA_t3,2))
        MACHST.append(round(M_t3, 3))
        
        Current_Pressure = Current_Pressure * PR
        Current_Temperature = Toi
        
        '''FOR PLOTS'''
        Pplot.append(Poi)
        Tplot.append(Toi)
        Rplot.append(Rsi)
        Hplot.append(H)
        Splot.append(S)
        LOADING_C.append(DH[STG] * 1000 / UU**2)
        FLOW_C.append(Vz_3 / UU) 
        Value_Reaction.append((HS[1 + 2*STG] - HS[0 + 2*STG]) / (HS[2 + 2*STG] - HS[0 + 2*STG]))
        
#####################################################
#####################################################
# PART B.v.4 - STATOR GEOMETRY
#####################################################
#####################################################
        
        '''STATOR GEOMETRY'''
        #Converting to Radians in order to use equations
        #Then converting back to degrees for Zetad
        
        #HUB
        GAMMA1_h_STATOR = ALPHA_2*(Pi/180)                           #STAGGER ANGLE
        GAMMA2_h_STATOR = (4*(S_Sol[STG]**0.5)*ALPHA_h3*(Pi/180) - ALPHA_2*(Pi/180)) / (4*(S_Sol[STG]**0.5) - 1)
        ZETA_h_STATOR   = (GAMMA1_h_STATOR + GAMMA2_h_STATOR) / 2     #ZETA
        ZETAd_h_STATOR  = ZETA_h_STATOR*(180/Pi)
        
        #PITCH
        GAMMA1_p_STATOR = ALPHA_2pd*(Pi/180)                           #STAGGER ANGLE
        GAMMA2_p_STATOR = (4*(S_Sol[STG]**0.5)*ALPHA_3d*(Pi/180) - ALPHA_2pd*(Pi/180)) / (4*(S_Sol[STG]**0.5) - 1)
        ZETA_p_STATOR   = (GAMMA1_p_STATOR + GAMMA2_p_STATOR) / 2     #ZETA
        ZETAd_p_STATOR  = ZETA_p_STATOR*(180/Pi)
        
        #TIP
        GAMMA1_t_STATOR = ALPHA_2t*(Pi/180)                           #STAGGER ANGLE
        GAMMA2_t_STATOR = (4*(S_Sol[STG]**0.5)*ALPHA_t3*(Pi/180) - ALPHA_2t*(Pi/180)) / (4*(S_Sol[STG]**0.5) - 1)
        ZETA_t_STATOR   = (GAMMA1_t_STATOR + GAMMA2_t_STATOR) / 2     #ZETA
        ZETAd_t_STATOR  = ZETA_t_STATOR*(180/Pi)
        
        #OVERALL
        WIDTH_STATOR    = (Height_2 + Height_3) / 2 / S_AR[STG]     #AXIAL WIDTH
        CHORD_STATOR    = WIDTH_STATOR / math.cos(ZETA_p_STATOR)      #CHORD LENGTH
        
        #REAL ASPECT RATIO
        REAL_AR_S       = (Height_2 + Height_3) / 2 / CHORD_STATOR   #BLADE HEIGHT / CHORD
        WIDTH_h_STATOR  = CHORD_STATOR*math.cos(ZETA_h_STATOR)        #WIDTH AT HUB
        WIDTH_p_STATOR  = CHORD_STATOR*math.cos(ZETA_p_STATOR)        #WIDTH AT PITCH
        WIDTH_t_STATOR  = CHORD_STATOR*math.cos(ZETA_t_STATOR)        #WIDTH AT TIP
        NUM_BLADES_S    = math.floor(Pi*(Radp_2 + Radp) / (CHORD_STATOR/S_Sol[STG]) + 1) #Number of Blades in ROTOR
        
        #INDEX VALUES INTO GEOMETRY VECTORS
        STAGGERANGLE1_Sh.append(round(GAMMA1_h_STATOR, 3))           #GAMMA ANGLE AT ROTOR INLET (HUB) RAD
        STAGGERANGLE1_Sp.append(round(GAMMA1_p_STATOR, 3))           #GAMMA ANGLE AT ROTOR INLET (PITCH) RAD
        STAGGERANGLE1_St.append(round(GAMMA1_t_STATOR, 3))           #GAMMA ANGLE AT ROTOR INLET (TIP) RAD 
        STAGGERANGLE2_Sh.append(round(GAMMA2_h_STATOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (HUB) RAD
        STAGGERANGLE2_Sp.append(round(GAMMA2_p_STATOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (PITCH) RAD
        STAGGERANGLE2_St.append(round(GAMMA2_t_STATOR, 3))           #GAMMA ANGLE AT ROTOR OUTLET (TIP) RAD
        ZETAANGLE_Sh.append(round(ZETA_h_STATOR, 3))                 #ZETA ANGLE AT ROTOR (HUB) RAD
        ZETAANGLE_Sp.append(round(ZETA_p_STATOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
        ZETAANGLE_St.append(round(ZETA_t_STATOR, 3))                 #ZETA ANGLE AT ROTOR (PITCH) RAD
        ZETAANGLE_ShD.append(round(ZETAd_h_STATOR, 3))               #ZETA ANGLE AT ROTOR (HUB) DEGREES
        ZETAANGLE_SpD.append(round(ZETAd_p_STATOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
        ZETAANGLE_StD.append(round(ZETAd_t_STATOR, 3))               #ZETA ANGLE AT ROTOR (PITCH) DEGREES
        AXIAL_WIDTH_S.append(round(WIDTH_STATOR, 6)*1000)            #millimeters
        CHORD_LENGTH_S.append(round(CHORD_STATOR, 6)*1000)           #millimeters
        RAR_S.append(round(REAL_AR_S, 6))                            #Real Aspect Ratio
        STATORWIDTH_h.append(round(WIDTH_h_STATOR, 6)*1000)          #Width of Rotor (HUB) (mm)
        STATORWIDTH_p.append(round(WIDTH_p_STATOR, 6)*1000)          #Width of Rotor (PITCH) (mm)
        STATORWIDTH_t.append(round(WIDTH_t_STATOR, 6)*1000)          #Width of Rotor (TIP) (mm)
        NBS.append(NUM_BLADES_S)                                     #NUMBER OF BLADES ON THE ROTOR
        BLADEROWS.append(NUM_BLADES_S)
        
        ''' Adds Geometry Info on Stage By Stage Basis '''
        FG = open(f'{Output}_Geometry.txt', "a") #Appends the Current File Created in Inlet_Calculations Function
        FG.write("**** STATOR GEOMETRY **** \n")
        FG.write("Chord Length            : " + f'{round(CHORD_STATOR, 6)*1000} mm' + "\n")
        FG.write("Axial Width             : " + f'{round(WIDTH_STATOR, 6)*1000} mm' + "\n")
        FG.write("Number of Blades        : " + f'{NUM_BLADES_S}' + "\n")
        FG.write("Real Aspect Ratio       : " + f'{REAL_AR_S}' + "\n")
        FG.write("                                            \n") 
        FG.write("Hub Stagger Angle 1     : " + f'{round(GAMMA1_h_STATOR, 3)} rad' + "\n")
        FG.write("Hub Stagger Angle 2     : " + f'{round(GAMMA2_h_STATOR, 3)} rad' + "\n")
        FG.write("Hub Zeta Angle          : " + f'{round(ZETA_h_STATOR, 3)} rad' + "\n")
        FG.write("Hub Zeta Angle          : " + f'{round(ZETAd_h_STATOR, 3)} deg' + "\n")
        FG.write("Hub Width               : " + f'{round(WIDTH_h_STATOR, 6)*1000} mm' + "\n")
        FG.write("                                            \n")  
        FG.write("Pitch Stagger Angle 1   : " + f'{round(GAMMA1_p_STATOR, 3)} rad' + "\n")
        FG.write("Pitch Stagger Angle 2   : " + f'{round(GAMMA2_p_STATOR, 3)} rad' + "\n")
        FG.write("Pitch Zeta Angle        : " + f'{round(ZETA_p_STATOR, 3)} rad' + "\n")
        FG.write("Pitch Zeta Angle        : " + f'{round(ZETAd_p_STATOR, 3)} deg' + "\n")
        FG.write("Pitch Width             : " + f'{round(WIDTH_p_STATOR, 6)*1000} mm' + "\n")
        FG.write("                                            \n") 
        FG.write("Tip Stagger Angle 1     : " + f'{round(GAMMA1_t_STATOR, 3)} rad' + "\n")
        FG.write("Tip Stagger Angle 2     : " + f'{round(GAMMA2_t_STATOR, 3)} rad' + "\n")
        FG.write("Tip Zeta Angle          : " + f'{round(ZETA_t_STATOR, 3)} rad' + "\n")
        FG.write("Tip Zeta Angle          : " + f'{round(ZETAd_t_STATOR, 3)} deg' + "\n")
        FG.write("Tip Width               : " + f'{round(WIDTH_t_STATOR, 6)*1000} mm' + "\n")
        FG.write("                                            \n") 

#####################################################
#####################################################
# PART B.v.5 - CROSS SECTION GEOMETRY
#####################################################
#####################################################
#X_RS AND Y_RS FOR WALL AND STACK FUNCTIONS LATER
#X_BLADES, Y_BLADES, X_HUB, Y_HUB, X_TIP, Y_TIP...> 
#FOR ROUGH CROSS SECTION PLOTTING
        
        ''' ROTOR AND STATOR CROSS SECTION GEOMETRY '''
        #COORDINATES FOR EACH STAGE
        X_RS = []      #X Coords
        Y_RS = []      #Y Coords
        
        #Got the per stage info above, but need to append the 
        #hub and tip arrays for each coordinate too
        #Also appending the blade coordinates to accurately map where the
        #Blades are

        X_RS.append(XRIF)
        Y_RS.append(round(Radh, 9)*1000)
        X_HUB.append(XRIF)
        Y_HUB.append(round(Radh, 9)*1000)
        X_BLADES.append(XRIF)
        Y_BLADES.append(round(Radh, 9)*1000)
        X_ROTOR.append(XRIF)
        Y_ROTOR.append(round(Radh, 9)*1000)
        
        X_RS.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR - WIDTH_t_ROTOR)/2, 9)*1000)        #{mm}
        Y_RS.append(round(Radt, 9)*1000)
        X_TIP.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR - WIDTH_t_ROTOR)/2, 9)*1000) 
        Y_TIP.append(round(Radt, 9)*1000)
        X_BLADES.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR - WIDTH_t_ROTOR)/2, 9)*1000)
        Y_BLADES.append(round(Radt, 9)*1000)
        X_ROTOR.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR - WIDTH_t_ROTOR)/2, 9)*1000)
        Y_ROTOR.append(round(Radt, 9)*1000)
        
        X_RS.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR + WIDTH_t_ROTOR)/2, 9)*1000)
        Y_RS.append(round(Radt, 9)*1000)
        X_TIP.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR + WIDTH_t_ROTOR)/2, 9)*1000)
        Y_TIP.append(round(Radt, 9)*1000)
        X_BLADES.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR + WIDTH_t_ROTOR)/2, 9)*1000)
        Y_BLADES.append(round(Radt, 9)*1000)
        X_ROTOR.append(round(X_RS[0]/1000 + (WIDTH_h_ROTOR + WIDTH_t_ROTOR)/2, 9)*1000)
        Y_ROTOR.append(round(Radt, 9)*1000)
        
        X_RS.append(round(X_RS[0]/1000 + WIDTH_h_ROTOR, 9)*1000)
        Y_RS.append(round(Radh, 9)*1000)
        X_HUB.append(round(X_RS[0]/1000 + WIDTH_h_ROTOR, 9)*1000)
        Y_HUB.append(round(Radh, 9)*1000)
        X_BLADES.append(round(X_RS[0]/1000 + WIDTH_h_ROTOR, 9)*1000)
        Y_BLADES.append(round(Radh, 9)*1000)
        X_ROTOR.append(round(X_RS[0]/1000 + WIDTH_h_ROTOR, 9)*1000)
        Y_ROTOR.append(round(Radh, 9)*1000)
        
        X_RS.append(round(X_RS[0]/1000 + (3*WIDTH_h_ROTOR + WIDTH_t_ROTOR)/4 + R_Coef_Space[STG]*WIDTH_p_ROTOR, 9)*1000)
        Y_RS.append(round(Radt_2, 9)*1000)
        X_TIP.append(round(X_RS[0]/1000 + (3*WIDTH_h_ROTOR + WIDTH_t_ROTOR)/4 + R_Coef_Space[STG]*WIDTH_p_ROTOR, 9)*1000)
        Y_TIP.append(round(Radt_2, 9)*1000)
        
        X_RS.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR - WIDTH_h_STATOR)/2, 9)*1000)
        Y_RS.append(round(Rad_2, 9)*1000)
        X_HUB.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR - WIDTH_h_STATOR)/2, 9)*1000)
        Y_HUB.append(round(Rad_2, 9)*1000)
        X_BLADES.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR - WIDTH_h_STATOR)/2, 9)*1000)
        Y_BLADES.append(round(Rad_2, 9)*1000)
        X_BLADES.append(round(X_RS[0]/1000 + (3*WIDTH_h_ROTOR + WIDTH_t_ROTOR)/4 + R_Coef_Space[STG]*WIDTH_p_ROTOR, 9)*1000)
        Y_BLADES.append(round(Radt_2, 9)*1000)
        X_STATOR.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR - WIDTH_h_STATOR)/2, 9)*1000)
        Y_STATOR.append(round(Rad_2, 9)*1000)
        X_STATOR.append(round(X_RS[0]/1000 + (3*WIDTH_h_ROTOR + WIDTH_t_ROTOR)/4 + R_Coef_Space[STG]*WIDTH_p_ROTOR, 9)*1000)
        Y_STATOR.append(round(Radt_2, 9)*1000)
        
        X_RS.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR + WIDTH_h_STATOR)/2, 9)*1000)
        Y_RS.append(round(Rad_2, 9)*1000)
        X_HUB.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR + WIDTH_h_STATOR)/2, 9)*1000)
        Y_HUB.append(round(Rad_2, 9)*1000)
        
        X_RS.append(round(X_RS[4]/1000 + WIDTH_t_STATOR, 9)*1000)
        Y_RS.append(round(Radt_2, 9)*1000)
        X_TIP.append(round(X_RS[4]/1000 + WIDTH_t_STATOR, 9)*1000)
        Y_TIP.append(round(Radt_2, 9)*1000)
        X_BLADES.append(round(X_RS[4]/1000 + WIDTH_t_STATOR, 9)*1000)
        Y_BLADES.append(round(Radt_2, 9)*1000)
        X_BLADES.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR + WIDTH_h_STATOR)/2, 9)*1000)
        Y_BLADES.append(round(Rad_2, 9)*1000)
        X_STATOR.append(round(X_RS[4]/1000 + WIDTH_t_STATOR, 9)*1000)
        Y_STATOR.append(round(Radt_2, 9)*1000)
        X_STATOR.append(round(X_RS[4]/1000 + (WIDTH_t_STATOR + WIDTH_h_STATOR)/2, 9)*1000)
        Y_STATOR.append(round(Rad_2, 9)*1000)
        
        X_RS.append(round(X_RS[4]/1000 + (3*WIDTH_t_STATOR + WIDTH_h_STATOR)/4 + S_Coef_Space[STG]*WIDTH_p_STATOR, 9)*1000)
        #X_RS.append(round(X_RS[6]/1000 + S_Coef_Space[STG]*WIDTH_p_STATOR, 9)*1000)
        Y_RS.append(round(Radt_3, 9)*1000)
        X_TIP.append(round(X_RS[4]/1000 + (3*WIDTH_t_STATOR + WIDTH_h_STATOR)/4 + S_Coef_Space[STG]*WIDTH_p_STATOR, 9)*1000)
        #X_TIP.append(round(X_RS[6]/1000 + S_Coef_Space[STG]*WIDTH_p_STATOR, 9)*1000)
        Y_TIP.append(round(Radt_3, 9)*1000)
        
        X_RS.append(round(X_RS[8]/1000, 9)*1000)
        Y_RS.append(round(Rad_h3, 9)*1000)
        X_HUB.append(round(X_RS[8]/1000, 9)*1000)
        Y_HUB.append(round(Rad_h3, 9)*1000)
        # X_BLADES.append(round(X_RS[8]/1000, 9)*1000)
        # Y_BLADES.append(round(Rad_h3, 9)*1000)
        
        XRIF = X_RS[8]
        
        #Add Coordinates to Entire Cross Section Vectors
        X_CROSS_SECTION.extend(X_RS)
        Y_CROSS_SECTION.extend(Y_RS)

        
        FG.write("**** ROTOR & STATOR CROSS SECTION GEOMETRY **** \n")
        FG.write("X COORDINATES      [mm] : " + f'{X_RS}' + "\n")
        FG.write("Y COORDINATES      [mm] : " + f'{Y_RS}' + "\n")
        FG.write("                                            \n") 
        FG.close
       
#####################################################
#####################################################
# PART B.v.6 - EXPORT BLADE INFORMATION TO .TINF FILE
#####################################################
##################################################### 
       
        '''EXPORTING COMPRESSOR BLADES PER ROW FILE'''
        if IGV_Q == 1: 
        #
            #
            if STG == 0:
            #
                #   
                FW = open(f'{Output}.tinf', "w+")
                FW.write(" BLADE_ROW  NB       R-V_THETA [m^2/sec] \n")
                R_vtheta1 = round(Vti_h * Radh, 6)
                FW.write(f'    {STG+1}       {NBS[STG]}         {R_vtheta1} \n')
                RVTHI.append(R_vtheta1)
                FW.close
                #
            #
            FW = open(f'{Output}.tinf', "a")
            R_vtheta2 = round(Rad_2 * Vt_2h, 6)
            FW.write(f'    {2*(STG+1)}       {NBR[STG]}         {R_vtheta2} \n')
            R_vtheta3 = round(Rad_h3 * Vt_h3, 6)
            FW.write(f'    {2*(STG+1)+1}       {NBS[STG+1]}         {R_vtheta3} \n')
            RVTHI.append(R_vtheta2)
            RVTHI.append(R_vtheta3)
            #
        #
        if IGV_Q == 0: 
        #
            #
            if STG == 0: 
            #
                # 
                FW = open(f'{Output}.tinf', "w+")
                FW.write(" BLADE_ROW  NB       R-V_THETA [m^2/sec] \n")
                R_vtheta1 = round(Vti_h * Radh, 4)
                FW.write(f'    {STG}       {0}         {R_vtheta1} \n')
                RVTHI.append(R_vtheta1)
                FW.close
                #
            #
            FW = open(f'{Output}.tinf', "a")
            R_vtheta2 = round(Rad_2 * Vt_2h, 4)
            FW.write(f'    {2*(STG)+1}       {NBR[STG]}         {R_vtheta2} \n')
            R_vtheta3 = round(Rad_h3 * Vt_h3, 4)
            FW.write(f'    {2*(STG)+2}       {NBS[STG+1]}         {R_vtheta3} \n')
            RVTHI.append(R_vtheta2)
            RVTHI.append(R_vtheta3)
            #
        #
        FW.close
        #
        #
        '''Update Stator Radii and Height (Height_3) into Height_1'''
        Height_1 = Height_3
        Radh     = Rad_h3
        Radt     = Radt_3
        Toi_p    = Toi
        Toi_t    = To_t3
        Vzi_p    = Vz_3
        Vzi_t    = Vz_t3
        Vti_p    = Vt_3
        Vti_t    = Vt_t3
        
        if STG >= (n - 1):
            break
            #
        #
    #
    j = 0
    DHH.append(0)
    for b in HH :
        DHH.append(round(HH[int(j)+1] - HH[int(j)],2))
        j += 1
        if j >= len(HH)-1:
            break    

#####################################################
#####################################################
# PART B.vi - FULL GEOMETRY VECTOR EXPORT
#####################################################
#####################################################
    
    '''EXPORTING COMPRESSOR GEOMETRY VECTORS'''
    FG = open(f'{Output}_GeometryVectors.txt', "w+")
    FG.write("**** GEOMETRY FILE, VECTOR DATA************* \n")
    FG.write("**** ROTOR GEOMETRY **** \n")
    FG.write("AXIAL WIDTH              [mm]  :" +f'{AXIAL_WIDTH_R}' + "\n")
    FG.write("CHORD LENGTH             [mm]  :" +f'{CHORD_LENGTH_R}' + "\n")
    FG.write("Number of Blades               :" +f'{NBR}' + "\n")
    FG.write("Real Aspect Ratio              :" +f'{RAR_R}' + "\n")
    FG.write("X Coordinates (Cross Sec)[mm]  :" +f'{X_CROSS_SECTION}' + "\n")
    FG.write("Y Coordinates (Cross Sec)[mm]  :" +f'{Y_CROSS_SECTION}' + "\n")
    FG.write("                                            \n") 
    FG.write("** HUB ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_h}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_h}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_h}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_hD}' + "\n")
    FG.write("Hub Width               [mm]   :" +f'{ROTORWIDTH_h}' + "\n")
    FG.write("                                            \n") 
    FG.write("** PITCH ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_p}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_p}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_p}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_pD}' + "\n")
    FG.write("Pitch Width             [mm]   :" +f'{ROTORWIDTH_p}' + "\n")
    FG.write("                                            \n") 
    FG.write("** TIP ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_t}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_t}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_t}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_tD}' + "\n")
    FG.write("Tip Width               [mm]   :" +f'{ROTORWIDTH_t}' + "\n")
    FG.write("                                            \n") 
    FG.write("**** STATOR GEOMETRY **** \n")
    FG.write("AXIAL WIDTH              [mm]  :" +f'{AXIAL_WIDTH_S}' + "\n")
    FG.write("CHORD LENGTH             [mm]  :" +f'{CHORD_LENGTH_S}' + "\n")
    FG.write("Number of Blades               :" +f'{NBS}' + "\n")
    FG.write("Real Aspect Ratio              :" +f'{RAR_S}' + "\n")
    FG.write("                                            \n") 
    FG.write("** HUB ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_Sh}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_Sh}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_Sh}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_ShD}' + "\n")
    FG.write("Hub Width               [mm]   :" +f'{STATORWIDTH_h}' + "\n")
    FG.write("                                            \n") 
    FG.write("** PITCH ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_Sp}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_Sp}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_Sp}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_SpD}' + "\n")
    FG.write("Pitch Width             [mm]   :" +f'{STATORWIDTH_p}' + "\n")
    FG.write("                                            \n") 
    FG.write("** TIP ** \n")
    FG.write("Stagger Angle 1         [rad]  :" +f'{STAGGERANGLE1_St}' + "\n")
    FG.write("Stagger Angle 2         [rad]  :" +f'{STAGGERANGLE1_St}' + "\n")
    FG.write("Zeta Angle              [rad]  :" +f'{ZETAANGLE_St}' + "\n")
    FG.write("Zeta Angle              [deg]  :" +f'{ZETAANGLE_StD}' + "\n")
    FG.write("Tip Width               [mm]   :" +f'{STATORWIDTH_t}' + "\n")
    FG.write("                                            \n")
    FG.close
        
#####################################################
#####################################################
# PART B.vii - FULL FLUID AND COMPRESSOR PROPERTY EXPORT
#####################################################
#####################################################
    
    '''EXPORTING COMPRESSOR AND FLUID PROPERTIES'''
    FW = open(f'{Output}.txt', "w+")
    FW.write("Enthalpy           [J/kg]   :" + f'{HH}' + "\n")
    FW.write("Delta Enthalpy     [J/kg]   :" + f'{DHH}' + "\n")
    FW.write("Entropy            [J/kgK]  :" + f'{SS}' + "\n")
    FW.write("Loading Coefficient [ ]     :" + f'{LOADING_C}' + "\n")
    FW.write("Flow Coefficient    [ ]     :" + f'{FLOW_C}' + "\n")
    FW.write("Diffusion Factor    [ ]     :" + f'{DIFFUSION}' + "\n")
    FW.write("Per Stage Efficiency [%]    :" + f'{AD_EFF_perstage}' + "\n")
    FW.write("                                            \n")
    FW.write("Pitch-Line Calculations \n")
    FW.write("Stage Values of Reaction    :" + f'{Value_Reaction}' + "\n")
    FW.write("Gamma Value         []      :" + f'{GAMMAVALUE}' + "\n")
    FW.write("Viscosity          [uPa/s]  :" + f'{VISCOSITY}' + "\n")
    FW.write("Pitch Line Radius   [m]     :" + f'{RAD}' + "\n")
    FW.write("Pitch Line Camber  [deg]    :" + f'{CAMBER}' + "\n")
    FW.write("Total Pressure     [kPa]    :" + f'{PT}' + "\n")
    FW.write("Static Pressure    [kPa]    :" + f'{PS}' + "\n")
    FW.write("Temperature        [K]      :" + f'{TT}' + "\n")
    FW.write("Static Temperature [K]      :" + f'{TS}' + "\n")
    FW.write("Static Density    [kg/m3]   :" + f'{RS}' + "\n")
    FW.write("Mach #             []       :" + f'{MACHSP}' + "\n")
    FW.write("Total Velocity     [m/s]    :" + f'{VV}' + "\n")
    FW.write("Velocity Axial     [m/s]    :" + f'{VZV}' + "\n")
    FW.write("Velocity Theta     [m/s]    :" + f'{VTV}' + "\n")
    FW.write("Relative Velocity  [m/s]    :" + f'{WW}' + "\n")
    FW.write("Relative Axial     [m/s]    :" + f'{WZW}' + "\n")
    FW.write("Relative Theta     [m/s]    :" + f'{WTW}' + "\n")
    FW.write("ALPHA              [deg]    :" + f'{ALP}' + "\n")
    FW.write("BETA               [deg]    :" + f'{BET}' + "\n")
    FW.write("                                            \n")
    FW.write("Hub Calculations \n")
    FW.write("Radius              [m]     :" + f'{RADh}' + "\n")
    FW.write("Total Pressure     [kPa]    :" + f'{PTh}' + "\n")
    FW.write("Static Pressure    [kPa]    :" + f'{PSh}' + "\n")
    FW.write("Temperature        [K]      :" + f'{TTh}' + "\n")
    FW.write("Static Temperature [K]      :" + f'{TSh}' + "\n")
    FW.write("Static Density    [kg/m3]   :" + f'{RSh}' + "\n")
    FW.write("Mach #             []       :" + f'{MACHSH}' + "\n")
    FW.write("Total Velocity     [m/s]    :" + f'{VVh}' + "\n")
    FW.write("Velocity Axial     [m/s]    :" + f'{VZVh}' + "\n")
    FW.write("Velocity Theta     [m/s]    :" + f'{VTVh}' + "\n")
    FW.write("Relative Velocity  [m/s]    :" + f'{WWh}' + "\n")
    FW.write("Relative Axial     [m/s]    :" + f'{WZWh}' + "\n")
    FW.write("Relative Theta     [m/s]    :" + f'{WTWh}' + "\n")
    FW.write("ALPHA              [deg]    :" + f'{ALPh}' + "\n")
    FW.write("BETA               [deg]    :" + f'{BETh}' + "\n")
    FW.write("                                            \n")
    FW.write("Tip Calculations \n")
    FW.write("Radius              [m]     :" + f'{RADt}' + "\n")
    FW.write("Total Pressure     [kPa]    :" + f'{PTt}' + "\n")
    FW.write("Static Pressure    [kPa]    :" + f'{PSt}' + "\n")
    FW.write("Temperature        [K]      :" + f'{TTt}' + "\n")
    FW.write("Static Temperature [K]      :" + f'{TSt}' + "\n")
    FW.write("Static Density    [kg/m3]   :" + f'{RSt}' + "\n")
    FW.write("Mach #             []       :" + f'{MACHST}' + "\n")
    FW.write("Total Velocity     [m/s]    :" + f'{VVt}' + "\n")
    FW.write("Velocity Axial     [m/s]    :" + f'{VZVt}' + "\n")
    FW.write("Velocity Theta     [m/s]    :" + f'{VTVt}' + "\n")
    FW.write("Relative Velocity  [m/s]    :" + f'{WWt}' + "\n")
    FW.write("Relative Axial     [m/s]    :" + f'{WZWt}' + "\n")
    FW.write("Relative Theta     [m/s]    :" + f'{WTWt}' + "\n")
    FW.write("ALPHA              [deg]    :" + f'{ALPt}' + "\n")
    FW.write("BETA               [deg]    :" + f'{BETt}' + "\n")
    FW.write("                                            \n")  
    
#####################################################
#####################################################
# PART B.viii - ISENTROPIC CALCULATIONS
#####################################################
#####################################################
    
    '''ISENTROPIC PROCESS'''
    P_isen = Poi
    S_isen = S_constant
    H_isen = HSP(S_isen, P_isen)
    T_isen = TSH(S_isen, H_isen)
    
    PR_F = Current_Pressure / Poi_initial
    TR_F = Current_Temperature / Toi_initial
    
    if Fluid == 0:
        AD_EFF   = (PR_F**((Gamma-1) / Gamma) - 1) / (TR_F - 1)
        POLY_EFF = (Gamma - 1) / Gamma * math.log(PR_F) / math.log(TR_F)
    elif Fluid == 1:
        AD_EFF   = (H_isen - H_initial) / (H - H_initial)
        POLY_EFF = 'NA'
        
#####################################################
#####################################################
# PART B.ix - FINAL COMPRESSOR DATA OUTPUT
#####################################################
#####################################################    
    
    print(' ')
    print('**** FINAL COMPRESSOR DATA ****')
    print(f'Total Pressure          = {round(Current_Pressure, 3)} kPa')
    print(f'Final Pressure Ratio    = {round(PR_F, 3)} ')
    print(f'Final Temperature Ratio = {round(TR_F, 3)} ')
    print(f'Adiabatic Efficiency    = {AD_EFF}')
    print(f'Polytropic Efficiency   = {POLY_EFF}')
    
    FW.write("**** FINAL COMPRESSOR DATA **** \n")
    FW.write("Total Pressure     [kPa]    :" + f'{round(Current_Pressure, 3)}' + "\n")
    FW.write("Final Pressure Ratio        :" + f'{round(PR_F, 3)}' + "\n")
    FW.write("Final Temperature Ratio     :" + f'{round(TR_F, 3)}' + "\n")
    FW.write("Adiabatic Efficiency        :" + f'{AD_EFF}' + "\n")
    FW.write("                                            \n")  
    FW.close
    #print(f'Enthalpy Values         = {Hplot} J/kg')    
    
#####################################################
#####################################################
# PART B.x - REFERENCE VALUES FOR .TSTK and .TWAL Files
#####################################################
#####################################################
    
    '''CREATING VALUES FOR T-AXI STACK AND WALL FILES'''
    if IGV_Q == 1:
        refl = Y_TIP[3]/1000    #{m}
    else:
        refl = Y_TIP[1]/1000    #{m}
        
    PT1 = Poi_initial
    TT1 = Toi_initial
    if Fluid == 1: 
        A0              = float(APT(PT1, TT1))
        viscosity_stack = float(VISCOSITY[0]/(10**6))
        Rhoin           = float(RPT(PT1, TT1))
        Re              = Rhoin*A0*refl/viscosity_stack
        mdot_norm       = float(MFR*A0 / (PT1*1000*GAMMAVALUE[0]*refl**2))
        rnorm           = refl
    elif Fluid == 0: 
        A0              = float((Gamma*Rgas*TT1)**0.5)
        viscosity_stack = float(1.458e-6*(TT1**1.5)/(TT1+110.4))
        Rhoin           = PT1*1000/(Rgas*TT1)
        Re              = Rhoin*A0*refl/viscosity_stack
        mdot_norm       = float(MFR*A0 / (PT1*1000*Gamma*refl**2))
        rnorm           = refl
        
    TOTAL_CLEARANCE     = refl*TCLC        #{mm}
    
#####################################################
#####################################################
# PART B.xi - CALL WALL FUNCTION (FILE)
#####################################################
#####################################################
    
    ROWS, RHUB, ZHUB, RTIP, ZTIP = WALLS(X_CROSS_SECTION, Y_CROSS_SECTION, refl, IGV_Q, Stage_Number, XRIF, Rad_h3, Radt_3, Output, EXTRAPOLATION)
   
#####################################################
#####################################################
# PART B.xii - CALL STACK FUNCTION (FILE)
#####################################################
##################################################### 
   
    STACK(Output, GAMMAVALUE, TT, mdot_norm, Re, refl, Gamma, Fluid, Stage_Number, IGV_Q, NBS, Y_BLADES, TCLC, X_BLADES, igv_Phi, igv_L, NBR, RAV, A0, R_Phi, Lambda, S_Phi, X_HUB, n, viscosity_stack, ROWS, RHUB, ZHUB, RTIP, ZTIP, RVTHI, Bleed, BLADEROWS)
    
#####################################################
#####################################################
# PART B.xiii - CALL PLOTS FUNCTION
#####################################################
#####################################################
    if Fluidtype == 5: 
            Fluidtypename = 'Custom Fluid'
            Output = f'{Fluidtypename}_' + day
            print(f'Fluid Type Name = {Fluidtypename} \n')
    PLOTS(n, Pplot, Tplot, Fluidtypename, Rplot, Hplot, Splot, HH, SS, RAD, PT, PS, TT, TS, RS, VV, VZV, VTV, ALP, BET, CAMBER, MACHSP, MACHSH, MACHST, LOADING_C, FLOW_C, DIFFUSION, Output, Poi_initial, Toi_initial, AD_EFF_perstage, Value_Reaction, RADh, RADt, GAMMAVALUE, VISCOSITY, SPECIFICHEAT, X_CROSS_SECTION, Y_CROSS_SECTION, STAGGERANGLE1_h, STAGGERANGLE1_p, STAGGERANGLE1_t, STAGGERANGLE2_h, STAGGERANGLE2_p, STAGGERANGLE2_t, ZETAANGLE_hD, ZETAANGLE_pD, ZETAANGLE_tD, AXIAL_WIDTH_R, CHORD_LENGTH_R, RAR_R, ROTORWIDTH_h, ROTORWIDTH_p, ROTORWIDTH_t, NBR, STAGGERANGLE1_Sh, STAGGERANGLE1_Sp, STAGGERANGLE1_St, STAGGERANGLE2_Sh, STAGGERANGLE2_Sp, STAGGERANGLE2_St, ZETAANGLE_ShD, ZETAANGLE_SpD, ZETAANGLE_StD, AXIAL_WIDTH_S, CHORD_LENGTH_S, RAR_S, STATORWIDTH_h, STATORWIDTH_p, STATORWIDTH_t, NBS, X_HUB, Y_HUB, X_TIP, Y_TIP, X_BLADES, Y_BLADES, X_ROTOR, Y_ROTOR, X_STATOR, Y_STATOR, IGV_Q, RHUB, ZHUB, RTIP, ZTIP, refl) 
#

#####################################################
#####################################################
# PART C - INLET CALCULATIONS FUNCTION
#####################################################
#####################################################

def Inlet_Calculations(STG, MFR, Rev, Poi, Toi, Gamma, Rgas, TCLC, Radius, Stage_Number, MACH1, alpha1, Q, igv_S, igv_AR, igv_Phi, igv_a, igv_M, igv_L, igv_R_Coef, igv_rad_tip, typer, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, Fluid, VPT, Output, Rad_tip, KKI):
#
    #
    Pi = 3.1415926

    print(' ')
    
    ''' Calculate INLET Properties '''
    
    CP     = Rgas * Gamma / (Gamma - 1)
     
#####################################################
#####################################################
# PART C.i - CALCULATIONS WITH AN IGV
#####################################################
#####################################################
    
    if Q == 1 : 
    #
        #       
#####################################################
#####################################################
# PART C.i.1 - IGV CALCULATIONS
#####################################################
#####################################################
        if Fluid == 1:
        #
            #
            '''OBTAIN ACTUAL STATIC PROPERTIES'''
            Tsi_igv, Psi_igv, Rsi_igv, Ai_igv, Vi_igv, Me, h_igv, S_igv, H_igv = STATIC_CONV(HPT, SPT, APT, RPT, TSH, PSH, Poi, Toi, igv_M)
            S_constant = S_igv
            #
        #
        elif Fluid == 0:
        #
            #
            '''PERFECT GAS ESTIMATES'''
            #Toi, Poi [IGV] = Toi, Poi [INLET]
            Tsi_igv    = Toi /  (1 + (((Gamma - 1) / 2) * igv_M**2))
            Psi_igv    = Poi / ((1 + (((Gamma - 1) / 2) * igv_M**2))**(Gamma / (Gamma - 1)))
            Rsi_igv    = Psi_igv / (Rgas /1000 * Tsi_igv)
            Ai_igv     = (Gamma * Rgas * Tsi_igv)**0.5
            Vi_igv     = igv_M * Ai_igv
            H_igv      = float(HPT(Poi, Toi))
            S_igv      = float(SPT(Poi, Toi))
            S_constant = S_igv
            #
        #
        ALPHA_igv  = igv_a * (Pi/180)
        Vzi_igv    = Vi_igv * math.cos(ALPHA_igv)
        Vti_igv    = Vi_igv * math.sin(ALPHA_igv)
        A_comp_igv = MFR / (Rsi_igv * Vzi_igv)
        
        if KKI == 0: 
            '''Radius Type'''
            Radh_igv, Radp_igv, Radt_igv, Hslice = Radius_Type(typer, Radius, A_comp_igv, Pi)
        elif KKI == 1: 
            Radh_igv, Radp_igv, Radt_igv, Hslice = Radius_Type(typer, igv_rad_tip, A_comp_igv, Pi)
        
        '''IGV Outlet Velocity'''
        Vzi_igv2     = Vzi_igv
        Vti_igv2     = Vti_igv * (Radp_igv / Radt_igv)
        ALPHA_igv2   = math.atan(Vti_igv2 / Vzi_igv2)
        ALPHA_igv2d  = ALPHA_igv2 * (180 / Pi)
        Vi_igv2      = Vzi_igv2 / (math.cos(ALPHA_igv2d))
        
        '''Different Velocity at different radii'''
        Vzh_igv2     = Vzi_igv
        Vzp_igv2     = Vzi_igv
        Vzt_igv2     = Vzi_igv
        Vth_igv2     = Vti_igv * (Radp_igv / Radh_igv)
        Vtp_igv2     = Vti_igv * (Radp_igv / Radp_igv)
        Vtt_igv2     = Vti_igv * (Radp_igv / Radt_igv)
        ALPHAh_igv   = math.atan(Vth_igv2 / Vzh_igv2)
        ALPHAp_igv   = math.atan(Vtp_igv2 / Vzp_igv2)
        ALPHAt_igv   = math.atan(Vtt_igv2 / Vzt_igv2)
        Vh_igv2      = Vzh_igv2 / math.cos(ALPHAh_igv)
        Vp_igv2      = Vzp_igv2 / math.cos(ALPHAp_igv)
        Vt_igv2      = Vzt_igv2 / math.cos(ALPHAt_igv)
        
        
        print('**** IGV CALCULATIONS ****')
        print(f'Speed of Sound        = {round(Ai_igv, 3)} m/s')
        print(f'Velocity Entering     = {round(Vi_igv, 3)} m/s')
        print(f'Mach Entering         = {round(igv_M, 3)} ')
        print(f'Static Density        = {round(Rsi_igv, 3)} kg/m3')
        print(f'Total Temperature     = {round(Toi, 3)} K')
        print(f'Static Temperature    = {round(Tsi_igv, 3)} K')
        print(f'Total Pressure        = {round(Poi, 3)} kPa')
        print(f'Static Pressure       = {round(Psi_igv, 3)} kPa')
        print(f'Velocity Exiting      = {round(Vi_igv2, 3)} m/s')
        print(' ')
            
#####################################################
#####################################################
# PART C.i.2 - ROTOR INLET CALCULATIONS
#####################################################
#####################################################
 
        if Fluid == 1: 
        #
            #
            '''OBTAIN ACTUAL STATIC PROPERTIES'''
            Poi    = Poi - (igv_Phi * (Poi - Psi_igv))
            Tsi, Psi, Rsi, Ai, Vi, Mi, h, S, H = STATIC_CONV(HPT, SPT, APT, RPT, TSH, PSH, Poi, Toi, MACH1[STG])
            #
        #
        elif Fluid == 0:
        #
            #
            '''PERFECT GAS ESTIMATES''' 
            H      = float(HPT(Poi, Toi))
            S      = float(SPT(Poi, Toi))
            Tsi    = Toi / (1 + (((Gamma - 1) / 2) * MACH1[STG]**2))
            Poi    = Poi - (igv_Phi * (Poi - Psi_igv))
            Psi    = Poi / ((Toi / Tsi)**(Gamma / (Gamma - 1)))
            Roi    = Poi / (Rgas / 1000 * Toi)
            Rsi    = Psi / (Rgas / 1000 * Tsi)
            Ai     = (Gamma * Rgas * Tsi)**0.5
            Vi     = MACH1[STG] * Ai
            h      = H - (Vi**2 / 2)
            Mi     = Vi / Ai
            #
        #
        ALPHA  = alpha1[STG] * (Pi/180)
        ALPHAd = alpha1[STG]
        Vzi    = Vi * math.cos(ALPHA)
        A_comp = MFR / (Rsi * Vzi * igv_L)      #{m^2}
        RAV    = Rev * (Pi / 30)              #Rad/s
        
        if KKI == 0: 
            '''Radius Type'''
            Radh, Radp, Radt, Height = Radius_Type(typer, Radius, A_comp, Pi)
        elif KKI == 1: 
            Radh, Radp, Radt, Height = Radius_Type(typer, Rad_tip[STG], A_comp, Pi)            
        
        #Constant Pitch...
        U      = RAV * Radp
        Vti    = Vi * math.sin(ALPHA)
        Wti    = Vti - U
        Wzi    = Vzi
        BETA   = math.atan(Wti / Wzi)
        BETAd  = BETA * (180 / Pi)
        Wi     = Wzi / math.cos(BETA)
            
        #HUB VELO TRIANGLES
        U_h, ALPHA_h, Vi_h, Vti_h, Vzi_h, BETA_h, Wi_h, Wti_h, Wzi_h, Poi_h, Toi_h, Tsi_h, Psi_h, M_h, Ai_h, Rsi_h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vzi, Vti, RAV, Radh, Radp, Toi, Poi, Gamma, CP, Rgas, Fluid)
        A_comp_h = MFR / (Rsi_h * Vzi_h * igv_L)
        
        #TIP VELO TRIANGLES
        U_t, ALPHA_t, Vi_t, Vti_t, Vzi_t, BETA_t, Wi_t, Wti_t, Wzi_t, Poi_t, Toi_t, Tsi_t, Psi_t, M_t, Ai_t, Rsi_t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vzi, Vti, RAV, Radt, Radp, Toi, Poi, Gamma, CP, Rgas, Fluid)
        A_comp_t = MFR / (Rsi_t * Vzi_t * igv_L)
    
        #Eventually for loop for radii types...
        if Fluid == 1: 
        #
            #
            '''OBTAIN ACTUAL RELATIVE TOTAL PROPERTIES'''
            H_rel, S_rel, TTR, PTR, MR = TOTAL_CONV(HPT, SPT, APT, RPT, TSH, PSH, h, S, Tsi, Wi, Ai)
            #
        #
        elif Fluid == 0:
        #
            #
            TTR    = Tsi + (Wi**2 / (2 * CP))
            PTR    = Psi * (((TTR / Tsi))**(Gamma / (Gamma - 1)))
            MR     = Wi / Ai
            #
        #
        Q      = 2
        
#####################################################
#####################################################
# PART C.i.3 - IGV GEOMETRY
#####################################################
#####################################################
        
        ''' IGV GEOMETRY '''
        #ALPHA (pitch) already in radians, Convert others to rad too
        #HUB
        GAMMA1_h_igv = ALPHAh_igv                                #STAGGER ANGLE
        GAMMA2_h_igv = (4*(igv_S**0.5)*ALPHA_h* (Pi/180) - ALPHAh_igv) / (4*(igv_S**0.5) - 1)
        ZETA_h_igv   = (GAMMA1_h_igv + GAMMA2_h_igv) / 2         #ZETA
        ZETAd_h_igv  = ZETA_h_igv*(180/Pi)
        
        #PITCH
        GAMMA1_p_igv = ALPHAp_igv                                #STAGGER ANGLE                             
        GAMMA2_p_igv = (4*(igv_S**0.5)*ALPHA - ALPHAp_igv) / (4*(igv_S**0.5) - 1)
        ZETA_p_igv   = (GAMMA1_p_igv + GAMMA2_p_igv) / 2         #ZETA
        ZETAd_p_igv  = ZETA_p_igv*(180/Pi)
        
        #TIP
        GAMMA1_t_igv = ALPHAt_igv                                #STAGGER ANGLE
        GAMMA2_t_igv = (4*(igv_S**0.5)*ALPHA_t* (Pi/180) - ALPHAt_igv) / (4*(igv_S**0.5) - 1)
        ZETA_t_igv   = (GAMMA1_t_igv + GAMMA2_t_igv) / 2         #ZETA
        ZETAd_t_igv  = ZETA_t_igv*(180/Pi)
        
        #OVERALL
        WIDTH_igv    = (Hslice + Height) / 2 / igv_AR            #AXIAL WIDTH
        CSTA_igv     = WIDTH_igv / math.cos(ZETA_p_igv)          #CHORD LENGTH
        
        #HUB
        WIDTH_h_igv  = CSTA_igv * math.cos(ZETA_h_igv)           #HUB WIDTH
        
        #PITCH
        WIDTH_p_igv  = CSTA_igv * math.cos(ZETA_p_igv)           #PITCH WIDTH
        
        #TIP
        WIDTH_t_igv  = CSTA_igv * math.cos(ZETA_t_igv)           #TIP WIDTH
        
        SNB = math.floor(Pi*(Radp_igv + Radp) / (CSTA_igv/igv_S) + 1) #NUM. OF BLADES
        
#####################################################
#####################################################
# PART C.i.4 - IGV CROSS SECTION GEOMETRY
#####################################################
#####################################################        

        ''' IGV CROSS SECTION GEOMETRY '''
        #COORDINATES
        X_igv = []      #X Coords
        Y_igv = []      #Y Coords
        X_igv_h = []    
        Y_igv_h = []
        X_igv_t = []
        Y_igv_t = []
        X_igv_blades = []
        Y_igv_blades = []
        
        X_igv.append(0*1000)
        Y_igv.append(round(Radt_igv, 9)*1000)
        X_igv_t.append(X_igv[0])
        Y_igv_t.append(Y_igv[0])
        
        
        X_igv.append(0*1000)
        Y_igv.append(round(Radh_igv, 9)*1000)
        X_igv_h.append(X_igv[1])
        Y_igv_h.append(Y_igv[1])
        X_igv_blades.append(X_igv[1])
        Y_igv_blades.append(Y_igv[1])
        X_igv_blades.append(X_igv[0])
        Y_igv_blades.append(Y_igv[0])
        
        X_igv.append(round(X_igv[1]/1000 + (WIDTH_t_igv - WIDTH_h_igv)/2, 9)*1000)
        Y_igv.append(round(Radh_igv, 9)*1000)
        X_igv_h.append(X_igv[2])
        Y_igv_h.append(Y_igv[2])
        
        X_igv.append(round(X_igv[1]/1000 + (WIDTH_t_igv + WIDTH_h_igv)/2, 9)*1000)
        Y_igv.append(round(Radh_igv, 9)*1000)
        X_igv_h.append(X_igv[3])
        Y_igv_h.append(Y_igv[3])
        
        X_igv.append(round(X_igv[0]/1000 + WIDTH_t_igv, 9)*1000)
        Y_igv.append(round(Radt_igv, 9)*1000)
        X_igv_t.append(X_igv[4])
        Y_igv_t.append(Y_igv[4])
        X_igv_blades.append(X_igv[4])
        Y_igv_blades.append(Y_igv[4])
        X_igv_blades.append(X_igv[3])
        Y_igv_blades.append(Y_igv[3])
        
        X_igv.append(round(X_igv[1]/1000 + (3*WIDTH_t_igv + WIDTH_h_igv)/4 + igv_R_Coef*WIDTH_p_igv, 9)*1000)
        Y_igv.append(round(Radh_igv, 9)*1000)
        X_igv_h.append(X_igv[5])
        Y_igv_h.append(Y_igv[5])
        
        X_igv.append(round(X_igv[5]/1000, 9)*1000)
        Y_igv.append(round(Radh, 9)*1000)
        X_igv_h.append(X_igv[6])
        Y_igv_h.append(Y_igv[6])
        
        X_igv.append(round(X_igv[5]/1000, 9)*1000)
        Y_igv.append(round(Radt, 9)*1000)
        X_igv_t.append(X_igv[7])
        Y_igv_t.append(Y_igv[7])
        
        XRIF = X_igv[5]
        
        FG = open(f'{Output}_Geometry.txt', "w+") #Replaces previous file or creates
        FG.write("**** GEOMETRY FILE, STAGE BY STAGE************* \n")
        FG.write("**** IGV GEOMETRY **** \n")
        FG.write("Chord Length            : " + f'{round(CSTA_igv, 6)*1000} mm' + "\n")
        FG.write("Axial Width             : " + f'{round(WIDTH_igv, 6)*1000} mm' + "\n")
        FG.write("                                            \n") 
        FG.write("Hub Stagger Angle 1     : " + f'{round(GAMMA1_h_igv, 3)} rad' + "\n")
        FG.write("Hub Stagger Angle 2     : " + f'{round(GAMMA2_h_igv, 3)} rad' + "\n")
        FG.write("Hub Zeta Angle          : " + f'{round(ZETA_h_igv, 3)} rad' + "\n")
        FG.write("Hub Zeta Angle          : " + f'{round(ZETAd_h_igv, 3)} deg' + "\n")
        FG.write("                                            \n")  
        FG.write("Pitch Stagger Angle 1   : " + f'{round(GAMMA1_p_igv, 3)} rad' + "\n")
        FG.write("Pitch Stagger Angle 2   : " + f'{round(GAMMA2_p_igv, 3)} rad' + "\n")
        FG.write("Pitch Zeta Angle        : " + f'{round(ZETA_p_igv, 3)} rad' + "\n")
        FG.write("Pitch Zeta Angle        : " + f'{round(ZETAd_p_igv, 3)} deg' + "\n")
        FG.write("                                            \n") 
        FG.write("Tip Stagger Angle 1     : " + f'{round(GAMMA1_t_igv, 3)} rad' + "\n")
        FG.write("Tip Stagger Angle 2     : " + f'{round(GAMMA2_t_igv, 3)} rad' + "\n")
        FG.write("Tip Zeta Angle          : " + f'{round(ZETA_t_igv, 3)} rad' + "\n")
        FG.write("Tip Zeta Angle          : " + f'{round(ZETAd_t_igv, 3)} deg' + "\n")
        FG.write("                                            \n") 
        FG.write("**** IGV CROSS SECTION GEOMETRY **** \n")
        FG.write("X COORDINATES           : " + f'{X_igv} mm' + "\n")
        FG.write("Y COORDINATES           : " + f'{Y_igv} mm' + "\n")
        FG.write("IGV END X COORDINATE    : " + f'{XRIF} mm' + "\n")
        FG.write("Number of Blades        : " + f'{SNB}' + "\n")
        FG.write("                                            \n")
        FG.close
        #
    #
    
#####################################################
#####################################################
# PART C.ii - CALCULATIONS WITHOUT AN IGV
#####################################################
#####################################################
    
    elif Q == 2 : 
    #
        #
        ''' 
            WITH NO IGV
        '''
        
#####################################################
#####################################################
# PART C.ii.1 - ROTOR INLET CALCULATIONS
#####################################################
#####################################################
        
        if Fluid == 1: 
        #
            #
            '''OBTAIN ACTUAL STATIC PROPERTIES'''
            Tsi, Psi, Rsi, Ai, Vi, Mi, h, S, H = STATIC_CONV(HPT, SPT, APT, RPT, TSH, PSH, Poi, Toi, MACH1[STG])
            S_constant = S
            #
        #
        elif Fluid == 0:
        #
            #
            '''PERFECT GAS ESTIMATES'''
            #Pitch 
            H      = float(HPT(Poi, Toi))
            S      = float(SPT(Poi, Toi))
            S_constant = S
            Tsi    = Toi /  (1 + (((Gamma - 1) / 2) * MACH1[STG]**2))
            Psi    = Poi / ((1 + (((Gamma - 1) / 2) * MACH1[STG]**2))**(Gamma / (Gamma - 1)))
            Roi    = Poi / (Rgas / 1000 * Toi)
            Rsi    = Psi / (Rgas / 1000 * Tsi)
            Ai     = (Gamma * Rgas * Tsi)**0.5
            Vi     = MACH1[STG] * Ai
            h      = H - (Vi**2 / 2)
            Mi     = MACH1[STG]
            #
        #
        ALPHA  = alpha1[STG] * (Pi/180)
        ALPHAd = alpha1[STG]
        Vzi    = Vi * math.cos(ALPHA)
        A_comp = MFR / (Rsi * Vzi)
        RAV    = Rev * (Pi / 30)
        
        if KKI == 0: 
            '''Radius Type'''
            Radh, Radp, Radt, Height = Radius_Type(typer, Radius, A_comp, Pi)
        elif KKI == 1: 
            Radh, Radp, Radt, Height = Radius_Type(typer, Rad_tip[STG], A_comp, Pi) 
            
        '''Pitch-Line'''
        U      = RAV * Radp
        Vti    = Vi * math.sin(ALPHA)
        Wti    = Vti - U
        Wzi    = Vzi
        BETA   = math.atan(Wti / Wzi)
        BETAd  = BETA * (180 / Pi)
        Wi     = Wzi / math.cos(BETA)
        
        '''CONSTANT HUB VELO TRIANGLES'''
        U_h, ALPHA_h, Vi_h, Vti_h, Vzi_h, BETA_h, Wi_h, Wti_h, Wzi_h, Poi_h, Toi_h, Tsi_h, Psi_h, M_h, Ai_h, Rsi_h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vzi, Vti, RAV, Radh, Radp, Toi, Poi, Gamma, CP, Rgas, Fluid)
        A_comp_h = MFR / (Rsi_h * Vzi_h)
        
        '''CONSTANT TIP VELO TRIANGLES'''
        U_t, ALPHA_t, Vi_t, Vti_t, Vzi_t, BETA_t, Wi_t, Wti_t, Wzi_t, Poi_t, Toi_t, Tsi_t, Psi_t, M_t, Ai_t, Rsi_t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vzi, Vti, RAV, Radt, Radp, Toi, Poi, Gamma, CP, Rgas, Fluid)
        A_comp_t = MFR / (Rsi_t * Vzi_t)
        
        Hslice = 0 #No IGV so setting height at IGV exit to 0
    
        #Eventually for loop for radii types...
        if Fluid == 1:
        #
            #
            '''OBTAIN ACTUAL RELATIVE TOTAL PROPERTIES'''
            H_rel, S_rel, TTR, PTR, MR = TOTAL_CONV(HPT, SPT, APT, RPT, TSH, PSH, h, S, Tsi, Wi, Ai)
            #
        #
        elif Fluid == 0:
        #
            #
            TTR    = Tsi + (Wi**2 / (2*CP))
            PTR    = Psi * (((TTR / Tsi))**(Gamma / (Gamma - 1)))
            MR     = Wi / Ai
            #
        #
        
#####################################################
#####################################################
# PART C.ii.2 - (no) IGV GEOMETRY CALCS
#####################################################
#####################################################
        
        ''' NO IGV FILE SO GEOMETRY STARTS AT x=0'''
        #COORDINATES
        X_igv = []      #X Coords
        Y_igv = []      #Y Coords
        X_igv_h = []    
        Y_igv_h = []
        X_igv_t = []
        Y_igv_t = []
        X_igv_blades = []
        Y_igv_blades = []
        XRIF = 0
        SNB  = 0
        '''Opening/Creating Geometry File even though there's no IGV in this case'''
        FG = open(f'{Output}_Geometry.txt', "w+") #Replaces previous file or creates
        FG.write("**** GEOMETRY FILE **** \n")
        FG.write("\n")
        FG.close
        #
    #
    else:
        print("How'd you get to this option!?")
        #
    #
    GammaValue = float(GPT(Psi, Tsi))
    Viscosity  = float(VPT(Psi, Tsi))
    CCPP       = float(CpPT(Psi, Tsi))
    
#####################################################
#####################################################
# PART C.iii - DISPLAY DATA
#####################################################
#####################################################
    
    #Display Data in Window
    print('**** INLET CALCULATIONS ****')
    print('**** PITCH-LINE ****')
    print(f'Radius                = {round(Radp, 3)} m')
    print(f'Static Temperature    = {round(Tsi, 3)} K')
    print(f'Total  Temperature    = {round(Toi, 3)} K')
    print(f'Static Pressure       = {round(Psi, 3)} kPa')
    print(f'Total Pressure        = {round(Poi, 3)} kPa')
    print(f'Total Enthalpy        = {round(H, 3)} J/kg')
    print(f'Total Entropy         = {round(S, 3)} J/kgK')
    print(f'Static Enthalpy       = {round(h, 3)} J/kg')
    print(f'Static Density        = {round(Rsi, 3)} kg/m3')
    print(f'Mach #                = {round(Mi, 3)} ')
    print(f'Total Velocity        = {round(Vi, 3)} m/s')
    print(f'Axial Velocity        = {round(Vzi, 3)} m/s')
    print(f'Velocity Theta        = {round(Vti, 3)} m/s')
    print(f'Abs Angle (Alpha)     = {round(ALPHAd, 3)} deg')
    print(f'Relative Angle (Beta) = {abs(round(BETAd, 3))} deg')
    print(f'Area Pitch            = {round(A_comp, 5)} m^2')
    
    print(' ')
    print(f'**** CONSTANT HUB ****')
    print(f'Radius                = {round(Radh, 3)} m')
    print(f'Total Temperature     = {round(Toi_h, 3)} K')
    print(f'Total Pressure        = {round(Poi_h, 3)} kPa')
    print(f'Static Temperature    = {round(Tsi_h, 3)} K')
    print(f'Static Pressure       = {round(Psi_h, 3)} kPa')
    print(f'Static Density        = {round(Rsi_h, 3)} kg/m3')
    print(f'Mach #                = {round(M_h, 3)} ')
    print(f'Velocity              = {round(Vi_h, 3)} m/s')
    print(f'Velocity Theta        = {round(Vti_h, 3)} m/s')
    print(f'Axial Velocity        = {round(Vzi_h, 3)} m/s')
    print(f'Abs Angle (Alpha)     = {round(ALPHA_h, 3)} deg')
    print(f'Relative Angle (Beta) = {abs(round(BETA_h, 3))} deg')
    print(f'Area Hub              = {round(A_comp_h, 5)} m^2')
    print(' ')
    print(f'**** TIP CALCULATIONS ****')
    print(f'Radius                = {round(Radt, 3)} m')
    print(f'Total Temperature     = {round(Toi_t, 3)} K')
    print(f'Total Pressure        = {round(Poi_t, 3)} kPa')
    print(f'Static Temperature    = {round(Tsi_t, 3)} K')
    print(f'Static Pressure       = {round(Psi_t, 3)} kPa')
    print(f'Static Density        = {round(Rsi_t, 3)} kg/m3')
    print(f'Mach #                = {round(M_t, 3)} ')
    print(f'Velocity              = {round(Vi_t, 3)} m/s')
    print(f'Velocity Theta        = {round(Vti_t, 3)} m/s')
    print(f'Axial Velocity        = {round(Vzi_t, 3)} m/s')
    print(f'Abs Angle (Alpha)     = {round(ALPHA_t, 3)} deg')
    print(f'Relative Angle (Beta) = {abs(round(BETA_t, 3))} deg')
    print(f'Area Tip              = {round(A_comp_t, 5)} m^2')
    print(' ')
    print('**** RELATIVE PROPERTIES ****')
    print(f'Relative Velocity     = {round(Wi, 3)} m/s')
    print(f'Relative Tot. Temp    = {round(TTR, 3)} K')
    print(f'Relative Tot. Press   = {round(PTR, 3)} kPa')
    #print(f'Relative Enthalpy     = {round(H_rel, 3)} J/kg')
    print(f'Relative Mach         = {round(MR, 3)}')

    ''' Return the values to go into the first Rotor '''
    return CP, Tsi_h, Psi, Rsi_h, Ai_h, Vi_h, ALPHA_h, Vzi_h, A_comp_h, Radt, Height, Radp, RAV, U_h, Vti_h, Wti_h, Wzi_h, BETA_h, Wi_h, TTR, PTR, MR, Q, Poi_h, H, S, Poi_h, S_constant, ALPHAd, BETAd, Toi_h, Poi_h, Radh, Psi_h, M_h, Vti_h, Vzi_h, ALPHA_h, BETA_h, Toi_t, Poi_t, Tsi_t, Psi_t, Rsi_t, M_t, Vi_t, Vti_t, Vzi_t, ALPHA_t, BETA_t, Tsi, Toi, Psi, Poi, Rsi, Mi, Vi, Vzi, Vti, ALPHAd, BETAd, Rsi, U, Wi, Wti, Wzi, h, GammaValue, Viscosity, Wi_t, Wti_t, Wzi_t, CCPP, Hslice, XRIF, X_igv, Y_igv, X_igv_h, Y_igv_h, X_igv_t, Y_igv_t, X_igv_blades, Y_igv_blades, SNB, H_igv, S_igv
    #
#

#####################################################
#####################################################
# PART D - ROTOR OUTLET CALCULATIONS FUNCTION
#####################################################
#####################################################

def Rotor_Outlet_Calculations(n, STG, Stage_Number, DH, DTT, VZR, R_Phi, Lambda, Radius, Radp, RAV, Toi, Vt_h, CP, Vzi, Gamma, Rgas, PTR, Psi, TTR, MFR, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, typer, H, S, Po_h, Fluid, R_Sol, Wi, Wti, Wzi, VPT, Wi_p, Wti_p, Wzi_p, Vt_p, Vt_t, Vz_p, Vz_t, Radh, Radt, To_p, To_t, KKI, Rad_tip):
#                                                                                                                                                                                                                                                                      
    #
    Pi = 3.1415926

    print('**** ROTOR/STATOR INTERFACE CALCS ****')
    
#####################################################
#####################################################
# PART D.i - CALCULATIONS (REAL GAS)
#####################################################
#####################################################
    
    #Eventually for loop for radii types...
    if Fluid == 1:
    #
        #
        if typer == 1: 
            '''CONSTANT HUB'''
            #VARIABLES   
            #Hub Case
            #Toi, Radius, Vt_h, Vz_2
            
            Rad_2    = Radius                                       #{m}
            Vt_2     = (((DH[STG] * 1000 / RAV) + (Radius * Vt_h))) / Rad_2  #{m/s}
            Vz_2     = VZR[STG] * Vzi                               #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            U_2      = RAV * Rad_2                                  #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{rad}
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            
            #DF       = 1 - (W_2/Wi) + (abs(Wti - Wt_2) / (2*R_Sol[STG]*Wi))
                        
            '''RELATIVE, PITCH-LINE'''
            WTR_2    = ((Rad_2 * Vt_2) / Radp) - (RAV * Radp)       #{m/s}
            WZR_2    = Wz_2                                         #{m/s}
            BETAR_2  = math.atan(WTR_2 / WZR_2)                     #{rad}
            BETAR_2d = BETAR_2 * (180/Pi)                           #{deg}
            WR_2     = WZR_2 / math.cos(BETAR_2)                    #{m/s}
            PTR_2    = PTR - (R_Phi[STG] * (PTR - Psi))             #{kPa}
            TTR_2    = TTR                                          #{K}
            
            #DF       = 1 - (WR_2/Wi) + (abs(Wti - W_2) / (2*R_Sol[STG]*Wi))
        
            '''STATIC CONVERSION, HUB'''
            H_2      = H + (DH[STG] * 1000)                         #{J/kg}
            h_2      = H_2 - (V_2**2 / 2)                           #{J/kg}
            #S_2      = S                                            #{J/kgK}
            S_2      = float(SPT(PTR_2, TTR_2))                     #{J/kgK}
            To_2     = float(TSH(S_2, H_2 ))                        #{K}
            Ts_2     = float(TSH(S_2, h_2))                         #{K}
            Ps_2     = float(PSH(S_2, h_2))                         #{kPa}
            Po_2     = float(PSH(S_2, H_2))                         #{kPa}   
            A_2      = float(APT(Ps_2, Ts_2))                       #{m/s}
            M_2      = V_2 / A_2                                    #{ }
            Rs_2     = float(RPT(Ps_2, Ts_2))                       #{kg/m3}
    
        if typer == 2: 
            '''CONSTANT PITCH'''
            #VARIABLES 
            #Pitch Case
            #To_p, Radp, Vt_p, Vz_p
            
            Rad_2    = Radp                                         #{m}
            Vt_2     = (((DH[STG] * 1000 / RAV) + (Radius * Vt_p))) / Rad_2  #{m/s}
            Vz_2     = VZR[STG] * Vz_p                               #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            U_2      = RAV * Rad_2                                  #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{rad}
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            
            #DF       = 1 - (W_2/Wi) + (abs(Wti - Wt_2) / (2*R_Sol[STG]*Wi))
                        
            '''RELATIVE, PITCH-LINE'''
            WTR_2    = ((Rad_2 * Vt_2) / Radp) - (RAV * Radp)       #{m/s}
            WZR_2    = Wz_2                                         #{m/s}
            BETAR_2  = math.atan(WTR_2 / WZR_2)                     #{rad}
            BETAR_2d = BETAR_2 * (180/Pi)                           #{deg}
            WR_2     = WZR_2 / math.cos(BETAR_2)                    #{m/s}
            PTR_2    = PTR - (R_Phi[STG] * (PTR - Psi))             #{kPa}
            TTR_2    = TTR                                          #{K}
            
            #DF       = 1 - (WR_2/Wi) + (abs(Wti - W_2) / (2*R_Sol[STG]*Wi))
        
            '''STATIC CONVERSION, PITCH'''
            H_2      = H + (DH[STG] * 1000)                         #{J/kg}
            h_2      = H_2 - (V_2**2 / 2)                           #{J/kg}
            #S_2      = S                                           #{J/kgK}
            S_2      = float(SPT(PTR_2, TTR_2))                     #{J/kgK}
            To_2     = float(TSH(S_2, H_2 ))                        #{K}
            Ts_2     = float(TSH(S_2, h_2))                         #{K}
            Ps_2     = float(PSH(S_2, h_2))                         #{kPa}
            Po_2     = float(PSH(S_2, H_2))                         #{kPa}   
            A_2      = float(APT(Ps_2, Ts_2))                       #{m/s}
            M_2      = V_2 / A_2                                    #{ }
            Rs_2     = float(RPT(Ps_2, Ts_2))                       #{kg/m3}
    
        if typer == 3: 
            '''CONSTANT TIP'''
            #VARIABLES   
            #Tip Case
            #To_t, Radt, Vt_t, Vz_t
            #Rad_tip (if Arbitrary Option)
            
            if KKI == 0: 
                Rad_2    = Radt                                   #{m}
            elif KKI == 1: 
                Rad_2    = Rad_tip[STG]
                
            Vt_2     = (((DH[STG] * 1000 / RAV) + (Radius * Vt_t))) / Rad_2  #{m/s}
            Vz_2     = VZR[STG] * Vz_t                              #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            U_2      = RAV * Rad_2                                  #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{rad}
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            
            #DF       = 1 - (W_2/Wi) + (abs(Wti - Wt_2) / (2*R_Sol[STG]*Wi))
                        
            '''RELATIVE, PITCH-LINE'''
            WTR_2    = ((Rad_2 * Vt_2) / Radp) - (RAV * Radp)       #{m/s}
            WZR_2    = Wz_2                                         #{m/s}
            BETAR_2  = math.atan(WTR_2 / WZR_2)                     #{rad}
            BETAR_2d = BETAR_2 * (180/Pi)                           #{deg}
            WR_2     = WZR_2 / math.cos(BETAR_2)                    #{m/s}
            PTR_2    = PTR - (R_Phi[STG] * (PTR - Psi))             #{kPa}
            TTR_2    = TTR                                          #{K}
            
            #DF       = 1 - (WR_2/Wi) + (abs(Wti - W_2) / (2*R_Sol[STG]*Wi))
        
            '''STATIC CONVERSION, TIP'''
            H_2      = H + (DH[STG] * 1000)                         #{J/kg}
            h_2      = H_2 - (V_2**2 / 2)                           #{J/kg}
            #S_2      = S                                            #{J/kgK}
            S_2      = float(SPT(PTR_2, TTR_2))                     #{J/kgK}
            To_2     = float(TSH(S_2, H_2 ))                        #{K}
            Ts_2     = float(TSH(S_2, h_2))                         #{K}
            Ps_2     = float(PSH(S_2, h_2))                         #{kPa}
            Po_2     = float(PSH(S_2, H_2))                         #{kPa}   
            A_2      = float(APT(Ps_2, Ts_2))                       #{m/s}
            M_2      = V_2 / A_2                                    #{ }
            Rs_2     = float(RPT(Ps_2, Ts_2))                       #{kg/m3}
    
        '''RELATIVE, STATIC CONVERSION, PITCH-LINE'''
        hR_2     = H_2 - (WR_2**2 / 2)                          #{J/kg}
        Ts_2p    = float(TSH(S_2, hR_2))                        #{K}
        Ps_2p    = float(PSH(S_2, hR_2))                        #{kPa}  
        AR_2p    = float(APT(Ps_2p, Ts_2p))                       #{m/s}
        MR_2p    = WR_2 / AR_2p                                 #{ }
        Rs_2p    = float(RPT(Ps_2p, Ts_2p))                     #{kg/m3}
    
        '''BACK TO TOTAL'''
        Po_2     = float(PSH(S_2, H_2))                         #{kPa}
        A_comp_2 = MFR / (Rs_2 * Vz_2 * Lambda[STG])            #{m2}
    
        '''OVERALL RELATIVE'''
        TTR_2 = TTR_2                                           #{K}
        PTR_2 = PTR_2                                           #{kPa}
        Po_2  = Po_2                                            #{kPa}
        S_2p  = float(SPT(Po_2, TTR_2))                         #{J/kgK}   
    
        '''Radius Type'''
        Radh_2, Radp_2, Radt_2, Height_2 = Radius_Type(typer, Rad_2, A_comp_2, Pi)
    
        '''Constant Pitch Again'''
        Vz_2p    = Vz_2                                         #{m/s}
        U_2p     = RAV * Radp_2                                 #{m/s}
        Vt_2p    = Vt_2 * Rad_2 / Radp_2                        #{m/s}
        ALPHA_2p = math.atan(Vt_2p / Vz_2p)                     #{rad}
        ALPHA_2pd= ALPHA_2p * (180/Pi)                          #{deg}  
        V_2p     = Vz_2p / math.cos(ALPHA_2p)                   #{m/s}
        Wt_2p    = Vt_2p - U_2p                                 #{m/s}
        Wz_2p    = Vz_2p                                        #{m/s}
        BETA_2p  = math.atan(Wt_2p / Wz_2p)                     #{rad}
        BETA_2pd = BETA_2p * (180/Pi)                           #{deg}
        W_2p     = Wz_2p / math.cos(BETA_2p)                    #{m/s}
        M_2p     = V_2p / AR_2p                                 #Mach #
        h_2p     = H_2 - (V_2p**2/2)                            #{J/kg}
        
        '''Rotation Speed in Mach # '''
        ROTS     = U_2p                                         #{ }
        
        '''END'''
        Rm_2     = Radp_2                                       #{m}
        Vtheta_2 = Vt_2                                         #{m/s} 
        #Dht_2    = CP * DTT[STG]                               #NA
    
        if typer == 1: 
            #CONSTANT TIP VELO TRIANGLES
            U_2t, ALPHA_2t, V_2t, Vt_2t, Vz_2t, BETA_2t, W_2t, Wt_2t, Wz_2t, Po_2t, To_2t, Ts_2t, Ps_2t, M_2t, A_2t, Rs_2t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radt_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2t = MFR / (Rs_2t * Vz_2t)
            
        elif typer == 2: 
            #CONSTANT HUB VELO TRIANGLES
            U_2h, ALPHA_2h, V_2h, Vt_2h, Vz_2h, BETA_2h, W_2h, Wt_2h, Wz_2h, Po_2h, To_2h, Ts_2h, Ps_2h, M_2h, A_2h, Rs_2h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radh_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2h = MFR / (Rs_2h * Vz_2h)
            
            #CONSTANT TIP VELO TRIANGLES
            U_2t, ALPHA_2t, V_2t, Vt_2t, Vz_2t, BETA_2t, W_2t, Wt_2t, Wz_2t, Po_2t, To_2t, Ts_2t, Ps_2t, M_2t, A_2t, Rs_2t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radt_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2t = MFR / (Rs_2t * Vz_2t)
            
        elif typer == 3: 
            #CONSTANT HUB VELO TRIANGLES
            U_2h, ALPHA_2h, V_2h, Vt_2h, Vz_2h, BETA_2h, W_2h, Wt_2h, Wz_2h, Po_2h, To_2h, Ts_2h, Ps_2h, M_2h, A_2h, Rs_2h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radh_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2h = MFR / (Rs_2h * Vz_2h)
            #        
        #
    #
    
#####################################################
#####################################################
# PART D.ii - CALCULATIONS (PERFECT GAS)
#####################################################
#####################################################

    elif Fluid == 0:
    #
        #
        '''PERFECT GAS ESTIMATE'''
        
        '''It depends which section we keep constant'''
        if typer == 1 : 
            '''Consant Hub...'''
            #VARIABLES   
            #Hub Case
            #Toi, Radius, Vt_h, Vz_2
            
            #print(f'Initial CP     = {CP} J/kgK')
            To_2 = Toi + DTT[STG]
            #CP   = float(CpPT(Po_h, To_2))
            #print(f'Adjusted CP    = {CP} J/kgK')
            Rad_2    = Radius                                       #{m}
            U_2      = RAV * (Rad_2)                                #{m/s}
            Vt_2     = Vt_h + (CP * DTT[STG] / U_2)                 #{m/s}
            Vz_2     = VZR[STG] * Vzi                               #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{ }
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            Ts_2     = To_2 - (V_2**2 /(2 * CP))                    #{K}
        
            if Ts_2 <= 0:
                print('Static Temp Negative')
                print('Velocity is too high')
                Ts_2 = 0.7 * To_2
            
            A_2      = (Gamma * Rgas * Ts_2)**0.5                   #{m/s}
            M_2      = V_2 / A_2                                    #{ }
            
        elif typer == 2: 
            '''Consant Pitch...'''
            #VARIABLES 
            #Pitch Case
            #To_p, Radp, Vt_p, Vz_p
            
            #print(f'Initial CP     = {CP} J/kgK')
            To_2 = To_p + DTT[STG]
            #CP   = float(CpPT(Po_h, To_2))
            #print(f'Adjusted CP    = {CP} J/kgK')
            Rad_2    = Radp                                         #{m}
            U_2      = RAV * (Rad_2)                                #{m/s}
            Vt_2     = Vt_p + (CP * DTT[STG] / U_2)                 #{m/s}
            Vz_2     = VZR[STG] * Vz_p                              #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{ }
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            Ts_2     = To_2 - (V_2**2 /(2 * CP))                    #{K}
        
            if Ts_2 <= 0:
                print('Static Temp Negative')
                print('Velocity is too high')
                Ts_2 = 0.7 * To_2
            
            A_2      = (Gamma * Rgas * Ts_2)**0.5                   #{m/s}
            M_2      = V_2 / A_2                                    #{ }            
    
        elif typer == 3: 
            '''Consant Tip...'''
            #VARIABLES   
            #Tip Case
            #To_t, Radt, Vt_t, Vz_t
            
            if KKI == 0: 
                Rad_2    = Radt                                   #{m}
            elif KKI == 1: 
                if STG == (n-1): #Last Stage, Tip Radius Remains Const
                    Rad_2    = Rad_tip[STG]
                else:            #Not Last Stage, Tip Radius is Average of this and next
                    Rad_2    = (Rad_tip[STG] + Rad_tip[STG+1])/2
            
            #print(f'Initial CP     = {CP} J/kgK')
            To_2 = To_t + DTT[STG]
            #CP   = float(CpPT(Po_h, To_2))
            #print(f'Adjusted CP    = {CP} J/kgK')        
            U_2      = RAV * (Rad_2)                                #{m/s}
            Vt_2     = Vt_t + (CP * DTT[STG] / U_2)                 #{m/s}
            Vz_2     = VZR[STG] * Vz_t                              #{m/s}
            ALPHA_2  = math.atan(Vt_2 / Vz_2)                       #{ }
            ALPHA_2d = ALPHA_2 * (180/Pi)                           #{deg}
            V_2      = Vt_2 / (math.sin(ALPHA_2))                   #{m/s}
            Wt_2     = Vt_2 - U_2                                   #{m/s}
            Wz_2     = Vz_2                                         #{m/s}
            BETA_2   = math.atan(Wt_2 / Wz_2)                       #{ }
            BETA_2d  = BETA_2 * (180/Pi)                            #{deg}
            W_2      = Wz_2 / (math.cos(BETA_2))                    #{m/s}
            Ts_2     = To_2 - (V_2**2 /(2 * CP))                    #{K}
        
            if Ts_2 <= 0:
                print('Static Temp Negative')
                print('Velocity is too high')
                Ts_2 = 0.7 * To_2
            
            A_2      = (Gamma * Rgas * Ts_2)**0.5                   #{m/s}
            M_2      = V_2 / A_2                                    #{ }  
            
        #RELATIVE, CONST PITCH...
        WTR_2    = ((Rad_2 * Vt_2) / Radp) - (RAV * Radp)           #{m/s}
        WZR_2    = Wz_2                                             #{m/s}
        BETAR_2  = math.atan(WTR_2 / WZR_2)                         #{rad}
        BETAR_2d = BETAR_2 * (180/Pi)                               #{deg}
        WR_2     = WZR_2 / math.cos(BETAR_2)                        #{m/s}
        PTR_2    = PTR - (R_Phi[STG] * (PTR - Psi))                 #{kPa}
        TTR_2    = TTR                                              #{K}
        Ts_2p    = TTR_2 - (WR_2**2 / (2 * CP))                     #{K}
        AR_2     = (Gamma * Rgas * Ts_2p)**0.5                      #{m/s}
        MR_2p     = WR_2 / AR_2                                     #Mach #
        Ps_2p    = PTR_2 / ((TTR_2 / Ts_2p)**(Gamma / (Gamma - 1))) #{kPa}
        RsR_2    = Ps_2p / (Rgas / 1000 * Ts_2p)                    #{kg/m3}
        PoR_2    = Ps_2p * ((To_2 / Ts_2p)**(Gamma / (Gamma - 1)))  #{kPa}
        A_comp_2 = MFR / (RsR_2 * Vz_2 * Lambda[STG])               #{m^2}
    
    
        #Hub, Pitch, and Tip all have same Relative and total Pressure
        TTR_2 = TTR_2                                               #{K}
        PTR_2 = PTR_2                                               #{kPa}
        Po_2  = PoR_2                                               #{kPa}
        
        #Based on Radius Type
        if typer == 1: 
            #HUB
            Ps_2  = PTR_2 / ((TTR_2 / Ts_2)**(Gamma / (Gamma - 1))) #{kPa}
            Rs_2  = Ps_2 / (Rgas / 1000 * Ts_2)                     #{kg/m3}
            Rs_2p = Ps_2p/ (Rgas / 1000 * Ts_2p)                    #{kg/m3}
        elif typer == 2: 
            #PITCH
            Ps_2  = PTR_2 / ((TTR_2 / Ts_2)**(Gamma / (Gamma - 1))) #{kPa}
            Rs_2  = Ps_2 / (Rgas / 1000 * Ts_2)                     #{kg/m3}
            Rs_2p = Ps_2p/ (Rgas / 1000 * Ts_2p)                    #{kg/m3}
        elif typer == 3: 
            #TIP
            Ps_2  = PTR_2 / ((TTR_2 / Ts_2)**(Gamma / (Gamma - 1))) #{kPa}
            Rs_2  = Ps_2 / (Rgas / 1000 * Ts_2)                     #{kg/m3}
            Rs_2p = Ps_2p/ (Rgas / 1000 * Ts_2p)                    #{kg/m3}
    
        '''Radius Type'''
        Radh_2, Radp_2, Radt_2, Height_2 = Radius_Type(typer, Rad_2, A_comp_2, Pi)
    
        #Constant Pitch Again (I, 2, 2)
        Vz_2p    = Vz_2                                         #{m/s}
        U_2p     = RAV * Radp_2                                 #{m/s}
        Vt_2p    = Vt_2 * Rad_2 / Radp_2                        #{m/s}
        ALPHA_2p = math.atan(Vt_2p / Vz_2p)                     #{rad}
        ALPHA_2pd= ALPHA_2p * (180/Pi)                          #{deg}  
        V_2p     = Vz_2p / math.cos(ALPHA_2p)                   #{m/s}
        Wt_2p    = Vt_2p - U_2p                                 #{m/s}
        Wz_2p    = Vz_2p                                        #{m/s}
        BETA_2p  = math.atan(Wt_2p / Wz_2p)                     #{rad}
        BETA_2pd = BETA_2p * (180/Pi)                           #{deg}
        W_2p     = Wz_2p / math.cos(BETA_2p)                    #{m/s}
        M_2p     = V_2p / AR_2                                  #Mach #
        Rm_2     = Radp_2                                       #{m}
        Vtheta_2 = Vt_2                                         #{m/s}
        Dht_2    = CP * DTT[STG]                                #NA
        
        '''Rotation Speed in Mach # '''
        ROTS     = U_2p                                         #{m/s}
        H_2      = float(HPT(Po_2, To_2))                       #{J/kg}
        h_2      = H_2 - (V_2**2 / 2)                           #{J/kg}
        S_2      = S                                            #{J/kgK}
        hR_2     = H_2 - (WR_2**2 / 2)                          #{J/kg}
        h_2p     = H_2 - (W_2p**2/2)                            #{J/kg}
        
        if typer == 1: 
            #CONSTANT TIP VELO TRIANGLES
            U_2t, ALPHA_2t, V_2t, Vt_2t, Vz_2t, BETA_2t, W_2t, Wt_2t, Wz_2t, Po_2t, To_2t, Ts_2t, Ps_2t, M_2t, A_2t, Rs_2t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radt_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2t = MFR / (Rs_2t * Vz_2t)
            
        elif typer == 2: 
            #CONSTANT HUB VELO TRIANGLES
            U_2h, ALPHA_2h, V_2h, Vt_2h, Vz_2h, BETA_2h, W_2h, Wt_2h, Wz_2h, Po_2h, To_2h, Ts_2h, Ps_2h, M_2h, A_2h, Rs_2h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radh_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2h = MFR / (Rs_2h * Vz_2h)
            
            #CONSTANT TIP VELO TRIANGLES
            U_2t, ALPHA_2t, V_2t, Vt_2t, Vz_2t, BETA_2t, W_2t, Wt_2t, Wz_2t, Po_2t, To_2t, Ts_2t, Ps_2t, M_2t, A_2t, Rs_2t = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radt_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2t = MFR / (Rs_2t * Vz_2t)
            
        elif typer == 3: 
            #CONSTANT HUB VELO TRIANGLES
            U_2h, ALPHA_2h, V_2h, Vt_2h, Vz_2h, BETA_2h, W_2h, Wt_2h, Wz_2h, Po_2h, To_2h, Ts_2h, Ps_2h, M_2h, A_2h, Rs_2h = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_2, Vt_2, RAV, Radh_2, Rad_2, To_2, Po_2, Gamma, CP, Rgas, Fluid)
            A_comp_2h = MFR / (Rs_2h * Vz_2h)
            
        #
    #
    
#####################################################
#####################################################
# PART D.iii - MISCELLANEOUS CALCULATIONS 
#####################################################
#####################################################
    
    DF       = 1 - (W_2p/Wi_p) + (abs(Wti_p - Wt_2p) / (2*R_Sol[STG]*Wi_p))
    GammaValue2 = float(GPT(Ps_2p, Ts_2p))
    Viscosity2 = float(VPT(Ps_2p, Ts_2p))
    CCPP       = float(CpPT(Ps_2p, Ts_2p))
    
#####################################################
#####################################################
# PART D.iv - DISPLAY DATA
#####################################################
#####################################################
 
        #NOW WILL HAVE TO SAVE INFO TO CORRESPONDING EXPORTS
        #PITCH IS SAME EVERY TIME
        #BUT TIP AND HUB HAS TO BE "UPDATED" DEPENDING ON 
        #WHICH RADIUS IS HELD CONSTANT
        
    if typer == 1: 
        #Display Data in Window
        print(' ')
        print('**** PITCH-LINE CALCULATIONS ****')
        print(f'Radius                = {round(Radp_2, 3)} m')
        print(f'Experimental BETA     = {round(BETAR_2d, 6)} deg')
        #print(f'Total Entropy PITCH   = {round(S_2p, 3)} J/kgK')
        print(f'Total Temperature     = {round(To_2, 3)} K')
        print(f'Static Temperature    = {round(Ts_2p, 3)} K')
        print(f'Static Pressure       = {round(Ps_2p, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2, 3)} kPa')
        print(f'Mach #                = {round(M_2p, 3)} ')
        print(f'Total Velocity        = {round(V_2p, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2p, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2p, 3)} m/s')
        print(f'Area                  = {round(A_comp_2, 5)} m^2')
        print(' ')
        print('**** CONSTANT HUB CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radh_2, 3)} m')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(h_2, 3)} J/kg')
        print(f'Total Entropy         = {round(S_2, 3)} J/kgK')
        #print(f'Total Entropy   HUB   = {round(S_2h, 3)} J/kgK')
        print(f'Static Temperature    = {round(Ts_2, 3)} K')
        print(f'Total Temperature     = {round(To_2, 3)} K')
        print(f'Static Pressure       = {round(Ps_2, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2, 3)} kPa')
        print(f'Mach #                = {round(M_2, 3)} ')
        print(f'Total Velocity        = {round(V_2, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2, 3)} m/s')
        print(' ')
        print('**** TIP CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radt_2, 3)} m')
        print(f'Static Temperature    = {round(Ts_2t, 3)} K')
        print(f'Total Temperature     = {round(To_2t, 3)} K')
        print(f'Static Pressure       = {round(Ps_2t, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2t, 3)} kPa')
        print(f'Mach #                = {round(M_2t, 3)} ')
        print(f'Total Velocity        = {round(V_2t, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2t, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2t, 3)} m/s')
        print(' ')
        print('**** RELATIVE PROPERTIES ****')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(hR_2, 3)} J/kg')
        print(f'Relative Velocity     = {round(W_2, 3)} m/s')
        print(f'Relative Axial Vel    = {round(Wz_2, 3)} m/s')
        print(f'Velocity Theta        = {abs(round(Wt_2, 3))} m/s')
        print(f'Mach #                = {round(MR_2p, 3)} ')
        print(f'Relative Tot. Temp    = {round(TTR_2, 3)} K')
        print(f'Relative Tot. Press   = {round(PTR_2, 3)} kPa')
        print(f'Relative Stat. Temp   = {round(Ts_2p, 3)} K')
        print(f'Relative Stat. Press  = {round(Ps_2p, 3)} kPa ')
        print(f'Radius Relative       = {round(Radp, 3)} m')
        
        ''' Return properties to go into Stator '''
        return To_2, Ts_2, Ps_2p, Po_2, A_2, M_2, Radt_2, Radp_2, Rad_2, V_2, H_2, S_2, V_2p, Vz_2p, Vt_2p, ALPHA_2pd, BETAR_2d, Ps_2, Vz_2, Vt_2, ALPHA_2d, BETA_2d, Ts_2t, To_2t, Ps_2t, Po_2t, M_2t, V_2t, Vz_2t, Vt_2t, Rs_2t, ALPHA_2t, BETA_2t, Ts_2p, M_2p, Rs_2, Rs_2p, DF, h_2p, GammaValue2, Viscosity2, W_2p, Wz_2p, Wt_2p, W_2, Wz_2, Wt_2, W_2t, Wz_2t, Wt_2t, CCPP, Height_2, ROTS, BETA_2pd
    
    elif typer == 2: 
        #Display Data in Window
        print(' ')
        print('**** PITCH-LINE CALCULATIONS ****')
        print(f'Radius                = {round(Radp_2, 3)} m')
        print(f'Experimental BETA     = {round(BETAR_2d, 6)} deg')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(h_2, 3)} J/kg')
        print(f'Total Entropy         = {round(S_2, 3)} J/kgK')
        #print(f'Total Entropy PITCH   = {round(S_2p, 3)} J/kgK')
        print(f'Total Temperature     = {round(To_2, 3)} K')
        print(f'Static Temperature    = {round(Ts_2p, 3)} K')
        print(f'Static Pressure       = {round(Ps_2p, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2, 3)} kPa')
        print(f'Mach #                = {round(M_2p, 3)} ')
        print(f'Total Velocity        = {round(V_2p, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2p, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2p, 3)} m/s')
        print(f'Area                  = {round(A_comp_2, 5)} m^2')
        print(' ')
        print('**** CONSTANT HUB CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radh_2, 3)} m')
        #print(f'Total Entropy   HUB   = {round(S_2h, 3)} J/kgK')
        print(f'Static Temperature    = {round(Ts_2h, 3)} K')
        print(f'Total Temperature     = {round(To_2h, 3)} K')
        print(f'Static Pressure       = {round(Ps_2h, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2h, 3)} kPa')
        print(f'Mach #                = {round(M_2h, 3)} ')
        print(f'Total Velocity        = {round(V_2h, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2h, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2h, 3)} m/s')
        print(' ')
        print('**** TIP CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radt_2, 3)} m')
        print(f'Static Temperature    = {round(Ts_2t, 3)} K')
        print(f'Total Temperature     = {round(To_2t, 3)} K')
        print(f'Static Pressure       = {round(Ps_2t, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2t, 3)} kPa')
        print(f'Mach #                = {round(M_2t, 3)} ')
        print(f'Total Velocity        = {round(V_2t, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2t, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2t, 3)} m/s')
        print(' ')
        print('**** RELATIVE PROPERTIES ****')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(hR_2, 3)} J/kg')
        print(f'Relative Velocity     = {round(W_2, 3)} m/s')
        print(f'Relative Axial Vel    = {round(Wz_2, 3)} m/s')
        print(f'Velocity Theta        = {abs(round(Wt_2, 3))} m/s')
        print(f'Mach #                = {round(MR_2p, 3)} ')
        print(f'Relative Tot. Temp    = {round(TTR_2, 3)} K')
        print(f'Relative Tot. Press   = {round(PTR_2, 3)} kPa')
        print(f'Relative Stat. Temp   = {round(Ts_2p, 3)} K')
        print(f'Relative Stat. Press  = {round(Ps_2p, 3)} kPa ')
        print(f'Radius Relative       = {round(Radp, 3)} m')
        
        ''' Return properties to go into Stator '''
        return To_2, Ts_2h, Ps_2p, Po_2, A_2, M_2h, Radt_2, Radp_2, Radh_2, V_2h, H_2, S_2, V_2p, Vz_2p, Vt_2p, ALPHA_2pd, BETAR_2d, Ps_2h, Vz_2h, Vt_2h, ALPHA_2h, BETA_2h, Ts_2t, To_2t, Ps_2t, Po_2t, M_2t, V_2t, Vz_2t, Vt_2t, Rs_2t, ALPHA_2t, BETA_2t, Ts_2p, M_2p, Rs_2h, Rs_2p, DF, h_2p, GammaValue2, Viscosity2, W_2p, Wz_2p, Wt_2p, W_2h, Wz_2h, Wt_2h, W_2t, Wz_2t, Wt_2t, CCPP, Height_2, ROTS, BETA_2pd
    
    elif typer == 3: 
        #Display Data in Window
        print(' ')
        print('**** PITCH-LINE CALCULATIONS ****')
        print(f'Radius                = {round(Radp_2, 3)} m')
        print(f'Experimental BETA     = {round(BETAR_2d, 6)} deg')
        #print(f'Total Entropy PITCH   = {round(S_2p, 3)} J/kgK')
        print(f'Total Temperature     = {round(To_2, 3)} K')
        print(f'Static Temperature    = {round(Ts_2p, 3)} K')
        print(f'Static Pressure       = {round(Ps_2p, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2, 3)} kPa')
        print(f'Mach #                = {round(M_2p, 3)} ')
        print(f'Total Velocity        = {round(V_2p, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2p, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2p, 3)} m/s')
        print(f'Area                  = {round(A_comp_2, 5)} m^2')
        print(' ')
        print('**** CONSTANT HUB CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radh_2, 3)} m')
        #print(f'Total Entropy   HUB   = {round(S_2h, 3)} J/kgK')
        print(f'Static Temperature    = {round(Ts_2h, 3)} K')
        print(f'Total Temperature     = {round(To_2h, 3)} K')
        print(f'Static Pressure       = {round(Ps_2h, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2h, 3)} kPa')
        print(f'Mach #                = {round(M_2h, 3)} ')
        print(f'Total Velocity        = {round(V_2h, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2h, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2h, 3)} m/s')
        print(' ')
        print('**** TIP CALCULATIONS ****')
        print(f'Radius at Rotor Exit  = {round(Radt_2, 3)} m')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(h_2, 3)} J/kg')
        print(f'Total Entropy         = {round(S_2, 3)} J/kgK')
        print(f'Static Temperature    = {round(Ts_2, 3)} K')
        print(f'Total Temperature     = {round(To_2, 3)} K')
        print(f'Static Pressure       = {round(Ps_2, 3)} kPa')
        print(f'Total Pressure        = {round(Po_2, 3)} kPa')
        print(f'Mach #                = {round(M_2, 3)} ')
        print(f'Total Velocity        = {round(V_2, 3)} m/s')
        print(f'Axial Velocity        = {round(Vz_2, 3)} m/s')
        print(f'Velocity Theta        = {round(Vt_2, 3)} m/s')
        print(' ')
        print('**** RELATIVE PROPERTIES ****')
        print(f'Total Enthalpy        = {round(H_2, 3)} J/kg')
        print(f'Static Enthalpy       = {round(hR_2, 3)} J/kg')
        print(f'Relative Velocity     = {round(W_2, 3)} m/s')
        print(f'Relative Axial Vel    = {round(Wz_2, 3)} m/s')
        print(f'Velocity Theta        = {abs(round(Wt_2, 3))} m/s')
        print(f'Mach #                = {round(MR_2p, 3)} ')
        print(f'Relative Tot. Temp    = {round(TTR_2, 3)} K')
        print(f'Relative Tot. Press   = {round(PTR_2, 3)} kPa')
        print(f'Relative Stat. Temp   = {round(Ts_2p, 3)} K')
        print(f'Relative Stat. Press  = {round(Ps_2p, 3)} kPa ')
        print(f'Radius Relative       = {round(Radp, 3)} m')
        
        ''' Return properties to go into Stator '''
        return To_2h, Ts_2h, Ps_2p, Po_2h, A_2, M_2h, Radt_2, Radp_2, Radh_2, V_2h, H_2, S_2, V_2p, Vz_2p, Vt_2p, ALPHA_2pd, BETAR_2d, Ps_2h, Vz_2h, Vt_2h, ALPHA_2h, BETA_2h, Ts_2, To_2, Ps_2, Po_2, M_2, V_2, Vz_2, Vt_2, Rs_2, ALPHA_2d, BETA_2d, Ts_2p, M_2p, Rs_2h, Rs_2p, DF, h_2p, GammaValue2, Viscosity2, W_2p, Wz_2p, Wt_2p, W_2h, Wz_2h, Wt_2h, W_2, Wz_2, Wt_2, CCPP, Height_2, ROTS, BETA_2pd
    #
#

#####################################################
#####################################################
# PART E - STATOR OUTLET CALCULATIONS FUNCTION
#####################################################
#####################################################

def Stator_Outlet_Calculations(n, alphaend, Mend, STG, alpha1, MACH1, Radius, To_2, H_2, Gamma, S_Phi, Po_2, Ps_2, Rgas, MFR, Lambda, RAV, CP, Bleed, Poi, HPT, SPT, APT, RPT, GPT, CpPT, TSH, PSH, Fluid, V_2p, Vt_2p, S_Sol, VPT, typer, KKI, Rad_tip):
#
    #
    Pi = 3.1415926

    print(' ')
    
#####################################################
#####################################################
# PART E.i - CHECK IF LAST STAGE
#####################################################
#####################################################
    
    if STG == (n - 1) :
        ALPHA_3d = alphaend                                     #{deg}
        M_3      = Mend                                         #{Mach}
        
        if KKI == 0: 
            Rad_3    = Radius                                   #{m}
        elif KKI == 1: 
            Rad_3    = Rad_tip[STG]
        
    else : 
        ALPHA_3d = alpha1[STG + 1] #Because it's the inlet for the next rotor
        M_3      = MACH1[STG + 1] 
        
        if KKI == 0: 
            Rad_3    = Radius                                   #{m}
        elif KKI == 1: 
            Rad_3    = Rad_tip[STG+1]
    
    ALPHA_3  = ALPHA_3d * (Pi / 180)    #{rad}
    
#####################################################
#####################################################
# PART E.ii - CALCULATIONS (REAL GAS)
#####################################################
#####################################################
    '''REAL GAS CALCS -- PITCH LINE'''
    if Fluid == 1:
        H_3      = H_2                                      #{J/kg}
        Po_3     = Po_2 - (S_Phi[STG] * (Po_2 - Ps_2))      #{kPa}
        To_3i    = To_2                                     #{K}
        resi     = 10 #residual
        nnn      = 0
        
        #We have the Same Enthalpy and New Pressure -- Find Temp that corresponds with their values (Real Gas)
        while abs(resi) >= 1:
            if resi >= 0: 
                To_3i = To_3i + 0.01
                Hi    = float(HPT(Po_3, To_3i))
                resi  = H_3 - Hi
            elif resi < 0:
                To_3i = To_3i - 0.001
                Hi    = float(HPT(Po_3, To_3i))
                resi  = H_3 - Hi 
            nnn  += 1
            if nnn >= 5000:
                break
    
        To_3     = To_3i    
        S_3      = float(SPT(Po_3, To_3))
    
    #
        #
        '''OBTAIN ACTUAL STATIC PROPERTIES'''
        Ts_3, Ps_3, Rs_3, A_3, V_3, M_3, h_3, S_3NA, H_3NA = STATIC_CONV(HPT, SPT, APT, RPT, TSH, PSH, Po_3, To_3, M_3)
        Vz_3     = V_3 * math.cos(ALPHA_3)
        Vt_3     = V_3 * math.sin(ALPHA_3)
        A_comp_3 = MFR / (Rs_3 * Vz_3 * Lambda[STG])
        #
    #
    
#####################################################
#####################################################
# PART E.iii - CALCULATIONS (PERFECT GAS)
#####################################################
#####################################################
    
    elif Fluid == 0:
    #
        #
        '''PERFECT GAS ESTIMATES -- PITCH LINE'''
        H_3      = H_2                                          #{J/kg}
        To_3     = To_2                                         #{K}
        Ts_3     = To_3 / (1 + (((Gamma - 1) / 2) * M_3**2))    #{K}
        Po_3     = Po_2 - (S_Phi[STG] * (Po_2 - Ps_2))          #{kpa}
        Ps_3     = Po_3 / ((To_3 / Ts_3)**(Gamma / (Gamma - 1)))  #{kPa}
        A_3      = (Gamma * Rgas * Ts_3)**0.5                   #{m/s}
        V_3      = M_3 * A_3                                    #{m/s}
        Vz_3     = V_3 * math.cos(ALPHA_3)                      #{m/s}
        Vt_3     = V_3 * math.sin(ALPHA_3)                      #{m/s}
        Rs_3     = Ps_3 / (Rgas / 1000 * Ts_3)                  #{kg/m3}
        Ro_3     = Po_3 / (Rgas / 1000 * To_3)                  #{kg/m3}
        A_comp_3 = MFR / (Rs_3 * Vz_3 * Lambda[STG])            #{m^2}
        S_3      = float(SPT(Po_3, To_3))                       #{J/kgK}
        h_3      = H_3 - (V_3**2 / 2)                           #{J/kg}
        #
    #
    
    DF = 1 - (V_3/V_2p) + (abs(Vt_2p - Vt_3) / (2*S_Sol[STG]*V_2p))  #{Diffusion Factor}
    
    '''Radius Type'''
    Radh_3, Radp_3, Radt_3, Height_3 = Radius_Type(typer, Rad_3, A_comp_3, Pi)
    
    #
    U_3      = RAV * Radp_3                 #{m/s}
    Wt_3     = Vt_3 - U_3                   #{m/s}
    BETA_3   = math.atan(Wt_3 / Vz_3)       #{rad}
    BETA_3d  = BETA_3 * (180 / Pi)          #deg
    W_3      = Vz_3 / math.cos(BETA_3)      #{m/s}
    Wz_3     = Vz_3                         #{m/s}
    
    if Fluid == 1:
    #
        #
        '''OBTAIN ACTUAL RELATIVE TOTAL PROPERTIES'''
        H_rel3, S_rel3, TTR_3, PTR_3, MR_3 = TOTAL_CONV(HPT, SPT, APT, RPT, TSH, PSH, h_3, S_3, Ts_3, W_3, A_3)
        #
    #
    elif Fluid == 0:
    #
        #
        TTR_3    = Ts_3 + (W_3**2 / (2 * CP))
        PTR_3    = Ps_3 * (TTR_3 / Ts_3)**(Gamma / (Gamma - 1))
    
    ######### NOT NEEDED ###################
    #Relative Mach, Total Temp, Total Press#
    MR_3 = W_3 / A_3                       #
    # TTR_3                                #
    # PTR_3                                #
    ########################################
    
#####################################################
#####################################################
# PART E.iv - MISCELLANEOUS CALCULATIONS
#####################################################
#####################################################
    
    MFR      = MFR * (1 - Bleed[STG] / 100)
    Rm_3     = Rad_3
    Vtheta_3 = Vt_3
    Dht_3    = 0.0
    
    PR       = Po_3 / Poi
    GammaValue = float(GPT(Ps_3, Ts_3))
    Viscosity = float(VPT(Ps_3, Ts_3))
    CCPP      = float(CpPT(Ps_3, Ts_3))
    
    #CONSTANT HUB VELO TRIANGLES
    U_h3, ALPHA_h3, V_h3, Vt_h3, Vz_h3, BETA_h3, W_h3, Wt_h3, Wz_h3, Po_h3, To_h3, Ts_h3, Ps_h3, M_h3, A_h3, Rs_h3 = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_3, Vt_3, RAV, Radh_3, Radp_3, To_3, Po_3, Gamma, CP, Rgas, Fluid)
    A_comp_h3 = MFR / (Rs_h3 * Vz_h3)
    TTR_h3    = Ts_h3 + (W_h3**2 / (2 * CP))
    PTR_h3    = Ps_h3 * (TTR_h3 / Ts_h3)**(Gamma / (Gamma - 1))
        
    #CONSTANT TIP VELO TRIANGLES
    U_t3, ALPHA_t3, V_t3, Vt_t3, Vz_t3, BETA_t3, W_t3, Wt_t3, Wz_t3, Po_t3, To_t3, Ts_t3, Ps_t3, M_t3, A_t3, Rs_t3 = VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, Vz_3, Vt_3, RAV, Radt_3, Radp_3, To_3, Po_3, Gamma, CP, Rgas, Fluid)
    A_comp_t3 = MFR / (Rs_t3 * Vz_t3)
    
#####################################################
#####################################################
# PART E.v - DISPLAY DATA
#####################################################
#####################################################    

    print('**** OUTLET OF STAGE ****')
    print('**** PITCH-LINE CALCULATIONS ***')
    print(f'Radius                = {round(Radp_3, 3)} m')
    print(f'Static Temperature    = {round(Ts_3, 3)} K ')
    print(f'Static Pressure       = {round(Ps_3, 3)} kPa ')
    print(f'Static Density        = {round(Rs_3, 3)} kg/m^3')
    print(f'Total Temperature     = {round(To_3, 3)} K ')
    print(f'Total Pressure        = {round(Po_3, 3)} kPa ')
    print(f'Mach #                = {round(M_3, 3)} ')
    print(f'Total Velocity        = {round(V_3, 3)} m/s')
    print(f'Axial Velocity        = {round(Vz_3, 3)} m/s')
    print(f'Velocity Theta        = {round(Vt_3, 3)} m/s')
    print(f'Relative Velocity     = {round(W_3, 3)} m/s')
    print(f'Area at Outlet        = {round(A_comp_3, 3)} m^2')
    print(f'New Mass Flow Rate    = {round(MFR, 3)} kg/s')
    print(f'Stage Pressure Ratio  = {round(PR, 3)} ')
    print(' ')
    print('**** CONSTANT HUB CALCULATIONS ****')
    print(f'Radius                = {round(Rad_3, 3)} m')
    print(f'Static Temperature    = {round(Ts_h3, 3)} K ')
    print(f'Static Pressure       = {round(Ps_h3, 3)} kPa ')
    print(f'Static Density        = {round(Rs_h3, 3)} kg/m^3')
    print(f'Total Temperature     = {round(To_h3, 3)} K ')
    print(f'Total Pressure        = {round(Po_h3, 3)} kPa ')
    print(f'Mach #                = {round(M_h3, 3)} ')
    print(f'Total Velocity        = {round(V_h3, 3)} m/s')
    print(f'Axial Velocity        = {round(Vz_h3, 3)} m/s')
    print(f'Velocity Theta        = {round(Vt_h3, 3)} m/s')
    print(f'Relative Velocity     = {round(W_h3, 3)} m/s')
    print(f'Area at Outlet        = {round(A_comp_h3, 3)} m^2')
    print(' ')
    print('**** TIP CALCULATIONS ****')
    print(f'Radius                = {round(Radt_3, 3)} m')
    print(f'Static Temperature    = {round(Ts_t3, 3)} K ')
    print(f'Static Pressure       = {round(Ps_t3, 3)} kPa ')
    print(f'Static Density        = {round(Rs_t3, 3)} kg/m^3')
    print(f'Total Temperature     = {round(To_t3, 3)} K ')
    print(f'Total Pressure        = {round(Po_t3, 3)} kPa ')
    print(f'Mach #                = {round(M_t3, 3)} ')
    print(f'Total Velocity        = {round(V_t3, 3)} m/s')
    print(f'Axial Velocity        = {round(Vz_t3, 3)} m/s')
    print(f'Velocity Theta        = {round(Vt_t3, 3)} m/s')
    print(f'Relative Velocity     = {round(W_t3, 3)} m/s')
    print(f'Area at Outlet        = {round(A_comp_t3, 3)} m^2')
    print(' ')
    print('**** RELATIVE PROPERTIES ****')
    print(f'Total Pressure        = {round(PTR_3, 3)} kPa')
    print(f'Total Temperature     = {round(TTR_3, 3)} K')
    print(f'Axial Velocity        = {round(Wz_3, 3)} m/s')
    print(f'Tangential Velocity   = {round(Wt_3, 3)} m/s')
    print(f'Mach #                = {round(MR_3, 3)} ')
    
    '''
        Return properties to go into next Rotor
    '''
    return Rad_3, Radp_3, To_3, Po_3, Ts_h3, Ps_3, Vt_h3, Vz_h3, PTR_3, TTR_3, MFR, H_3, S_3, PR, Po_h3, Rs_3, V_3, Vz_3, Vt_3, ALPHA_3d, BETA_3d, Radh_3, Ts_h3, Ps_h3, Rs_h3, To_h3, Po_h3, M_h3, V_h3, Vz_h3, Vt_h3, ALPHA_h3, BETA_h3, Radt_3, Ts_t3, Ps_t3, Rs_t3, To_t3, Po_t3, M_t3, V_t3, Vz_t3, Vt_t3, ALPHA_t3, BETA_t3, Rs_3, Ts_3, M_3, U_3, W_h3, Wt_h3, Wz_h3, DF, h_3, GammaValue, Viscosity, W_3, Wt_3, Wz_3, W_t3, Wt_t3, Wz_t3, CCPP, Height_3
    #
#

#####################################################
#####################################################
# PART F - RADIUS TYPE FUNCTION
#####################################################
##################################################### 

def Radius_Type(typer, Radius, Area, Pi): 
    if typer == 1 :
    #
        #
        '''Constant Hub'''
        Radh   = Radius
        Radt   = ((Area / Pi) + (Radius**2))**0.5
        Height = Radt - Radius
        Radp   = Radius + (Height / 2)
        #
    #
    elif typer == 2 : 
    #
        #
        '''Constant Pitch'''
            
        Radp   = Radius
        Height = Area / (2 * Pi * Radp)
        Radh   = Radp - (Height / 2)
        if Radh < 0 :
        #
            #
            print('ERROR - Hub Radius is less than zero for this flow rate, cannot continue.')
            print('Set a larger radius or run with constant hub radius.')
            sys.exit()
            #
        #
        Radt  = Radp + (Height / 2)
    #
    elif typer == 3 : 
    #
        #
        '''Constant Tip'''
            
        Radt   = Radius
        rterm      = Radt**2 - (Area / Pi)
        if rterm < 0 :
        #
            #
            print('ERROR - Tip radius is too small for this flow and condition, cannot continue.')
            print('Set a larger radius or run with constant hub radius.')
            sys.exit()
            #
        #
        Radh   = rterm**0.5
        Height = Radt - Radh
        Radp   = Radt - (Height / 2)
        #
    #        
    return Radh, Radp, Radt, Height
    #
#

#####################################################
#####################################################
# PART G - VELOCITY TRIANGLES FUNCTION
#####################################################
##################################################### 

def VELO_TRIANGLES(HPT, SPT, RPT, APT, TSH, PSH, VZ, VT, RAV, R, RI, To, Po, Gamma, CP, Rgas, Fluid):
    Pi = 3.1415926
    
    '''
        FREE VORTEX
    '''
    VZ_2    = VZ                                #{m/s}
    VT_2    = VT * RI / R                       #{m/s}
    U       = RAV * R                           #{m/s}
    ALPHA   = math.atan(VT_2 / VZ_2)            #{rad}
    ALPHA_d = ALPHA * (180/Pi)                  #{deg}
    V_2     = VZ / math.cos(ALPHA)              #{m/s}
    WT_2    = VT_2 - U                          #{m/s}
    WZ_2    = VZ_2                              #{m/s}
    BETA    = math.atan(WT_2 / WZ_2)            #{rad}
    BETA_d  = BETA * (180/Pi)                   #{deg}
    W       = WZ_2 / math.cos(BETA)             #{m/s}
    To_2    = To                                #{K}
    Po_2    = Po                                #{kPa}
    
    if Fluid == 1: 
        '''STATIC CONVERSION'''
        H_2     = float(HPT(Po_2, To_2))        #{J/kg}
        S_2     = float(SPT(Po_2, To_2))        #{J/kgK}
        h_2     = H_2 - (V_2**2 / 2)            #{J/kg}
        Ts_2    = float(TSH(S_2, h_2))          #{K}
        Ps_2    = float(PSH(S_2, h_2))          #{kPa}
        A_2     = float(APT(Ps_2, Ts_2))        #{m/s}
        Rs_2    = float(RPT(Ps_2, Ts_2))        #{kg/m3}
        M_2     = V_2 / A_2                     #Mach #
        
    if Fluid == 0: 
        Ts_2    = To_2 - V_2**2/(2*CP)              #{K}
        A_2     = (Gamma * Rgas * Ts_2)**0.5        #{m/s}
        M_2     = V_2 / A_2                         #{ }
        Ps_2    = Po_2 / ((1 + (((Gamma - 1) / 2) * M_2**2))**(Gamma / (Gamma - 1))) #{kPa}
        Rs_2    = Ps_2 / (Rgas / 1000 * Ts_2)       #{kg/m3}
        H_2     = float(HPT(Po_2, To_2))
        S_2     = float(SPT(Po_2, To_2))
        h_2     = H_2 - (V_2**2 / 2)
        
    return U, ALPHA_d, V_2, VT_2, VZ_2, BETA_d, W, WT_2, WZ_2, Po_2, To_2, Ts_2, Ps_2, M_2, A_2, Rs_2

#####################################################
#####################################################
# PART H - STATIC PROPERTY CONVERSION FUNCTION
#####################################################
##################################################### 

def STATIC_CONV(HPT, SPT, APT, RPT, TSH, PSH, Po, To, M):
    
    '''TAKES TOTAL PRESSURE, TOTAL TEMPERATURE, AND THE DESIRED
        MACH... THEN ITERATES TO FIND THE VELOCITY AND CONSEQUENTLY
        STATIC PROPERTIES...
        
        USING BI-SECTION METHOD'''
    
    '''REAL sCO2 PROPERTIES'''
    H = float(HPT(Po, To))
    S = float(SPT(Po, To))
    
    '''FIRST ESTIMATE V'''
    A          = float(APT(Po, To))
    V_0        = M * A
    residual   = 1.0
    n          = 0
    a          = 0
    b          = 2 * V_0
    
    '''LOWER BOUND -- Mach(Ho, So, a)'''
    h          = H - (a**2 / 2)
    Ts         = float(TSH(S, h))
    Ps         = float(PSH(S, h))
    Rs         = float(RPT(Ps, Ts))
    As         = float(APT(Ps, Ts))
    MachA      = a / As
    residualA  = M - MachA
    SignA      = numpy.sign(residualA)
    
    '''UPPER BOUND -- Mach(Ho, So, b)'''
    h          = H - (b**2 / 2)
    Ts         = float(TSH(S, h))
    Ps         = float(PSH(S, h))
    Rs         = float(RPT(Ps, Ts))
    As         = float(APT(Ps, Ts))
    MachB      = b / As
    residualB  = M - MachB
    SignB      = numpy.sign(residualB)
    
    while residual >= 1e-10 :
    #
        #
        n += 1
        m  = (a + b) / 2
        
        '''Mach(Ho, So, m)'''
        h          = H - (m**2 / 2)  
        Ts         = float(TSH(S, h))
        Ps         = float(PSH(S, h))
        Rs         = float(RPT(Ps, Ts))
        As         = float(APT(Ps, Ts))
        MachM      = m / As
        residualM  = M - MachM
        SignM      = numpy.sign(residualM)
        
        if SignM != SignA : 
            b     = m
            MachB = MachM
        else:
            a     = m
            MachA = MachM
        residual  = abs(residualM)
        if n > 1000:
            print('UNABLE TO CONVERGE')
            break
            #
        #
    #
    Ve = m
    Me = MachM
    return Ts, Ps, Rs, As, Ve, Me, h, S, H
    #
#

#####################################################
#####################################################
# PART I - TOTAL PROPERTY CONVERSION FUNCTION
#####################################################
##################################################### 

def TOTAL_CONV(HPT, SPT, APT, RPT, TSH, PSH, h, s, Ti, V, A):
    
    '''
        TAKES STATIC ENTHALPY, STATIC TEMPERATURE, AND THE VELOCITY... 
        THEN OBTAINS THE TOTAL PROPERTIES...
    '''
    
    H      = h + (V**2 / 2)
    S      = s 
    To     = float(TSH(S, H))
    Po     = float(PSH(S, H))
    M      = V / A
    
    return H, S, To, Po, M
    #
#

#####################################################
#####################################################
# PART J - WALLS FUNCTION
#####################################################
##################################################### 

def WALLS(XCOORD, YCOORD, refl, IGV_Q, Stage_Number, XRIF, Radh, Radt, Output, EXTRAPOLATION): 
    XNEW_IGV = []
    YNEW_IGV = []
    XNEW = []
    YNEW = []
    ZHUB = []
    RHUB = []
    ZTIP = []
    RTIP = []
    n = int(max(Stage_Number)) - 1
    
    #Finding the Extra Coordinates for IGV
    #(3 New Points for IGV)
    if IGV_Q == 1: 
        #9
        Xi, Yi = EXTRAPOLATION(XCOORD[0], XCOORD[2], XCOORD[1], XCOORD[6], YCOORD[0], YCOORD[2], YCOORD[1], YCOORD[6])
        XNEW_IGV.append(Xi)  #{mm}    
        YNEW_IGV.append(Yi)  #{mm}
        
        #10
        Xi, Yi = EXTRAPOLATION(XCOORD[3], XCOORD[4], XCOORD[1], XCOORD[6], YCOORD[3], YCOORD[4], YCOORD[1], YCOORD[6])
        XNEW_IGV.append(Xi)  #{mm}    
        YNEW_IGV.append(Yi)  #{mm}
        
        #11
        Xi, Yi = EXTRAPOLATION(XCOORD[3], XCOORD[4], XCOORD[0], XCOORD[8 + 1], YCOORD[3], YCOORD[4], YCOORD[0], YCOORD[8 + 1])
        XNEW_IGV.append(Xi)  #{mm}    
        YNEW_IGV.append(Yi)  #{mm}
        
    elif IGV_Q == 0:
        #No need for the extra IGV coordinates
        indent = 1
    
    #Finding the Extra Coordinates for Each Stage
    #(5 New Points for each stage)
    for j in Stage_Number:
        if IGV_Q == 1: 
            tt = 8
        elif IGV_Q == 0:
            tt = 0
        i = int(j) - 1

        if i < n-1: #Every stage before the last
            
            #Update the 5th Y coordinate for every stage
            YCOORD[10*i + 4 + tt] = YCOORD[10*i + 1 + tt] + (XCOORD[10*i + 4 + tt] - XCOORD[10*i + 1 + tt]) / (XCOORD[10*(i+1) + 1 + tt] - XCOORD[10*i + 1 + tt]) * (YCOORD[10*(i+1) + 1 + tt] - YCOORD[10*i + 1 + tt])
        
            #Getting new points 11-15 for each stage
            #11
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 2  + tt], XCOORD[10*i + 3  + tt], XCOORD[10*i + 1  + tt], XCOORD[10*(i+1) + 1  + tt], YCOORD[10*i + 2  + tt], YCOORD[10*i + 3  + tt], YCOORD[10*i + 1  + tt], YCOORD[10*(i+1) + 1  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #12
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 2  + tt], XCOORD[10*i + 3  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9  + tt], YCOORD[10*i + 2  + tt], YCOORD[10*i + 3  + tt], YCOORD[10*i + 0  + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #13
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 4  + tt], XCOORD[10*i + 5  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9  + tt], YCOORD[10*i + 4  + tt], YCOORD[10*i + 5  + tt], YCOORD[10*i + 0  + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #14
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 6  + tt], XCOORD[10*i + 7  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9 + tt], YCOORD[10*i + 6  + tt], YCOORD[10*i + 7  + tt], YCOORD[10*i + 0  + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #15
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 4  + tt], XCOORD[10*(i+1) + 1 + tt], XCOORD[10*i + 7  + tt], XCOORD[10*i + 6  + tt], YCOORD[10*i + 4  + tt], YCOORD[10*(i+1) + 1 + tt], YCOORD[10*i + 7 + tt], YCOORD[10*i + 6 + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
        else: 
            #Update the 5th Y coordinate for every stage
            YCOORD[10*i + 4 + tt] = YCOORD[10*i + 1 + tt] + (XCOORD[10*i + 4 + tt] - XCOORD[10*i + 1 + tt]) / (XCOORD[10*i + 8 + tt] - XCOORD[10*i + 1 + tt]) * (YCOORD[10*i + 8 + tt] - YCOORD[10*i + 1 + tt])
        
            #Getting new points 11-15 for each stage
            #11
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 2  + tt], XCOORD[10*i + 3  + tt], XCOORD[10*i + 1  + tt], XCOORD[10*i + 8 + tt], YCOORD[10*i + 2  + tt], YCOORD[10*i + 3  + tt], YCOORD[10*i + 1  + tt], YCOORD[10*i + 8 + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #12
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 2  + tt], XCOORD[10*i + 3  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9  + tt], YCOORD[10*i + 2  + tt], YCOORD[10*i + 3  + tt], YCOORD[10*i + 0  + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #13
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 4  + tt], XCOORD[10*i + 5  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9  + tt], YCOORD[10*i + 4  + tt], YCOORD[10*i + 5  + tt], YCOORD[10*i + 0  + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #14
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 6  + tt], XCOORD[10*i + 7  + tt], XCOORD[10*i + 0  + tt], XCOORD[10*i + 9 + tt], YCOORD[10*i + 6  + tt], YCOORD[10*i + 7  + tt], YCOORD[10*i + 0 + tt], YCOORD[10*i + 9  + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
            #15
            Xi, Yi = EXTRAPOLATION(XCOORD[10*i + 4  + tt], XCOORD[10*i + 8 + tt], XCOORD[10*i + 7  + tt], XCOORD[10*i + 6  + tt], YCOORD[10*i + 4  + tt], YCOORD[10*i + 8 + tt], YCOORD[10*i + 7 + tt], YCOORD[10*i + 6 + tt])
            XNEW.append(Xi)  #{mm}    
            YNEW.append(Yi)  #{mm}
            
        if i == n:
            break
            
    if IGV_Q == 1: 
        #IGV Coords
        NPTS = 0
        ZHUB.append(XCOORD[0] / 1000 / refl)
        RHUB.append(YCOORD[1] / 1000 / refl)
        ZTIP.append(XCOORD[0] / 1000 / refl)
        RTIP.append(YCOORD[0] / 1000 / refl)

        NPTS += 1
        ZHUB.append(XNEW_IGV[0] / 1000 / refl)
        RHUB.append(YNEW_IGV[0] / 1000 / refl)
        ZTIP.append(XCOORD[0] / 1000 / refl)
        RTIP.append(YCOORD[0] / 1000 / refl)        

        NPTS += 1
        ZHUB.append(XNEW_IGV[1] / 1000 / refl)
        RHUB.append(YNEW_IGV[1] / 1000 / refl)
        ZTIP.append(XNEW_IGV[2] / 1000 / refl)
        RTIP.append(YNEW_IGV[2] / 1000 / refl)
        
    elif IGV_Q == 0: 
        NPTS = 0
        #NON Existent IGV Coords
        ZHUB.append(-XNEW[1]/2 / 1000 / refl)
        RHUB.append(YCOORD[0 + tt] / 1000 / refl)
        ZTIP.append(-XNEW[1]/2 / 1000 / refl)
        RTIP.append(YCOORD[1 + tt] / 1000 / refl)        
    
    for j in Stage_Number:
        if IGV_Q == 1: 
            tt = 8
        elif IGV_Q == 0:
            tt = 0
        i = int(j) - 1
        
        NPTS += 1
        ZHUB.append(XCOORD[10*i + 0  + tt] / 1000 / refl)
        RHUB.append(YCOORD[10*i + 0  + tt] / 1000 / refl)
        ZTIP.append(XCOORD[10*i + 1  + tt] / 1000 / refl)
        RTIP.append(YCOORD[10*i + 1  + tt] / 1000 / refl)

        NPTS += 1
        ZHUB.append(XNEW[5*i + 1] / 1000 / refl)
        RHUB.append(YNEW[5*i + 1] / 1000 / refl)
        ZTIP.append(XNEW[5*i + 0] / 1000 / refl)
        RTIP.append(YNEW[5*i + 0] / 1000 / refl)        

        NPTS += 1
        ZHUB.append(XNEW[5*i + 2] / 1000 / refl)
        RHUB.append(YNEW[5*i + 2] / 1000 / refl)
        ZTIP.append(XCOORD[10*i + 4 + tt] / 1000 / refl)
        RTIP.append(YCOORD[10*i + 4 + tt] / 1000 / refl) 
        
        NPTS += 1
        ZHUB.append(XNEW[5*i + 3] / 1000 / refl)
        RHUB.append(YNEW[5*i + 3] / 1000 / refl)
        ZTIP.append(XNEW[5*i + 4] / 1000 / refl)
        RTIP.append(YNEW[5*i + 4] / 1000 / refl) 
        if i == n: 
            break
        
    NPTS += 1    
    ZHUB.append(XRIF / 1000 / refl)
    RHUB.append(Radh / refl)
    ZTIP.append(XRIF / 1000 / refl)
    RTIP.append(Radt / refl)
    
    if IGV_Q == 1: 
        ii = (i+1)*2 + 1
    if IGV_Q == 0: 
        ii = (i+1)*2
        
    ''' Write to Wall File '''
    FW = open(f'{Output}.twal', "w+")
    print(f' There are {NPTS+1} points defining the flow path with {ii} blade rows. \n')
    FW.write(f'PY-C-DES output, extension {Output} \n')
    
    #Write out HUB values
    FW.write(f'\t{round(ZHUB[0] - .15, 6)} \t{round(ZHUB[NPTS] + .15, 6)} \t0.000000 \t0.000000 \n')
    FW.write(f'\t{round(ZHUB[NPTS] + 0.30, 6)} \t{round(RHUB[NPTS], 6)} \n')
    FW.write(f'\t{round(ZHUB[NPTS] + 0.20, 6)} \t{round(RHUB[NPTS], 6)} \n')    
    FW.write(f'\t{round(ZHUB[NPTS] + 0.10, 6)} \t{round(RHUB[NPTS], 6)} \n')
    FW.write(f'\t{round(ZHUB[NPTS] + 0.05, 6)} \t{round(RHUB[NPTS], 6)} \n')
    
    for i in range(0, NPTS+1, 1): #Going 'backwards' through the Hub Vectors
        j = int(NPTS - i) 
        FW.write(f'\t{round(ZHUB[j], 6)} \t{round(RHUB[j], 6)} \n')

    FW.write(f'\t{round(ZHUB[0] - 0.05, 6)} \t{round(RHUB[0], 6)} \n')
    FW.write(f'\t{round(ZHUB[0] - 0.10, 6)} \t{round(RHUB[0], 6)} \n')    
    FW.write(f'\t{round(ZHUB[0] - 0.20, 6)} \t{round(RHUB[0], 6)} \n')
    FW.write(f'\t{round(ZHUB[0] - 0.30, 6)} \t{round(RHUB[0], 6)} \n')
    FW.write(f'\t{999.000} \t{999.000} \n')
    
    #Write out TIP values
    FW.write(f'\t{round(ZTIP[0] - 0.30, 6)} \t{round(RTIP[0], 6)} \n')
    FW.write(f'\t{round(ZTIP[0] - 0.20, 6)} \t{round(RTIP[0], 6)} \n')    
    FW.write(f'\t{round(ZTIP[0] - 0.10, 6)} \t{round(RTIP[0], 6)} \n')
    FW.write(f'\t{round(ZTIP[0] - 0.05, 6)} \t{round(RTIP[0], 6)} \n')
    
    for i in range(0, NPTS+1, 1): #Going 'forwards' through the Tip Vectors
        j = int(i)
        FW.write(f'\t{round(ZTIP[j], 6)} \t{round(RTIP[j], 6)} \n')
        
    FW.write(f'\t{round(ZTIP[NPTS] + 0.05, 6)} \t{round(RTIP[NPTS], 6)} \n')
    FW.write(f'\t{round(ZTIP[NPTS] + 0.10, 6)} \t{round(RTIP[NPTS], 6)} \n')    
    FW.write(f'\t{round(ZTIP[NPTS] + 0.20, 6)} \t{round(RTIP[NPTS], 6)} \n')
    FW.write(f'\t{round(ZTIP[NPTS] + 0.30, 6)} \t{round(RTIP[NPTS], 6)} \n')
    
    FW.close
    
    return ii, RHUB, ZHUB, RTIP, ZTIP

#####################################################
#####################################################
# PART K - STACK FUNCTION
#####################################################
##################################################### 

def STACK(Output, GAMMAVALUE, TT, mdot_norm, Re, refl, Gamma, Fluid, Stage_Number, IGV_Q, NBS, Y_BLADES, TCLC, X_BLADES, igv_Phi, igv_L, NBR, RAV, A0, R_Phi, Lambda, S_Phi, X_HUB, n, viscosity_stack, ROWS, RHUB, ZHUB, RTIP, ZTIP, RVTHI, Bleed, BLADEROWS):
    
    FW = open(f'{Output}.tstk', "w+")
    FW.write(f'PyC_DES output, extension{Output}-sc \n')
    FW.write("          Gamma,TT_in(deg R),type(0-comp,1-turb),ma, Re,BL-switch,visc-flag,refl \n")
    if Fluid == 1:
        FW.write(f'      {GAMMAVALUE[0]}    {round(TT[0]*1.8, 8)}   0   {round(mdot_norm, 6)}  {round(Re, 6)}   1   1     {round(refl, 6)} \n')
    elif Fluid == 0:
        FW.write(f'      {Gamma}    {round(TT[0]*1.8, 8)}   0   {round(mdot_norm, 6)}  {round(Re, 6)}   1   1     {round(refl, 6)} \n')
    FW.write("The overall information for each blade row is: \n")
    FW.write("   Number_of_blades  rotation_speed  Hub_cl Tip_cl \n")
    FW.write("Geometry, loss, blockage and angular momentum are input starting at the hub \n")
    FW.write("    r_le    r_te   x_le   x_te   Loss  Blk  rV_thetav \n")

    #GET NORMALIZED ROTATION SPEED & CLEARANCE VALUES
    Omegan = numpy.zeros(ROWS)
    HCL    = numpy.zeros(ROWS)
    TCL    = numpy.zeros(ROWS)
    Omegan[0] = 0
    HCL[0] = 0
    TCL[0] = 0
    for j in range(1, int(n)+1, 1):
        i = int(j)
        if IGV_Q == 1: 
            Omegan[2*i-1] = RAV * refl / A0
            Omegan[2*i]   = 0
            HCL[2*i-1]    = 0
            
            if i < n:     
                HCL[2*i] = TCLC
            else:    
                HCL[2*i] = 0
                
            TCL[2*i-1] = TCLC
            TCL[2*i] = 0
            
        elif IGV_Q == 0:
            Omegan[2*i-1] = 0
            Omegan[2*i-2] = (RAV * refl / A0)
            
            if i < n: 
                HCL[2*i-1] = TCLC
                HCL[2*i-2] = 0
                TCL[2*i-1] = 0
                TCL[2*i-2] = TCLC
                
    blkin = igv_L
    
    #WITH FREE VORTEX
    r_le_hub = []
    r_te_hub = []
    r_le_tip = []
    r_te_tip = []
    
    z_le_hub = []
    z_te_hub = []
    z_le_tip = []
    z_te_tip = []
    
    loss_hub = []
    loss_tip = []
    
    block_hub= []
    block_tip= []
    
    rvth_hub = []
    rvth_tip = []
    
    for j in range(0, ROWS, 1): 
        i = int(j)
        r_le_hub.append(RHUB[2*(i+1)-1] - 0.004)
        r_te_hub.append(RHUB[2*(i+1)] - 0.004)
        r_le_tip.append(RTIP[2*(i+1)-1] + 0.004)
        r_te_tip.append(RTIP[2*(i+1)] + 0.004)
        
        z_le_hub.append(ZHUB[2*(i+1)-1])
        z_te_hub.append(ZHUB[2*(i+1)])
        z_le_tip.append(ZTIP[2*(i+1)-1])
        z_te_tip.append(ZTIP[2*(i+1)])
        
        loss_hub.append(0)
        loss_tip.append(0)
        
        if i == 0:
            block_hub.append(1 - blkin)
            block_tip.append(1 - blkin)
        else: 
            block_hub.append(0)
            block_tip.append(0)
        if IGV_Q == 1: 
            rvth_hub.append(RVTHI[i] / A0 / refl)
            rvth_tip.append(rvth_hub[i])
        elif IGV_Q == 0: 
            rvth_hub.append(RVTHI[i+1] / A0 / refl)
            rvth_tip.append(rvth_hub[i])
        
    # for j in range(0, int(n)-1, 1): 
    #     i = int(j)
    #     if IGV_Q == 1: 
    #         block_hub[2*(i+1)] = -Bleed[i]/100
    #         block_tip[2*(i+1)] = -Bleed[i]/100
            
    #     elif IGV_Q == 0: 
    #         block_hub[2*(i+1)-1] = -Bleed[i]/100
    #         block_tip[2*(i+1)-1] = -Bleed[i]/100
 
    for j in range(1, ROWS+1, 1): 
    #
        #
        i= int(j) 
        if IGV_Q == 1: 
            FW.write(f'   {BLADEROWS[i-1]}     {round(Omegan[i-1], 6)}     {round(HCL[i-1], 6)}     {round(TCL[i-1], 6)}    \n')
        if IGV_Q == 0: 
            FW.write(f'   {BLADEROWS[i]}     {round(Omegan[i-1], 6)}     {round(HCL[i-1], 6)}     {round(TCL[i-1], 6)}    \n')
        FW.write(f'   {round(r_le_hub[i-1],6)}   {round(r_te_hub[i-1],6)}   {round(z_le_hub[i-1],6)}   {round(z_te_hub[i-1],6)}   {round(loss_hub[i-1],6)}   {round(block_hub[i-1],6)}   {round(rvth_hub[i-1], 6)}    \n')
        FW.write(f'   {round(r_le_tip[i-1],6)}   {round(r_te_tip[i-1],6)}   {round(z_le_tip[i-1],6)}   {round(z_te_tip[i-1],6)}   {round(loss_tip[i-1],6)}   {round(block_tip[i-1],6)}   {round(rvth_tip[i-1], 6)}    \n')
        
        if i < ROWS: 
            FW.write(f'   {999.00000000000000}   \n')
    
    FW.close
    '''Add some reference values at the bottom of .tinf file'''
    FW = open(f'{Output}.tinf', "a")
    FW.write(' \n')
    FW.write(f' The reference length is :  {refl} m \n')
    FW.write(f' The reference velocity is :   {A0}     m/s \n')
    FW.write(f' Fluid Viscosity is :   {viscosity_stack} N-s/m^2 \n')
    FW.write(f' The reference reynolds number is :   {Re} \n')
    FW.write(f' The normalized mass flow is :  {mdot_norm}  \n')
    FW.close
    
    #No Forced Vortex Option In This Code

#####################################################
#####################################################
# PART L - EXTRAPOLATION FUNCTION (FOR COORDINATES IN WALL AND STACK FILES)
#####################################################
##################################################### 

def EXTRAPOLATION(XA, XB, XC, XD, YA, YB, YC, YD): 
    #Takes 4 points from the coordinates and interpolates to 
    #find two "new" points for the geometry
    
    #****
    #THIS ISNT REALLY "EXTRAPOLATION" MORE SO INTERPOLATION
    #DIFFERENT NAME USED TO AVOID MIX UP WITH SCIPY INTERP FUNCTION
    #****
    
    Xint = (XA*(YB - YA) + YC*(XB - XA) - ((YD - YC)/(XD - XC))*XC*(XB - XA) - YA*(XB - XA)) / ((YB-YA) - ((YD - YC)*(XB - XA))/(XD - XC))
    
    Yint = YC + ((Xint - XC) / (XD - XC) * (YD - YC))
    
    return Xint, Yint

#####################################################
#####################################################
# PART M - OBTAIN sCO2 PROPERTIES FUNCTION
#####################################################
##################################################### 

def sCO2_Properties():
    
    #Get Current Working Directory
    cwd = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd))
    #Change to Folder Directory with CO2 Files
    os.chdir(cwd + '/CO2')
    cwd2 = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd2))
    
    '''
        READ THE ENTHALPY FILE...
    '''
    
    HPT         = open('CO2_HPT.csv', 'r')
    lines4      = HPT.readlines()
    Pressure    = []
    Temperature = []
    HPT.close()
    
    Le       = lines4[7]
    Le2      = Le.split()
    LeT      = lines4[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    Enthalpy = [[0 for i in range(LZ)] for j in range(LZ)]
    Entropy  = [[0 for i in range(LZ)] for j in range(LZ)]
    SoS      = [[0 for i in range(LZ)] for j in range(LZ)]
    Density  = [[0 for i in range(LZ)] for j in range(LZ)]
    Gamma    = [[0 for i in range(LZ)] for j in range(LZ)]
    Cp       = [[0 for i in range(LZ)] for j in range(LZ)]
    ENTH     = [[0 for i in range(LZ)] for j in range(LZ)]
    Viscosity= [[0 for i in range(LZ)] for j in range(LZ)]
    
    i = 1
    while i < LZ+1 :
    #
        #
        Pressure.append(float(Le2[i]))
        Temperature.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeH  = []
        LeH2 = []
        LeH  = lines4[i+11]
        LeH2 = LeH.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeH2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Enthalpy[i][b] = (float(LeH2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTROPY FILE
    '''
    
    SPT          = open('CO2_SPT.csv', 'r')
    lines5       = SPT.readlines()
    SPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines5[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Entropy[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE SoS FILE
    '''
    
    APT          = open('CO2_APT.csv', 'r')
    lines6       = APT.readlines()
    APT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines6[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                SoS[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                SoS[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE DENSITY FILE
    '''
    
    RPT          = open('CO2_RPT.csv', 'r')
    lines7       = RPT.readlines()
    RPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines7[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Density[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Density[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE GAMMA FILE
    '''
    
    GPT          = open('CO2_GPT.csv', 'r')
    lines8       = GPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Gamma[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Gamma[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE VISCOSITY FILE
    '''
    
    VPT          = open('CO2_VPT.csv', 'r')
    lines8       = VPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Viscosity[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Viscosity[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE CP FILE
    '''
    
    CpPT         = open('CO2_CpPT.csv', 'r')
    lines9       = CpPT.readlines()
    CpPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines9[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Cp[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Cp[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE TEMPERATURE FILE
    '''
    
    TSH         = open('CO2_TSH.csv', 'r')
    lines8      = TSH.readlines()
    S           = []
    H           = []
    TSH.close()
    
    LeS       = lines8[7]
    LeS2      = LeS.split()
    LeH      = lines8[9]
    LeH2     = LeH.split()
    TempS    = [[0 for i in range(LZ)] for j in range(LZ)]
    PresS    = [[0 for i in range(LZ)] for j in range(LZ)]

    i = 1
    while i < LZ+1 :
    #
        #
        S.append(float(LeS2[i]))
        H.append(float(LeH2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeT  = []
        LeT2 = []
        LeT  = lines8[i+11]
        LeT2 = LeT.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                TempS[i][b-1] = (float(LeT2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                TempS[i][b] = (float(LeT2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE PRESSURE FILE
    '''
    
    PSH          = open('CO2_PSH.csv', 'r')
    lines9       = PSH.readlines()
    PSH.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeP  = []
        LeP2 = []
        LeP  = lines9[i+11]
        LeP2 = LeP.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                PresS[i][b-1] = (float(LeP2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                PresS[i][b] = (float(LeP2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTHALPY-SP FILE
    '''
    
    HSP         = open('CO2_HSP.csv', 'r')
    lines10     = HSP.readlines()
    ENTR        = []
    PRE         = []
    HSP.close()
    
    Le       = lines10[7]
    Le2      = Le.split()
    LeT      = lines10[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    
    i = 1
    while i < LZ+1 :
    #
        #
        ENTR.append(float(Le2[i]))
        PRE.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines10[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                ENTH[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                ENTH[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    # '''
    #     READ THE GAMMA FILE
    # '''
    
    # GSH          = open('GSH_Tcdes.csv', 'r')
    # lines10       = GSH.readlines()
    # GSH.close()
    
    # i = 0
    # while i < LZ :
    # #
    #     #
    #     LeG  = []
    #     LeG2 = []
    #     LeG  = lines10[i+2]
    #     LeG2 = LeG.split()
    #     b = 0
    #     while b < LZ:
    #     #
    #         #
    #         Gamma[i][b] = (float(LeG2[b]))
    #         b += 1
    #         #
    #     #
    #     i += 1
    #     #
    #
    HPT = interpolate.interp2d(Pressure, Temperature, Enthalpy, kind = 'cubic')
    SPT = interpolate.interp2d(Pressure, Temperature, Entropy, kind = 'cubic')
    APT = interpolate.interp2d(Pressure, Temperature, SoS, kind = 'cubic')
    RPT = interpolate.interp2d(Pressure, Temperature, Density, kind = 'cubic')
    GPT = interpolate.interp2d(Pressure, Temperature, Gamma, kind = 'cubic')
    VPT = interpolate.interp2d(Pressure, Temperature, Viscosity, kind = 'cubic')
    CpPT= interpolate.interp2d(Pressure, Temperature, Cp, kind = 'cubic')
    HSP = interpolate.interp2d(ENTR,  PRE,    ENTH, kind = 'cubic')
    TSH = interpolate.interp2d(S, H, TempS, kind = 'linear')
    PSH = interpolate.interp2d(S, H, PresS, kind = 'linear')
    #GSH = interpolate.interp2d(S, H, Gamma, kind = 'linear')
    
    #Return to Original Directory
    os.chdir(cwd)
    
    return HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT #, GSH  
  
#####################################################
#####################################################
# PART N - OBTAIN AIR PROPERTIES FUNCTION
#####################################################
#####################################################   
  
def AIR_Properties(): 
    
    #Get Current Working Directory
    cwd = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd))
    #Change to Folder Directory with AIR Files
    os.chdir(cwd + '/AIR')
    cwd2 = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd2))
    
    '''
        READ THE AIR FILE...
    '''
    
    HPT         = open('AIR_HPT.csv', 'r')
    lines4      = HPT.readlines()
    Pressure    = []
    Temperature = []
    HPT.close()
    
    Le       = lines4[7]
    Le2      = Le.split()
    LeT      = lines4[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    Enthalpy = [[0 for i in range(LZ)] for j in range(LZ)]
    Entropy  = [[0 for i in range(LZ)] for j in range(LZ)]
    SoS      = [[0 for i in range(LZ)] for j in range(LZ)]
    Density  = [[0 for i in range(LZ)] for j in range(LZ)]
    Gamma    = [[0 for i in range(LZ)] for j in range(LZ)]
    Cp       = [[0 for i in range(LZ)] for j in range(LZ)]
    ENTH     = [[0 for i in range(LZ)] for j in range(LZ)]
    Viscosity= [[0 for i in range(LZ)] for j in range(LZ)]
    
    i = 1
    while i < LZ+1 :
    #
        #
        Pressure.append(float(Le2[i]))
        Temperature.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeH  = []
        LeH2 = []
        LeH  = lines4[i+11]
        LeH2 = LeH.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeH2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Enthalpy[i][b] = (float(LeH2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTROPY FILE
    '''
    
    SPT          = open('AIR_SPT.csv', 'r')
    lines5       = SPT.readlines()
    SPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines5[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Entropy[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE SoS FILE
    '''
    
    APT          = open('AIR_APT.csv', 'r')
    lines6       = APT.readlines()
    APT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines6[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                SoS[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                SoS[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE DENSITY FILE
    '''
    
    RPT          = open('AIR_RPT.csv', 'r')
    lines7       = RPT.readlines()
    RPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines7[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Density[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Density[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE GAMMA FILE
    '''
    
    GPT          = open('AIR_GPT.csv', 'r')
    lines8       = GPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Gamma[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Gamma[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE VISCOSITY FILE
    '''
    
    VPT          = open('AIR_VPT.csv', 'r')
    lines8       = VPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Viscosity[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Viscosity[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE CP FILE
    '''
    
    CpPT         = open('AIR_CpPT.csv', 'r')
    lines9       = CpPT.readlines()
    CpPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines9[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Cp[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Cp[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE TEMPERATURE FILE
    '''
    
    TSH         = open('AIR_TSH.csv', 'r')
    lines8      = TSH.readlines()
    S           = []
    H           = []
    TSH.close()
    
    LeS       = lines8[7]
    LeS2      = LeS.split()
    LeH      = lines8[9]
    LeH2     = LeH.split()
    TempS    = [[0 for i in range(LZ)] for j in range(LZ)]
    PresS    = [[0 for i in range(LZ)] for j in range(LZ)]

    i = 1
    while i < LZ+1 :
    #
        #
        S.append(float(LeS2[i]))
        H.append(float(LeH2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeT  = []
        LeT2 = []
        LeT  = lines8[i+11]
        LeT2 = LeT.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                TempS[i][b-1] = (float(LeT2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                TempS[i][b] = (float(LeT2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE PRESSURE FILE
    '''
    
    PSH          = open('AIR_PSH.csv', 'r')
    lines9       = PSH.readlines()
    PSH.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeP  = []
        LeP2 = []
        LeP  = lines9[i+11]
        LeP2 = LeP.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                PresS[i][b-1] = (float(LeP2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                PresS[i][b] = (float(LeP2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTHALPY-SP FILE
    '''
    
    HSP         = open('AIR_HSP.csv', 'r')
    lines10     = HSP.readlines()
    ENTR        = []
    PRE         = []
    HSP.close()
    
    Le       = lines10[7]
    Le2      = Le.split()
    LeT      = lines10[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    
    i = 1
    while i < LZ+1 :
    #
        #
        ENTR.append(float(Le2[i]))
        PRE.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines10[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                ENTH[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                ENTH[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    # '''
    #     READ THE GAMMA FILE
    # '''
    
    # GSH          = open('GSH_Tcdes.csv', 'r')
    # lines10       = GSH.readlines()
    # GSH.close()
    
    # i = 0
    # while i < LZ :
    # #
    #     #
    #     LeG  = []
    #     LeG2 = []
    #     LeG  = lines10[i+2]
    #     LeG2 = LeG.split()
    #     b = 0
    #     while b < LZ:
    #     #
    #         #
    #         Gamma[i][b] = (float(LeG2[b]))
    #         b += 1
    #         #
    #     #
    #     i += 1
    #     #
    #
    HPT = interpolate.interp2d(Pressure, Temperature, Enthalpy, kind = 'cubic')
    SPT = interpolate.interp2d(Pressure, Temperature, Entropy, kind = 'cubic')
    APT = interpolate.interp2d(Pressure, Temperature, SoS, kind = 'cubic')
    RPT = interpolate.interp2d(Pressure, Temperature, Density, kind = 'cubic')
    GPT = interpolate.interp2d(Pressure, Temperature, Gamma, kind = 'cubic')
    VPT = interpolate.interp2d(Pressure, Temperature, Viscosity, kind = 'cubic')
    CpPT= interpolate.interp2d(Pressure, Temperature, Cp, kind = 'cubic')
    HSP = interpolate.interp2d(ENTR,  PRE,    ENTH, kind = 'cubic')
    TSH = interpolate.interp2d(S, H, TempS, kind = 'linear')
    PSH = interpolate.interp2d(S, H, PresS, kind = 'linear')
    #GSH = interpolate.interp2d(S, H, Gamma, kind = 'linear')
    
    #Return to Original Directory
    os.chdir(cwd)
    
    return HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT #, GSH  

#####################################################
#####################################################
# PART O - OBTAIN CH4 PROPERTIES FUNCTION
#####################################################
##################################################### 

def CH4_Properties(): 
    
    #Get Current Working Directory
    cwd = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd))
    #Change to Folder Directory with CH4 Files
    os.chdir(cwd + '/CH4')
    cwd2 = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd2))
    
    '''
        READ THE ENTHALPY FILE...
    '''
    
    HPT         = open('CH4_HPT.csv', 'r')
    lines4      = HPT.readlines()
    Pressure    = []
    Temperature = []
    HPT.close()
    
    Le       = lines4[7]
    Le2      = Le.split()
    LeT      = lines4[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    Enthalpy = [[0 for i in range(LZ)] for j in range(LZ)]
    Entropy  = [[0 for i in range(LZ)] for j in range(LZ)]
    SoS      = [[0 for i in range(LZ)] for j in range(LZ)]
    Density  = [[0 for i in range(LZ)] for j in range(LZ)]
    Gamma    = [[0 for i in range(LZ)] for j in range(LZ)]
    Cp       = [[0 for i in range(LZ)] for j in range(LZ)]
    ENTH     = [[0 for i in range(LZ)] for j in range(LZ)]
    Viscosity= [[0 for i in range(LZ)] for j in range(LZ)]
    
    i = 1
    while i < LZ+1 :
    #
        #
        Pressure.append(float(Le2[i]))
        Temperature.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeH  = []
        LeH2 = []
        LeH  = lines4[i+11]
        LeH2 = LeH.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeH2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Enthalpy[i][b] = (float(LeH2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTROPY FILE
    '''
    
    SPT          = open('CH4_SPT.csv', 'r')
    lines5       = SPT.readlines()
    SPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines5[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Entropy[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE SoS FILE
    '''
    
    APT          = open('CH4_APT.csv', 'r')
    lines6       = APT.readlines()
    APT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines6[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                SoS[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                SoS[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE DENSITY FILE
    '''
    
    RPT          = open('CH4_RPT.csv', 'r')
    lines7       = RPT.readlines()
    RPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines7[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Density[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Density[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE GAMMA FILE
    '''
    
    GPT          = open('CH4_GPT.csv', 'r')
    lines8       = GPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Gamma[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Gamma[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE VISCOSITY FILE
    '''
    
    VPT          = open('CH4_VPT.csv', 'r')
    lines8       = VPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Viscosity[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Viscosity[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE CP FILE
    '''
    
    CpPT         = open('CH4_CpPT.csv', 'r')
    lines9       = CpPT.readlines()
    CpPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines9[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Cp[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Cp[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE TEMPERATURE FILE
    '''
    
    TSH         = open('CH4_TSH.csv', 'r')
    lines8      = TSH.readlines()
    S           = []
    H           = []
    TSH.close()
    
    LeS       = lines8[7]
    LeS2      = LeS.split()
    LeH      = lines8[9]
    LeH2     = LeH.split()
    TempS    = [[0 for i in range(LZ)] for j in range(LZ)]
    PresS    = [[0 for i in range(LZ)] for j in range(LZ)]

    i = 1
    while i < LZ+1 :
    #
        #
        S.append(float(LeS2[i]))
        H.append(float(LeH2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeT  = []
        LeT2 = []
        LeT  = lines8[i+11]
        LeT2 = LeT.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                TempS[i][b-1] = (float(LeT2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                TempS[i][b] = (float(LeT2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE PRESSURE FILE
    '''
    
    PSH          = open('CH4_PSH.csv', 'r')
    lines9       = PSH.readlines()
    PSH.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeP  = []
        LeP2 = []
        LeP  = lines9[i+11]
        LeP2 = LeP.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                PresS[i][b-1] = (float(LeP2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                PresS[i][b] = (float(LeP2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTHALPY-SP FILE
    '''
    
    HSP         = open('CH4_HSP.csv', 'r')
    lines10     = HSP.readlines()
    ENTR        = []
    PRE         = []
    HSP.close()
    
    Le       = lines10[7]
    Le2      = Le.split()
    LeT      = lines10[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    
    i = 1
    while i < LZ+1 :
    #
        #
        ENTR.append(float(Le2[i]))
        PRE.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines10[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                ENTH[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                ENTH[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    # '''
    #     READ THE GAMMA FILE
    # '''
    
    # GSH          = open('GSH_Tcdes.csv', 'r')
    # lines10       = GSH.readlines()
    # GSH.close()
    
    # i = 0
    # while i < LZ :
    # #
    #     #
    #     LeG  = []
    #     LeG2 = []
    #     LeG  = lines10[i+2]
    #     LeG2 = LeG.split()
    #     b = 0
    #     while b < LZ:
    #     #
    #         #
    #         Gamma[i][b] = (float(LeG2[b]))
    #         b += 1
    #         #
    #     #
    #     i += 1
    #     #
    #
    HPT = interpolate.interp2d(Pressure, Temperature, Enthalpy, kind = 'cubic')
    SPT = interpolate.interp2d(Pressure, Temperature, Entropy, kind = 'cubic')
    APT = interpolate.interp2d(Pressure, Temperature, SoS, kind = 'cubic')
    RPT = interpolate.interp2d(Pressure, Temperature, Density, kind = 'cubic')
    GPT = interpolate.interp2d(Pressure, Temperature, Gamma, kind = 'cubic')
    VPT = interpolate.interp2d(Pressure, Temperature, Viscosity, kind = 'cubic')
    CpPT= interpolate.interp2d(Pressure, Temperature, Cp, kind = 'cubic')
    HSP = interpolate.interp2d(ENTR,  PRE,    ENTH, kind = 'cubic')
    TSH = interpolate.interp2d(S, H, TempS, kind = 'linear')
    PSH = interpolate.interp2d(S, H, PresS, kind = 'linear')
    #GSH = interpolate.interp2d(S, H, Gamma, kind = 'linear')
    
    #Return to Original Directory
    os.chdir(cwd)
    
    return HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT #, GSH  
    
def OTHER_Properties(Fluidtypename): 
    
    #Get Current Working Directory
    cwd = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd))
    #Change to Folder Directory with {Fluidtypename} Files
    os.chdir(cwd + '/' + f'{Fluidtypename}')
    cwd2 = os.getcwd()
    print('Current Working Directory: {0}'.format(cwd2))
    
    '''
        READ THE ENTHALPY FILE...
    '''
    
    HPT         = open(f'{Fluidtypename}' + '_HPT.csv', 'r')
    lines4      = HPT.readlines()
    Pressure    = []
    Temperature = []
    HPT.close()
    
    Le       = lines4[7]
    Le2      = Le.split()
    LeT      = lines4[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    Enthalpy = [[0 for i in range(LZ)] for j in range(LZ)]
    Entropy  = [[0 for i in range(LZ)] for j in range(LZ)]
    SoS      = [[0 for i in range(LZ)] for j in range(LZ)]
    Density  = [[0 for i in range(LZ)] for j in range(LZ)]
    Gamma    = [[0 for i in range(LZ)] for j in range(LZ)]
    Cp       = [[0 for i in range(LZ)] for j in range(LZ)]
    ENTH     = [[0 for i in range(LZ)] for j in range(LZ)]
    Viscosity= [[0 for i in range(LZ)] for j in range(LZ)]
    
    i = 1
    while i < LZ+1 :
    #
        #
        Pressure.append(float(Le2[i]))
        Temperature.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeH  = []
        LeH2 = []
        LeH  = lines4[i+11]
        LeH2 = LeH.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeH2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Enthalpy[i][b] = (float(LeH2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTROPY FILE
    '''
    
    SPT          = open(f'{Fluidtypename}' + '_SPT.csv', 'r')
    lines5       = SPT.readlines()
    SPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines5[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Enthalpy[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Entropy[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE SoS FILE
    '''
    
    APT          = open(f'{Fluidtypename}' + '_APT.csv', 'r')
    lines6       = APT.readlines()
    APT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines6[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                SoS[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                SoS[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE DENSITY FILE
    '''
    
    RPT          = open(f'{Fluidtypename}' + '_RPT.csv', 'r')
    lines7       = RPT.readlines()
    RPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines7[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Density[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Density[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE GAMMA FILE
    '''
    
    GPT          = open(f'{Fluidtypename}' + '_GPT.csv', 'r')
    lines8       = GPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Gamma[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Gamma[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE VISCOSITY FILE
    '''
    
    VPT          = open(f'{Fluidtypename}' + '_VPT.csv', 'r')
    lines8       = VPT.readlines()
    GPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines8[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Viscosity[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Viscosity[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE CP FILE
    '''
    
    CpPT         = open(f'{Fluidtypename}' + '_CpPT.csv', 'r')
    lines9       = CpPT.readlines()
    CpPT.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines9[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                Cp[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                Cp[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE TEMPERATURE FILE
    '''
    
    TSH         = open(f'{Fluidtypename}' + '_TSH.csv', 'r')
    lines8      = TSH.readlines()
    S           = []
    H           = []
    TSH.close()
    
    LeS       = lines8[7]
    LeS2      = LeS.split()
    LeH      = lines8[9]
    LeH2     = LeH.split()
    TempS    = [[0 for i in range(LZ)] for j in range(LZ)]
    PresS    = [[0 for i in range(LZ)] for j in range(LZ)]

    i = 1
    while i < LZ+1 :
    #
        #
        S.append(float(LeS2[i]))
        H.append(float(LeH2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeT  = []
        LeT2 = []
        LeT  = lines8[i+11]
        LeT2 = LeT.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                TempS[i][b-1] = (float(LeT2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                TempS[i][b] = (float(LeT2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE PRESSURE FILE
    '''
    
    PSH          = open(f'{Fluidtypename}' + '_PSH.csv', 'r')
    lines9       = PSH.readlines()
    PSH.close()
    
    i = 0
    while i < LZ :
    #
        #
        LeP  = []
        LeP2 = []
        LeP  = lines9[i+11]
        LeP2 = LeP.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                PresS[i][b-1] = (float(LeP2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                PresS[i][b] = (float(LeP2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    
    '''
        READ THE ENTHALPY-SP FILE
    '''
    
    HSP         = open(f'{Fluidtypename}' + '_HSP.csv', 'r')
    lines10     = HSP.readlines()
    ENTR        = []
    PRE         = []
    HSP.close()
    
    Le       = lines10[7]
    Le2      = Le.split()
    LeT      = lines10[9]
    LeT2     = LeT.split()
    LZ       = len(Le2)-1
    
    i = 1
    while i < LZ+1 :
    #
        #
        ENTR.append(float(Le2[i]))
        PRE.append(float(LeT2[i]))
        i += 1
        #
    #
    
    i = 0
    while i < LZ :
    #
        #
        LeS  = []
        LeS2 = []
        LeS  = lines10[i+11]
        LeS2 = LeS.split()
        if i == 0 :
        #
            #
            b = 1
            while b < LZ+1:
            #
                #
                ENTH[i][b-1] = (float(LeS2[b]))
                b += 1
            #
        else: 
        #
            #
            b = 0
            while b < LZ:
            #
                #
                ENTH[i][b] = (float(LeS2[b]))
                b += 1
                #
            #
        #
        i += 1
        #
    #
    # '''
    #     READ THE GAMMA FILE
    # '''
    
    # GSH          = open('GSH_Tcdes.csv', 'r')
    # lines10       = GSH.readlines()
    # GSH.close()
    
    # i = 0
    # while i < LZ :
    # #
    #     #
    #     LeG  = []
    #     LeG2 = []
    #     LeG  = lines10[i+2]
    #     LeG2 = LeG.split()
    #     b = 0
    #     while b < LZ:
    #     #
    #         #
    #         Gamma[i][b] = (float(LeG2[b]))
    #         b += 1
    #         #
    #     #
    #     i += 1
    #     #
    #
    HPT = interpolate.interp2d(Pressure, Temperature, Enthalpy, kind = 'cubic')
    SPT = interpolate.interp2d(Pressure, Temperature, Entropy, kind = 'cubic')
    APT = interpolate.interp2d(Pressure, Temperature, SoS, kind = 'cubic')
    RPT = interpolate.interp2d(Pressure, Temperature, Density, kind = 'cubic')
    GPT = interpolate.interp2d(Pressure, Temperature, Gamma, kind = 'cubic')
    VPT = interpolate.interp2d(Pressure, Temperature, Viscosity, kind = 'cubic')
    CpPT= interpolate.interp2d(Pressure, Temperature, Cp, kind = 'cubic')
    HSP = interpolate.interp2d(ENTR,  PRE,    ENTH, kind = 'cubic')
    TSH = interpolate.interp2d(S, H, TempS, kind = 'linear')
    PSH = interpolate.interp2d(S, H, PresS, kind = 'linear')
    #GSH = interpolate.interp2d(S, H, Gamma, kind = 'linear')
    
    #Return to Original Directory
    os.chdir(cwd)
    
    return HPT, SPT, APT, RPT, GPT, CpPT, HSP, TSH, PSH, VPT #, GSH  
    
#####################################################
#####################################################
# PART Q - PLOTTING FUNCTION
#####################################################
##################################################### 

def PLOTS(n, Pplot, Tplot, Fluidtypename, Rplot, Hplot, Splot, HH, SS, RAD, PT, PS, TT, TS, RS, VV, VZV, VTV, ALP, BET, CAMBER, MACHSP, MACHSH, MACHST, LOADING_C, FLOW_C, DIFFUSION, FILE, Poi_initial, Toi_initial, AD_EFF_perstage, Value_Reaction, RADh, RADt, GAMMAVALUE, VISCOSITY, SPECIFICHEAT, X_CROSS_SECTION, Y_CROSS_SECTION, STAGGERANGLE1_h, STAGGERANGLE1_p, STAGGERANGLE1_t, STAGGERANGLE2_h, STAGGERANGLE2_p, STAGGERANGLE2_t, ZETAANGLE_hD, ZETAANGLE_pD, ZETAANGLE_tD, AXIAL_WIDTH_R, CHORD_LENGTH_R, RAR_R, ROTORWIDTH_h, ROTORWIDTH_p, ROTORWIDTH_t, NBR, STAGGERANGLE1_Sh, STAGGERANGLE1_Sp, STAGGERANGLE1_St, STAGGERANGLE2_Sh, STAGGERANGLE2_Sp, STAGGERANGLE2_St, ZETAANGLE_ShD, ZETAANGLE_SpD, ZETAANGLE_StD, AXIAL_WIDTH_S, CHORD_LENGTH_S, RAR_S, STATORWIDTH_h, STATORWIDTH_p, STATORWIDTH_t, NBS, X_HUB, Y_HUB, X_TIP, Y_TIP, X_BLADES, Y_BLADES, X_ROTOR, Y_ROTOR, X_STATOR, Y_STATOR, IGV_Q, RHUB, ZHUB, RTIP, ZTIP, refl): 

    '''PLOTS'''
    
    '''Display Plots in Window? YES = 1, NO = 0'''
    DISPLAY = 0
    
    #plt.plot([1, 2, 3, 4], [10, 20, 25, 30], color='lightblue', linewidth=3)
    Nplot    = int(n) + 1
    NNplot   = len(HH)
    STGplot  = numpy.linspace(0, 1, Nplot) * (Nplot - 1)
    Fullplot = numpy.linspace(0, 1, NNplot)* NNplot
    
    f1 = plt.figure()
    plt.scatter(STGplot, Pplot, color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Total Pressure [kPa]')
    if DISPLAY == 0 :
        plt.close(f1)
    
    f2 = plt.figure()
    plt.scatter(STGplot, Tplot, color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Total Temperature [K]')
    if DISPLAY == 0 :
        plt.close(f2)
    
    f3 = plt.figure() 
    plt.scatter(STGplot, Rplot, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Static Density [kg/m3]')
    if DISPLAY == 0 :
        plt.close(f3)
    
    f4 = plt.figure()
    plt.scatter(Splot, Hplot, color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Total Entropy [J/kgK]')
    plt.ylabel('Total Enthalpy [J/kg]')
    if DISPLAY == 0 :
        plt.close(f4)

    f5 = plt.figure()
    plt.scatter(STGplot, Splot, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Entropy [J/kgK]')
    if DISPLAY == 0 :
        plt.close(f5)
    
    f6 = plt.figure()
    plt.plot(STGplot[1:], LOADING_C[1:])
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Pitch-Line Loading Coefficient []')
    if DISPLAY == 0 :
        plt.close(f6)
    
    f7 = plt.figure()
    plt.plot(STGplot[1:], FLOW_C[1:])
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Pitch-Line Flow Coefficient []')
    if DISPLAY == 0 :
        plt.close(f7)
    
    f8 = plt.figure()
    plt.scatter(Fullplot, SS, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Entropy [J/kgK]')
    if DISPLAY == 0 :
        plt.close(f8)
    
    f9 = plt.figure()
    plt.scatter(Fullplot, RAD, color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch Radius [m]')
    if DISPLAY == 0 :
        plt.close(f9)
    
    f10 = plt.figure()
    plt.scatter(PT, TT, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Pitch-Line Total Pressure [kPa]')
    plt.ylabel('Pitch-Line Total Temperature [K]')
    if DISPLAY == 0 :
        plt.close(f10)
    
    f11 = plt.figure()
    plt.scatter(PS, TS, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Pitch-Line Static Pressure [kPa]')
    plt.ylabel('Pitch-Line Static Temperature [K]')
    if DISPLAY == 0 :
        plt.close(f11)
    
    f12 = plt.figure()
    plt.scatter(Fullplot, RS, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Density [kg/m3]')
    if DISPLAY == 0 :
        plt.close(f12)
    
    f13 = plt.figure()
    plt.plot(Fullplot, VV )
    plt.plot(Fullplot, VZV )
    plt.plot(Fullplot, VTV )
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Velocity Magnitude [m/s]')
    plt.legend(['Total', 'Axial', 'Theta'])
    if DISPLAY == 0 :
        plt.close(f13)
    
    f14 = plt.figure()
    plt.plot(Fullplot, ALP, 'b' )
    plt.plot(Fullplot, BET , 'r')
    plt.plot(Fullplot, CAMBER, 'g')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Angle of Flow [deg]')
    plt.legend(['ALPHA', 'BETA', 'CAMBER'])
    if DISPLAY == 0 :
        plt.close(f14)
    
    res1 =  [abs(ele) for ele in ALP]
    res2 =  [abs(ele) for ele in BET] 
    res3 =  [abs(ele) for ele in CAMBER] 
    U    =  []
    RotorALPHA   = []
    RotorBETA    = []
    StatorALPHA  = []
    StatorBETA   = []
    RotorAngleTurning = []
    StatorAngleTurning = []
    PRplot = []
    TRplot = []
    PRploti = []
    TRploti = []
    
    for W in Fullplot :
        Ww = int(W)
        RotorALPHA.append(res1[1+2*Ww])
        RotorBETA.append(res2[1+2*Ww])
        StatorALPHA.append(res1[2*Ww])
        StatorBETA.append(res2[2*Ww])
        RotorAngleTurning.append(res3[1+2*Ww])
        StatorAngleTurning.append(res3[2*Ww])
        U.append(Ww+1)
        if Ww >= max(Fullplot)/2-2:
            StatorAngleTurning.append(res3[2*Ww+2])
            break
    
    f15 = plt.figure()
    plt.plot(U, RotorALPHA, 'b' )
    plt.plot(U, StatorALPHA, 'r' )
    plt.title("Stator Alpha Flow Angle " + f'[ {Fluidtypename} ]')
    plt.xlabel('Stage #')
    plt.ylabel('Pitch-Line Alpha [deg]')
    plt.legend(['Leading Edge', 'Trailing Edge'])
    if DISPLAY == 0 :
        plt.close(f15)
    
    f16 = plt.figure()
    plt.plot(U, RotorBETA, 'r' )
    plt.plot(U, StatorBETA, 'b' )
    plt.title("Rotor Beta Flow Angle " + f'[ {Fluidtypename} ]')
    plt.xlabel('Stage #')
    plt.ylabel('Pitch-Line Beta [deg]')
    plt.legend(['Trailing Edge', 'Leading Edge'])
    if DISPLAY == 0 :
        plt.close(f16)
    
    f17 = plt.figure()
    plt.plot(U, RotorAngleTurning, 'b')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Abs. Pitch-Line Angle Turning of Rotor [deg]')
    if DISPLAY == 0 :
        plt.close(f17)
    
    f18 = plt.figure()
    plt.plot(U, StatorAngleTurning[1:], 'r')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Abs. Pitch-Line Angle Turning of Stator [deg]')
    if DISPLAY == 0 :
        plt.close(f18)
    
    f19 = plt.figure()
    plt.plot(Fullplot, MACHSP )
    plt.plot(Fullplot, MACHSH )
    plt.plot(Fullplot, MACHST)
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Mach # []')
    plt.legend(['Pitch-Line', 'Hub', 'Tip'])
    if DISPLAY == 0 :
        plt.close(f19)
    
    f20 = plt.figure()
    plt.scatter(Fullplot[1:], DIFFUSION[1:])
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Diffusion Factor []')
    if DISPLAY == 0 :
        plt.close(f20)
    
    f21 = plt.figure()
    plt.scatter(FLOW_C[1:], LOADING_C[1:], color='darkblue', marker='*')
    plt.scatter(0.597,0.436, color='orange')
    plt.title(Fluidtypename)
    plt.xlabel('Flow Coefficent')
    plt.ylabel('Loading Coefficent')
    count = 1
    for a,b in zip(FLOW_C[1:], LOADING_C[1:]) :
        label = count
        count +=1
        plt.annotate(label,(a,b), textcoords = "offset points", xytext = (0,8), ha = 'center')
    
    if DISPLAY == 0 :
        plt.close(f21)
    
    f22 = plt.figure()
    plt.scatter(FLOW_C[1:], LOADING_C[1:], color='darkblue', marker='.')
    plt.scatter(0.597,0.436, color='orange', marker='.')
    plt.title(Fluidtypename)
    plt.xlabel('Flow Coefficent')
    plt.ylabel('Loading Coefficent')
    plt.xlim(xmin = 0, xmax = 1)
    plt.ylim(ymin = 0, ymax = 1)
    if DISPLAY == 0 :
        plt.close(f22)
            
    PRplot.append(Pplot[0] / Poi_initial)
    TRplot.append(Tplot[0] / Toi_initial)
    Ww = 0
    for W in STGplot :
        Wt = int(W)
        PRplot.append(Pplot[Ww+1]/Pplot[Ww])
        TRplot.append(Tplot[Ww+1]/Tplot[Ww])
        PRploti.append(Pplot[Ww]/Poi_initial)
        TRploti.append(Tplot[Ww]/Toi_initial)
        Ww += 1
        if Wt >= max(STGplot)-1:
            break
      
    f23 = plt.figure()
    plt.scatter(STGplot[1:], PRploti, color='darkgreen', marker='*')
    plt.xlabel('Stage #')
    plt.ylabel("Pressure Ratio (From Initial " + f'{Poi_initial}' + " kPa )")
    if DISPLAY == 0 :
        plt.close(f23)
        
    f24 = plt.figure()
    plt.scatter(STGplot[1:], TRploti, color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel("Temperature Ratio (From Initial " + f'{Toi_initial}' + " K )")
    if DISPLAY == 0 :
        plt.close(f24)
        
    f25 = plt.figure()   
    plt.scatter(STGplot[1:], PRplot[1:], color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Pressure Ratio per Stage')
    if DISPLAY == 0 :
        plt.close(f25)
    
    f26 = plt.figure()   
    plt.scatter(STGplot[1:], TRplot[1:], color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Temperature Ratio per Stage')
    if DISPLAY == 0 :
        plt.close(f26)
        
    f27 = plt.figure()   
    plt.scatter(STGplot[1:], AD_EFF_perstage, color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Efficiency [ ]')
    if DISPLAY == 0 :
        plt.close(f27)
        
    f28 = plt.figure()
    plt.scatter(STGplot[1:], Value_Reaction, color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Stage Value of Reaction (Pitch-Line)')
    if DISPLAY == 0 :
        plt.close(f28)
    
    f29 = plt.figure()
    plt.scatter(Fullplot, RADh, color='darkred', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Hub Radius [m]')
    if DISPLAY == 0 :
        plt.close(f29)
    
    f30 = plt.figure()
    plt.scatter(Fullplot, RADt, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Tip Radius [m]')
    if DISPLAY == 0 :
        plt.close(f30)
    
    f31 = plt.figure()
    plt.scatter(Fullplot, RAD, color='darkgreen', marker='*')
    plt.scatter(Fullplot, RADh, color='darkred', marker='*')
    plt.scatter(Fullplot, RADt, color='darkblue', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Radius [m]')
    plt.legend(['Pitch-Line', 'Hub', 'Tip'])
    if DISPLAY == 0 :
        plt.close(f31)
    
    f32 = plt.figure()
    plt.scatter(Fullplot, GAMMAVALUE, color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Gamma [CP/CV]')
    if DISPLAY == 0 :
        plt.close(f32)
    
    f33 = plt.figure()
    plt.scatter(Fullplot, VISCOSITY, color='darkgreen', marker='*')
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Viscosity [uPa/s]')
    if DISPLAY == 0 :
        plt.close(f33)
    
    f34 = plt.figure()
    plt.plot(PS, GAMMAVALUE)
    plt.title(Fluidtypename)
    plt.xlabel('Static Pressure [kPa]')
    plt.ylabel('Pitch-Line Gamma [CP/CV]')
    if DISPLAY == 0 :
        plt.close(f34)
    
    f35 = plt.figure()
    plt.plot(TS, GAMMAVALUE)
    plt.title(Fluidtypename)
    plt.xlabel('Static Temperature [K]')
    plt.ylabel('Pitch-Line Gamma [CP/CV]')
    if DISPLAY == 0 :
        plt.close(f35)
    
    f36 = plt.figure()
    plt.plot(Fullplot, SPECIFICHEAT)
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('Pitch-Line Specific Heat [J/(kgK)]')
    if DISPLAY == 0 :
        plt.close(f36)
    
    f37 = plt.figure()
    plt.scatter(Fullplot[1::2], DIFFUSION[1::2])
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('ROTOR -- Pitch-Line Diffusion Factor')
    if DISPLAY == 0 :
        plt.close(f37)
    
    f38 = plt.figure()
    plt.scatter(Fullplot[2::2], DIFFUSION[2::2])
    plt.title(Fluidtypename)
    plt.xlabel('Station #')
    plt.ylabel('STATOR -- Pitch-Line Diffusion Factor')
    if DISPLAY == 0 :
        plt.close(f38)
    
    
    ''' Plot Full Compressor Rough Cross Section Geometry '''    
    f39 = plt.figure()
    plt.scatter(X_CROSS_SECTION, Y_CROSS_SECTION, color='black', marker='*')
    plt.plot(X_HUB, Y_HUB,                        color='blue', linewidth=5)
    plt.plot(X_TIP, Y_TIP,                        color='blue',linewidth=5)
    plt.plot(X_BLADES, Y_BLADES,                  color='grey', linewidth=1.5)
    plt.plot(X_BLADES, Y_BLADES,                  color='black', linewidth=1.5)
    plt.title("Compressor Rough Cross Section Coordinates " + Fluidtypename)
    plt.xlabel('X-Axis [mm]')
    plt.ylabel('Y-Axis [mm]')
    plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
    plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
    plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
    plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
    if DISPLAY == 0 :
        plt.close(f39)
    
    ''' Plot Stage By Stage Cross Section Geometry '''  
    ''' Having an IGV changes the # of points graphed to get specific stages '''
    
    if IGV_Q == 1:
        i  = 0    #First index for Cross Section Coordinates
        ii = 29   #Second index for Cross Section Coordinates
        v  = 0    #First index for HUB/TIP Array Coordinates
        vv = 16   #Second index for HUB/TIP Array Coordinates
        x  = 0    #First index for BLADES Array Coordinates
        xx = 24   #Second index for BLADES Array Coordinates
        STAGE = 3 #Using this variable to graph only to how many stages we have
    elif IGV_Q == 0:
        i  = 0      #First index for Cross Section Coordinates
        ii = 29-4   #Second index for Cross Section Coordinates
        v  = 0      #First index for HUB/TIP Array Coordinates
        vv = 16-4   #Second index for HUB/TIP Array Coordinates
        x  = 0      #First index for BLADES Array Coordinates
        xx = 24-4   #Second index for BLADES Array Coordinates
        STAGE = 3   #Using this variable to graph only to how many stages we have
    
    if STAGE <= n:
        ''' First and Second Stages '''
        f40 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("1st-2nd Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f40)
    
    if STAGE <= n: 
        STAGE += 1
        ''' Second to Third Stages '''
    
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 8
        xx += 3
       
        f41 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("2nd-3rd Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f41)
    
    if STAGE <= n: 
        STAGE += 1
        ''' Third to Fourth Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 8
        xx += 8
    
        f42 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("3rd-4th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f42)
    
    if STAGE <= n: 
        STAGE += 1
        ''' Fourth to Fifth Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 8
        
        f43 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("4th-5th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f43)
    
    if STAGE <= n: 
        STAGE += 1
        ''' Fifth to Sixth Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 7
    
        f44 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("5th-6th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f44)        

    if STAGE <= n: 
        STAGE += 1
        ''' Sixth to Seventh Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 8
    
        f45 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("6th-7th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f45)  
        
    if STAGE <= n: 
        STAGE += 1
        ''' Seventh to Eighth Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 8
    
        f46 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("7th-8th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f46) 
    
    if STAGE <= n: 
        STAGE += 1
        ''' Eighth to Ninth Stages '''
        i  += 0
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 8
    
        f47 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("8th-9th Stage Rough ross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f47) 

    if STAGE <= n: 
        STAGE += 1
        ''' Ninth to Tenth Stages '''
        i  += 14
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 7
    
        f48 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("8th-9th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f48) 
    
    if STAGE <= n: 
        STAGE += 1
        ''' Tenth to Eleventh Stages '''
        i  += 14
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 7
    
        f49 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("8th-9th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx])
        if DISPLAY == 0 :
            plt.close(f49) 
    
    if STAGE <= n: 
        STAGE += 1
        ''' Eleventh to Twelvth Stages '''
        i  += 14
        ii += 15
        v  += 0
        vv += 5
        x  += 7
        xx += 7
    
        f50 = plt.figure()
        plt.scatter(X_CROSS_SECTION[i:ii], Y_CROSS_SECTION[i:ii], color='black', marker='*')
        plt.plot(X_HUB[v:vv+2], Y_HUB[v:vv+2],                        color='blue', linewidth=5)
        plt.plot(X_TIP[v:vv], Y_TIP[v:vv],                        color='blue',linewidth=5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='grey', linewidth=1.5)
        plt.plot(X_BLADES[x:xx], Y_BLADES[x:xx],                  color='black', linewidth=1.5)
        plt.title("11th-12th Stage Rough Cross Section Coordinates " + Fluidtypename)
        plt.xlabel('X-Axis [mm]')
        plt.ylabel('Y_Axis [mm]')
        plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
        plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
        plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
        plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
        plt.xlim(X_BLADES[x], X_BLADES[xx+2])
        if DISPLAY == 0 :
            plt.close(f50)             

    
    ''' ADDING WALL FILE COORDS TO ROUGHT GEOMETRY ''' 
    RHUB = [element * refl*1000 for element in RHUB]
    ZHUB = [element * refl*1000 for element in ZHUB]
    RTIP = [element * refl*1000 for element in RTIP]
    ZTIP = [element * refl*1000 for element in ZTIP]
    f52 = plt.figure()
    plt.plot(ZHUB, RHUB,      color='blue', linewidth=5)
    plt.plot(ZTIP, RTIP,      color='blue',linewidth=5)
    plt.plot(X_BLADES, Y_BLADES,                  color='grey', linewidth=1.5)
    plt.plot(X_BLADES, Y_BLADES,                  color='black', linewidth=1.5)
    plt.title("Compressor Wall File + Rough Blades " + Fluidtypename)
    plt.xlabel('X-Axis [mm]')
    plt.ylabel('Y-Axis [mm]')
    plt.legend(['Hub', 'Tip', 'Stators', 'Rotors'], prop={'size': 5})
    plt.fill_between(X_BLADES, Y_BLADES, color = 'grey')
    plt.fill_between(X_ROTOR, Y_ROTOR, color = 'black')
    plt.fill_between(X_STATOR, Y_STATOR, color = 'grey')
    plt.gca().set_aspect('equal', adjustable='box')
    plt.ylim(min(Y_BLADES)/1.15, max(Y_BLADES)*1.15)
    if DISPLAY == 0 :
        plt.close(f52)

    f53 = plt.figure()
    plt.plot(STGplot[1:], STAGGERANGLE1_h, color='red', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_h, color='red', marker='o')
    plt.plot(STGplot[1:], STAGGERANGLE1_p, color='blue', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_p, color='blue', marker='o')
    plt.plot(STGplot[1:], STAGGERANGLE1_t, color='black', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_t, color='black', marker='o')
    plt.title(" Stagger Angles for ROTORS " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Stagger Angle [rad]')
    plt.legend(['HUB -- Stagger 1', 'HUB -- Stagger 2', 'PITCH -- Stagger 1', 'PITCH -- Stagger 2', 'TIP -- Stagger 1', 'TIP -- Stagger 2'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f53)
    
    f54 = plt.figure()
    plt.plot(STGplot[1:], STAGGERANGLE1_Sh, color='red', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_Sh, color='red', marker='o')
    plt.plot(STGplot[1:], STAGGERANGLE1_Sp, color='blue', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_Sp, color='blue', marker='o')
    plt.plot(STGplot[1:], STAGGERANGLE1_St, color='black', marker='*')
    plt.plot(STGplot[1:], STAGGERANGLE2_St, color='black', marker='o')
    plt.title(" Stagger Angles for STATORS " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Stagger Angle [rad]')
    plt.legend(['HUB -- Stagger 1', 'HUB -- Stagger 2', 'PITCH -- Stagger 1', 'PITCH -- Stagger 2', 'TIP -- Stagger 1', 'TIP -- Stagger 2'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f54)
    
    f55 = plt.figure()
    plt.plot(STGplot[1:], ZETAANGLE_hD, color='red', marker='*')
    plt.plot(STGplot[1:], ZETAANGLE_pD, color='blue', marker='o')
    plt.plot(STGplot[1:], ZETAANGLE_tD, color='black', marker='P')
    plt.title("Zeta Angles for ROTORS " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Zeta Angle [deg]')
    plt.legend(['Hub', 'Pitch', 'Tip'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f55)
    
    f56 = plt.figure()
    plt.plot(STGplot[1:], ZETAANGLE_ShD, color='red', marker='o')
    plt.plot(STGplot[1:], ZETAANGLE_SpD, color='blue', marker='*')
    plt.plot(STGplot[1:], ZETAANGLE_StD, color='black', marker='P')
    plt.title("Zeta Angles for STATORS " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Zeta Angle [deg]')
    plt.legend(['Hub', 'Pitch', 'Tip'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f56)
    
    f57 = plt.figure()
    #plt.plot(STGplot[1:], AXIAL_WIDTH_R, marker='o')
    plt.plot(STGplot[1:], ROTORWIDTH_h, color='red', marker='*')
    plt.plot(STGplot[1:], ROTORWIDTH_p, color='blue', marker='P')
    plt.plot(STGplot[1:], ROTORWIDTH_t, color='black', marker='v')
    plt.title("Rotor Axial Width " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Axial Width [mm]')
    plt.legend(['Hub', 'Pitch', 'Tip'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f57)
    
    f58 = plt.figure()
    #plt.plot(STGplot[1:], AXIAL_WIDTH_S, marker='o')
    plt.plot(STGplot[1:], STATORWIDTH_h, color='red', marker='*')
    plt.plot(STGplot[1:], STATORWIDTH_p, color='blue', marker='P')
    plt.plot(STGplot[1:], STATORWIDTH_t, color='black', marker='v')
    plt.title("Stator Axial Width " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Axial Width [mm]')
    plt.legend(['Hub', 'Pitch', 'Tip'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f58)
    
    f59 = plt.figure()
    plt.plot(STGplot[1:], CHORD_LENGTH_R, color = 'red', marker='o')
    plt.plot(STGplot[1:], CHORD_LENGTH_S, color = 'black',  marker='*')
    plt.title("Chord Length " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Chord Length [mm]')
    plt.legend(['Rotors', 'Stators'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f59)
    
    f60 = plt.figure()
    plt.plot(STGplot[1:], RAR_R, color = 'red', marker='o')
    plt.plot(STGplot[1:], RAR_S, color = 'black',  marker='*')
    plt.title("Real Aspect Ratio " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Real Aspect Ratio')
    plt.legend(['Rotors', 'Stators'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f60)
    
    f61 = plt.figure()
    plt.plot(STGplot[1:], NBR, color = 'red', marker='o')
    plt.plot(STGplot[0:]+0.5, NBS, color = 'black',  marker='*')
    plt.title("Total Blade Count Per Stage " + Fluidtypename)
    plt.xlabel('Stage #')
    plt.ylabel('Blade Count')
    plt.legend(['Rotors', 'Stators'], prop={'size': 7})
    if DISPLAY == 0 :
        plt.close(f61)
    
    
    # f33,ax37 = plt.subplots()
    # img = mpimg.imread('SmithDiagram.png')
    # ax37.imshow(img, extent=[0,1,0,1])
    # plt.scatter(FLOW_C[1:], LOADING_C[1:], color='darkblue', marker='.')
    # plt.scatter(0.597,0.436, color='orange', marker='.')
    # plt.title(Fluidtypename)
    # plt.xlabel('Flow Coefficent')
    # plt.ylabel('Loading Coefficent')
    # plt.xlim(xmin = 0, xmax = 1)
    # plt.ylim(ymin = 0, ymax = 1)
    # plt.show()
     
#####################################################
#####################################################
# PART Q.i - SAVE FIGURES TO PDF
#####################################################
##################################################### 
    
    '''SAVE'''
    PDF = PdfPages(f'{FILE}.pdf')
    PDF.savefig(f1)
    PDF.savefig(f2)
    PDF.savefig(f3)
    PDF.savefig(f4)
    PDF.savefig(f5)
    PDF.savefig(f6)
    PDF.savefig(f7)
    PDF.savefig(f8)
    PDF.savefig(f9)
    PDF.savefig(f10)
    PDF.savefig(f11)
    PDF.savefig(f12)
    PDF.savefig(f13)
    PDF.savefig(f14)
    PDF.savefig(f15)
    PDF.savefig(f16)
    PDF.savefig(f17)
    PDF.savefig(f18)
    PDF.savefig(f19)
    PDF.savefig(f20)
    PDF.savefig(f21)
    PDF.savefig(f22)
    PDF.savefig(f23)
    PDF.savefig(f24)
    PDF.savefig(f25)
    PDF.savefig(f26)
    PDF.savefig(f27)
    PDF.savefig(f28)
    PDF.savefig(f29)
    PDF.savefig(f30)
    PDF.savefig(f31)
    PDF.savefig(f32)
    PDF.savefig(f33)
    PDF.savefig(f34)
    PDF.savefig(f35)
    PDF.savefig(f36)
    PDF.savefig(f37)
    PDF.savefig(f38)
    PDF.savefig(f39)
    if n > 2:
        PDF.savefig(f40)
    if n >= 3:
        PDF.savefig(f41)
    if n >= 4:
        PDF.savefig(f42)
    if n >= 5:
        PDF.savefig(f43)
    if n >= 6:
        PDF.savefig(f44)
    if n >= 7:
        PDF.savefig(f45)
    if n >= 8:
        PDF.savefig(f46)
    if n >= 9:
        PDF.savefig(f47)
    if n >= 10:
        PDF.savefig(f48)
    if n >= 11:
        PDF.savefig(f49)
    if n >= 12:
        PDF.savefig(f50)
    
    # PDF.savefig(f51)
    PDF.savefig(f52)
    PDF.savefig(f53)
    PDF.savefig(f54)
    PDF.savefig(f55)
    PDF.savefig(f56)
    PDF.savefig(f57)
    PDF.savefig(f58)
    PDF.savefig(f59)
    PDF.savefig(f60)
    PDF.savefig(f61)
    PDF.close()
    
#####################################################
#####################################################
# ACTUALLY RUNS THE ENTIRE CODE -->      
if __name__ == '__main__': main(sys.argv[1])
#####################################################
#####################################################
