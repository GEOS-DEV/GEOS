import numpy as np
import matplotlib.pyplot as plt

def main():
    # Input geometric parameters
    r_casing_in  = 0.1 # meter
    r_casing_out = 0.106
    r_hole       = 0.133

    # Material properties
    G_casing = 80.8e9 # Pa
    K_casing = 175e9
    lambda_casing = K_casing - 2.0/3.0 * G_casing

    G_cement = 6.45e9
    K_cement = 10.3e9
    lambda_cement = K_cement - 2.0/3.0 * G_cement

    G_rock = 4.16667e9
    K_rock = 5.5556e9

    # Loading
    # Pressure applied on the inner face of the casing
    P0 = -10.0 # MPa

    # Analytical results
    # Rigidity of the casing-cement-rock system
    rigidity = np.array([ [r_casing_out, 1.0/r_casing_out, -r_casing_out, -1.0/r_casing_out, 0.0],
                          [2.0*(lambda_casing+G_casing), -2.0*G_casing/r_casing_out/r_casing_out, -2.0*(lambda_cement+G_cement), 2.0*G_cement/r_casing_out/r_casing_out, 0.0],
                          [0.0, 0.0, r_hole, 1.0/r_hole, -1.0/r_hole],
                          [0.0, 0.0, 2.0*(lambda_cement+G_cement), -2.0*G_cement/r_hole/r_hole, 2.0*G_rock/r_hole/r_hole],
                          [2.0*(lambda_casing+G_casing), -2.0*G_casing/r_casing_in/r_casing_in, 0.0, 0.0, 0.0] ])

    # Vector of force
    force = np.array([ 0.0, 0.0, 0.0, 0.0, P0 ])

    # Compute the coefficients describing the closed-form solutions of stress/strain
    vectorCoefficientAB = np.dot( np.linalg.inv( rigidity ), force )
    coeffA_cement = vectorCoefficientAB[2]
    coeffB_cement = vectorCoefficientAB[3]

    # Radial coordinate
    r_anal = np.arange( r_casing_out, r_hole, 0.01 * ( r_hole - r_casing_out ) )

    # Radial and hoop (tangent) stresses
    tmpVal1 = ( 2.0 * lambda_cement + 2.0 * G_cement ) * coeffA_cement
    tmpVal2 = 2.0 * G_cement * coeffB_cement / r_anal / r_anal

    sig_rr_anal = tmpVal1 - tmpVal2
    sig_tt_anal = tmpVal1 + tmpVal2

    # GEOSX results
    # Stresses are extracted along x-axis 
    # where radial stress coincides with stress_11 
    # and hoop stress coincides with stress_22

    # Radial stress, coordinates
    r_geosx, sig_rr_geosx = [], []

    for line in open('stress_11.curve', 'r' ):
        if not ( line.strip().startswith( "#" ) or line.strip()=='' ):
            values = [ float(s) for s in line.split() ]
            rVal = values[0]
            sigVal = values[1]*1e-6 # convert to MPa

            r_geosx.append( rVal )
            sig_rr_geosx.append( sigVal )

    # Tangent stress
    sig_tt_geosx = []
    for line in open( 'stress_22.curve', 'r' ):
        if not ( line.strip().startswith( "#" ) or line.strip()=='' ):
            values = [ float(s) for s in line.split() ]
            sigVal = values[1]*1e-6 # convert to MPa

            sig_tt_geosx.append( sigVal )

    # Plots
    fig = plt.figure( figsize=[ 13, 5 ] )

    plt.subplot( 121 )
    plt.plot( r_geosx, sig_rr_geosx, 'ko', label='GEOSX result' )
    plt.plot( r_anal, sig_rr_anal,  'k', linewidth=2, label='Analytic' )
    plt.xlim( r_casing_out, r_hole )
    plt.xlabel( "r (m)" )
    plt.ylabel( "Radial stress (MPa)" )

    plt.subplot( 122 )
    plt.plot( r_geosx, sig_tt_geosx, 'ko', label='GEOSX result' )
    plt.plot( r_anal, sig_tt_anal,  'k', linewidth=2, label='Analytic' )
    plt.xlim( r_casing_out, r_hole )
    plt.xlabel( "r (m)" )
    plt.ylabel( "Hoop stress (MPa)" )

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
