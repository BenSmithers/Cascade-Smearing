try:
    import nuSQuIDS as nsq
except ImportError:
    import nuSQUIDSpy as nsq

import numpy as np
from math import pi

con = nsq.Const()

"""
Here, we will implement the cross sections for neutrino-elecron scattering 
I'll be using nuSQuIDS units for consistency with the rest of my code 

"""

PARANOID = True

M_E = 0.51099908*con.MeV
M_MU = 105.658389*con.MeV
M_Z = 90.188*con.GeV
M_W = 80.22*con.GeV
WEIN_SIN = 0.226 # sin^2(theta_w)
GAM_W = 2.08*con.GeV

# chiral couplings 
LE = 2*WEIN_SIN - 1
RE = 2*WEIN_SIN

BR_HAD = 0.6832
BR_MU = 0.1050
BR_E = 0.1046
BR_TAU = 1.0 - BR_HAD - BR_E - BR_MU


def get_diff_xs(_enu:float,_bjork_y:float, flav_in:nsq.NeutrinoCrossSections_NeutrinoFlavor, nutype_in:nsq.NeutrinoCrossSections_NeutrinoType)->float:
    """
        
    """
    scale = 1.0/(con.cm**2)

    if np.min(_bjork_y)<0.0 or np.max(_bjork_y)>1.0 or np.min(_enu)<0.0:
        print("Min bjork y {}".format(np.min(_bjork_y)))
        print("Max bjork y {}".format(np.max(_bjork_y)))
        print("Min energy y {}".format(np.min(_enu)))
        raise ValueError("Invalid arguments")

    enu, bjork_y = np.meshgrid(_enu, _bjork_y)

    term1 = (con.GF**2)*M_E*enu/(2*pi)

    if flav_in == nsq.NeutrinoCrossSections_NeutrinoFlavor.electron:
        if nutype_in==nsq.NeutrinoCrossSections_NeutrinoType.neutrino:
            """
                nue + e -> nue + e
            """
            classic_deno = 1 + 2*M_E*enu*bjork_y/(M_Z**2)

            sub_term_1 = (RE*(1-bjork_y)/classic_deno)**2
            sub_term_2 = LE/classic_deno
            sub_term_3 = 2/(1 + 2*M_E*enu*(1-bjork_y)/(M_W**2))

            return term1*(sub_term_1 + (sub_term_2 + sub_term_3)**2)*scale

        elif nutype_in==nsq.NeutrinoCrossSections_NeutrinoType.antineutrino:
            """
                nue_bar + e -> nue_bar + e
            """
            classic_deno = 1 + 2*M_E*enu*bjork_y/(M_Z**2)
            
            sub_term_1 = (RE**2)/(classic_deno**2)
            sub_term_2 = LE/classic_deno
            sub_term_3 = 2/(1- 2*M_E*enu/(M_W**2) + (0+1.0j)*GAM_W/M_W)

            scary = sub_term_2 + sub_term_3
            scary = (scary*scary.conjugate()).real
            scary *= (1-bjork_y**2)

            final = term1 * (sub_term_1 + scary)

            """
                nue_bar + e -> hadrons 
            """

            nuterm = term1*4*((1-bjork_y)**2)*((1 - (M_MU**2 - M_E**2)/(2*M_E*enu))**2)/((1 - 2*M_E*enu/(M_W**2))**2 + (GAM_W/M_W)**2 )
            nuterm *= BR_HAD/BR_MU

            """
            The final states of these two channels are different, so we just add the cross sections together directly 
            """

            return (final + nuterm)*scale
        else:
            raise ValueError("What's {}".format(nutype_in))

    elif flav_in == nsq.NeutrinoCrossSections_NeutrinoFlavor.muon:
        classic_deno = 1 + 2*M_E*enu*bjork_y/(M_Z**2)

        if nutype_in==nsq.NeutrinoCrossSections_NeutrinoType.neutrino:
            """
            numu_bar + e -> numu_bar + e 
            """
            value = (RE*RE*(1-bjork_y)**2 + LE**2)/(classic_deno**2)

            return term1*value*scale

        elif nutype_in==nsq.NeutrinoCrossSections_NeutrinoType.antineutrino:
            """
            numu + e -> numu + e

            there's also numue -> mu nue, but I leave that out because it's a track 
            """
            value = (LE*LE*(1-bjork_y)**2 + RE**2)/(classic_deno**2)

            return term1*value*scale

        else:
            raise ValueError("What's {}".format(nutype_in))

    elif flav_in == nsq.NeutrinoCrossSections_NeutrinoFlavor.tau:
        return  np.zeros( shape = np.shape(bjork_y)) 
    else:
        if PARANOID:
            raise NotImplementedError("Unknown neutrino flavor: {}".format(flav_in))
        return 0.0*bjork_y