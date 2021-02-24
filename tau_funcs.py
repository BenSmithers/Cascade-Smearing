"""
Another class is introduced here to handle the weird way that taus decay.

So normally the energy deposited is (hadronic only) for NC events, and (all) for CC events. 
That's not the case with taus since they decay very quickly 
So instead, for taus, it's (hadronic only) or NC events and (hadronic with some leptonic) for CC events. This figures out how much "some" is 
"""

# system tools
import os
import pickle 

# numeric tools
import numpy as np
from scipy import integrate 

from cascade.utils import config, savefile, bilinear_interp, get_loc

# tau file verison should be updated with changes to this code! 
tau_file_version = "4.1.6"
tau_file_name = config["tau_data"]

full_path = os.path.join(config["datapath"], tau_file_name)

class TauData:
    """
    This pre-calculates some of the expectation values for the tau decay
    
    My reasoning for this system is that the integration is very costly. Since we're moving into a regime of higher dimensionality (3!), costly integration would cause tremendous slowdowns. 
    So, instead, we calculate a lot of integrals beforehand. Then, we can just interpolate between neighboring points for any arbitrary point and quickly get an accurate (enough) result 

    TODO: maybe do this with splines? 
    """
    def __init__(self):
        """
        Try to load the data in from a file, if it's not there, it's illformatted, or it's outdated, then just make a new one 
        """
        if os.path.exists(full_path):
            # load it, verify version number 
            file_object = open(full_path,'rb')
            tau_file = pickle.load(file_object)
            file_object.close()

            if not isinstance(tau_file, dict):
                self._gen_tau_data()
                return

            try:
                if tau_file["version"]==tau_file_version:
                    self._total_energies    = tau_file["total_energies"]
                    self._depo_energies      = tau_file["depo_energies"]
                    self._xs_scale          = tau_file["xs_scale"]
                else:
                    self._gen_tau_data()
            except KeyError: # somehow a different dict got saved here? 
                self._gen_tau_data()        
        else:
            self._gen_tau_data()
 
         
    def _gen_tau_data(self):
        print("Generating Tau Data")
        self._nodes = 30
        self._total_energies = np.logspace(1, 10, self._nodes) #represents E_tau
        self._depo_energies = np.logspace(1,10,self._nodes)

        self._xs_scale = np.zeros( (2,self._nodes, self._nodes)) 

        for e_tot_i in range(self._nodes):
            for e_tau_i in range(self._nodes):
                if self._depo_energies[e_tau_i]>=self._total_energies[e_tot_i]:
                    continue

                self._xs_scale[0][e_tot_i][e_tau_i] = self._get_mpv_xs(self._total_energies[e_tot_i],self._depo_energies[e_tau_i], 1)
                self._xs_scale[1][e_tot_i][e_tau_i] = self._get_mpv_xs(self._total_energies[e_tot_i],self._depo_energies[e_tau_i], -1)
                

            #boostedMichel( E_e, E_tau )
            #vals = integrate.quad( lambda z: TauDecayToAll(self._energy_nodes[e_i], z*self._energy_nodes[e_i],-1)*z*self._energy_nodes[e_i], 0.0, 0.999)
            #self._expected_michell[0][e_i] = vals[0]
        
        savefile(full_path, version=tau_file_version, total_energies=self._total_energies, depo_energies=self._depo_energies, xs_scale=self._xs_scale) 

    def _get_mpv_xs(self, E_total, E_depo, P):
        """
        Here, we integrate over the possible E_tau values to find the median one, then evaluate the cross section scale factor there 
        """
        if (P!=-1) and (P!=1):
            raise ValueError("Unepected P {}".format(P))
       
        # Consider the interaction chain that we care about
        """
                  init shower 
                 /
        nu_tau  <            sec_nu_tau
                 \         /
                  tau-----<
                           \ 
                            < secondary_shower 
        """
        # we have the primary tau's energy, but we need to integrate over all possible sec_nu_tau energies
        minimum = 0.
        maximum = E_total
        BR_l = 0.18
        def function(secondary_shower):
            """
            \int x*f(x) dx 

            since we want the mean value 
            """
            E_charged_tau = E_total - (E_depo - secondary_shower)
            e_sec_nu = E_charged_tau - secondary_shower 
            val = TauDecayToAllHadrons(E_charged_tau, e_sec_nu, P)+BR_l*TauLeptonDecay(e_sec_nu, E_charged_tau, P)
            return( val if val>0 else 0. )
       
        return( integrate.quad(function, minimum, maximum)[0]/maximum )


    @property
    def version(self):
        return(self._version)

    def __call__(self, E_total, E_depo, P):
        if not (isinstance(E_depo, float) or isinstance(E_depo, int)):
            raise TypeError("Expected {} for E_depo, not {}".format(float, type(E_depo)))

        if E_depo>E_total:
            return 0.

        p0 = (E_total, E_depo)
        
        total_e_lower, total_e_upper = get_loc(E_total, self._total_energies)
        had_e_lower, had_e_upper = get_loc(E_depo, self._depo_energies)

        p1 = (self._total_energies[total_e_lower], self._depo_energies[had_e_lower])
        p2 = (self._total_energies[total_e_upper], self._depo_energies[had_e_upper])

        if P==1:
            p_i = 0
        elif P==-1:
            p_i = 1
        else:
            raise ValueError("P not recognized: {}, use 1 or -1".format(P))

        q11 = self._xs_scale[p_i][total_e_lower][had_e_lower]
        q12 = self._xs_scale[p_i][total_e_lower][had_e_upper]
        q21 = self._xs_scale[p_i][total_e_upper][had_e_lower]
        q22 = self._xs_scale[p_i][total_e_upper][had_e_upper]

        return(bilinear_interp( p0, p1, p2, q11, q12, q21, q22))

def args_are_floats(*args):
    """
    Simply makes sure that all the args you pass it are floats or ints 

    Returns NOTHING
    Raisese TypeError if something isn't an int or a float 
    """
    for arg in args:
        if not (isinstance(arg, float), isinstance(arg, int)):
            raise TypeError("Found an arg that's a {}, not a {}: {}".format(float, type(arg), arg))


"""
Much of this code is modified from 
https://github.com/IceCubeOpenSource/TauRunner/blob/master/python/Casino.py

I've added the type checking and some functionality to pre-calculate the expected tau decay stuff 
"""

#branching ratios of the charged tau 
RPion = 0.07856**2 
RRho = 0.43335**2 
RA1 = 0.70913**2
BrLepton = 0.18
BrPion = 0.12
BrRho = 0.26
BrA1 = 0.13
BrHad = 0.13

def TauLeptonDecay( Etau, Enu, P):
    if not Enu<= Etau:
        raise ValueError("{} should be <= {}".format(Enu, Etau))
    z = Enu/Etau
    z2 = z*z
    z3 = z*z*z

    factor1 = (5./3) - (3.*z2)  + (4./3)*z3
    factor2 = (-1./3) + (3.*z2) - (8./3)*z3

    if P==1:
        return ( factor1 + (P*factor2)) 
    elif P==-1:
        return ( factor1 - (P*factor2))
    else:
        raise TypeError("Invalid P {}".format(P))

#    z = Enu/Etau
#    g0 = (5./3.) - 3.*z**2 + (4./3.)*z**3
#    g1 = (1./3.) - 3.*z**2 + (8./3.)*z**3
#    return(g0+P*g1)

    # For a given Etau, we need the expected E_e 

def TauDecayToPion(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RPion - z)  > 0.0):
        g0 = 1./(1. - RPion)
        g1 = -(2.*z - 1. - RPion)/(1. - RPion)**2
    return(g0+P*g1)

def TauDecayToRho(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RRho - z) > 0.0):
        g0 = 1./(1. - RRho)
        g1 = -((2.*z-1.+RRho)/(1.-RRho))*((1.-2.*RRho)/(1.+2.*RRho))
    return(g0+P*g1)

def TauDecayToA1(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0 = 0.
    g1 = 0.
    if((1. - RA1 - z) > 0.0):
        g0 = (1./(1.-RA1))
        g1 = -((2.*z-1.+RA1)/(1.-RA1))*((1.-2.*RA1)/(1.+2.*RA1))
    return(g0 + P*g1)

def TauDecayToHadrons(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    z = Enu/Etau
    g0=0.
    g1=0.
    if((0.3 - z) > 0.):
        g0 = 1./0.3
    return(g0+P*g1)

def TauDecayToAllHadrons(Etau, Enu, P):
    args_are_floats(Etau, Enu, P)
    """
    This adds up the differential decay rates 

    Enu is energy of outgoing tau neutrino
    ETau is energy of the intermediate tau lepton
    P is a polarization quantity (-1 for TauMinus)
    """
    #Etau is the energy of the tau lepton, Enu is the energy of the nu_tau after the tau decays
    decay_spectra = 0
    decay_spectra+=BrPion*TauDecayToPion(Etau, Enu, P)
    decay_spectra+=BrRho*TauDecayToRho(Etau, Enu, P)
    decay_spectra+=BrA1*TauDecayToA1(Etau, Enu, P)
    decay_spectra+=BrHad*TauDecayToHadrons(Etau, Enu, P)
    return decay_spectra

def TauDecayToAll(Etau, Enu, P=1):
    decay_spectra = 0
    decay_spectra += TauDecayToAllHadrons(Etau, Enu, P)
    decay_spectra += 2*BrLepton*TauLeptonDecay(Etau, Enu, P)
    return decay_spectra

# integrate.quad( function, min, max)

if False:
    # print( integrate.quad( lambda z:TauDecayToAll(1.0, z), 0, 1) )

    from cascade.utils import bhist
    zs = np.logspace(-8,0, 1000)
    cens = bhist([zs]).centers

    dz = np.array([TauDecayToAll(1.0, z) for z in cens])
    wids = bhist([zs]).widths

    import matplotlib 
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt

    plt.plot(cens,dz*wids)
    plt.show()
