# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-

import cpinterface

class GAMESSUS(cpinterface.MolecularCalculator):
    def __init__(self):
        self.methods = ['semiempirical:am1', 'semiempirical:pm3',
                        'semiempirical:rm1', 'semiempirical:mndo']
    
    def create_geometry(self, molecule,
                        options={'coordinates' : 'cartesian'}):
        """Create input geometry for a subsequent calculation.

        options:
         coordinates: 'cartesian' or 'zmatrix'

        @param molecule: molecular system data to convert to input geometry
        @type molecule : cinfony molecule
        @return: a GAMESS-US input with geometry specifications
        @rtype : str
        """
        
        csystem = self.check_coordinates(options.get('coordinates'))

        if csystem == 'cartesian':
            geometry = molecule.write('inp')

            #change now-obsolete CART coordinate designation to
            #synonymous PRINAXIS
            geometry=geometry.replace("COORD=CART", "COORD=PRINAXIS")

        elif csystem == 'zmatrix':
            pass

        #add title
        geometry = geometry.replace(" $DATA\n\n",
                                    " $DATA\n\n" + molecule.title + '\n')

        return geometry

    def make_energy_job(self, molecule, method, options={}):
        """Create an input specification for a single point energy calculation.

        @param molecule: molecular system for energy calculation
        @type molecule : cinfony molecule
        @param method: calculation method
        @type method : str
        @return: a GAMESS-US input for single point energy calculation
        @rtype : str
        """
        
        self.check_method(method)
        if method.startswith('semiempirical'):
            return self.make_semiempirical_job(molecule, method, 'ENERGY',
                                               options=options)

    def make_semiempirical_job(self, molecule, method, runtyp, options={}):
        """Create a semiempirical input specification for a calculation.
        GAMESS-US supports MNDO, RM1, AM1, and PM3 methods.
        """

        geometry = self.create_geometry(molecule, options=options)
        semethod = method.split('semiempirical:')[-1].upper()

        control_options = {'runtyp' : runtyp, 'mult' : molecule.spin,
                           'icharg' : molecule.charge}
        contrl = "SCFTYP=RHF RUNTYP={runtyp} ICHARG={icharg} MULT={mult} ".format(**control_options)
        job = geometry.replace("$CONTRL ", "$CONTRL " + contrl)
        basis = " $BASIS GBASIS={0} $END".format(semethod)
        joblines=job.split('\n')
        joblines.insert(1, basis)
        job = "\n".join(joblines)

        return job

