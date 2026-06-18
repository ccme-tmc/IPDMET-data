from pyscf import gto, scf, ao2mo, lib
from embed_sim import rdiis, ssdmet, sacasscf_mixer, siso, myavas
from pyscf.tools import molden, cubegen
import basis_set_exchange as bse
from functools import reduce
import numpy as np
import os

mol = gto.M(atom='structure.xyz',
        basis={
            'Dy':gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='Dy',fmt='nwchem')),
            'C1':gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='C',fmt='nwchem')),
            'C':gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='C',fmt='nwchem')),
            'H':gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='H',fmt='nwchem'))
            },
        symmetry=0,spin=5,charge=1,verbose=4)
title = 'DyCp2'
mf = scf.ROHF(mol).x2c().density_fit()
mf.chkfile = 'diis.chk'
mf.init_guess = 'chk'
#mf.level_shift = 0.2
mf.conv_tol = 1e-6
mf.diis = rdiis.RDIIS(rdiis_prop='dS',imp_idx=mol.search_ao_label(['Dy.*f']),power=0.)
mf.max_cycle = 0
mf.conv_check = False
mf.kernel()

ncas, nelec, mo = myavas.avas(mf, 'Dy 4f', minao=gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='Dy',fmt='nwchem')), openshell_option=3, threshold=0.5)
mycas = sacasscf_mixer.sacasscf_mixer(mf, ncas, nelec, statelis=[0,0,0,0,0,21])
mycas.kernel(mo)

ecorr = sacasscf_mixer.sacasscf_nevpt2(mycas)
mycas.fcisolver.e_states = mycas.fcisolver.e_states + ecorr
Ha2cm = 219474.63
np.savetxt(title+'_opt.txt',(mycas.fcisolver.e_states-np.min(mycas.fcisolver.e_states))*Ha2cm,fmt='%.6f')

mysiso = siso.SISO(title, mycas, verbose=5).density_fit()
mysiso.kernel()

from embed_sim.psh import PSH
mypsh = PSH(mysiso,Nstate=int(16),lmax=int(12),TFO=1,TFB=0)
mypsh.gen_ESO()
mypsh.get_E_J()
mypsh.get_g_tensor()
