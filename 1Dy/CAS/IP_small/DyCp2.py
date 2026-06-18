import sys, os

from pyscf import gto, lib, scf,data
from pyscf.lib import chkfile
import yaml, getopt
import numpy as np
import basis_set_exchange as bse


workdir = "./"



mol = gto.Mole()
mol.atom = '''
    Dy    9.556044    2.129450   14.226799
     C1    7.814595    3.997071   14.041436
     C1    7.200727    2.965486   14.838974
     C1    6.938404    1.876663   13.965859
     H    6.531568    1.059770   14.229882
     C1    7.361681    2.171775   12.655193
     C1    7.921230    3.469282   12.716848
     H    8.312761    3.926100   11.980707
     C    8.503646    5.313898   14.413355
     C    9.878548    4.970582   15.031895
     H   10.418109    4.481665   14.375368
     H    9.749413    4.415174   15.827782
     H   10.339512    5.798664   15.281140
     C    7.687535    6.218609   15.334204
     H    7.510011    5.752330   16.177664
     H    6.837283    6.443968   14.900887
     H    8.189838    7.040229   15.514196
     C    8.787666    6.155878   13.152411
     H    9.391098    5.661021   12.558016
     H    9.207887    7.002180   13.411719
     H    7.946631    6.338322   12.683971
     C    6.723208    2.906036   16.310740
     C    6.296226    1.497291   16.700558
     H    5.572673    1.199676   16.108768
     H    5.980750    1.495641   17.627630
     H    7.059214    0.888400   16.613129
     C    5.484590    3.792792   16.483772
     H    5.710639    4.717319   16.251671
     H    5.183826    3.753919   17.415895
     H    4.769248    3.473009   15.895940
     C    7.815513    3.303655   17.317108
     H    8.569194    2.680836   17.244694
     H    7.450916    3.270548   18.225904
     H    8.124358    4.213174   17.121820
     C    6.950428    1.380720   11.412148
     C1   11.309863    0.269903   14.059336
     C1   11.901112    1.272952   14.910574
     C1   12.181399    2.391313   14.077235
     H   12.582384    3.198662   14.375984
     C1   11.789215    2.143232   12.756625
     C1   11.242659    0.830959   12.754637
     H   10.883315    0.393475   11.991466
     C   10.609173   -1.066309   14.369601
     C    9.222604   -0.715259   14.958306
     H    8.705729   -0.207899   14.298180
     H    9.338013   -0.174634   15.767936
     H    8.744044   -1.540739   15.183129
     C   10.356125   -1.873578   13.098710
     H    9.946109   -2.732191   13.332961
     H   11.207101   -2.032895   12.638725
     H    9.754502   -1.374060   12.508433
     C   11.388838   -1.989050   15.300392
     H   11.514218   -1.548549   16.167759
     H   12.263964   -2.187113   14.906337
     H   10.890915   -2.823034   15.427303
     C   12.344826    1.283796   16.390293
     C   13.586277    0.384404   16.537472
     H   14.269858    0.662100   15.891605
     H   13.943606    0.465759   17.445988
     H   13.336493   -0.547943   16.366686
     C   12.791193    2.697941   16.809946
     H   12.036365    3.316899   16.724663
     H   13.524064    2.993955   16.231184
     H   13.096076    2.680385   17.741016
     C   11.241718    0.894893   17.376776
     H   10.517681    1.554932   17.333874
     H   11.608333    0.871062   18.284595
     H   10.890656    0.010540   17.141989
     C   12.213807    2.961852   11.525514
     C   13.524995    2.371619   11.032273
     H   13.378725    1.447438   10.740426
     H   14.183485    2.384972   11.757515
     H   13.859285    2.902245   10.278193
     C   12.475589    4.435982   11.907376
     H   11.655479    4.828820   12.273767
     H   12.750122    4.936535   11.111351
     H   13.186252    4.477696   12.581167
     C   11.159631    2.944490   10.455500
     H   11.017158    2.024059   10.151839
     H   11.450239    3.497238    9.701022
     H   10.321317    3.301810   10.818330
     C    6.745838   -0.112074   11.766167
     H    7.576732   -0.477195   12.134665
     H    6.499476   -0.608765   10.958526
     H    6.029816   -0.193875   12.431066
     C    5.590468    1.907575   11.008406
     H    4.970463    1.823686   11.762766
     H    5.251590    1.390881   10.247941
     H    5.667800    2.851174   10.753294
     C    7.952284    1.464883   10.296391
     H    7.658968    0.900753    9.551498
     H    8.824540    1.155279   10.618130
     H    8.026235    2.394094    9.994340
    '''


#Setting molecule Properties
mol.charge = 1
mol.spin = 5
mol.basis = {'Dy': gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='Dy',fmt='nwchem')),'C1' : gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='C', fmt='nwchem')),
                                               'C' : gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='C', fmt='nwchem')),
                                               'H' : gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='H', fmt='nwchem'))
                                               }
mol.verbose = 4
mol.build()
title = 'DyCp2'
#1DM calculation by pyscf
mf = scf.rohf.ROHF(mol).x2c().density_fit()
mf.max_cycle = 10000
mf.verbose = 4
mf.conv_tol = 1e-6
#chk_fname = title + '_rohf.chk'
from pyscf.lib import chkfile

mf.chkfile='DyCp2.chk'
#mf.init_guess = 'atom'
mf.init_guess = 'chk'
scfdat = chkfile.load(mf.chkfile,'scf')
mf.e_tot = scfdat['e_tot']
mf.mo_coeff = scfdat['mo_coeff']
mf.mo_occ = scfdat['mo_occ']
mf.mo_energy = scfdat['mo_energy']
mf.level_shift = 2
#mf.kernel()


print()
print("Enter to AO-DMET procedure.")
title = 'Dy'

imp_inds = mol.search_ao_label(['Dy.*'])
thres = 1e-12

from embed_sim import ssdmet, sacasscf_mixer, siso
mydmet = ssdmet.SSDMET(mf, title=title, imp_idx=imp_inds, threshold=thres, es_natorb=False).density_fit()
conv_tol = mf.conv_tol
mydmet.build(conv_tol, restore_imp = True)
#es_mf = mydmet.ROHF()
#es_mf.kernel(mydmet.es_dm)
es_mf = mydmet.es_mf



title = 'Dy'
ncas, nelec, es_mo = mydmet.avas(['Dy 4f'], minao=mol._basis['Dy'], threshold=0.5, openshell_option=2)
ncas=7
nelec=9

es_cas = sacasscf_mixer.sacasscf_mixer(es_mf, ncas, nelec, statelis = [0,0,0,0,0,21])
es_cas.kernel(es_mo)
Ha2cm = 219474.63
np.savetxt(title+'_cas_NO_SOC.txt',(es_cas.fcisolver.e_states-np.min(es_cas.fcisolver.e_states))*Ha2cm,fmt='%.6f')
#===================PT2========================
#es_ecorr = sacasscf_mixer.sacasscf_nevpt2(es_cas, method='SC')
#es_cas.fcisolver.e_states = es_cas.fcisolver.e_states + es_ecorr
#===================PT2========================
#Ha2cm = 219474.63
#np.savetxt(mydmet.title+'_opt.txt',(es_cas.fcisolver.e_states-np.min(es_cas.fcisolver.e_states))*Ha2cm,fmt='%.6f')
total_cas = mydmet.total_cas(es_cas)

mysiso = siso.SISO(title, total_cas, verbose=6).density_fit()
mysiso.kernel()

