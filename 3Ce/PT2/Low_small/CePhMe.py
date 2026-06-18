import sys, os
from embed_sim import rdiis

from pyscf import gto, lib, scf,data
from pyscf.lib import chkfile
import yaml, getopt
import numpy as np
import basis_set_exchange as bse


workdir = "./"



mol = gto.Mole()
mol.atom = '''
       Ce    0.000000    0.000000    0.000000
        N   -1.922661    0.723329   -1.565825
        N    0.334909   -2.026738   -1.565825
        N    1.587752    1.303409   -1.565825
        N   -3.091228    0.118782   -1.179810
        N    1.442746   -2.736473   -1.179810
        N    1.648482    2.617691   -1.179810
        N   -3.035207   -0.019312    1.331579
        N    1.534329   -2.618910    1.331579
        N    1.500878    2.638223    1.331579
        N   -1.878574    0.650064    1.660681
        N    0.376316   -1.951925    1.660681
        N    1.502259    1.301861    1.660681
        C   -1.212817    2.383973   -3.262449
        C   -1.458172   -2.242316   -3.262449
        C    2.670989   -0.141656   -3.262449
        H   -0.999904    3.152180   -2.692533
        H   -2.229916   -2.442032   -2.692533
        H    3.229820   -0.710148   -2.692533
        H   -0.401688    1.855150   -3.414279
        H   -1.405764   -1.275447   -3.414279
        H    1.807452   -0.579703   -3.414279
        C   -2.240547    1.536986   -2.591667
        C   -0.210796   -2.708863   -2.591667
        C    2.451343    1.171877   -2.591667
        C   -3.596961    1.450924   -2.865452
        C    0.541943   -3.840521   -2.865452
        C    3.055017    2.389597   -2.865452
        H   -4.077906    1.911597   -3.543084
        H    0.383460   -4.487370   -3.543084
        H    3.694446    2.575770   -3.543084
        C   -4.095629    0.563579   -1.954233
        C    1.559741   -3.828708   -1.954233
        C    2.535888    3.265130   -1.954233
        H   -5.006200    0.302874   -1.879908
        H    2.240803   -4.486933   -1.879908
        H    2.765397    4.184059   -1.879908
        B   -3.084165   -0.839861    0.024164
        B    2.269423   -2.251034    0.024164
        B    0.814742    3.090895    0.024164
        H   -2.117435   -1.472816   -0.047726
        H    2.334213   -1.097345   -0.047726
        H   -0.216777    2.570160   -0.047726
        H   -3.991301   -1.492502    0.020170
        H    3.288195   -2.710316    0.020170
        H    0.703105    4.202819    0.020170
        H    3.120901   -0.000000   -4.121686
        H   -1.560450    2.702779   -4.121686
        H   -1.560450   -2.702779   -4.121686
        C   -3.980827    0.244079    2.243597
        C    1.779035   -3.569537    2.243597
        C    2.201792    3.325458    2.243597
        C   -2.143903    1.324502    2.792963
        C   -0.075101   -2.518925    2.792963
        C    2.219004    1.194424    2.792963
        H   -4.866131   -0.100468    2.238024
        C   -3.458613    1.095141    3.188362
        C    0.780885   -3.542818    3.188362
        H    2.520073   -4.163960    2.238024
        C    2.677727    2.447676    3.188362
        H    2.346058    4.264426    2.238024
        C   -1.119367    2.147862    3.502087
        C   -1.300420   -2.043331    3.502087
        C    2.419787   -0.104531    3.502087
        H   -3.901220    1.451220    3.949808
        H    0.693817   -4.104165    3.949808
        H    3.207403    2.652946    3.949808
        H   -0.406673    1.564845    3.837559
        H   -0.739681    2.803387    2.880390
        H   -1.539460    2.615836    4.253786
        H   -1.495650   -2.641129    4.253786
        H   -2.057964   -2.042276    2.880390
        H   -1.151859   -1.134611    3.837559
        H    2.797645   -0.761111    2.880390
        H    1.558532   -0.430233    3.837559
        H    3.035110    0.025293    4.253786
    '''


#Setting molecule Properties
mol.charge = 0
mol.spin = 1
mol.basis = {'Ce': gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='Ce',fmt='nwchem')),'C' : gto.basis.parse(bse.get_basis('ANO-R0',elements='C', fmt='nwchem')),
                                               'N' : gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='N', fmt='nwchem')),
                                               'B' : gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='B', fmt='nwchem')),
                                               'H' : gto.basis.parse(bse.get_basis('ANO-R0',elements='H', fmt='nwchem'))}
mol.verbose = 4
mol.build()
title = 'CePhMe'
#1DM calculation by pyscf
mf = scf.rohf.ROHF(mol).x2c().density_fit()
mf.max_cycle = 10000
mf.verbose = 4
mf.conv_tol = 1e-7
#chk_fname = title + '_rhf.chk'
from pyscf.lib import chkfile

mf.chkfile='CePhMe_rohf.chk'
mf.init_guess = 'chk'
scfdat = chkfile.load(mf.chkfile,'scf')
mf.e_tot = scfdat['e_tot']
mf.mo_coeff = scfdat['mo_coeff']
mf.mo_occ = scfdat['mo_occ']
mf.mo_energy = scfdat['mo_energy']
#mf.kernel()


print()
print("Enter to AO-DMET procedure.")
title = 'Ce'

imp_inds = mol.search_ao_label(['Ce.*'])
thres = 1e-12

from embed_sim import ssdmet, sacasscf_mixer, siso
mydmet = ssdmet.SSDMET(mf, title=title, imp_idx=imp_inds, threshold=thres, es_natorb=False).density_fit()
conv_tol = mf.conv_tol
mydmet.build(conv_tol)
#es_mf = mydmet.ROHF()
#es_mf.kernel(mydmet.es_dm)
es_mf = mydmet.es_mf

'''
from pyscf.tools import molden
with open(mydmet.title+' imp_rohf_orbs.molden', 'w')as f1:
    molden.header(MF.mol, f1)
    molden.orbital_coeff(MF.mol, f1, mydmet.es_orb @ es_mf.mo_coeff, ene=es_mf.mo_energy, occ=es_mf.mo_occ)


'''
#ncas, nelec, es_mo = mydmet.avas(['Ce 4f','Ce 5d'], minao=CLUS_MOL._basis['Ce'], threshold=0.5, openshell_option=3)

title = 'Ce'
ncas, nelec, es_mo = mydmet.avas(['Ce 4f','Ce 5d'], minao=mol._basis['Ce'], threshold=0.5, openshell_option=2)
ncas=12
nelec=1

es_cas = sacasscf_mixer.sacasscf_mixer(es_mf, ncas, nelec)
es_cas.kernel(es_mo)
Ha2cm = 219474.63
np.savetxt(title+'_cas_NO_SOC.txt',(es_cas.fcisolver.e_states-np.min(es_cas.fcisolver.e_states))*Ha2cm,fmt='%.6f')
#===================PT2========================
es_ecorr = sacasscf_mixer.sacasscf_nevpt2(es_cas, method='SC')
es_cas.fcisolver.e_states = es_cas.fcisolver.e_states + es_ecorr
#===================PT2========================
Ha2cm = 219474.63
np.savetxt(mydmet.title+'_opt.txt',(es_cas.fcisolver.e_states-np.min(es_cas.fcisolver.e_states))*Ha2cm,fmt='%.6f')
total_cas = mydmet.total_cas(es_cas)

mysiso = siso.SISO(title, total_cas, verbose=6).density_fit()
mysiso.kernel()

