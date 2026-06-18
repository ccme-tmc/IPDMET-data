import sys, os

from pyscf import gto, lib, scf,data
from pyscf.lib import chkfile
import yaml, getopt
import numpy as np
import basis_set_exchange as bse


workdir = "./"



mol = gto.Mole()
mol.atom = '''
Dy         4.13440       14.71800        6.08340
N          3.61030       16.73800        5.33580
N          3.51710       12.72680        6.80140
C1          6.22680       12.86060        5.87140
C1          6.64740       13.78080        6.68970
H          6.11570       13.97670        7.45090
Si         3.91070       16.68940        3.64460
Si         3.36790       17.98660        6.51650
Si         4.82870       11.66070        6.53700
C          3.60960       17.02380        8.13020
H          3.02880       16.21450        8.05260
C          4.15800       14.78000        3.37560
H          5.01910       14.51570        3.78570
H          3.43920       14.30080        3.85900
C          4.15410       14.30980        1.91460
H          4.28470       13.33870        1.88290
H          3.29540       14.53680        1.49970
H          4.87950       14.75240        1.42620
C          5.07600       16.48870        8.23730
H          5.32550       16.04900        7.39870
H          5.68510       17.23580        8.41020
H          5.13360       15.84330        8.97410
C          4.72800       19.32020        6.30450
H          5.48900       18.82850        5.88050
C          6.92190       12.48760        4.55170
H          6.42190       11.77070        4.10970
H          6.95330       13.27330        3.96540
H          7.83510       12.18340        4.73720
C          7.91900       14.54690        6.50280
H          8.02780       15.18330        7.24100
H          8.67720       13.92630        6.49320
H          7.88530       15.03510        5.65390
C          5.53040       17.52850        3.13860
H          5.51270       18.45620        3.51230
C          2.51540       17.22450        2.49130
H          2.77010       16.93040        1.56900
C          1.62410       18.67980        6.42980
H          1.51300       19.02260        5.49690
C          5.56960       10.73230        8.00260
H          4.95160        9.96520        8.17590
C          4.48880       10.48100        5.09190
H          5.39840       10.46870        4.67720
C          3.18150       17.72710        9.44760
H          2.26220       18.05570        9.35780
H          3.23000       17.08920       10.18860
H          3.78160       18.48160        9.62480
C          5.31330       19.88780        7.55810
H          5.53470       19.15790        8.17360
H          6.12520       20.39050        7.34110
H          4.66130       20.48490        7.98070
C          4.34560       20.37630        5.27200
H          3.97130       19.93780        4.48020
H          3.67810       20.98170        5.65720
H          5.14200       20.88730        5.01760
C          6.74450       16.80690        3.76540
H          6.60990       16.72080        4.73300
H          6.83570       15.91560        3.36850
H          7.55770       17.32590        3.59430
C          5.76420       17.65210        1.62060
H          5.00200       18.10850        1.20730
H          6.58180       18.16790        1.45570
H          5.85980       16.75760        1.22990
C          2.34760       18.75270        2.44340
H          3.20800       19.17100        2.23320
H          1.69210       18.98870        1.75380
H          2.03280       19.07390        3.31510
C          1.17070       16.55760        2.81950
H          1.28640       15.58540        2.84810
H          0.85460       16.87470        3.69210
H          0.51390       16.78950        2.13080
C          0.56760       17.59740        6.60310
H          0.77940       16.83970        6.01770
H          0.55730       17.29730        7.53510
H         -0.31160       17.95790        6.36470
C          1.28060       19.86140        7.35070
H          1.95500       20.56390        7.24720
H          0.39870       20.21640        7.11230
H          1.26760       19.55570        8.28260
C          5.67080       11.48430        9.31540
H          4.80080       11.88160        9.53160
H          5.93250       10.86460       10.02750
H          6.34220       12.19320        9.23500
C          6.91400       10.10600        7.64700
H          6.83560        9.62370        6.79820
H          7.58950       10.81070        7.56220
H          7.17980        9.48210        8.35470
C          4.23290        9.03180        5.13070
H          4.72540        8.62830        5.87110
H          3.26770        8.87540        5.25860
H          4.51280        8.62510        4.28490
C          3.69620       11.07490        3.92720
H          3.85590       12.04140        3.88180
H          3.98560       10.65780        3.08810
H          2.74040       10.90820        4.06180
Si         1.97570       12.40450        7.51020
C         -0.33200       14.29160        7.35070
H         -0.21140       14.45520        8.31010
H         -0.66290       15.10770        6.92030
H         -0.97930       13.56820        7.22170
C          1.98020       12.64770        9.41570
H          2.60770       11.96220        9.78430
C          1.01650       13.89230        6.72620
H          1.61490       14.68970        6.81030
C          2.54090       14.01990        9.83970
H          3.39380       14.17970        9.38530
H          1.90300       14.72300        9.59370
H          2.67940       14.03000       10.80950
C          0.63710       12.45920       10.16560
H          0.24570       11.59520        9.92260
H          0.79680       12.48670       11.13200
H          0.02040       13.17960        9.91620
C          0.85330       13.62880        5.20360
H          1.71870       13.37550        4.82120
H          0.20990       12.90190        5.06650
H          0.52630       14.44140        4.76510
C          1.89210        9.71890        7.95470
H          2.86830        9.77660        7.90040
H          1.60710        9.85940        8.88210
H          1.59930        8.83320        7.65410
C          1.26510       10.79520        7.06120
H          1.53810       10.59810        6.11820
C         -0.27240       10.73440        7.12050
H         -0.64900       11.42990        6.54140
H         -0.57700        9.85460        6.81440
H         -0.56920       10.88080        8.04240
    '''


#Setting molecule Properties
mol.charge = 1
mol.spin = 5
mol.basis = {'Dy': gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='Dy',fmt='nwchem')),
             'N' : gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='N', fmt='nwchem')),
             'C1' : gto.basis.parse(bse.get_basis('ANO-RCC-VTZP',elements='C', fmt='nwchem')),
             'Si' : gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='Si', fmt='nwchem')),
             'C' : gto.basis.parse(bse.get_basis('ANO-RCC-VDZP',elements='C', fmt='nwchem')),
             'H' : gto.basis.parse(bse.get_basis('ANO-R0',elements='H', fmt='nwchem'))
                                               }
mol.verbose = 4
mol.build()
title = '2Dy'
#1DM calculation by pyscf
mf = scf.rohf.ROHF(mol).x2c().density_fit()
mf.max_cycle = 10000
mf.verbose = 4
mf.conv_tol = 1e-7
#chk_fname = title + '_rohf.chk'
from pyscf.lib import chkfile

mf.chkfile='2Dy.chk'
#mf.init_guess = 'atom'
mf.init_guess = 'chk'
scfdat = chkfile.load(mf.chkfile,'scf')
mf.e_tot = scfdat['e_tot']
mf.mo_coeff = scfdat['mo_coeff']
mf.mo_occ = scfdat['mo_occ']
mf.mo_energy = scfdat['mo_energy']
mf.level_shift = 0.2
#mf.kernel()



print()
print("Enter to LO(IC)-DMET procedure.")
title = 'Dy'

imp_inds = mol.search_ao_label(['Dy.*'])
thres = 1e-12

from embed_sim import aodmet, sacasscf_mixer, siso
mydmet = aodmet.AODMET(mf, title=title, imp_idx=imp_inds, threshold=thres, es_natorb=False).density_fit()
conv_tol = mf.conv_tol
mydmet.build(conv_tol)
#es_mf = mydmet.ROHF()
#es_mf.kernel(mydmet.es_dm)
es_mf = mydmet.es_mf



title = 'Dy'
ncas, nelec, es_mo = mydmet.avas(['Dy 4f'], minao=mol._basis['Dy'], threshold=0.5, openshell_option=2)
print("This is ncas", ncas)
print("This is nelec", nelec)

es_cas = sacasscf_mixer.sacasscf_mixer(es_mf, ncas, nelec, statelis = [0,0,0,0,0,21])
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