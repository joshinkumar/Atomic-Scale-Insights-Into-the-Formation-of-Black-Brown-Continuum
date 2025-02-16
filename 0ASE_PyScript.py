from ase.calculators.vasp import Vasp
from ase.io.vasp import read_vasp
import os
import shutil
import subprocess

# Create directory function
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: Creating directory. ' + directory)

# Read initial POSCAR
model = read_vasp('POSCAR')

############### Relaxing (geometry optimization) at 300 K ## Added Step (12/18/2024)
# Full relax=> ibrion=2, isif=3, nbands&ediff&nedos=> same for relax & optics;  lreal=False for high accuracy; tebeg=0, teend=0, 
createFolder('5_Relaxing')


calc = Vasp(directory = '5_Relaxing', xc = 'PBE', kpts = (1,1,1), system = 'c120_15ps', prec = 'Normal', ediff = 1e-5,
            ismear = 0, sigma = 0.05, algo = 'Normal', maxmix = 40, isym = 0, nelm=40, nelmin = 3, nelmdl=-10,
            ibrion = 0, nsw = 500, potim = 2, nwrite = 0, nblock = 10, lcharg = False,
            lwave = False, lorbit=0, tebeg = 0, teend = 8000, mdalgo = 2, smass = 0, isif = 0,
            ncore = 8, kpar=1, encut=400, lreal='Auto') #npar=8 ##, maxmem=1792, ncshmem=4
model.calc = calc

calc.set(directory = '5_Relaxing', tebeg = 8000, teend = 8000, potim = 2, nsw = 5000, sigma=0.1)
model.calc = calc
calc.set(directory = f'5_Relaxing', tebeg = 8000, teend = 600, potim = 2, nsw = 7500)
model.calc = calc
calc.set(directory = '5_Relaxing', tebeg = 300, teend = 300, potim = 2, nsw = 500)
model.calc = calc
calc.set(directory='5_Relaxing', ibrion=2, isif=3, nsw=1000, potim=1, lorbit=11, lwave=True, prec='Accurate', lreal=False, lscalapack=False, ediff=1e-6, ediffg=1e-3)
model.calc = calc

print("5_Relaxing: ", model.get_potential_energy())
os.rename('5_Relaxing/XDATCAR',f'5_Relaxing/XDATCAR_5_Relaxing')
subprocess.Popen(f'cd 5_Relaxing && nohup python $Plot_MD &', shell=True)

# Read relaxed POSCAR (new model reading not required)
#model = read_vasp('5_Relaxing/CONTCAR')



################ Relaxing (geometry optimization) at 300 K
## Full relax=> ibrion=2, isif=3, nbands&ediff&nedos=> same for relax & optics;  lreal=False for high accuracy; tebeg=0, teend=0, 
#createFolder('5_Relaxing')
#
#calc = Vasp(directory = '1_Pre_Heating', xc = 'PBE', kpts = (1,1,1), system = 'c120_15ps', prec = 'Normal', ediff = 1e-5,
#            ismear = 0, sigma = 0.05, algo = 'Normal', maxmix = 40, isym = 0, nelm=40, nelmin = 3, nelmdl=-10,
#            ibrion = 0, nsw = 500, potim = 2, nwrite = 0, nblock = 10, lcharg = False,
#            lwave = False, lorbit=0, tebeg = 0, teend = 8000, mdalgo = 2, smass = 0, isif = 0,
#            ncore = 8, kpar=1, encut=400, lreal='Auto') #npar=8 ##, maxmem=1792, ncshmem=4
#
#model.calc = calc
#
#calc.set(directory='5_Relaxing', ibrion=2, isif=2, nsw=1000, potim=1, lwave=True, prec='Accurate', lreal=False, lscalapack=False, ediff=1e-6, ediffg=1e-3)
#
#print("5_Relaxing: ", model.get_potential_energy())
#os.rename('5_Relaxing/XDATCAR',f'5_Relaxing/XDATCAR_5_Relaxing')
#subprocess.Popen(f'cd 5_Relaxing && nohup python $Plot_MD &', shell=True)




############### Optics Script Instructions:
# Depending upon the material being simulated (POSCAR file), change the following tags everywhere in the script:
# 1. nbands (keep 3 times of the original OUTCAR value) 2. ncore 3. kpar

# Global variables:
ENCUT = 400
NBANDS = 1024 # OUTCAR value = 312 | 288 DNW -> NBANDS needs to be factor of 64 | *
NCORE = 8
KPAR = 1
#EMIN = 0.4132823 # Ahuja(2004) would have commented EMIN and EMAX
#EMAX = 4.1328230
NEDOS = 8000
OMEGAMIN = 0.4132823 # 1.0332058 -> 1200nm | 0.4132823->3000nm
OMEGAMAX = 4.1328230 #for 300nm
OMEGATL = 30
LORBIT = 11

############### 1_Standard_DFT: Whatever we may choose to do afterwards in terms of dielectric response calculations, we have to start with a standard DFT (or hybrid functional) calculation
createFolder('6_Standard_DFT')

calc = Vasp(directory = '6_Standard_DFT', xc='PBE', ismear=0, sigma=0.01, ediff=1e-8, kpts=(4, 4, 4), gamma=True, 
            ncore=NCORE, kpar=KPAR, nbands=NBANDS, lorbit=LORBIT, encut=ENCUT, lscalapack=False)
model.calc = calc

print("6_Standard_DFT: ", model.get_potential_energy())
subprocess.Popen(f'cd 6_Standard_DFT && nohup python $Plot_MD &', shell=True)

############### 2_IP_LOPTICS: To compute the frequency dependent dielectric function in the independent-particle (IP) picture we restart from the WAVECAR of the previous run, with the following INCAR. 
createFolder('7_IP_LOPTICS')

# COPY WAVECAR FROM PREVIOUS RUN HERE.
shutil.copy('6_Standard_DFT/WAVECAR', '7_IP_LOPTICS/WAVECAR')

calc = Vasp(directory = '7_IP_LOPTICS', xc='PBE', algo='Normal', loptics=True, cshift=0.1, ismear=0, sigma=0.01, ediff=1e-8, kpts=(4, 4, 4), gamma=True,
            nedos=NEDOS, ncore=NCORE, kpar=KPAR, nbands=NBANDS, encut=ENCUT, lscalapack=False)
model.calc = calc

print("7_IP_LOPTICS: ", model.get_potential_energy())
subprocess.Popen(f'cd 7_IP_LOPTICS && nohup python $Plot_MD &', shell=True)
subprocess.Popen(f'cd 7_IP_LOPTICS && echo -e "71\n711\n1\n" | vaspkit && nohup python $Plot_Optics &', shell=True)

################ 3_LOCAL_RPA_CHI: To determine the frequency dependent dielectric function including local field effects one needs the WAVECAR and WAVEDER files from the previous calculation (ALGO=Exact and LOPTICS=.TRUE., and sufficient virtual orbitals). Per default, for ALGO=CHI, local field effects are included at the level of the RPA (LRPA=.TRUE.), i.e., limited to Hartree contributions only.
#createFolder('8_LOCAL_RPA_CHI')
#
## COPY WAVECAR and WAVEDER FROM PREVIOUS RUN HERE.
#shutil.copy('7_IP_LOPTICS/WAVECAR', '8_LOCAL_RPA_CHI/WAVECAR')
#shutil.copy('7_IP_LOPTICS/WAVEDER', '8_LOCAL_RPA_CHI/WAVEDER')
#
#calc = Vasp(directory = '8_LOCAL_RPA_CHI', xc='PBE', algo='CHI', ismear=0, sigma=0.01, ediff=1e-8, lwave=False, lcharg=False, kpts=(4, 4, 4), gamma=True,
#            omegamin=OMEGAMIN, omegamax=OMEGAMAX, omegatl=OMEGATL, nbands=NBANDS, encut=ENCUT, lscalapack=False) #maxmem=(256*1024/64) - 200
#            # Don't use KPAR and NCORE (especially) here!
#            # Don't use nedos=NEDOS, emax=EMAX, emin=EMIN  here!
#model.calc = calc
#
#
#print("8_LOCAL_RPA_CHI: ", model.get_potential_energy())
#subprocess.Popen(f'cd 8_LOCAL_RPA_CHI && nohup python $Plot_MD &', shell=True)
#subprocess.Popen(f'cd 8_LOCAL_RPA_CHI && echo -e "71\n711\n1\n" | vaspkit && nohup python $Plot_Optics &', shell=True)