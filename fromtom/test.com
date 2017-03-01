%mem=2GB
%nprocs=2
# RHF/GEN
  SCF(Tight,NoVarAcc,NoIncFock)

junk

0,1
H        0.000     0.000      1.000
H        1.732     0.000      0.000
H        0.000     0.000     -1.000
H       -1.732     0.000      0.000

@H.gbs

