%chk=P4model.chk
#T qcisd(tq)/sto-3g scan test

linear D4 model pes scan over varying h2-h3 bond length using ump4 with mbs

0 1
H
H 1 a
H 2 R 1 deltaPlus2pi
H 3 a 2 deltaPlus2pi 1 0.0
  Variables:
R 1.0 60 0.05
a 1.6
deltaPlus2pi 180.0
