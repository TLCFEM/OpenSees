wipe
model BasicBuilder -ndm 2 -ndf 2

node 1 0 0
node 2 1 0

uniaxialMaterial Elastic 1 100

recorder Node -file U1.txt -dT 0.01 -node 2 -dof 1 disp

element truss 1 1 2 1 1

mass 2 100 0

fix 1 1 1
fix 2 0 1

timeSeries Path 1 -dt .02 -values {1. 1. 0.}

pattern UniformExcitation 1 1 -accel 1  -fact 1

constraints Plain
numberer Plain
algorithm Newton
system LeeSparse
integrator LeeNewmarkFullKC 0.5 0.25 -type4 1 .05 1 1 1 1 2.
analysis Transient

test NormDispIncr 1E-12 20 1

analyze 2000 0.01