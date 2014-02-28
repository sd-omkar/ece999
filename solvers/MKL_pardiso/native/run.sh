#!/bin/bash

list1=(2cubes_sphere 2D_54019_highK 31770-lhs 88950-lhs a2nnsnsl a5esindl ABACUS_shell_ud af23560 ancfBigDan Andrews apache1 ASIC_100k ASIC_100ks av41092 bauru5727 bayer01 bcircuit blockqp1 boyd1)

list2=(c-59 c-61 c-62 cant case39 case39_A_01 cfd1 cfd2 circuit_4 ckt11752_tr_0 cond-mat-2005 cont-201 dawson5 dixmaanl Dubcova2 ecl32 epb3 ex11 F2 filter3D finan512 g7jac140 Ga3As3H12 GaAsH6 garon2)

list3=(gas_sensor gearbox gridgena H2O hcircuit ibm_matrix_2 jan99jac120 juba40k k1_san laminar_duct3D lhr10c lhr71 Lin lung2 mark3jac100 mark3jac140 matrix_9 minsurfo ncvxbqp1 oilpan olesnik0)

list4=(OPF_10000 pdb1HYS poisson3Db qa8fm rail_79841 rajat26 rma10 shallow_water1 shallow_water2 sparsine stomach t3dh_a thermal1 TSOPF_FS_b162_c4 TSOPF_FS_b39_c19 vanbody venkat25 vibrobox xenon1)

for mat in "${list1[@]}"; do
  echo Copying $mat.mtx
  scp $mat.mtx mic0:~/matrices/
  echo Runnig $mat.mtx
  micnativeloadex ./pardiso.native -a "/home/odeshmukh/matrices/${mat}.mtx" > $mat.log
done
