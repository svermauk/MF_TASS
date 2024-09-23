gfortran -o reweight.x ../spline/bspline-fortran/src/bspline_oo_module.f90 ../spline/bspline-fortran/src/bspline_sub_module.f90 ../spline/bspline-fortran/src/bspline_kinds_module.f90 ../spline/bspline-fortran/src/bspline_module.f90 reweight_mf.f90

ln -s ../bmtass_0.000/tass/2D_analysis_40ns/COLVAR COLVAR_1
ln -s ../bmtass_0.000/tass/2D_analysis_40ns/HILLS HILLS_1

./reweight.x -T0 300.d0 -T 888.d0 -dT 6300.d0 -tmin 2000 -tmax XXXX -grid -2.983d0 3.140d0 .157d0 -3.1415926d0 3.1415926d0 6.34664646464646459E-002 -pfrqMD 10 -dtMTD 500 -nr 40  -read_vbias -read_ct
