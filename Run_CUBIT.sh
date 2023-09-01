

yes | cp -rf Make_geometry.jou $CUBITDIR
cd $CUBITDIR
./coreform_cubit -batch -initfile Make_geometry.jou -nographics 
yes | cp -rf geometry.exo $cwd
cd $cwd


