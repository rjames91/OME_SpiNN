DIR=$(pwd)/../
cd $DIR/spinnaker_tools
source $(pwd)/setup	
cd $DIR/sPyNNaker/neural_modelling
source $(pwd)/setup.bash
#source $(pwd)/setup
cd $DIR/DRNL_SpiNN/c_models
make clean
make || exit $?
cd $DIR/IHC_AN_SpiNN/c_models
cd ./AN_group_node
make clean
make || exit $?
cd ../IHC_AN_float_MAP_14Copy
make clean
make || exit $?
cd $DIR/OME_SpiNN/c_models
cd ./MCack
make clean
make || exit $?
cd ../OME
make clean
make || exit $?
cd $DIR/OME_SpiNN/model_binaries
cp ./SpiNNakEar_* $DIR/sPyNNaker/spynnaker/pyNN/model_binaries
