
cwd=$(pwd)

##Save Current Path

PATHVAR=$PATH
PYLITHDIR="$HOME/pylith/pylith-2.2.2-linux-x86_64"
CUBITDIR="$HOME/Coreform-Cubit-2022.4/bin"/


echo "Starting Simulation ##########################################"
echo " "
echo " "
echo " "
echo "                                         _.oo.         "
echo "                 _.u[[/;:,.         .odMMMMMM          "
echo "              .o888UU[[[/;:-.  .o@P^    MMM^           "
echo "             oN88888UU[[[/;::-.        dP^             "
echo "            dNMMNN888UU[[[/;:--.   .o@P^               "
echo "           ,MMMMMMN888UU[[/;::-. o@^                   "
echo "           NNMMMNN888UU[[[/~.o@P^                      "
echo "           888888888UU[[[/o@^-..                       "
echo "          oI8888UU[[[/o@P^:--..                        "
echo "       .@^  YUU[[[/o@^;::---..                         "
echo "     oMP     ^/o@P^;:::---..                           "
echo "  .dMMM    .o@^ ^;::---...                             "
echo " dMMMMMMM@^`       `^^^^                               "
echo "YMMMUP^                                                "
echo "^^                                                     " 
echo " "
echo " "
echo " "
source shell/fix_path.sh

python Specify_Params.py

cd $cwd
