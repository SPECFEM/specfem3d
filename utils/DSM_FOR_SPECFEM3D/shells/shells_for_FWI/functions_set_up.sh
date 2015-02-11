function create_link ()
{
#
## $1 : tractions ref
## $2 : numero source
#
cd $TRACTION/$EARTHQUAKE$2
rm *vel.bin *tract.bin *tract.indx

ls $1/*vel.bin > liste.txt
awk '{print "ln -s " $0} '  liste.txt  > commande_linux.sh
bash ./commande_linux.sh

ls $1/*tract.bin > liste.txt
awk '{print "ln -s " $0} '  liste.txt  > commande_linux.sh
bash ./commande_linux.sh

rm liste.txt
rm commande_linux.sh

cd ..
}
