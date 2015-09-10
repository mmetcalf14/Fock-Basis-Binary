#/bin/bash
echo "to " $1
rsync -avzhe ssh Fock\ Basis\ Binary/* edge0:~/GitRepo/ed/
