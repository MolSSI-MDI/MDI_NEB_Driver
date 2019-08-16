#location of required codes
DRIVER_LOC=$(cat ../locations/MDI_NEB_Driver)
LAMMPS_LOC=$(cat ../locations/LAMMPS)
#NENGINES=$1
#SPRING_CONST=$2
#ENERGY_THRESHOLD=$3
#FORCE_THRESHOLD=$4
NENGINES=5
SPRING_CONST=0.00001
ENERGY_THRESHOLD=0.00001
FORCE_THRESHOLD=0.00001

#remove old files
if [ -d work ]; then 
  rm -r work
fi

#create work directory
cp -r data work
cd work

#launch QE
#${QE_LOC} -mdi "-role ENGINE -name QM -method TCP -port 8021 -hostname localhost" -in qe.in > qe.out &

for i in $(seq 1 ${NENGINES})
do
	echo "Iteration $i"
	#launch LAMMPS
	${LAMMPS_LOC} -mdi "-role ENGINE -name MM$i -method TCP -port 8021 -hostname localhost" -in $i.in > $i.out &
done
#launch drive
${DRIVER_LOC} -mdi "-role DRIVER -name driver -method TCP -port 8021" -engines ${NENGINES} -spring ${SPRING_CONST} -energy_threshold ${ENERGY_THRESHOLD} -force_threshold ${FORCE_THRESHOLD}  &

wait
