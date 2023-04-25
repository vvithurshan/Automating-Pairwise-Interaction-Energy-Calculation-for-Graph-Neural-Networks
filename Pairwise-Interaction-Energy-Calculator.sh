mkdir CHARMM-GRINN
cd CHARMM-GRINN
I=1
for pdb in ../PDB_CATNAP/*.pdb; do
	pdb_name=$(basename "$pdb" .pdb)
	mkdir $pdb_name
	cp ../PDB_CATNAP/$pdb ./$pdb_name
	cd $pdb_name

	##psf

## define this 

dir=/home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter
namd2=/home/vvithurshan/0014/vvarenthirarajah/Softwares/NAMD_2.13_Linux-x86_64-multicore/namd2
grinn=/home/vvithurshan/0014/vvarenthirarajah/Softwares/grinn_linux_v110_hf1/grinn
namd3=/home/vvithurshan/0014/vvarenthirarajah/Documents/NAMD_3.0alpha12_Linux-x86_64-multicore-CUDA/namd3
	

para="
parameters $dir/par_all36_carb.prm
parameters $dir/par_all36_cgenff.prm
parameters $dir/par_all36_lipid.prm
parameters $dir/par_all36m_prot.prm
parameters $dir/par_all36_na.prm
parameters $dir/par_interface.prm
parameters $dir/toppar_water_ions_namd.str"


cat <<EOF > psf_maker.tcl
mol load pdb ${pdb_name}.pdb

set A [atomselect top "chain A and protein"]
set B [atomselect top "chain B and protein"]
set C [atomselect top "chain C and protein"]

\$A writepdb A.pdb
\$B writepdb B.pdb
\$C writepdb C.pdb

resetpsf
package require psfgen
topology $dir/top_all36_prot.rtf
pdbalias residue HIS HSD
segment A {pdb A.pdb}
segment B {pdb B.pdb}
segment C {pdb C.pdb}

coordpdb A.pdb A
coordpdb B.pdb B
coordpdb C.pdb C

guesscoord
writepsf PSF_${pdb_name}.psf
writepdb pdb_${pdb_name}.pdb

EOF


	vmd -dispdev text -eofexit < psf_maker.tcl > psf_maker.log

cat <<EOF > solvate.tcl
package require solvate
solvate PSF_${pdb_name}.psf pdb_${pdb_name}.pdb -o solvate -s WT -x 15 -y 15 -z 15 +x 15 +y 15 +z 15 -b 2.4
EOF

	vmd -dispdev text -eofexit < solvate.tcl > solvate_c.log

cat <<EOF > remove_water.tcl
mol load psf solvate.psf pdb solvate.pdb
set prot [atomselect top "protein and not water"]
\$prot writepdb PDB_${pdb_name}.pdb
EOF

	vmd -dispdev text -eofexit < remove_water.tcl > remove_water.log

# 	cat <<EOF > wrap.tcl

# # Load the PSF and PDB files
# mol load psf PSF_${pdb_name}.psf pdb unwrapped_PDB_${pdb_name}.pdb

# pbc wrap -center com -centersel "protein" -compound residue -all
# set sel [atomselect top "protein"]
# \$sel moveby [vecscale -1 [$sel center]]

# set outfile [open PDB_${pdb_name}.pdb w]
# \$sel writepdb \$outfile
# \$sel delete
# close \$outfile

# EOF
# 	vmd -dispdev text -eofexit < wrap.tcl > wrap.log

	##restraint

	cat <<EOF > restraint.tcl
mol load psf PSF_${pdb_name}.psf pdb PDB_${pdb_name}.pdb
mkdir restraint
set all [atomselect top all]
\$all set beta 0
set protein_1 [atomselect top "backbone"]
\$protein_1 set beta 1
\$all writepdb protein_1.pdb
mv protein_1.pdb ./restraint

EOF

	vmd -dispdev text -eofexit < restraint.tcl > restraint.log

cat <<EOF > boxsize.tcl
mol load psf PSF_${pdb_name}.psf pdb PDB_${pdb_name}.pdb
proc get_cell {{molid top}} {
    set outfile [open "boxsize.txt" w]
    set all [atomselect \$molid all]
    set minmax [measure minmax \$all]
    set vec [vecsub [lindex \$minmax 1] [lindex \$minmax 0]]
    puts \$outfile "cellBasisVector1 [lindex \$vec 0] 0 0"
    puts \$outfile "cellBasisVector2 0 [lindex \$vec 1] 0"
    puts \$outfile "cellBasisVector3 0 0 [lindex \$vec 2]"
    set center [measure center \$all]
    puts \$outfile "cellOrigin \$center"
    \$all delete
    close \$outfile
}
get_cell 
EOF
[[ -f ./boxsize.txt ]] && echo "boxsize.txt file created" >>../CHARMM-GRINN.log|| echo "boxsize.txt file not created" >>../CHARMM-GRINN.log

	vmd -dispdev text -eofexit < boxsize.tcl > boxsize.log


cat <<EOF > modify_box.py
with open("boxsize.txt") as f, open("boxsize_modified.txt", 'w') as ff:
    cnt = 0
    for line in f:
        if cnt < 3:
            parts = line.split()
            value = float(parts[cnt+1]) + 30
            parts[cnt+1] = str(value)
            out = ' '.join(parts)
            ff.write(out + '\n')
        else:
            ff.write(line)
        cnt += 1
EOF

python3 modify_box.py

	boxsize_contents=$(cat boxsize_modified.txt)

	### minimization
	cat <<EOF > Minimization.conf


############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

# Simulation conditions
temperature 0


CUDASOAintegrate off
# Harmonic constraints

constraints on
consref ./restraint/protein_1.pdb
conskfile ./restraint/protein_1.pdb
constraintScaling 10
consexp 2
conskcol B


# Output Parameters
set output Minimization
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 500
XSTFreq 500
restartfreq 500
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on

$para


exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on



#boundary
$boxsize_contents


# Script

minimize 1000


###########
EOF

cat <<EOF > Minimization2.conf


############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Minimization
binCoordinates \$input.restart.coor
#binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate off

# Harmonic constraints

constraints on
consref ./restraint/protein_1.pdb
conskfile ./restraint/protein_1.pdb
constraintScaling 5
consexp 2
conskcol B


# Output Parameters
set output Minimization2
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters

usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpiston off
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on

set dir dirr

$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on


minimize 5000


###############

EOF


cat <<EOF > Annealing.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Minimization2
binCoordinates \$input.restart.coor
#binVelocities \$inputname.restart.coor
extendedSystem \$input.restart.xsc

# Simulation conditions
temperature 0

CUDASOAintegrate on

# Harmonic constraints

constraints on
consref ./restraint/protein_1.pdb
conskfile ./restraint/protein_1.pdb
constraintScaling 2
consexp 2
conskcol B

# Output Parameters
set output Annealing

binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 10000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 0
langevinHydrogen    off
langevindamping 1

# Barostat Parameters


langevinpiston off
usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on

# Script
set Temp 300
set barostat 0
set nSteps  500
for {set t 0} {\$t <= \$Temp} {incr t} {run \$nSteps;langevintemp \$t;if {\$barostat} {langevinpistontemp \$t}}


####################

EOF


## Equilibration
cat <<EOF > Equilibration.conf

############################################################################
##START HERE###
##Simulation Template##
# Simulation conditions
coordinates PDB_${pdb_name}.pdb
structure PSF_${pdb_name}.psf

set input Annealing

binCoordinates \$input.restart.coor
binVelocities \$input.restart.vel
extendedSystem \$input.restart.xsc

# Simulation conditions

CUDASOAintegrate on

# Harmonic constraints

constraints on
consref ./restraint/protein_1.pdb
conskfile ./restraint/protein_1.pdb
constraintScaling 1
consexp 2
conskcol B

# Output Parameters
set output MD_0
binaryoutput no
outputname \$output
outputenergies 100
outputtiming 100
outputpressure 100
binaryrestart yes
dcdfile \$output.dcd
dcdfreq 1000
XSTFreq 1000
restartfreq 1000
restartname \$output.restart


# Thermostat Parameters
langevin on
langevintemp 300
langevinHydrogen    off
langevindamping 1

# Barostat Parameters


langevinpiston off
usegrouppressure yes
useflexiblecell no
useConstantArea no
langevinpistontarget 1
langevinpistonperiod 200
langevinpistondecay 100
langevinpistontemp 300

# Integrator Parameters

timestep 2
firstTimestep 0
fullElectFrequency 2
nonbondedfreq 1
stepspercycle 10

# Force Field Parameters

paratypecharmm on
#set dir /home/vithurshan/2021/toppar_c36_jul20


$para

exclude scaled1-4
1-4scaling 1.0
rigidbonds all
cutoff 12.0
pairlistdist 14.0
switching on
switchdist 10.0
PME on
PMEGridspacing 1
wrapAll on
wrapWater on


run 500000 
## for 1 ns


###################
EOF

[[ -f ./MD_0.dcd ]] && echo "MD_0.dcd file created" >>../CHARMM-GRINN.log|| echo "MD_0.dcd file not created" >>../CHARMM-GRINN.log

cat <<EOF > bash_for_gpu.sh

ppn=4
pe=0-3


echo "Minimization" >>../CHARMM-GRINN.log;
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Minimization.conf > Minimization1.log

echo "Minimization2" >>../CHARMM-GRINN.log;
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Minimization2.conf > Minimization2.log

for (( min=1; min<=2; min++))
do
	mkdir min_\$min
	mv Minimization\$min.log ./min_\$min
	cd min_\$min
	source /home/vvithurshan/0014/vvarenthirarajah/Documents/2021_Reseach/scripts/energy.sh
	cd ../../
done

echo "Annealing" >>../CHARMM-GRINN.log;
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe   Annealing.conf > Annealing.log

echo "Equilibration" >>../CHARMM-GRINN.log;
\$namd3 +p\${ppn} +idlepoll +setcpuaffinity   +devices 0  +pemap \$pe  Equilibration.conf > Equilibration.log


echo "Done" >>../CHARMM-GRINN.log;

EOF

source bash_for_gpu.sh


## align frames

cat <<EOF > align.tcl

mol load psf PSF_${pdb_name}.psf dcd MD_0.dcd

# Get the number of frames
set numframes [molinfo top get numframes]

# Get the selection of backbone atoms
set sel [atomselect top "backbone"]

# Loop over all frames
for {set i 0} {\$i < \$numframes} {incr i} {

  # Set the current frame
  animate goto \$i

  # Align the selection to the first frame
  measure fit \$sel \$sel frame 0

}

# Write the aligned coordinates to a new DCD file
animate write dcd MD_0_aligned.dcd sel \$sel

EOF

	# vmd -dispdev text -eofexit < align.tcl > align.log


cat <<EOF > grinn.sh

##grinn
\$grinn -calc --pdb PDB_${pdb_name}.pdb --top PSF_${pdb_name}.psf --traj MD_0_aligned.dcd --exe \$namd2 --stride 20  --outfolder grinn_output --numcores 200 --sel1 'chain A B' --sel2 'chain C' --parameterfile /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_all36_carb.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_all36_cgenff.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_all36_lipid.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_all36_na.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_all36m_prot.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/par_interface.prm /home/vvithurshan/0014/vvarenthirarajah/Documents/Parameter/toppar_water_ions_namd.str

EOF
source grinn.sh

cd ..
done



