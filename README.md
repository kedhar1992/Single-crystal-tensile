# Single-crystal-tensile
#Twin identification

#Kedharnath & Bharat
#Tensile simulation of single crystal Zr

#Basic simulation conditions
units metal
dimension 3
boundary p p p

atom_style atomic
neighbor 2.0 bin
neigh_modify delay 10 check yes

####################################################################################
#Sample part

#Create simulation box
region whole block -1000 1000 -1000 1000 -1000 1000 units box
create_box 3 whole

#Inputing equilibrium single crystal
read_dump dump.zr 0 x y z box yes add yes

#setting the atomic mass (Zr-91.224 amu)
mass 1 91.224

pair_style eam/fs
pair_coeff * * Zr_3.eam.fs Zr Zr Zr



#Print the initial details of the system
variable xlen equal lx
variable ylen equal ly
variable zlen equal lz
print "lx: ${xlen}"
print "ly: ${ylen}"
print "lz: ${zlen}"

# EQUILIBRATION
reset_timestep	0
timestep 0.001
# atoms are given a random velocity based on a temperature of 100K.
thermo 100
thermo_style custom step lx ly lz press pxx pyy pzz pe ke temp
velocity all create 100 123456 mom yes rot no
# temperature and pressure are set to 100 and 0
fix 1 all npt temp 100 100 1 aniso 0 0 1 drag 1

run 25000

dump 1 all atom  100 dump.eq
unfix 1

            reset_timestep 0 
            thermo 100
            thermo_style custom step pe lx ly lz press pxx pyy pzz 
            min_style cg
            minimize 1e-15 1e-15 10000 100000


#Compute 
compute pot all pe
compute peratom all pe/atom
compute kin all ke
compute temper all temp

#Start deformation
reset_timestep 0
timestep 0.0001
#Perform deformation for every timestep along x direction with engg strain rate 0.1/fs
fix 1 all deform 100 z erate 0.01 units box remap x
#Reset the pressure in y and z direction to zero for every 10 timestep
fix 2 all npt temp 100 100 1 x 0 0 1 y 0 0 1

compute 1 all pressure temper
variable p1 equal "-pxx/10000"
variable p2 equal "-pyy/10000"
variable p3 equal "-pzz/10000"
variable p4 equal "-pxy/10000"
variable p5 equal "-pxz/10000"
variable p6 equal "-pyz/10000"
variable potential equal "pe*1.6021*(10^-16)"
variable kinetic equal "ke*1.6021*(10^-16)"

#Calculate the strain in x, y, z directions
variable tmp equal "lx"
variable L1 equal ${tmp}
print "Initial Length x, L1: ${L1}"
variable strainx equal "(lx - v_L1)/v_L1"

fix press1 all print 100 "${strainx} ${p1} ${p2} ${p3} ${p4} ${p5} ${p6} ${potential} ${kinetic}" file Zr_press.txt screen no

dump 2 all custom 100 dump.* id type xs ys zs fx fy fz

#Output thermodynamic info for every 100 timestep
timestep 0.0001
thermo 100
thermo_style custom step lx ly lz press pxx pyy pzz pe ke temp
#Run the simulation/deformation for 50 ps
run 200000
print "Deformation successful"

