#!/usr/bin/env python
# Perform fragment docking with SILCS-generated maps, fylin 2015

import sys,os
from operator import itemgetter, attrgetter, methodcaller
import molcal

#declare some variables
molcal_data="/home/fylin/sw/molcal/data"
cgenff_parfile=molcal_data+'/'+'par_all36_cgenff.prm'
cgenff_rulefile=molcal_data+'/'+'cgenff.rules'
silcs_rulefile=molcal_data+'/'+'silcs_classification_rules_oct14.dat'

#create charmm force-field object, and read a parfile
cff=molcal.ffcharmm()
cff.read_par_charmm(cgenff_parfile)

#create cgenff interface object and initialize with rulefile and parfile
cgenff=molcal.cgenff_interface()
cgenff.init(cgenff_rulefile,cgenff_parfile)

#create a molecule
mol=molcal.mcmol()

#populate molecule data structures from mol2file
#True argument generates rotors (rotatable bonds) needed for MC simulation
lig="min.<f>"
mol.init_from_mol2("/home/fylin/RingSystems/new_min_frag/minimized_mol2/min.<f>.mol2",True)

#assign CGenFF toppar for mol, first argument is charmm ff object created above
cgenff.assign_toppar(cff,mol)

#assign silcs classes to atoms in molecule GENN,NCLA,etc..
mol.assign_silcs_classes_cgenff(silcs_rulefile,"Silcs")

#create a simulation object
sim=molcal.simulation()

#assign force-field
sim.force_field=cff

#set simulation params
sim.dielectric=4.0
sim.rdielectric=True #variable dielectric
sim.silcs_gfe_cap=3.0
sim.switch_silcs_lgfe_norm=True #normalizes LGFE
sim.switch_silcs_grid=True #always should be true for silcs mc
sim.switch_silcs_excl_grid=True #apply silcs exclusion grid

#add the molecule to simulation
sim.add_molecule(mol,False,True)

FragMaps=[]
one=1.0
FragMaps.append({'name':"BENC", 'file':"1r2b.benc.gfe.map", 'scale':one})
FragMaps.append({'name':"PRPC", 'file':"1r2b.prpc.gfe.map", 'scale':one})
FragMaps.append({'name':"MEOO", 'file':"1r2b.meoo.gfe.map", 'scale':one})
FragMaps.append({'name':"FORN", 'file':"1r2b.forn.gfe.map", 'scale':one})
FragMaps.append({'name':"FORO", 'file':"1r2b.foro.gfe.map", 'scale':one})
FragMaps.append({'name':"MAMN", 'file':"1r2b.mamn.gfe.map", 'scale':one})
FragMaps.append({'name':"ACEO", 'file':"1r2b.aceo.gfe.map", 'scale':one})
FragMaps.append({'name':"AALO", 'file':"1r2b.aalo.gfe.map", 'scale':one})
FragMaps.append({'name':"IMIN", 'file':"1r2b.imin.gfe.map", 'scale':one})
FragMaps.append({'name':"IMIH", 'file':"1r2b.imin2.gfe.map", 'scale':one})
FragMaps.append({'name':"GENN", 'file':"1r2b.genhpb.gfe.map", 'scale':one})
FragMaps.append({'name':"GEND", 'file':"1r2b.gendon.gfe.map", 'scale':one})
FragMaps.append({'name':"GENA", 'file':"1r2b.genacc.gfe.map", 'scale':one})
FragMaps.append({'name':"NCLA", 'file':"1r2b.ncla.map", 'scale':one})

FragMap_Dir="/home/fylin/RingSystems/bcl6_map"
#read and add SILCS FragMaps to the simulation
for gr in FragMaps:
	g=molcal.grid()
	g.read_grid_autodock(FragMap_Dir + "/"+ gr['file'])
	g.category = "Silcs"
	g.atomtype = gr['name']
	g.ScaleFactor=gr['scale']
	sim.add_silcs_grid(g) #this line actually adds the map

#read and add SILCS exclusion map
g=molcal.grid()
g.read_grid_autodock(FragMap_Dir + "/"+ "1r2b.excl.map")
g.category = "Silcs";
g.atomtype = "EXCL"
sim.add_silcs_excl_grid(g) #this line adds exclusion map

#assigns 
sim.assign_silcs_grids(0)
sim.finalize()

#print out info about the simulation you created and modified
sim.query()

#sim.run_md(10000,1000,"md.pdb","md.log","MINI")
center=molcal.vec3(46.85,30.21,35.13)
radius=10.0
radius1=10.0
radius2=2.0

sim.set_center_v3(center)
sim.sim_radius=radius
#sim.sim_radius=radius1
#sim.sim_radius=radius2
sim.switch_confinement=True	

sim.prnlev=3
numruns=1000
move_range1=[1,180,180]
move_range2=[0.2,9,9]

final_conf_pdb="final_"+str(lig)+".pdb"
if(os.path.exists(final_conf_pdb)):
	os.remove(final_conf_pdb)

for irun in range(0,numruns):
	for i in range(0,1000):
		mol.place_in_sphere_rand_orie(center,radius)
		lgfe=sim.get_silcs_gfe_mol(0)
		if(lgfe<0.0):
			break
	sim.sim_radius=radius1
	mol.set_move_probabilities(0.2,0.2,0.6)
	mol.set_move_ranges(move_range1[0],move_range1[1],move_range1[2])
	mol.query()
	pdb="mc."+str(irun)+".pdb"
	log="mc."+str(irun)+".log"
	sim.run_mc(10000,1000,pdb,log,"MC")
	nframes=10000/1000

	mol.populate_conformations_from_pdb(pdb,True)
	min_lgfe=1000000.0
	min_id=-1
	for i in range(0,nframes):
		mol.get_conformation_into_main(i)
		lgfe=sim.get_silcs_gfe_mol(0)
		if(lgfe<min_lgfe):
			min_lgfe=lgfe
			min_id=i
			
	mol.get_conformation_into_main(min_id)
	
	sim.sim_radius=radius2
	mol.set_move_probabilities(0.2,0.2,0.6)
	mol.set_move_ranges(move_range2[0],move_range2[1],move_range2[2])

	pdb="mcmin."+str(irun)+".pdb"
	log="mcmin."+str(irun)+".log"
	sim.run_mc(40000,1000,pdb,log,"MINI")

	lgfe=sim.get_silcs_gfe_mol(0)
	remark="LGFE "+str(lgfe)
	mol.append_pdb(final_conf_pdb,remark)



