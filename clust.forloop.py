#!/usr/bin/env python

import sys,os
from operator import itemgetter, attrgetter, methodcaller
import molcal

#declare some variables
molcal_data="/home/prabhu/sw/molcal/data"
cgenff_parfile=molcal_data+'/'+'par_all36_cgenff.prm'
cgenff_rulefile=molcal_data+'/'+'cgenff.rules'
silcs_rulefile=molcal_data+'/'+'silcs_classification_rules_oct14.dat'

#create charmm force-field object, and read a parfile
cff=molcal.ffcharmm()
cff.read_par_charmm(cgenff_parfile)

#create cgenff interface object and initialize with rulefile and parfile
cgenff=molcal.cgenff_interface()
cgenff.init(cgenff_rulefile,cgenff_parfile)



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


for frag in range(1,101):
	#create a molecule
	mol=molcal.mcmol()

	#populate molecule data structures from mol2file
	#True argument generates rotors (rotatable bonds) needed for MC simulation
	lig="min."+str(frag)
        ligmol2=lig+".mol2"
        mini_dir = "/home/fylin/RingSystems/new_min_frag/minimized_mol2/"
	mol.init_from_mol2(mini_dir+"/"+str(ligmol2),True)

	#assign CGenFF toppar for mol, first argument is charmm ff object created above
	cgenff.assign_toppar(cff,mol)

	#assign silcs classes to atoms in molecule GENN,NCLA,etc..
	mol.assign_silcs_classes_cgenff(silcs_rulefile,"Silcs")

	#add the molecule to simulation
	sim.add_molecule(mol,False,True)

	#assigns 
	sim.assign_silcs_grids(0)
	sim.finalize()
	center=molcal.vec3(46.85,30.21,35.13)
	radius=5.0
	sim.set_center_v3(center)
	sim.sim_radius=radius
	sim.switch_confinement=True	

	sim.prnlev=3

	numruns=1000    ##original 1000
        final_name="final_"+str(lig)+".pdb"
        sim_dir="../new_bcl6_sim10/"
	final_conf_pdb=sim_dir+"fsim_"+str(frag)+"/"+str(final_name)
	#now cluster
	mol.populate_conformations_from_pdb(final_conf_pdb,True)
	cluster_assignment=molcal.int_vector()
	mol.cluster_conformations_com(5,cluster_assignment)

	unique_cluster_ids=[]
	unique_clusters=[]
	for cl in cluster_assignment:
		if cl not in unique_cluster_ids:
			unique_cluster_ids.append(cl)
			unique_clusters.append({'id':cl,'n':0})

	for cl in cluster_assignment:
		for ucl in unique_clusters:
			if(cl==ucl['id']):
				ucl['n']=ucl['n']+1

	unique_min_lgfe=[]
	distinct_cluster_representatives=[]
	for ucl in unique_clusters:
		min_lgfe=1000000.0
		min_id=-1
		for run in range(0,numruns):
			cl=cluster_assignment[run]
			if(ucl['id']==cl):
				mol.get_conformation_into_main(run)
				lgfe=sim.get_silcs_gfe_mol(0)
				if(lgfe<min_lgfe):
					min_lgfe=lgfe
					min_id=run
		#print ucl,cluster_assignment[min_id],min_lgfe,min_id
		distinct_cluster_representatives.append({'cl_id':ucl['id'],'cl_n':ucl['n'],'conf_id':min_id,'min_lgfe':min_lgfe})

	distinct_cluster_representatives=sorted(distinct_cluster_representatives,key=itemgetter('min_lgfe'))

	selected_conf_pdb="selected_"+str(lig)+".pdb"
	if(os.path.exists(selected_conf_pdb)):
		os.remove(selected_conf_pdb)

	for cl in distinct_cluster_representatives:
		mol.get_conformation_into_main(cl['conf_id'])
		remark="LGFE "+str(cl['min_lgfe'])+" N "+str(cl['cl_n'])
		mol.append_pdb(selected_conf_pdb,remark)
	sim.delete_molecule(0)
