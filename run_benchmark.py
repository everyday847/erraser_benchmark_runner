"""
ERRASER benchmark runner, v. 0.1
The idea is, you hand it a folder with a PDB and an MTZ. It figures out your state,
it cleans up your messes, it does analysis.
"""

import argparse
from logging import info, warning, error, basicConfig, INFO
import os
import glob

def erraser2(d: str) -> None:
	"""
	cd into a directory, ensure it's suited for erraser2 to be run there, and
	do so
	"""
	try:
		os.chdir(d)
	except:
		error("Directory {} does not exist; panicking.".format(d))
		return

	#info('in directory {}'.format(d))

	the_pdb: str = glob.glob('*.pdb')[0]
	the_mtz: str = the_pdb.replace('pdb', 'mtz')

	###
	# Phase 1 -- we have at least a pdb and mtz. make maps.
	if len(glob.glob('./*ccp4')) == 0:
		#the_pdb, the_mtz = glob.glob('*.pdb')[0], glob.glob('*.mtz')[0]
		os.system('phenix.maps {} {}'.format(the_pdb, the_mtz))
		if len(glob.glob('./*ccp4')) == 0:
			error("Map creation with Phenix requires manual intervention \
				for {}. Please intervene".format(d))
			return

	try:
		os.mkdir('./output')
	except:
		pass

	###
	# Phase 1.5 -- 'phenix.ready_set'
	if not os.path.exists('./output/phenix_rd1.pdb'):
		#print('phenix.ready_set {}'.format(the_pdb))
		os.system('phenix.ready_set {}'.format(the_pdb))
		os.rename('{}.updated.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_updated.pdb')

	the_cifs = list(glob.glob('[A-Z0-9][A-Z0-9][A-Z0-9].cif')) # how do we select only ligand cifs?
	# if there is a H2U cif, replace with our gold standard
	if 'H2U.cif' in the_cifs:
		the_cifs = [c for c in the_cifs if c != 'H2U.cif']
		warning('replacing automatically generated H2U parameters')
		# provide from somewhere.
		the_cifs.append('h2u_two_plane.cif')

	###
	# Phase 2 -- refine using phenix against the mtz.
	if not os.path.exists('./output/phenix_rd1.pdb'):
		try:
			os.system('phenix.refine {} {} {} --overwrite '.format('./output/phenix_updated.pdb', the_mtz, " ".join(the_cifs)))
			os.rename('phenix_updated_refine_001.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_rd1.pdb')
		except FileNotFoundError:
			warning("didn't find the phenix.refine output file, assuming bad rfree flags")
			os.system('phenix.refine {} {} {} refinement.input.xray_data.r_free_flags.generate=True --overwrite'.format('./output/phenix_updated.pdb', the_mtz, " ".join(the_cifs)))
			os.rename('phenix_updated_refine_001.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_rd1.pdb')
	
	###
	# Phase 2.5 -- cull ions and stuff from PDB. How should we do this?
	# Useful approach -- just get rid of most common stuff that will make
	# our lives the worst for the least gain.
	os.system("cat ./output/phenix_rd1.pdb | grep -v MG | grep -v SR > ./output/phenix_rd1_culled.pdb")

	os.system("pdb2fasta.py ./output/phenix_rd1_culled.pdb > {}".format(the_pdb.replace('pdb', 'fasta')))
	###
	# Phase 3 -- run erraser2 on a culled portion of the phenix pdb.
	#  
	if not os.path.exists('./output/erraser_rd1.pdb'):
		try:
			exe: str = '/Users/amw579/jd3_dev/main/source/bin/erraser2'
			score_flags: str = ' -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior'
			weight_adj: str = ' -set_weights elec_dens_fast 10.0 cart_bonded 5.0 linear_chainbreak 10.0 chainbreak 10.0 fa_rep 1.5 fa_intra_rep 0.5 rna_torsion 10 suiteness_bonus 5 rna_sugar_close 10 atom_pair_constraint 10 '
			mutes: str = '-mute  core.scoring.electron_density.xray_scattering core.pack.rotamer_set.RotamerSet_ core.scoring.CartesianBondedEnergy -render_density true'
			mapfile: str = '{}_2mFo-DFc_map.ccp4'.format(the_pdb.replace('.pdb', ''))
			fasta: str = the_pdb.replace('pdb', 'fasta')
			#-sampler_num_pose_kept 10 for quick?
			os.system('{} -s ./output/phenix_rd1_culled.pdb {} {} {} -rmsd_screen 3.0 -edensity:mapfile {} -fasta {} -allow_virtual_side_chains false -sampler_num_pose_kept 10'.format(
				exe, score_flags, weight_adj, mutes, mapfile, fasta
			))
			os.system('mv ./output/phenix_rd1_culledFINISHED.pdb ./output/erraser_rd1.pdb')
		except:
			quit()
	
	os.chdir('..')

def main(directories: str) -> None:
	for d in directories:
		info("Running ERRASER2 on contents of directory {}.".format(d))
		erraser2(d, execution_fn=write_to_all_sbatch)	

def test_args():
	pass

if __name__ == '__main__':

	basicConfig(level=INFO)

	parser = argparse.ArgumentParser(description='Run ERRASER2 on a working directory')
	parser.add_argument('directory', type=str, nargs='+',
						help='dir(s) to run')
	args = parser.parse_args()
	main(args.directory)