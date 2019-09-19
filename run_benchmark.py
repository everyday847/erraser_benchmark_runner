"""
ERRASER benchmark runner, v. 0.1
The idea is, you hand it a folder with a PDB and an MTZ. It figures out your state,
it cleans up your messes, it does analysis.
"""

import argparse
from logging import info, warning, error, basicConfig, INFO
import os
import glob
import functools # for partial
import shutil # for rmtree
from typing import Callable, List
from helpers import sbatch_header

def erraser2(d: str, execution_fn: Callable, nstruct: int) -> None:
    """
    cd into a directory, ensure it's suited for erraser2 to be run there, and
    do so
    """

    #'/home/groups/rhiju/amw579/Rosetta/main/source/bin/erraser2.default.linuxgccrelease' #'$Gsource/bin/erraser2'
    exe: str = '{}/main/source/bin/erraser2'.format(os.environ['ROSETTA'])


    try:
        os.chdir(d)
    except:
        error("Directory {} does not exist; panicking.".format(d))
        return

    #info('in directory {}'.format(d))

    the_pdb: str = glob.glob('*.pdb')[0]
    the_mtz: str = the_pdb.replace('pdb', 'mtz')

    info("About to run phase one: mapmaking.")
    ###
    # Phase 1 -- we have at least a pdb and mtz. make maps.
    if len(glob.glob('./*ccp4')) == 0:
        info("Map not yet generated, making it.")
        #the_pdb, the_mtz = glob.glob('*.pdb')[0], glob.glob('*.mtz')[0]
        # We don't use execution fn here b/c you do this ONCE.
        os.system('phenix.maps {} {} > phenix_maps.log'.format(the_pdb, the_mtz))
        if len(glob.glob('./*ccp4')) == 0:
            error("Map creation with Phenix requires manual intervention \
                for {}. Please intervene".format(d))
            return

    try:
        os.mkdir('./output')
    except:
        pass

    info("About to run phase one point five: ready-set.")
    ###
    # Phase 1.5 -- 'phenix.ready_set'
    if not os.path.exists('./output/phenix_updated.pdb'):
        info("Ready-set PDB not yet created; making it.")
        #print('phenix.ready_set {}'.format(the_pdb))
        # Ditto: do once.
        os.system('phenix.ready_set {}'.format(the_pdb))
        try:
            shutil.move('{}.updated.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_updated.pdb')
        except:
            error('moving {}.updated.pdb to output/phenix_updated.pdb failed, likely fatal error'.format(the_pdb.replace('.pdb', '')))

    the_cifs = list(glob.glob('[A-Z0-9][A-Z0-9][A-Z0-9].cif')) # how do we select only ligand cifs?
    # if there is a H2U cif, replace with our gold standard
    if 'H2U.cif' in the_cifs:
        the_cifs = [c for c in the_cifs if c != 'H2U.cif']
        warning('replacing automatically generated H2U parameters')
        # provide from somewhere.
        #os.system('cp {}/h2u_two_plane.cif ./output/'.format(os.path.dirname(os.path.realpath(__file__))))
        os.system('cp {}/h2u_two_plane.cif ./'.format(os.path.dirname(os.path.realpath(__file__))))
        the_cifs.append('h2u_two_plane.cif')


    ###
    # Phase 1.75: initial molprobity analysis.
    info("About to run phase one point seven five: molprobity.")
    if not os.path.exists("{}_start.molprobity".format(the_pdb.replace(".pdb", ""))):
        info("Starting molprobity not preexisting; generating it.")
        try:
            #os.system('phenix.molprobity ./output/phenix_updated.pdb {} {} >molprobity.log'.format(the_mtz, " ".join(the_cifs)))
            os.system('phenix.molprobity ./output/phenix_updated.pdb {} {} '.format(the_mtz, " ".join(the_cifs)))
            #os.system('mv molprobity.out {}_start.molprobity'.format(the_pdb.replace(".pdb", "")))
            shutil.move('molprobity.out', '{}_start.molprobity'.format(the_pdb.replace(".pdb", "")))
        except:
            warning("didn't find the phenix.molprobity output file, assuming bad rfree flags")
            #os.system('phenix.molprobity input.xray_data.r_free_flags.generate=True ./output/phenix_updated.pdb {} {} >molprobity.log'.format(the_mtz, " ".join(the_cifs)))
            os.system('phenix.molprobity input.xray_data.r_free_flags.generate=True ./output/phenix_updated.pdb {} {} '.format(the_mtz, " ".join(the_cifs)))
            shutil.move('molprobity.out', '{}_start.molprobity'.format(the_pdb.replace(".pdb", "")))


    ###
    # Phase 2 -- refine using phenix against the mtz.
    #if not os.path.exists('./output/phenix_rd1.pdb'):
    #    try:
    #        # We also do the initial phenix refinement once.
    #        os.system('phenix.refine {} {} {} --overwrite '.format('./output/phenix_updated.pdb', the_mtz, " ".join(the_cifs)))
    #        #os.rename('phenix_updated_refine_001.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_rd1.pdb')
    #        os.system('cp phenix_updated_refine_001.pdb ./output/phenix_rd1.pdb')
    #    except FileNotFoundError:
    #        warning("didn't find the phenix.refine output file, assuming bad rfree flags")
    #        os.system('phenix.refine {} {} {} refinement.input.xray_data.r_free_flags.generate=True --overwrite'.format('./output/phenix_updated.pdb', the_mtz, " ".join(the_cifs)))
    #        #os.rename('phenix_updated_refine_001.pdb'.format(the_pdb.replace('.pdb', '')), './output/phenix_rd1.pdb')
    #        os.system('cp phenix_updated_refine_001.pdb ./output/phenix_rd1.pdb')
   
    
    ###
    # Phase 2.05 -- cull ions and stuff from PDB. How should we do this?
    # Useful approach -- just get rid of most common stuff that will make
    # our lives the worst for the least gain.

    # AMW TODO: figure out a system where we aren't constantly overwriting this if we fix it in one instance.
    info("Culling waters and ions to produce a reduced updated file")
    os.system("cat ./output/phenix_updated.pdb | grep -v HOH | grep -v MG |grep -v SO4| grep -v 'K     K' | grep -v SR > ./output/phenix_updated_culled.pdb")
    info("Running pdb2fasta to generate a fasta.")
    warning("Update to new version that can't fail!!")
    os.system("if [[ ! -e ./output/{} ]] ; then\n\tpdb2fasta.py ./output/phenix_updated_culled.pdb > ./output/{}\nfi\n".format(the_pdb.replace('pdb', 'fasta'), the_pdb.replace('pdb', 'fasta')))
    
    ###
    # Phase 2.1 -- optionally, IF we are setting up jobs, copy this critical file to the runs dirs.
    if nstruct != 0:
        try:
            os.mkdir('output/runs/')
        except:
            pass
        try:
            os.mkdir('sbatch_files/')
        except:
            pass
    
        # write sbatch headers
        # not using execution_fn b/c we know it doesn't apply
        for i in range(nstruct):
            try:
                os.mkdir('output/runs/{}'.format(i))
            except:
                pass

            with open('sbatch_files/job{}.sbatch'.format(i), 'w') as g:
                g.write(sbatch_header.format(the_pdb=the_pdb, dirnum=i, pwd=os.getcwd()))
            
            shutil.copyfile('./output/phenix_updated_culled.pdb', './output/runs/{}/phenix_updated_culled.pdb'.format(i))
            shutil.copyfile('./output/{}'.format(the_pdb.replace('pdb', 'fasta')), './output/runs/{}/{}'.format(i, the_pdb.replace('pdb', 'fasta')))


    ###
    # Phase 3 -- run erraser2 on a culled portion of the phenix pdb.
    #  
    try:
        score_flags: str = ' -score:weights stepwise/rna/rna_res_level_energy4.wts -restore_talaris_behavior'
        weight_adj: str = ' -set_weights elec_dens_fast 10.0 cart_bonded 5.0 linear_chainbreak 10.0 chainbreak 10.0 fa_rep 1.5 fa_intra_rep 0.5 rna_torsion 10 suiteness_bonus 5 rna_sugar_close 10 atom_pair_constraint 10 '
        mutes: str = '-mute  core.scoring.electron_density.xray_scattering core.pack.rotamer_set.RotamerSet_ core.scoring.CartesianBondedEnergy -render_density true'
        mapfile: str = '{}_2mFo-DFc_map.ccp4'.format(the_pdb.replace('.pdb', ''))
        fasta: str = the_pdb.replace('pdb', 'fasta')
        #-sampler_num_pose_kept 10 for quick?
        info("About to run ERRASER2")
        #execution_fn('if [[ ! -e ./output/erraser_rd1.pdb ]]; then \n\tif [[ ! -e ./output/phenix_updated_culled.pdbFINISHED_1.pdb ]] ; then\n\t\t{} -s ./output/phenix_updated_culled.pdb {} {} {} -rmsd_screen 3.0 -edensity:mapfile {} -fasta ./output/{} -allow_virtual_side_chains false -sampler_num_pose_kept 100 -missing_density_to_jump true\n\tfi\n\n\tmv ./output/phenix_updated_culled.pdbFINISHED_1.pdb ./output/erraser_rd1_unmerged.pdb\nfi\n'.format(
        execution_fn('if [ ! -e ./output/erraser_rd1.pdb ]; then \n\tif [ ! -e ./output/phenix_updated_culled.pdbFINISHED_1.pdb ] ; then\n\t\t{} -s ./output/phenix_updated_culled.pdb {} {} {} -rmsd_screen 3.0 -edensity:mapfile {} -fasta ./output/{} -allow_virtual_side_chains false -sampler_num_pose_kept 100 -missing_density_to_jump true\n\tfi\n\n\tmv ./output/phenix_updated_culled.pdbFINISHED_1.pdb ./output/erraser_rd1_unmerged.pdb\nfi\n'.format(
            exe, score_flags, weight_adj, mutes, mapfile, fasta
        ))
    except:
        error("Disaster!")
        quit()
    
    ###
    # Phase 3.5 -- merge culled ions back.
    #  
    
    # Note that we first have to strip Hs
    execution_fn("cat ./output/erraser_rd1_unmerged.pdb | grep '^ATOM\|^HETATM' |grep -v 'H  $' > ./output/erraser_rd1.pdb")
    execution_fn("cat ./output/phenix_updated.pdb | grep \"HOH\|MG\|SO4\|K     K\|SR\" >> ./output/erraser_rd1.pdb")
    execution_fn("cat ./output/erraser_rd1_unmerged.pdb | grep -v '^ATOM\|^HETATM' >> ./output/erraser_rd1.pdb")


    ###
    # Phase 4 -- final phenix refinement
    #try:
    #    # We also do the final phenix refinement once.
    #    execution_fn('if [[ ! -e ./output/phenix_rd2.pdb ]] ; then \n\tcp {} ./output/\n\tcp *.cif ./output/\n\tpushd ./output/\n\tphenix.refine {} {} {} --overwrite\n\tmv erraser_rd1_refine_001.pdb phenix_rd2.pdb\n\tpopd\nfi'.format(the_mtz, 'erraser_rd1.pdb', the_mtz, " ".join(the_cifs)))
    #except:
    #    quit()

          
    ###
    # Phase 5 -- phenix.molprobity and analysis
    try:
        #execution_fn('cd ./output/ && phenix.molprobity phenix_rd2.pdb {} {} && mv molprobity.out {}_end.molprobity && cd ..'.format(the_mtz, " ".join(the_cifs), the_pdb.replace(".pdb", "")))
        execution_fn('cd ./output/ && phenix.molprobity erraser_rd1.pdb {} {} && mv molprobity.out {}_end.molprobity && cd ..'.format(the_mtz, " ".join(the_cifs), the_pdb.replace(".pdb", "")))
    except:
        #execution_fn('cd ./output/ && phenix.molprobity phenix_rd2.pdb {} {} && mv molprobity.out {}_end.molprobity && cd ..'.format(the_mtz, " ".join(the_cifs), the_pdb.replace(".pdb", "")))
        execution_fn('cd ./output/ && phenix.molprobity input.xray_data.r_free_flags.generate erraser_rd1.pdb {} {} && mv molprobity.out {}_end.molprobity && cd ..'.format(the_mtz, " ".join(the_cifs), the_pdb.replace(".pdb", "")))
        #quit()
    
    os.chdir('..')

def write_to_all_sbatch(text, nstruct):
    try:
        os.mkdir('sbatch_files')
    except:
        pass

    
    for i in range(nstruct):
        with open('sbatch_files/job{}.sbatch'.format(i), 'a') as f:
            f.write("{}\n".format(text.replace('output', 'output/runs/{}/'.format(i))))
    

def main(directories: List[str], run_test: bool, nstruct: int) -> None:
    for d in directories:
        info("Running ERRASER2 on contents of directory {}.".format(d))
        if run_test:
            erraser2(d, nstruct=0, execution_fn=os.system)
        else:
            try:
                shutil.rmtree('{}/sbatch_files'.format(d))
            except:
                pass
            erraser2(d, execution_fn=functools.partial(write_to_all_sbatch, nstruct=nstruct), nstruct=nstruct)

def test_args():
    pass

if __name__ == '__main__':

    basicConfig(level=INFO)

    parser = argparse.ArgumentParser(description='Run ERRASER2 on a working directory')
    parser.add_argument('directory', type=str, nargs='+',
                        help='dir(s) to run')
    parser.add_argument('--run_test', action="store_true",
                        help='just run a test or write a bunch of sbatch scripts?')
    parser.add_argument('--nstruct', type=int, nargs='?', default=0,
                        help='how many trajectories to set up? only meaningful if --run_test')
    args = parser.parse_args()
    main(args.directory, args.run_test, args.nstruct)
