#!/usr/bin/env python3
import os
import glob
import re
import math
import subprocess
from shutil import copyfile, copy2

def thingToList(thing):
    l = []
    for i in range(len(thing)):
        l.append(thing[i])
    return l

def l_vec_longer(level_vector):
    yIsX = os.environ.get('ADAPTATION_PARAM_Y_EQUALS_X',
                          'False').lower() in ['true', '1']
    if yIsX and len(level_vector) < 5:
        level_vector = [level_vector[0]] + thingToList(level_vector)
    assert (len(level_vector) < 5)
    return level_vector

def l_vec_to_string(l):
    return "prob_" + str("_".join([str(l_i) for l_i in l]))

def get_prob_path(prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'), level_vector=None):
    level_vector = l_vec_longer(level_vector)
    return os.path.join(prob_prepath, l_vec_to_string(level_vector))

def check_folder_exists(prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'), level_vector=None):
    level_vector = l_vec_longer(level_vector)
    return os.path.isfile(os.path.join(get_prob_path(prob_prepath, level_vector), 'parameters'))

def check_finished(prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'), level_vector=None):
    level_vector = l_vec_longer(level_vector)
    # if the simulation has run, but not long enough
    if os.path.isfile(os.path.join(get_prob_path(prob_prepath, level_vector), 'GENE.finished')):
        return True
    return False

def restart_sim(prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'), gene_path=os.environ.get('ADAPTATION_GENE3D_PATH'), level_vector=None, newSimTime=None):
    assert (level_vector)
    level_vector = l_vec_longer(level_vector)
    assert(os.environ.get('ADAPTATION_SUBMIT_COMMAND') is not None)
    assert (check_finished(prob_prepath, level_vector))
    # if the simulation has run, but not long enough
    prob_dir=get_prob_path(prob_prepath, level_vector)
    # rename output files
    assign_index = int(1)
    while os.path.isfile(os.path.join(get_prob_path(prob_prepath, level_vector),"out/nrg_"+str(assign_index))):
        assign_index += 1
    subprocess.run(gene_path+"/tools/runassign " + str(assign_index),
                   shell=True, cwd=prob_dir+"/out")
    # delete GENE.finished, adapt parameters
    with open(prob_dir + "/parameters", 'r') as pfile:
        oldpdata = pfile.read()
    newpdata = oldpdata.replace(
        "read_checkpoint = F", "read_checkpoint = T")
    newpdata = newpdata.replace(
        "read_checkpoint  = F", "read_checkpoint  = T")
    if newSimTime:
        newpdata, numSubstitutions = re.subn(r'simtimelim.*$',r'simtimelim = '+ str(newSimTime),newpdata, flags=re.MULTILINE)
        assert(numSubstitutions > 0)

    with open(prob_dir + '/parameters', 'w') as pfile:
        pfile.write(newpdata)
    qname = "submit.sh"
    os.remove(prob_dir + '/GENE.finished')
    print("resubmitting for " + str(prob_dir))
    # submit
    subprocess.run(os.environ.get('ADAPTATION_SUBMIT_COMMAND')+ " " + qname, shell=True, cwd=prob_dir)


def dispatch_to_run(prob_prepath=os.environ.get('ADAPTATION_PROB_PATH'), level_vector=None, gene_directory=os.environ.get('ADAPTATION_GENE3D_PATH')):
    level_vector = l_vec_longer(level_vector)
    assert not (check_folder_exists(prob_prepath, level_vector))
    assert(os.environ.get('ADAPTATION_SUBMIT_COMMAND') is not None)
    with open(os.environ.get('ADAPTATION_PARAMETER_FILE_TEMPLATE'), 'r') as pfile:
        pdata = pfile.read()
    with open(os.environ.get('ADAPTATION_SUBMIT_SCRIPT_TEMPLATE'), 'r') as qfile:
        qdata = qfile.read()

    # create the folders
    # make prob dir
    prob_dir=get_prob_path(prob_prepath, level_vector)
    os.mkdir(prob_dir)
    os.mkdir(prob_dir+"/out")
    print("generated " + prob_dir)
    subprocess.run("ln -s " + gene_directory + "/bin/gene3d",
                   shell=True, cwd=prob_dir)
    subprocess.run("ln -s " + gene_directory + "/gvec",
                   shell=True, cwd=prob_dir)

    if level_vector[2] > 5:
        p = [0, 0, 5, 0, 1, 1]  # todo adapt
    else:
        p = [1, 0, 4, 0, 1, 1]
    l_par = sum(p)
    l = level_vector
    nprocspernode = 48

    # adapt the "standard" parameter file
    # generate parameters files from template
    newpdata = pdata.replace("$n_procs_x", str(2**p[0]))
    newpdata = newpdata.replace("$n_procs_y", str(2**p[1]))
    newpdata = newpdata.replace("$n_procs_z", str(2**p[2]))
    newpdata = newpdata.replace("$n_procs_v", str(2**p[3]))
    newpdata = newpdata.replace("$n_procs_w", str(2**p[4]))
    newpdata = newpdata.replace("$n_procs_s", str(2**p[5]))

    newpdata = newpdata.replace("$nx0", str(2**l[0]))
    newpdata = newpdata.replace("$ny0", str(2**l[1]))
    newpdata = newpdata.replace("$nz0", str(2**l[2]))
    newpdata = newpdata.replace("$nv0", str(2**l[3]))
    newpdata = newpdata.replace("$nw0", str(2**l[4]))

    newpdata = newpdata.replace("$timelim", str(71500))
    newpdata = newpdata.replace("$simtimelim", str(1000))

    # Write the file out again, link to prob folder
    with open(prob_dir + '/parameters', 'w') as file:
        file.write(newpdata)

    # adapt the submit script
    probname = "G-3D-" + l_vec_to_string(level_vector)[5:]
    newqdata = qdata.replace("$probname", probname)
    newqdata =  newqdata.replace("$nnodes", str(
        int(math.ceil(2**int(l_par)/nprocspernode))))
    newqdata = newqdata.replace("$nprocesses", str(int(2**l_par)))
    newqdata = newqdata.replace("$p_level", str(int(l_par)))
    newqdata = newqdata.replace("$nprocspernode", str(nprocspernode))
    #newqdata = newqdata.replace("$samples_suffix", samples_suffix)
    # Write the file out again
    qname = "submit.sh"
    with open(os.path.join(prob_dir, qname), 'w') as file:
        file.write(newqdata)
    # Geometry file not needed any more
    #geometryfile = template_directory + "/iterdb_w7x_eim00"  # "/W7X_*dat"
    #geometryfile = glob.glob(geometryfile)[0]
    ## print(geometryfile)
    #copy2(geometryfile, prob_dir)

    # submit
    subprocess.run(os.environ.get('ADAPTATION_SUBMIT_COMMAND')+ " " + qname, shell=True, cwd=prob_dir)


# if called directly, run tests
if __name__ == "__main__":
    #prob_prepath="/hppfs/scratch/02/di39qun2/gene3d-flw-simulations/"
    prob_prepath=os.environ.get('ADAPTATION_PROB_PATH')
    #gene_path="/hppfs/work/pn34mi/di68xux2/myGene3d/"
    gene_path=os.environ.get('ADAPTATION_GENE3D_PATH')
    #lvs = [[6, 5, 5, 5, 3],[5, 6, 5, 5, 3],[5, 5, 5, 6, 3],[5, 5, 5, 5, 4]]
    lvs = [[7,5,5,5,3], [6,6,5,5,3]]
    lvs = [[7,7,7,7,7]]
    for lv in lvs:
    #lv = [5, 5, 5, 5, 3]
        print (get_prob_path(level_vector=lv))
        dispatch_to_run(prob_prepath, lv)
        #restart_sim(prob_prepath, gene_path, lv)
        print(check_folder_exists(prob_prepath, lv), check_finished(prob_prepath, lv))
