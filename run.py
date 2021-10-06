import numpy as np
import numpy as np
import os
import shutil
import glob
import subprocess
import time
import sys
import getopt


def import_aux():
    with open('run.aux', 'r') as f:
        skip            = f.readline()
        output          = str(f.readline().split()[0])
        skip            = f.readline()
        run_string      = f.readline().split()[0]
        rescale         = f.readline().split()[0]
        skip            = f.readline()
        restart_from    = f.readline().split()[0]
        input_string    = f.readline().split()[0]
        skip            = f.readline()
        low_ring_size   = int(f.readline().split()[0])
        upper_ring_size = int(f.readline().split()[0])
        skip            = f.readline()
        MC_T            = int(f.readline().split()[0])
        seed_start      = int(f.readline().split()[0])
        seed_end        = int(f.readline().split()[0])
        skip            = f.readline()
        alpha_target    = float(f.readline().split()[0])
        alpha_lenience  = float(f.readline().split()[0])
        p6_target       = float(f.readline().split()[0])
        p6_lenience     = float(f.readline().split()[0])
        skip            = f.readline()
        restart         = str(f.readline().split()[0])
        skip            = f.readline()
        upper_A         = float(f.readline().split()[0])
        lower_A         = float(f.readline().split()[0])
        step_A          = float(f.readline().split()[0])
        skip            = f.readline()
        repulsion       = str(f.readline().split()[0])
        repulsion_except= int(f.readline().split()[0])
        repulsion_min   = float(f.readline().split()[0])
        repulsion_max   = float(f.readline().split()[0])
        repulsion_step  = float(f.readline().split()[0])
    options = {}

    # Ouput Prefix
    if not os.path.isdir(output):
        os.mkdir(output)
    rundir = output

    #Ring Size String eg. 567
    ring_size_string = ""
    i = low_ring_size

    while i <= upper_ring_size:
        ring_size_string=ring_size_string+"{:}".format(i)
        i+=1

    # Update rundir to include Ring size string
    if not os.path.isdir(rundir+'/'+ring_size_string+'/'):
        os.mkdir(rundir+'/'+ring_size_string+'/')
    rundir = rundir+'/'+ring_size_string+'/'

    options['output_prefix'] = rundir

    if run_string=='harmonic':
        options["run"] = 0
        print('Harmonic potential')

    elif run_string=='static':
        options["run"] = 1
        print('Static Ions')
        print('Beware : running this without harmonic potential will slow down calculations')

    elif run_string=='dynamic':
        options["run"] = 2
        print('Beware : running this without harmonic potential will slow down calculations')
    else:
        print('unrecognised run type')


    if input_string=='netmc':

        if restart_from=='no':
            # start with harmlin, using factors from the aux file
            options['input_location'] = ring_size_string+'/T{:}/t{:}_s'.format(MC_T, MC_T)
            options['input_sublocations'] = []

            for i in range(seed_start, seed_end+1):
                for j in range(repulsion_min, repulsion_max+repulsion_step, repulsion_step):
                    for k in range(1+int((upper_A-lower_A)/step_A)):
                        area = (lower_A+i*step_A)
                        options['file_sublocations'].append('{:}_same/harmlin/intercept_{:}/area_{:}'.format(i, int(j*100), int(1000*area)))
                if repulsion_except==1:
                    for k in range(1+int((upper_A-lower_A)/step_A)):
                        area = (lower_A+i*step_A)
                        options['file_sublocations'].append('{:}_same/harmlin/no_repulsion/area_{:}'.format(i, int(1000*area)))



        else:
            # use the string from aux file to describe the runtype its finding, and take files from there
            options['input_location'] = rundir+'/T{:}/t{:}_s'.format(MC_T, MC_T)
            options['input_sublocations'] = []
            for i in range(seed_start, seed_end+1):
                for files in os.listdir(options['input_location']+'{:}_same/{:}/'.format(i, restart_from)):
                    options['file_sublocations'].append('{:}_same/{:}/{:}'.format(i, restart_from, files))


        options['output_location'] = options['output_prefix'] + '/T{:}/t{:}_s'.format(MC_T, MC_T)
        if not os.path.isdir(options['output_prefix']+'/T{:}'.format(MC_T)):
            os.mkdir(options['output_prefix']+'/T{:}'.format(MC_T))
        for folders in options['file_sublocations']:
            if not os.path.isdir(options['output_location']+folders):
                os.mkdir(options['output_location']+folders)


            for i in range(1+int((upper_A-lower_A)/step_A)):
                area = (lower_A+i*step_A)
                if not os.path.isdir(options['output_location'] + folders+'/area_{:}'.format(int(1000*area))):
                    os.mkdir(options['output_location'] + folders+'/area_{:}'.format(int(1000*area)))
                shutil.copytree(options['input_location']+folders, options['output_location']+folders+'/area_{:}'.format(int(1000*area)))
                scale_crds(options['output_location'] + folders+'/area_{:}'.format(int(1000*area)))




    elif input_string=='alpha':

        options['input_location'] = ring_size_string+'/'
        options['input_sublocations'] = []
        with open('p6_alpha_details.out') as f:
            array = np.genfromtxt(f)
        p6_allowed      = np.argwhere((array[:,2]>=(p6_target-0.5*p6_lenience)) & (array[:,2]<=(p6_target+0.5*p6_lenience)))
        alpha_allowed   = np.argwhere((array[:,3]>=(alpha_target-0.5*alpha_lenience)) & (array[:,3]<=(alpha_target+0.5*alpha_lenience)))


        options['output_location'] = options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target, p6_target)
        if not os.path.isdir(options['output_prefix']+'/alpha_{:}'.format(alpha_target)):
            os.mkdir(options['output_prefix']+'/alpha_{:}'.format(alpha_target))
        if not os.path.isdir(options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target,p6_target)):
            os.mkdir(options['output_prefix']+'/alpha_{:}/p6_{:}/'.format(alpha_target,p6_target))

        for vals in p6_allowed:
            if vals in alpha_allowed:
                options['file_sublocations'].append('T{:}/t{:}_s{:}_same'.format(array[vals,0], array[vals,0], array[vals,1]))

                if not os.path.isdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/'.format(alpha_target, p6_target, array[vals,0], array[vals,1])):
                    os.mkdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/'.format(alpha_target, p6_target, array[vals,0], array[vals,1]))

                for i in range(1+int((upper_A-lower_A)/step_A)):
                    area = (lower_A+i*step_A)
                    if not os.path.isdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area))):
                        os.mkdir(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))


                    shutil.copytree(options['input_location'] + 'T{:}/t{:}_s{:}_same'.format(array[vals,0], array[vals,0], array[vals,1]), options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))
                    scale_crds(options['output_prefix'] + '/alpha_{:}/p6_{:}/t{:}_s{:}/{:}'.format(alpha_target, p6_target, array[vals,0], array[vals,1], int(1000*area)))



print('Number of arguments:', len(sys.argv), 'arguments.')
print('Argument List:', str(sys.argv))

def make_extract_restart(working_dir, area_value):
    with open(working_dir+'/area_{:}/extract_restart.sh'.format(area_value), 'w') as filename:
        filename.write('#!/bin/bash\n')
        filename.write('\n')

        filename.write('if [ -d "output" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output\n')
        filename.write('fi\n')

        filename.write('if [ -d "output_restart" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output_restart\n')
        filename.write('fi\n')

        filename.write('\n')
        filename.write('mkdir output_restart\n')
        filename.write('\n')
        filename.write('\n')
        filename.write('mv fort* output_restart/\n')
        filename.write('mv *.out* output_restart/\n')
        filename.write('mv header.hdr output_restart/\n')
        filename.write('#mv testout.rst output_restart/\n')

def make_extract(working_dir, area_value):
    with open(working_dir+'/area_{:}/extract.sh'.format(area_value), 'w') as filename:
        filename.write('#!/bin/bash\n')
        filename.write('\n')

        filename.write('if [ -d "output" ]\n')
        filename.write('then\n')
        filename.write('    rm -r output\n')
        filename.write('fi\n')

        filename.write('\n')
        filename.write('mkdir output_restart\n')
        filename.write('\n')
        filename.write('\n')
        filename.write('mv fort* output_restart/\n')
        filename.write('mv *.out* output_restart/\n')
        filename.write('mv header.hdr output_restart/\n')
        filename.write('mv testout.rst output_restart/\n')


def make_subscript(working_dir, area_value):
    with open(working_dir+'/area_{:}/sub_script'.format(area_value), 'w') as filename:
        filename.write('!/bin/csh \n')
        filename.write('#$ -j y\n')
        filename.write('#$ -o . \n')
        filename.write('#$ -l s_rt=1600:00:00\n')
        filename.write('#$ -l h_rt=1600:00:00\n')
#        filename.write('#$ -pe smp 16')
        filename.write('#$ -N area_{:}\n'.format(area_value))
        filename.write('')
        filename.write('set fromdir = {:}\n'.format(working_dir+'/area_{:}/'.format(area_value)))
        filename.write('set todir =   {:}\n'.format(working_dir+'/area_{:}/'.format(area_value)))
        filename.write('')
        filename.write('setenv WORK /tmp/$USER/code_cycle/job_data/$JOB_ID\n')
        filename.write('')
        filename.write('mkdir -p $WORK\n')
        filename.write('cd $fromdir\n')
        filename.write('cp -r * $WORK\n')
        filename.write('cd $WORK\n')
        filename.write('{:}\n'.format(working_dir+'/area_{:}/dip_10000.x'.format(area_value)))
        filename.write('cp * $todir\n')
        filename.write('rm -Rf $WORK\n')

def scale_crds(dir):
    shutil.move(dir+'/crys.crds', dir+'/unscaled_crys.crds')
    with open(dir+'/unscaled_cryrs.crds', 'r') as filein:
        crds = np.genfromtxt(filein, skip_header=5, skip_footer=6)
        dim = np.genfromtxt(filein, skip_header=8+crds.shape[0])
    crds = np.multiply(crds, 1.0/0.52917721090380)
    dim = np.multiply(dim, 1.0/0.52917721090380)
    with open(dir+'/crys.crds', 'w') as fileout:
        fileout.write('T\nF\nF\nF\nF\n')
        for i in range(crds.shape[0]):
            fileout.write('{:<26}{:<26}{:<26}\n'.format(crds[i,0], crds[i,1], crds[i,2]))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(1.0,0.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,1.0,0.0))
        fileout.write('{:<26}{:<26}{:<26}\n'.format(0.0,0.0,1.0))
        fileout.write('{:<26}'.format(dim[0]))
        fileout.write('{:<26}'.format(dim[1]))
        fileout.write('{:<26}'.format(dim[2]))
    return


def main(options):
    if (options["restart"]==False):
        folders = glob.glob('area_*', recursive=True)
        for folder in folders:
            try:
                shutil.rmtree(folder)
            except OSError as e:
                print("Error : %s" % (e.strerror))


    working_dir = os.getcwd()


#    with open('test_opt_pe.dat', 'r') as f:
#        array = np.genfromtxt(f)

    with open("crys.crds", "r") as f:
        array = np.genfromtxt(f, skip_header=5, skip_footer=6)



    for i in range(array.shape[0]):
        area_i = array[i,0]
        area_label = str(int(area_i*1000))
        print(area_label)

        pbx = array[i,1]
        pby = array[i,2]
        pbz = array[i,3]
        print(pbx,'  ', pby, '  ', pbz)

        factor = 1.0/0.52917721090380
        print(pbx*factor,'  ', pby*factor, '  ', pbz*factor)
        if (options["restart"]==True):
            os.chdir("area_{}".format(area_label))
            make_extract_restart(working_dir, area_label)
            subprocess.call(['./extract_restart.sh'])

            if (options["server"]==False):
                subprocess.call(['nq', './dip_11000.x'])
            else:
                make_subscript(working_dir, area_label)
                subprocess.call(['qsub', 'sub_script'])
            os.chdir('../')

        if (options["restart"]==False):
            os.mkdir("area_{}".format(area_label))

            shutil.copy('test_area_{}_sample_0_crys.crds'.format(i), 'area_{}/'.format(area_label))
            with open('test_area_{}_sample_0_crys.crds'.format(i), 'r') as f:
                crds = np.genfromtxt(f)

            print(max(crds[:,0]))
            print(max(crds[:,1]))
            print('\n')
            print(pbx*factor)
            print(pby*factor)

            with open("area_{}/crys.crds".format(area_label), 'w') as f:
                f.write('T\nF\nF\nF\nF\n')
                for i in range(crds.shape[0]):
                    f.write('{:<20}{:<20}{:<20}\n'.format(crds[i,0], crds[i,1], crds[i,2]))

            for files in ['cf.inpt','a.x','a.f','dip_10000.x','crystal_cell.inpt','dip_11000.x','tube_19_0_15.crds','ts_sio2.inpt','runtime.inpt','q','lj.inpt','lightcf_licl.inpt','filenames.inpt','extract.sh','dynmat.inpt','test_opt_pe.dat','split_to_areas.py']:
                shutil.copy(files, 'area_{}/'.format(area_label))

            with open('area_{}/crystal_cell.inpt'.format(area_label), 'w') as g:


                g.write('   1.0d0  0.0d0  0.0d0\n')
                g.write('   0.0d0  1.0d0  0.0d0     Cell matrix.\n')
                g.write('   0.0d0  0.0d0  1.0d0\n')
                g.write('{:<25} L_x\n'.format(pbx*factor))
                g.write('{:<25} L_y\n'.format(pby*factor))
                g.write('{:<25} L_z\n'.format(200.0))
                g.write('           20       Number of unit cells in x direction (starting from a crystal only)\n')
                g.write('           12       Number of unit cells in y direction (starting from a crystal only)\n')
                g.write('           1        Number of unit cells in z direction (starting from a crystal only)\n')
                g.write('           2        Number of type I particles in unit cell.\n')
                g.write('x_tet.mat\n')
                g.write('           2        Number of type II particles in unit cell.\n')
                g.write('m_tet.mat\n')
                g.write('   4.6d0            Unit cell length - x direction.\n')
                g.write('   8.0d0            Unit cell length - x direction.\n')
                g.write('   100.0d0          Unit cell length - x direction.\n')
            with open('area_{:}/runtime.inpt'.format(area_label), 'w') as filename:
                filename.write('{:}         Number of steps in the run.\n'.format(options["steps"]))
                filename.write('0.0001d0    Translational temperature.\n')
                filename.write('151.315     Relative molecular mass (gmol^-1)\n')
                filename.write('800         Number of ionic molecular units.\n')
                filename.write('2           Number of ionic species.\n')
                filename.write('1600        Number of species of type 1.\n')
                filename.write('-1.38257d0  Permanent charge on type 1 ions.\n')
                filename.write('16.0d0      Atomic mass of type 1.\n')
                filename.write('.true.      Polarizable ?\n')
                filename.write('800         Number of species of type 2.\n')
                filename.write('2.76514d0   Permanent charge on type 2 ions.\n')
                filename.write('28.08d0     Atomic mass of type 2.\n')
                filename.write('.false.     Polarizable ?\n')
                filename.write('25.0d0      Timestep (a.u.).\n')
                filename.write('rim         Type of run (rim,dippim,quadpim)\n')
                filename.write('ts          Type of run (epp,cim,aim)\n')
                filename.write('.false.     Steepest descent minimisation ? (if .true. require time step + annealing steps + T/F for anneal radii only).\n')
                filename.write('.false.     Conjugate gradient minimisation?\n')
                filename.write('.true.      Restart ? (if .true. then next line is \'inpt\' or \'rst\' if cell lengths are to be read from the \'crystal_cell\' or restart files).\n')
                filename.write('inpt\n')
                filename.write('.false.     Set up velocities?\n')
                filename.write('.false.     Velocity rescale prior to main run ?\n')
                filename.write('.false.     Random displacement of ions.\n')
                filename.write('.true.      Move ions?\n')
                filename.write('.false.     Do a dynamical matrix calculation?\n')
                filename.write('.true.      Relax input structure? (.true. will \'kill\' the velocities at a local k.e. maximum - steepest descent minimisation).\n')
                filename.write('1           Number of steps inbetween periodic output (energies).\n')
                filename.write('100         Number of steps inbetween periodic output (velocities etc).\n')
                filename.write('1           Number of steps inbetween periodic output (frictions etc).\n')
                filename.write('1           Number of steps inbetween periodic output (pressure etc).\n')
                filename.write('10          Number of steps inbetween rdf call in main loop.\n')
                filename.write('10          Number of steps inbetween S(k) call in main loop\n')
                filename.write('1           Number of ions to monitor.\n')
                filename.write('1           Ion number to monitor. (1)\n')
                filename.write('5.60d0      eta = <x>/boxlen.\n')
                filename.write('45.0d0      rcut (au).\n')
                filename.write('1.0d-7      convergence parameter.\n')
                filename.write('.false.     Nose-Hoover thermostat? (if true then enter a relaxation time)\n')
                filename.write('.false.     Periodic rescale of temperature?\n')
                filename.write('.false.     Periodic re-annealing?\n')
                filename.write('.false.     Isotropic barostat?\n')
                filename.write('.false.     Anisotropic barostat?\n')
                filename.write('.false.     Orthorhombic cell?\n')
                filename.write('.false.     T ramping ?\n')
                filename.write('.false.     p ramping?\n')
                filename.write('.false.     EFG calc. ?\n')
                filename.write('.false.     Box ramping?\n')
                filename.write('.false.     Box ramping?\n')





            os.chdir('area_{}'.format(area_label))
            subprocess.call(['nq', './dip_11000.x'])
            os.chdir('../')
            #    os.chdir('area_{}'.format(area_label))
            #    subprocess.call(['./dip_11000.x'])
            #    time.sleep(1000)
            #    os.chdir('../')
options = get_options()
main(options)



