#!/home/jgarcia/anaconda3/bin/python

import qetools as qt
import numpy as np
import matplotlib.pyplot as plt
import os.path
import subprocess, os

import sys

num_proc =sys.argv[1];
inputf =sys.argv[2];
print("Using num_proc",num_proc," to converge ",inputf);

my_env = os.environ.copy()
my_env["OMP_NUM_THREADS"] = "1";

QEnl    = qt.read_QEnamelist(inputf);
atm_spec= qt.read_atomic_species(inputf);
atm_pos = qt.read_atomic_positions(inputf);
kpoints = qt.read_kpoints(inputf);

#CONTROL PARAMETERS
etol=1000*float(QEnl["&CONTROL"]["etot_conv_thr"]);
ftol=float(QEnl["&CONTROL"]["forc_conv_thr"]);


print("Evaluating the WFC convergence");
wfc_conv= float(QEnl["&SYSTEM"]["ecutwfc"]);
wfc_min = wfc_conv;
wfc_max = 300;
dwfc    = 10;
wfcoffs = np.arange(wfc_min, wfc_max, dwfc );

energies=[];
magns  = [];
forces = [];
filename = inputf.replace(".in","");
for wfc in wfcoffs:
    QEnl["&SYSTEM"]["ecutwfc"] = wfc; 
    QEnl["&SYSTEM"]["ecutrho"] = 10*wfc; 
    ifilename = filename+".wfc"+str(wfc)+".in";
    ofilename = filename+".wfc"+str(wfc)+".out"
    qt.write_QEnamelist(QEnl,ifilename);
    qt.write_atomic_specties(atm_spec,ifilename);
    qt.write_atomic_positions(atm_pos, ifilename);
    qt.write_kpoints(kpoints, ifilename);

    try:
        energy        = qt.read_finalenergy(ofilename);
        magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
        force         = qt.read_finalocurrence("Total force",ofilename);
        if energy==None or force == None or magnetization==None:
            raise cls()   
    except:
        ofile = open(ofilename, "w")
        subprocess.run(["mpirun", "-np", str(num_proc), "pw.x", "-inp", ifilename], stdout=ofile, env=my_env);
        ofile.close();
        energy        = qt.read_finalenergy(ofilename);
        magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
        force         = qt.read_finalocurrence("Total force",ofilename);
        
    print("wfc",wfc,"energy=",energy," magnetization=",magnetization," force=",force);
    energies.append(energy); magns.append(magnetization); forces.append(force);
    
    DEN = qt.final_diff(energies);
    DFN = qt.final_diff(forces  ); 
    if DEN != 0 and (DEN < etol) and (DFN < ftol) :
        print("converged wf=",wfc," led to a threshold smaller than ", DEN,"<",etol);
        print("converged wf=",wfc," led to a threshold smaller than ", DFN,"<",ftol);
        wfc_conv = wfc;
        break;
    
plt.ylabel("rel. chg. Enr ");plt.xlabel("wfc");
qt.plot_conv(wfcoffs,energies, "wfcconv.energy_"+filename+".pdf"); plt.show();
plt.ylabel("rel. chg. Forc ");plt.xlabel("wfc");
qt.plot_conv(wfcoffs,forces, "wfcconv.force_"+filename+".pdf"); plt.show();

print("Evaluating the KPOINTS convergence");

kshift= (0,1);
conv_kps = [];
for ks in kshift:
    kpoints = qt.read_kpoints(inputf);
    kp_conv= np.max(kpoints[:2]);
    kpmin  = kp_conv; kpmax  = 100; dk=3;
    kps    = np.arange(kpmin,kpmax, dk );
    kps_pershift = [];

    
    energies=[];
    magns  = [];
    forces = [];
    filename = inputf.replace(".in","")+".kshift"+str(ks);
    nmax   = 0 ;
    for n,kp in enumerate(kps):
        QEnl["&SYSTEM"]["ecutwfc"] = wfc_conv;
        QEnl["&SYSTEM"]["ecutrho"] = 10*wfc_conv; 
        ifilename = filename+".kp"+str(kp)+".in";
        ofilename = filename+".kp"+str(kp)+".out";
        qt.write_QEnamelist(QEnl,ifilename);
        qt.write_atomic_specties(atm_spec,ifilename);
        qt.write_atomic_positions(atm_pos, ifilename);
        #Change kpoint and shift
        kpoints= qt.mod_kpoints(kpoints,n*dk,ks);
        qt.write_kpoints(kpoints, ifilename);

        try:
            energy        = qt.read_finalenergy(ofilename);
            magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
            force         = qt.read_finalocurrence("Total force",ofilename);
            if energy==None or force == None or magnetization==None:
                raise cls()   
        except:
            ofile = open(ofilename, "w")
            subprocess.run(["mpirun", "-np", str(num_proc), "pw.x", "-inp", ifilename], stdout=ofile, env=my_env);
            ofile.close();
            energy        = qt.read_finalenergy(ofilename);
            magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
            force         = qt.read_finalocurrence("Total force",ofilename);


        print("kp",kp,"energy=",energy," magnetization=",magnetization," force=",force);
        energies.append(energy); magns.append(magnetization); forces.append(force);

        DEN = qt.final_diff(energies);
        DFN = qt.final_diff(forces  ); 
        if DEN != 0 and (DEN < etol) and (DFN < ftol) :
            conv_kps.append([n,ks,kp]);
            print("converged kp=",kp," led to a threshold smaller than ", DEN,"<",etol);
            print("converged kp=",kp," led to a threshold smaller than ", DFN,"<",ftol);
            break;
 
    plt.ylabel("rel. chg. Enr ");plt.xlabel("nkp");
    qt.plot_conv(kps,energies, "kpconv.energy_"+filename+".pdf");plt.show();

    plt.ylabel("rel. chg. Forc ");plt.xlabel("nkp");
    qt.plot_conv(kps,forces, "kpconv.force_"+filename+".pdf");   plt.show();
    

print(" Evaluating the lattice constant convergence" );
lata_conv= 0;
lata  =  float(QEnl['&SYSTEM']['celldm(1)']);
latas = np.linspace(lata*0.7,lata*1.3, 20);
kpoints = qt.best_kpoints( kpoints, conv_kps );
energies=[];
magns  = [];
forces = [];
filename = inputf.replace(".in","");
for lata in latas:
    QEnl["&SYSTEM"]["ecutwfc"] = wfc_conv;
    QEnl["&SYSTEM"]["ecutrho"] = 10*wfc_conv; 
    qt.best_kpoints( kpoints, conv_kps );

    QEnl['&SYSTEM']['celldm(1)'] = lata;
    ifilename = filename+".lata"+str(lata)+".in";
    ofilename = filename+".lata"+str(lata)+".out";
    qt.write_QEnamelist(QEnl,ifilename);
    qt.write_atomic_specties(atm_spec,ifilename);
    qt.write_atomic_positions(atm_pos, ifilename);
    qt.write_kpoints(kpoints, ifilename);

    try:
        energy        = qt.read_finalenergy(ofilename);
        magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
        force         = qt.read_finalocurrence("Total force",ofilename);
        if energy==None or force == None or magnetization==None:
            raise cls()   
    except:
        ofile = open(ofilename, "w")
        subprocess.run(["mpirun", "-np", str(num_proc), "pw.x", "-inp", ifilename], stdout=ofile, env=my_env);
        ofile.close();
        energy        = qt.read_finalenergy(ofilename);
        magnetization = qt.read_finalocurrence("absolute magnetization",ofilename);
        force         = qt.read_finalocurrence("Total force",ofilename);
 
    print("lata",lata,"energy=",energy," magnetization=",magnetization," force=",force);
    energies.append(energy); magns.append(magnetization); forces.append(force);

latamin = latas[np.argmin(energies)];
print("lat_const_minimum",latamin);
strain = latas/latamin-1;
plt.plot(strain, energies,'-o'); plt.ylabel("Energy(eV)");plt.xlabel("a (bohr)");
plt.savefig("envslat."+inputf+".pdf"); plt.show();
plt.plot(strain, magns,'-o');  plt.ylabel("Magentization ");plt.xlabel("a (bohr)");
plt.savefig("magslat."+inputf+".pdf"); plt.show();


ifilename = filename+".conv.in";
print("Printing convergeed parameters in ", ifilename)
QEnl["&SYSTEM"]["ecutwfc"] = wfc_conv;
QEnl["&SYSTEM"]["ecutrho"] = 10*wfc_conv; 
qt.best_kpoints( kpoints, conv_kps );
QEnl['&SYSTEM']['celldm(1)'] = latamin;
qt.write_QEnamelist(QEnl,ifilename);
qt.write_atomic_specties(atm_spec,ifilename);
qt.write_atomic_positions(atm_pos, ifilename);
qt.write_kpoints(kpoints, ifilename);

print("FINISHED")

