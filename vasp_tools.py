import os
import copy

"""
Functions for setting up inputs for VASP (and eventually other DFT codes.)
"""

# TODO: also make a function that, given an existing calc, tries to move atoms between cells to match that calc closely.
# TODO: Make it easy to change constraints, and to view selective dynamics. Not sure if ASE will properly print partially fixed atoms, so I might need to use pymatgen.


def make_job_script(
    computer_name="Cypress",
    executable="vasp_std",
    system_name="job",
    job_script_string=None,
    job_minutes=60,
    num_job_nodes=1,
    job_cores=None,
    job_partition="debug",
    job_qos="normal",
    job_queue=None,
    job_out_file="job.out",
    job_error_file="job.err",
    job_script_name="jobScript.sh",
    type_of_run="Relaxation",
    magnetic_continuation=False,
    job_memory=None,
    job_mail_condition="end",
    job_mail_address="",
    job_architecture=None,
    job_allocation=None,
    continuation=False,
):
    # TODO:Add an input that's just a string that gets added (probably add it below keyword block and above job run block)
    # clear the file:
    with open(job_script_name, "w+") as file:
        file.write("")
    jobManager = {"Cori": "SLURM", "Cypress": "SLURM", "QueenBee": "PBS"}
    moduleLines = {
        "Cori": "module load vasp/5.4.1_vtst",
        "Cypress": "module load intel-psxe",
        "QueenBee": """module purge 
module load intel/19.0.5 intel-mpi
export TASKS_PER_HOST=20 # number of MPI tasks per host
export THREADS_HOST=1    # number of OpenMP threads spawned by each task on the host. 1 for pure MPI.
export TASKS_PER_MIC=30  # number of MPI tasks per MIC. Uses coprocessor. Can vary.
export THREADS_MIC=1     # number of OpenMP threads spawned by each task on the MIC. 1 for pure MPI.
cd $PBS_O_WORKDIR        # go to where your PBS job is submitted if necessary


""",
    }
    executionCommand = {"Cori": "srun", "Cypress": "mpirun", "QueenBee": "mpirun"}
    defaultCores = {
        "Cori": 32,
        "Cypress": 20,
        "QueenBee": 20,
    }  # Use this as the number of cores (per node) unless job_cores is set.
    if job_cores != None:
        coresWrite = job_cores
    else:
        coresWrite = defaultCores[computer_name]
    if job_script_string != None:
        jobScriptOut = job_script_string
        with open(job_script_name, "w+") as file:
            file.write(job_script_string)
    else:
        if (
            jobManager[computer_name] == "SLURM"
        ):  # TODO: put general slurm stuff here, and computer-specific stuff elsewhere (maybe in a dictionary)
            with open(job_script_name, "a") as file:
                file.write("#!/bin/bash")
                file.write("\n#SBATCH --job-name=" + system_name)
                file.write("\n#SBATCH-N " + str(num_job_nodes))
                file.write("\n#SBATCH --ntasks-per-node=" + str(coresWrite))
                file.write("\n#SBATCH --cpus-per-task=1")
                file.write("\n#SBATCH-t " + str(job_minutes))
                file.write("\n#SBATCH-o " + job_out_file)
                file.write("\n#SBATCH-e " + job_error_file)
                file.write("\n#SBATCH --mail-type=" + job_mail_condition)
                file.write("\n#SBATCH --mail-user=" + job_mail_address)
                if job_partition != "":
                    file.write("\n#SBATCH-p " + job_partition)
                if job_qos != None:
                    file.write("\n#SBATCH --qos=" + job_qos)
                if job_architecture != None:
                    file.write("\n#SBATCH-C " + job_architecture)
                if job_memory != None:
                    file.write("\n#SBATCH --mem=" + str(job_memory))
                file.write("\n\n")
        elif (
            jobManager[computer_name] == "PBS"
        ):  # TODO: put general slurm stuff here, and computer-specific stuff elsewhere
            # In preparation for getting time in correct format:
            s = job_minutes * 60
            hours, remainder = divmod(s, 3600)
            minutes, seconds = divmod(remainder, 60)
            with open(job_script_name, "a") as file:
                file.write("#!/bin/bash")
                file.write("\n#PBS -N " + system_name)
                file.write("\n#PBS -V")
                file.write(
                    "\n#PBS -l nodes="
                    + str(num_job_nodes)
                    + ":"
                    + "ppn="
                    + str(coresWrite)
                )
                file.write(
                    "\n#PBS -l walltime="
                    + "{:02}:{:02}:{:02}".format(int(hours), int(minutes), int(seconds))
                )
                file.write("\n#PBS -o " + job_out_file)
                file.write("\n#PBS -e " + job_error_file)
                if job_queue != None:
                    file.write("\n#PBS -q " + job_queue)
                if job_allocation != None:
                    file.write("\n#PBS -A " + job_allocation)
                file.write("\n\n")
        runString = (
            "\n"
            + executionCommand[computer_name]
            + " -n "
            + str(coresWrite * num_job_nodes)
            + " "
            + executable
            + "\n\n"
        )
        with open(job_script_name, "a") as file:
            with open(job_script_name, "a") as file:
                file.write("\n" + moduleLines[computer_name])
                file.write(runString)
        ### More exotic runs:
        if magnetic_continuation:
            with open(job_script_name, "a") as file:
                file.write("\ncontPrep.sh #WILLDELETE")
                file.write("\nmv INCAR_mag INCAR #WILLDELETE")
                file.write(runString + " #WILLDELETE")
                file.write("\nsed -i '/#WILLDELETE/d' ./jobScript.sh #WILLDELETE")
        if type_of_run.lower() == "lattice_constant":
            with open(job_script_name, "a") as file:
                # file.write(runString)
                file.write("\ncontPrep.sh")
                file.write(runString + "\n")
        if type_of_run.lower() == "dos":
            with open(job_script_name, "a") as file:
                # file.write(runString)
                file.write("\nmv INCAR INCAR_Charge \nmv KPOINTS KPOINTS_Charge")
                file.write("\nmv INCAR_DOS INCAR \nmv KPOINTS_DOS KPOINTS")
                file.write(runString + "\n")
        #        else:
        #            with open(job_script_name,'a') as file:
        #                file.write(runString+"\n")

        if continuation == True:
            with open(job_script_name, "r") as file:
                data = file.readlines()
                c = 0
                for a in data:
                    c = c + 1
                    if "mpirun" in a:
                        reduce_time = str(job_minutes - 3)
                        data[c - 1] = "timeout " + reduce_time + "m " + data[c - 1]
            with open(job_script_name, "w") as file:
                file.writelines(data)
                file.write("continuation.sh")


def make_vasp_input(
    input_atoms,
    calculation_directory,
    meta_parameters,
    kpoint_parameters,
    potcar_parameters,
    job_parameters,
    dimer_parameters,
):
    try:
        os.chdir(calculation_directory)
    except FileNotFoundError:
        os.mkdir(calculation_directory)
        os.chdir(calculation_directory)

    os.chdir(calculation_directory)
    make_incar(**meta_parameters)
    geometry_guess_quality = meta_parameters["geometry_guess_quality"]
    meta_parameters.setdefault("magnetic_continuation", False)
    if geometry_guess_quality.lower() == "bad":
        # print("Making faster calc in subdirectory.")
        os.chdir("FastRelax")
        input_atoms.write("POSCAR")
        os.chdir("..")
    else:
        input_atoms.write("POSCAR")
    make_kpoints(
        geometry_guess_quality=geometry_guess_quality,
        **kpoint_parameters,
        cell=input_atoms.cell
    )
    makePOTCAR(
        input_atoms, potcar_parameters, geometry_guess_quality=geometry_guess_quality
    )
    make_job_script(**job_parameters)
    type_of_run = meta_parameters["type_of_run"]
    # For DOS run, need two version of some files to run charge calc and then DOS calc
    if type_of_run.lower() == "dos":
        kpoint_parameters_dos = copy.deepcopy(kpoint_parameters)
        parameters_overwrite = {
            "NEDOS": 5000,
            "NPAR": 1,
            "LORBIT": 10,
            "LVTOT": True,
            "NBANDS": 9 * len(input_atoms.get_atomic_numbers()),
            "SIGMA": 0.1,
            "EDIFF": 10**-8,
            "ICHARG": 11,
        }  # TODO: Make this so it depends on the parameter_set. Will probably require some refactoring.
        kpoint_parameters_dos["factor"] = [2.75, 2.75, 1]
        kpoint_parameters_dos["kpoints_name"] = "KPOINTS_DOS"
        make_kpoints(**kpoint_parameters_dos)
        meta_parametersDOS = copy.deepcopy(meta_parameters)
        meta_parametersDOS["overwrite_parameters"].update(parameters_overwrite)
        meta_parametersDOS["incar_name"] = "INCAR_DOS"
        make_incar(**meta_parametersDOS)
    elif type_of_run.lower() == "dimer":
        makeMODECAR(dimer_parameters["bond_break"], input_atoms)
    if meta_parameters[
        "magnetic_continuation"
    ]:  # Could apply to multiple types of runs (relaxation, dimer, MD)
        os.rename("INCAR", "INCAR_mag")
        parameters_overwrite = {"ISPIN": 1, "NSW": 0, "ICHARG": 2}
        meta_parameters_magnetic_continuation = copy.deepcopy(meta_parameters)
        meta_parameters_magnetic_continuation["overwrite_parameters"].update(
            parameters_overwrite
        )
        make_incar(**meta_parameters_magnetic_continuation)


def makeMODECAR(bond_break, input_atoms, file_name="MODECAR"):
    import numpy as np

    modecar = np.zeros((len(input_atoms), 3))
    modecar[bond_break[0]] = (
        input_atoms.get_positions()[bond_break[0]]
        - input_atoms.get_positions()[bond_break[1]]
    )
    modecar[bond_break[1]] = -(
        input_atoms.get_positions()[bond_break[0]]
        - input_atoms.get_positions()[bond_break[1]]
    )
    modecar = modecar / np.linalg.norm(modecar)
    np.savetxt(file_name, modecar)


def makePOTCAR(input_atoms, potcar_parameters, geometry_guess_quality="Good"):
    # TODO: Set default PP types for each element for POTCAR.
    # TODO: make it so the code still runs if sym_potcar_map isn't given
    # If I ever go back to using ASE, will probably want to use something like this: potcar_parameters = {"setups": "recommended", "xc":"PW91", "pp": "potpaw_GGA" }
    import pymatgen
    import pymatgen.io.vasp.inputs
    from itertools import groupby
    import os

    if "coreLevelIndex" not in potcar_parameters.keys():
        forPOTCAR = [x[0] for x in groupby(input_atoms.get_chemical_symbols())]
    else:  ### Handles the case we need to select out an atom for Janak-Slater or final state core levels.
        coreLevelIndex = potcar_parameters["coreLevelIndex"]
        symbols = input_atoms.get_chemical_symbols()
        symbols[coreLevelIndex] = symbols[coreLevelIndex] + "_corelevel"
        forPOTCARTemp = [x[0] for x in groupby(symbols)]
        # coreLevelSpecies = [idx for idx, s in enumerate(forPOTCARTemp) if '_corelevel' in s][0] + 1
        forPOTCAR = [sym.replace("_corelevel", "") for sym in forPOTCARTemp]
    if "replaceNamesDict" in potcar_parameters.keys():
        replaceNamesDict = potcar_parameters["replaceNamesDict"]
        for origName in replaceNamesDict.keys():
            forPOTCAR = [
                replaceNamesDict[name] if name == origName else name
                for name in forPOTCAR
            ]
    try:  # Takes care of case if sym_potcar_map not set
        symPotcarMap = potcar_parameters["sym_potcar_map"]
    except KeyError:
        symPotcarMap = None
    potcarpmg = pymatgen.io.vasp.inputs.Potcar(
        symbols=forPOTCAR,
        functional=potcar_parameters["xc"],
        sym_potcar_map=symPotcarMap,
    )
    potcarpmg.write_file("POTCAR")
    if geometry_guess_quality.lower() == "bad":
        # TODO: Rewrite so it just copies the POTCAR to this folder.
        os.chdir("FastRelax")
        ### Next two lines shouldn't be needed, but just commenting them for now. Delete soon if things still work.
        # forPOTCAR = [x[0] for x in groupby(input_atoms.get_chemical_symbols())]
        # potcarpmg = pymatgen.io.vasp.inputs.Potcar(symbols = forPOTCAR, functional=potcar_parameters["xc"],sym_potcar_map = symPotcarMap)
        potcarpmg.write_file("POTCAR")


def make_kpoints(
    explicit_kpoints=[],
    gamma=True,
    kpoint_density=0,
    vacuum_directions=None,
    cell=None,
    kpoints_name="KPOINTS",
    factor=[1, 1, 1],
    geometry_guess_quality="Good",
    overwrite=True,
):
    import numpy as np
    from ase.calculators.vasp import Vasp
    import os

    # from ase.io.vasp import write_vasp
    """make_kpoints(explicit_kpoints=(7,7,1), gamma=True) or 
    make_kpoints(kpoint_density = 10, vacuum_directions=(False,False,True), cell = atoms.cell)
    Cell only needed if kPointImplicit used"""
    # Maybe TODO: allow control over k-point scheme.
    # Maybe TODO: look for conflicts in keywords.
    # Maybe TODO: Implement faster k-point grids
    if overwrite == False:
        if os.path.isfile(kpoints_name):
            print(kpoints_name + " already exists! Delete or set overwrite = True")
            return
    if overwrite == True:
        try:
            os.remove(kpoints_name)
        except OSError:
            pass
    if explicit_kpoints != []:
        explicit_kpoints = explicit_kpoints
    elif kpoint_density != 0:
        explicit_kpoints = [1, 1, 1]
        for i in range(3):
            if vacuum_directions[i] == True:
                explicit_kpoints[i] = 1
            else:
                explicit_kpoints[i] = kpoint_density / np.linalg.norm(cell[i])
    kpoints_scaled = list(np.array(explicit_kpoints) * factor)
    kpoints = {"kpts": kpoints_scaled, "gamma": gamma}
    calculate_kpoints = Vasp(**kpoints)

    try:  # This block is so that an existing KPOINTS file won't get overwritten if the new KPOINTS is supposed to have a different name.
        os.rename("KPOINTS", "KPOINTS_tmpfile")
    except (FileNotFoundError, IOError):
        pass
    calculate_kpoints.write_kpoints()
    os.rename("KPOINTS", kpoints_name)
    try:  # This block is so that an existing INCAR file won't get overwritten if the new INCAR is supposed to have a different name.
        os.rename("KPOINTS_tmpfile", "KPOINTS")
    except (FileNotFoundError, IOError):
        pass

    if (
        geometry_guess_quality.lower() == "bad"
    ):  # TODO: Pull this out and put it somewhere else (not in this function), and use the "factor" argument.
        # Make faster, lower level calc in subdirectory
        os.chdir("FastRelax")
        kpointsFastRaw = [
            int(np.rint(num * fac * 0.5)) for num, fac in zip(explicit_kpoints, factor)
        ]
        kpointsFastRaw = [max(1, kpoints) for kpoints in kpointsFastRaw]
        kpointsFast = {"kpts": kpointsFastRaw, "gamma": gamma}
        calculate_kpointsFast = Vasp(**kpointsFast)
        calculate_kpointsFast.write_kpoints()
        os.chdir("..")


def potcarSummary(file="POTCAR"):
    # file = 'POTCAR'
    with open(file) as potcar:
        for line in potcar:
            if "TITEL" in line:
                print(line)


def continueParams(continuation="False"):
    """A function that returns a good parameter set for relaxation, based on how good the guess is."""
    possibleSettings = ["False", "CHGCAR", "WAVECAR", "CHGCAR+WAVECAR"]
    if continuation not in possibleSettings:
        print("continuation must be one of these: ", possibleSettings)
    elif continuation == "False":
        stringCont = """ISTART = 0
ICHARG = 2
"""
    elif continuation == "CHGCAR":
        stringCont = """
TODO
"""
    elif continuation == "WAVECAR":
        stringCont = """ISTART = 1
ICHARG = 2
"""
    return stringCont


def make_incar(
    parameter_set,
    indcar_string="",
    geometry_guess_quality="Good",
    spin_polarized=True,
    magnetic_moments=None,
    magnetic_continuation=False,
    copy_incar_from=None,
    system_name=None,
    num_geometry_steps=0,
    incar_name="INCAR",
    verbose=True,
    continuation="False",
    type_of_run="Relaxation",
    write_directory="",
    overload_parameters=None,
    overwrite_parameters={},
    overwrite_file=True,
):
    """Makes an INCAR file. (Need to describe options.)
    magnetic_continuation: Run a spin-nonpolarized run at the beginning, then continue from that charge density. This often helps with magnetic convergence to the correct state.


    #### Examples: ####

    #Using predefined parameter set (will return ASE calculator object, in case you want to access calculator params):
    calc = make_incar(parameter_set='PW91_1',num_geometry_steps=300)

    #Manual INCAR, from a string:

    indcar_string='''SYSTEM=Name
    ENCUT=400
    GGA=91
    VOSKOWN=1'''
    make_incar(parameter_set='String',indcar_string=indcar_string)

    #Copy an existing INCAR:
    copy_incar_from = '/Users/Matt/Dropbox/Data/CUBackupData/VASPonShared/Done/MethylAds/CarbonAg/FHollow/INCAR'
    make_incar(parameter_set='Copy',copy_incar_from=copy_incar_from)

    """
    from pymatgen.io.vasp import Incar
    import copy
    import os

    if overwrite_file == False:
        if os.path.isfile(incar_name):
            print(incar_name + " already exists! Delete or set overwrite = True")
            return
    if overwrite_file == True:
        try:
            os.remove(incar_name)
        except OSError:
            pass
    if parameter_set == "String":
        import os

        if indcar_string == "":
            print(
                "WARNING: parameter_set is 'String', but indcar_string is empty, so I'm writing an empty INCAR file."
            )
        # This take a strings, and writes that to a file called "INCAR" in the write_directory.
        with open(os.path.join(write_directory, incar_name), "w+") as file:
            file.write(indcar_string)
        return
    elif parameter_set == "Copy":
        import os
        import shutil

        shutil.copy(copy_incar_from, os.path.join(write_directory, incar_name))
        return

        # Pre-defined parameter sets. For now, just export a preset text file for each; eventually may want something more sophisticated to easily do low-level initial calcs.
        # TODO: Make a separate function that processes the parameter sets, accounts for type_of_run, etc.
    elif parameter_set == "PW91_1":
        # A parameter set that's been used quite a bit for adsorption on transition metals.
        if type_of_run.lower() == "dos":
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 396.0,
                "PREC": "Normal",
                "GGA": "91",
                "VOSKOWN": 1,
                "NELMIN": 6,
                "EDIFF": 10**-6,
                "ISMEAR": -5,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": True,
                "LCHARG": True,
                "ICHARG": 2,
                "LORBIT": False,
            }
        elif type_of_run.lower() == "vibrational":
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 396.0,
                "PREC": "Normal",
                "GGA": "91",
                "VOSKOWN": 1,
                "NELMIN": 6,
                "EDIFF": 2 * 10**-6,
                "ISMEAR": 2,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": False,
                "LCHARG": False,
                "LORBIT": False,
                "IBRION": 7,
                "NSW": 1,
            }
        else:
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 396.0,
                "PREC": "Normal",
                "GGA": "91",
                "VOSKOWN": 1,
                "NELMIN": 6,
                "EDIFF": 10**-5,
                "NWRITE": 1,
                "ISMEAR": 2,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "EDIFFG": -0.03,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": False,
                "LCHARG": False,
                "LORBIT": False,
            }

    elif parameter_set == "PBE+TS_1":
        # A parameter set that's been used quite a bit for adsorption on transition metals.
        if type_of_run.lower() == "dos":
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 400.0,
                "PREC": "Accurate",
                "GGA": "PE",
                "NELMIN": 5,
                "EDIFF": 10**-6,
                "ISMEAR": -5,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "IVDW": 2,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": True,
                "LCHARG": True,
                "ICHARG": 2,
                "LORBIT": False,
            }
        elif type_of_run.lower() == "vibrational":
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 400.0,
                "PREC": "Accurate",
                "GGA": "PE",
                "NELMIN": 5,
                "EDIFF": 2 * 10**-6,
                "ISMEAR": 2,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "IVDW": 2,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": False,
                "LCHARG": False,
                "LORBIT": False,
                "IBRION": 7,
                "NSW": 1,
            }
        else:
            base_params = {
                "LREAL": "Auto",
                "ENCUT": 400.0,
                "PREC": "Accurate",
                "GGA": "PE",
                "NELMIN": 5,
                "EDIFF": 10**-5,
                "NWRITE": 1,
                "ISMEAR": 2,
                "SIGMA": 0.2,
                "AMIX": 0.2,
                "BMIX": 0.00001,
                "AMIX_MAG": 0.8,
                "BMIX_MAG": 0.00001,
                "EDIFFG": -0.03,
                "IVDW": 2,
                "ALGO": "VeryFast",
                "LWAVE": False,
                "LAECHG": False,
                "LVTOT": False,
                "LCHARG": False,
                "LORBIT": False,
            }
    elif parameter_set == "PBE_1":
        # Similar to PBE+TS_1, but without TS.
        base_params = {
            "LREAL": "Auto",
            "ENCUT": 400.0,
            "PREC": "Normal",
            "GGA": "PE",
            "NELMIN": 5,
            "EDIFF": 10**-5,
            "NWRITE": 1,
            "ISMEAR": 2,
            "SIGMA": 0.2,
            "AMIX": 0.2,
            "BMIX": 0.00001,
            "AMIX_MAG": 0.8,
            "BMIX_MAG": 0.00001,
            "EDIFFG": -0.03,
            "ALGO": "VeryFast",
            "LWAVE": False,
            "LAECHG": False,
            "LVTOT": False,
            "LCHARG": False,
            "LORBIT": False,
        }
    # elif parameter_set == 'PBE_Rodrigo':
    elif parameter_set == "PBE_Rodrigo":
        base_params = {
            "LREAL": "False",
            "ENCUT": 450.0,
            "GGA": "PE",
            "NELMIN": 5,
            "NELMDL": -7,
            "EDIFF": 10**-6,
            "NWRITE": 0,
            "ISMEAR": 0,
            "SIGMA": 0.03,
            "AMIX": 0.10,
            "BMIX": 3.00,
            "EDIFFG": -0.02,
            "ALGO": "Fast",
            "LWAVE": False,
            "LCHARG": False,
            "LVDW": True,
            "VDW_VERSION": 2,
            "VDW_RADIUS": 40,
            "VDW_SCALING": 0.75,
        }
    if spin_polarized == True:
        base_params.update({"ISPIN": 2})
    if magnetic_moments is not None:
        base_params.update({"MAGMOM": magnetic_moments})
    if geometry_guess_quality.lower() not in ["verygood", "good", "bad"]:
        print(
            "ERROR: geometry_guess_quality is not one of: ['VeryGood','Good','Bad'] (case shouldn't matter)"
        )

    if type_of_run.lower() == "lattice_constant":
        base_params.update({"ISIF": 3})

    if type_of_run.lower() == "dimer":
        base_params.update(
            {
                "IBRION": 3,
                "POTIM": 0.0,
                "IOPT": 2,
                "ICHAIN": 2,
                "DDR": 0.01,
                "DROTMAX": 4,
                "DFNMIN": 0.01,
                "DFNMAX": 1.0,
            }
        )

    if type_of_run.lower() == "charge":  ## Lightly tested
        base_params.update(
            {"PREC": "Accurate", "LVTOT": True, "LCHARG": True, "LAECHG": True}
        )

    if (
        geometry_guess_quality.lower() == "verygood"
        and type_of_run.lower() == "relaxation"
    ):
        base_params.update({"IBRION": 1, "POTIM": 0.1})
    elif (
        geometry_guess_quality.lower() == "good" and type_of_run.lower() == "relaxation"
    ):
        base_params.update({"IBRION": 3, "POTIM": 0.1})
    elif geometry_guess_quality.lower() == "bad":
        fastParams = copy.deepcopy(base_params)
        fastParams.update({"ENCUT": base_params["ENCUT"] * 0.66})
        fastParams.update({"EDIFF": 10 ** (-4)})
        fastParams.update({"EDIFFG": base_params["EDIFFG"] * 1.5})
        fastParams.update({"NELMIN": 1})
        if type_of_run.lower() == "relaxation":
            base_params.update({"IBRION": 1, "POTIM": 0.1})
            fastParams.update({"IBRION": 3, "POTIM": 0.02})
        elif type_of_run.lower() == "dimer":
            pass  # Could adjust dimer relaxation params here.
        try:
            os.mkdir("FastRelax")
        except FileExistsError:
            pass
        os.chdir("FastRelax")
        fastParams.update({"NSW": num_geometry_steps})
        fastParams.update(overwrite_parameters)
        fastIncar = Incar(fastParams)
        fastIncar.write_file("INCAR")
        os.chdir("..")

    # elif type_of_run == "MolecularDynamics":
    #    calc.set(ibrion=3)
    #     calc.set(smass=0)
    # Need to deal with temp (TEBEG, TEEND). Maybe just require these as overwrite_parameters, or just run at room temp the whole time by default

    base_params.update({"NSW": num_geometry_steps})
    base_params.update(overwrite_parameters)
    if magnetic_continuation:
        base_params.update({"MAXMIX": 35, "ISTART": 0, "ICHARG": 1, "LCHARG": True})
    ### Warnings here
    if type_of_run.lower() == "vibrational" and "NCORE" in base_params:
        print("WARNING!! This is a vibrational run but it has NCORE set. Job may fail.")
    #     calc.set(nsw=num_geometry_steps)
    #     calc.set(**overwrite_parameters)

    #     initializeForINCAR(calc, input_atoms)
    #     input_atoms.set_calculator(calc)
    try:  # This block is so that an existing INCAR file won't get overwritten if the new INCAR is supposed to have a different name.
        os.rename("INCAR", "INCAR_tmpfile")
    except (FileNotFoundError, IOError):
        pass
    #     calc.write_incar(input_atoms)
    incar = Incar(base_params)
    incar.write_file(incar_name)
    #     os.rename("INCAR",incar_name)
    try:  # This block is so that an existing INCAR file won't get overwritten if the new INCAR is supposed to have a different name.
        os.rename("INCAR_tmpfile", "INCAR")
    except (FileNotFoundError, IOError):
        pass


#     return calc


# May want to make this much more sophisticated in the future
def print_output_summary(outcarName="OUTCAR", printConvergenceState=True):
    jobFinished = False  # Initialize
    geomConverged = False
    with open(outcarName, "r") as f:
        for line in f.readlines():
            if (
                "reached required accuracy - stopping structural energy minimisation"
                in line
            ):
                geomConverged = True
            if "Elapsed time (sec):" in line:
                jobFinished = True
    if printConvergenceState:
        if jobFinished:
            print("Calculation finished.")
        else:
            print("WARNING: Calculation did not finish.")
        if geomConverged:
            print("The geometry converged.")
        else:
            print("WARNING: The geometry did not converge.")


def print_input_summary(calculation_directory):
    """Prints a summary of the input files in calculation_directory"""
    import os

    os.chdir(calculation_directory)
    if "FastRelax" in os.listdir():
        print("Subdirectory FastRelax contains files:")
        os.chdir("FastRelax")
        print(os.listdir())
        print("\n")
        os.chdir("..")
    print("Contents of files in " + str(calculation_directory))
    for file in os.listdir():
        if file in ["INCAR", "KPOINTS", "POSCAR", "MODECAR", "jobScript.sh"]:
            print(file + " file:")
            with open(file, "r") as openFile:
                print(openFile.read())
        if file == "POTCAR":
            print(file + " file:")
            potcarSummary(file)
