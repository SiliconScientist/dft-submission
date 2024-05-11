import toml
from config import Config
from vasp_tools import makeAllVASPInput, print_input_summary, print_output_summary
from server_interface import upsync, downsync
import ase
from ase.io import read
from ase.visualize import view
from pathlib import Path
import os

import time


def main():
    config = Config(**toml.load("config.toml"))
    structure_names = os.listdir(config.paths.structures)
    calculation_directory_list = [
        config.paths.local_personal_directory
        + config.paths.project_name
        + config.paths.sub_project_name
        + structure_name
        for structure_name in structure_names
    ]
    local_project_directory = (
        config.paths.local_personal_directory
        + config.paths.project_name
        + config.paths.sub_project_name
    )
    remote_project_directory = (
        config.paths.remote_personal_directory
        + "/"
        + config.paths.project_name
        + config.paths.sub_project_name
    )
    if config.submit == True:
        input_atoms_list = [
            read(
                filename=config.paths.structures + "/" + structure_name,
                index=0,
                format="vasp",
            )
            for structure_name in structure_names
        ]
        for input_atoms, calculation_directory in zip(
            input_atoms_list, calculation_directory_list
        ):
            makeAllVASPInput(
                input_atoms,
                calculation_directory,
                config.meta_parameters.dict(),
                config.kpoint_parameters.dict(),
                config.potcar_parameters.dict(),
                config.job_parameters.dict(),
                config.dimer_parameters,
            )
            print_input_summary(calculation_directory)

        time.sleep(3)  # Maybe not necessary
        upsync(
            calculation_directory_list,
            config.paths.remote_host_name,
            remote_project_directory,
            options="rv",
        )
    else:
        remote_project_directory_list = [
            remote_project_directory / structure_name
            for structure_name in structure_names
        ]
        downsync(
            remote_project_directory_list,
            config.paths.remote_host_name,
            local_project_directory,
            options="rv",
        )

        for structure_name in structure_names:
            os.chdir(os.path.join(local_project_directory, structure_name))

            try:
                print_output_summary()
                atoms = ase.io.read("OUTCAR")
                print("Final energy: ", atoms.get_potential_energy())
                atoms = ase.io.read("CONTCAR")
                atoms.wrap(center=(0.3, 0.3, 0.8))
                atoms.set_constraint()  # Needed for repeat to work correctly
                view(atoms.repeat((2, 2, 1)))
            except FileNotFoundError:
                print("No OUTCAR found.")
                continue
            except ValueError:
                print("ValueError")
                continue
            print()

        # TODO: Look at how slapdash this is. Make it legit!
        if config.meta_parameters.type_of_run == "lattice_constant":
            atomsAll1 = ase.io.read("OUTCAR-1@:")
            atomsIn = atomsAll1[0]
            atomsAll2 = ase.io.read("OUTCAR-1@:")
            atomsOut = atomsAll2[-1]
            print(atomsIn.cell)
            print(atomsOut.cell)
            print("lattice constant, if fcc:", atomsOut.cell[0, 1] * 2)


if __name__ == "__main__":
    main()
