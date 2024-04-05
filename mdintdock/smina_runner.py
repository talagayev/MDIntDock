import subprocess
import os

class SminaRunner:
    def __init__(self, output_dir="output", num_modes=9, exhaustiveness=8):
        self.output_dir = output_dir
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness

    def run_smina(self, receptor_file, sdf_file):
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Write ligand data to a temporary file
        ligand_filename = os.path.basename(sdf_file)
        ligand_basename = os.path.splitext(ligand_filename)[0]
        ligand_output_dir = os.path.join(self.output_dir, f"smina_{ligand_basename}")
        os.makedirs(ligand_output_dir, exist_ok=True)

        # Run smina for this ligand
        command = [
            "smina",
            "--receptor", receptor_file,
            "--ligand", sdf_file,
            "--out", f'{ligand_output_dir}/output.pdbqt',
            "--num_modes", str(self.num_modes),
            "--exhaustiveness", str(self.exhaustiveness),
            "--log", f'{ligand_output_dir}/conf.txt'
        ]
        try:
            subprocess.run(command, check=True)
            print(f"Ligand {ligand_basename} docking completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running smina for ligand {ligand_basename}: {e}")

    def run_ref_docking(self, receptor_file, sdf_file):
        # Create output directory for ref docking if it doesn't exist
        ref_docking_output_dir = os.path.join(self.output_dir, "ref_docking")
        os.makedirs(ref_docking_output_dir, exist_ok=True)

        # Get ligand number from the ligand file name
        ligand_number = int(os.path.basename(sdf_file).split("_")[1].split(".")[0])

        # Run reference docking for the ligand
        output_file = os.path.join(ref_docking_output_dir, f"lig_{ligand_number}_ref_dock.sdf")
        command = [
            "smina",
            "--receptor", receptor_file,
            "--ligand", sdf_file,
            "--out", output_file,
            "--num_modes", "25",
            "--autobox_ligand", sdf_file,
            "--exhaustiveness", "50",
            "--log", f'{output_file}.log'
        ]
        try:
            subprocess.run(command, check=True)
            print(f"Reference docking for ligand {os.path.basename(sdf_file)} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running reference docking for ligand {os.path.basename(sdf_file)}: {e}")



