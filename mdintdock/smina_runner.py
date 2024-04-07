import subprocess
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import SDWriter

class SminaRunner:
    def __init__(self, output_dir="output", num_modes=10, exhaustiveness=5):
        self.output_dir = output_dir
        self.num_modes = num_modes
        self.exhaustiveness = exhaustiveness

    def run_docking(self, receptor_file, sdf_file, ligand_name, ligand_output_dir, reference_ligand=None):
        # Create output directory if it doesn't exist
        os.makedirs(self.output_dir, exist_ok=True)

        # Read ligands from the SDF file
        ligands = self._read_sdf_ligands(sdf_file, output_dir=ligand_output_dir)

        # Run docking for each ligand
        for ligand_name, ligand_data in ligands.items():
            # Run smina for this ligand-frame pair
            output_folder = os.path.join(self.output_dir, f"docking_{ligand_name}")
            os.makedirs(output_folder, exist_ok=True)

            # Prepare output file name
            frame_number = os.path.basename(receptor_file).split("_")[1].split(".")[0]
            output_file = os.path.join(output_folder, f"{ligand_name}_frame_{frame_number}.sdf")

            # Prepare log file name
            log_file = os.path.join(output_folder, f"{ligand_name}_frame_{frame_number}_log.txt")

            # Prepare command for smina
            command = [
                "smina",
                "--receptor", receptor_file,
                "--ligand", ligand_data,
                "--out", output_file,
                "--num_modes", str(self.num_modes),
                "--exhaustiveness", str(self.exhaustiveness),
            ]

            if reference_ligand:
                command.extend(["--autobox_ligand", reference_ligand])

            # Run smina
            try:
                subprocess.run(command, check=True)
                print(f"Docking for ligand {ligand_name} and frame {frame_number} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running smina for ligand {ligand_name} and frame {frame_number}: {e}")

    def _read_sdf_ligands(self, sdf_file, output_dir):
        """Read ligands from an SDF file and write each ligand to a separate SDF file."""
        ligands = {}

        try:
            suppl = Chem.SDMolSupplier(sdf_file)
            for idx, mol in enumerate(suppl):
                if mol is not None:
                    ligand_name = mol.GetProp("_Name")
                    output_file = os.path.join(output_dir, f"{ligand_name}.sdf")
                    writer = Chem.SDWriter(output_file)
                    writer.write(mol)
                    writer.close()
                    ligands[ligand_name] = output_file
        except FileNotFoundError:
            print(f"Error: SDF file '{sdf_file}' not found.")
            return None

        if not ligands:
            print(f"Error: No ligands found in SDF file '{sdf_file}'.")
            return None

        return ligands



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
            "--num_modes", "10",
            "--autobox_ligand", sdf_file,
            "--exhaustiveness", "10",
        ]
        try:
            subprocess.run(command, check=True)
            print(f"Reference docking for ligand {os.path.basename(sdf_file)} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running reference docking for ligand {os.path.basename(sdf_file)}: {e}")
