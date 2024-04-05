import os
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import SDWriter

class FrameWriter:
    def __init__(self, topology_file, trajectory_file, lig_resname, output_directory=None):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.lig_resname = lig_resname

        if output_directory is None:
            # Extract filename without extension from topology_file
            folder_name = os.path.splitext(os.path.basename(topology_file))[0]
            self.output_directory = os.path.join(folder_name + f"_{lig_resname}")
        else:
            self.output_directory = output_directory

        # Create output directory if it doesn't exist
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

        # Create subdirectories
        self.protein_pdb_directory = os.path.join(self.output_directory, "frames_protein_pdb")
        self.refligand_sdf_directory = os.path.join(self.output_directory, "frames_refligand_sdf")

        if not os.path.exists(self.protein_pdb_directory):
            os.makedirs(self.protein_pdb_directory)
        if not os.path.exists(self.refligand_sdf_directory):
            os.makedirs(self.refligand_sdf_directory)

    def write_protein_frame(self, frame_number):
        u = mda.Universe(self.topology_file, self.trajectory_file)
        protein = u.select_atoms("protein")
        frame = protein.positions
        with mda.Writer(os.path.join(self.protein_pdb_directory, f"frame_{frame_number}.pdb"), protein.n_atoms) as pdb:
            for atom in protein.atoms:
                atom.position = frame[atom.index]
            pdb.write(protein)

    def write_ligand_frame(self, frame_number):
        u = mda.Universe(self.topology_file, self.trajectory_file)
        ligand = u.select_atoms(f"resname {self.lig_resname}")
        mol = ligand.atoms.convert_to("RDKIT")
        with Chem.SDWriter(os.path.join(self.refligand_sdf_directory, f'ligand_{frame_number}.sdf')) as w:
            w.write(mol)

