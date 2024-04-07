import os
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import SDWriter

class FrameWriter:
    def __init__(self, topology_file, trajectory_file, lig_resname):
        self.topology_file = topology_file
        self.trajectory_file = trajectory_file
        self.lig_resname = lig_resname
        self.protein_pdb_directory = "protein_frames"
        self.refligand_sdf_directory = "ligand_frames"

        # Create directories if they don't exist
        os.makedirs(self.protein_pdb_directory, exist_ok=True)
        os.makedirs(self.refligand_sdf_directory, exist_ok=True)

    def write_frames(self):
        u = mda.Universe(self.topology_file, self.trajectory_file)

        # Create output directories if they don't exist
        os.makedirs(self.protein_pdb_directory, exist_ok=True)
        os.makedirs(self.refligand_sdf_directory, exist_ok=True)

        # Loop through each frame in the trajectory
        for ts in u.trajectory:
            frame_number = ts.frame + 1  # Frame number starts from 1

            # Write protein frame
            protein = u.select_atoms(f"(protein and around 10.0 resname {self.lig_resname})")
            with mda.Writer(os.path.join(self.protein_pdb_directory, f"frame_{frame_number}.pdb"), protein.n_atoms) as pdb:
                pdb.write(protein)

            # Write ligand frame
            ligand = u.select_atoms(f"resname {self.lig_resname}")
            mol = ligand.atoms.convert_to("RDKIT")
            with Chem.SDWriter(os.path.join(self.refligand_sdf_directory, f'ligand_{frame_number}.sdf')) as w:
                w.write(mol)

