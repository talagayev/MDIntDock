import MDAnalysis as mda
import os
import prolif as plf
from rdkit import Chem
from rdkit.Chem import KekulizeException

def create_selections(topology_file, trajectory_file, ligand_resname):
    """
    Create selections for the ligand and protein based on the given topology file,
    trajectory file, and ligand residue name.

    Parameters:
        - topology_file (str): Path to the topology file (PDB format).
        - trajectory_file (str): Path to the trajectory file (DCD format).
        - ligand_resname (str): Residue name of the ligand.

    Returns:
        - ligand_selection (MDAnalysis.core.groups.AtomGroup): Selection for the ligand.
        - protein_selection (MDAnalysis.core.groups.AtomGroup): Selection for the protein.
    """
    # Load topology and trajectory
    u = mda.Universe(topology_file, trajectory_file)

    # Create selections for the ligand and protein
    ligand_selection = u.select_atoms(f"resname {ligand_resname}")
    protein_selection = u.select_atoms(
        f"(protein or resname WAT) and byres around 20.0 resname {ligand_resname}",
        ligand=ligand_selection,
    )

    return ligand_selection, protein_selection

def process_poses(frame_file, ligand_file, ligand_name):
    try:
        # Load the universe
        u = mda.Universe(frame_file)

        # Create protein molecule
        protein_mol = plf.Molecule.from_mda(u)

        # Load ligand molecule
        pose_iterable = plf.sdf_supplier(ligand_file)

        # Process poses with ligand and protein
        fp = plf.Fingerprint()
        fp.run_from_iterable(pose_iterable, protein_mol)

        # Convert fingerprint results to DataFrame
        df = fp.to_dataframe()

        # Add 'Frame' column and fill with frame number
        df.insert(0, 'Frame', int(os.path.basename(frame_file).split("_")[1].split(".")[0]))

        # Add 'ligand' column and fill with ligand name
        df.insert(1, 'ligand', ligand_name)

        return df

    except KekulizeException as e:
        print(f"Failed to kekulize molecule: {e}. Applying RDKit Path.")

        # Use RDKit to read Molecule from PDB
        rdkit_prot = Chem.MolFromPDBFile(frame_file, removeHs=False)

        # Create protein molecule
        protein_mol = plf.Molecule(rdkit_prot)

        # Load ligand molecule
        pose_iterable = plf.sdf_supplier(ligand_file)

        # Process poses with ligand and protein
        fp = plf.Fingerprint()
        fp.run_from_iterable(pose_iterable, protein_mol)

        # Convert fingerprint results to DataFrame
        df = fp.to_dataframe()

        # Add 'Frame' column and fill with frame number
        df.insert(0, 'Frame', int(os.path.basename(frame_file).split("_")[1].split(".")[0]))

        # Add 'ligand' column and fill with ligand name
        df.insert(1, 'ligand', ligand_name)

        return df
    

