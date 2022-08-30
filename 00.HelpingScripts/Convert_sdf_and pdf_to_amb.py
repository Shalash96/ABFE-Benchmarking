import BioSimSpace as BSS

import parmed

from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField as OFF_ForceField

try:
    from openmm import XmlSerializer, app, unit
    from openmm.app import HBonds, NoCutoff, PDBFile, ForceField
except ImportError:
    from simtk import unit
    from simtk.openmm import XmlSerializer, app
    from simtk.openmm.app import HBonds, NoCutoff, PDBFile




def reparam(lig_no):
    # Load in the ligand.
    ligand_molecule = Molecule(f"structures/cyclod/ligand{lig_no}.sdf")

    # Specify the "Sage" forcefield.
    force_field = OFF_ForceField("openff_unconstrained-2.0.0.offxml")

    # Parametrize the ligand molecule by creating a Topology object from it.
    ligand_system = force_field.create_openmm_system(ligand_molecule.to_topology())

    # Read in the coordinates of the ligand from the PDB file.
    ligand_pdbfile = PDBFile(f"structures/cyclod/ligand{lig_no}.pdb")

    # Convert the ligand system to a ParmEd object.
    ligand_parmed_structure = parmed.openmm.load_topology(ligand_pdbfile.topology,
                                                        ligand_system,
                                                        ligand_pdbfile.positions)

    # Write the ligand structure to AMBER format.
    ligand_parmed_structure.save(f"tmp/ligand{lig_no}.rst7", overwrite=True)
    ligand_parmed_structure.save(f"tmp/ligand{lig_no}.parm7", overwrite=True)

    # Parse the protein PDB file.
    protein_pdbfile = PDBFile("structures/cyclod/protein_water.pdb")

    # Load the AMBER protein force field through OpenMM.
    omm_forcefield = app.ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")

    # Parameterize the protein.
    protein_system = omm_forcefield.createSystem(protein_pdbfile.topology,
                                                nonbondedCutoff=1*unit.nanometer,
                                                nonbondedMethod=app.NoCutoff,
                                                constraints=None,
                                                rigidWater=False)

    # Convert the protein System into a ParmEd Structure.
    protein_parmed_structure = parmed.openmm.load_topology(
                                            protein_pdbfile.topology,
                                            protein_system,
                                            xyz=protein_pdbfile.positions)

    # Write the protein structure to AMBER format.
    protein_parmed_structure.save("tmp/protein.rst7", overwrite=True)
    protein_parmed_structure.save("tmp/protein.parm7", overwrite=True)

    # Load the ligand and protein with BioSimSpace.
    ligand = BSS.IO.readMolecules(f"tmp/ligand{lig_no}.*7")[0]
    protein = BSS.IO.readMolecules("tmp/protein.*7")[0]

    # Combine to form a complex.
    complx = (ligand + protein).toSystem()

    #Save system
    BSS.IO.saveMolecules(f"reparam_cyclod/lig{lig_no}_complex", complx, ["rst7", "prm7"])





for i in range(1,4):
    reparam(i)
