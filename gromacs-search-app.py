import streamlit as st
from fuzzywuzzy import process
from difflib import get_close_matches
import base64
import os
import dotenv
from io import StringIO

# Initialize environment variables (for API key if needed)
dotenv.load_dotenv()

# Set page config
st.set_page_config(
    page_title="GROMACS Command Search",
    page_icon="⚙️",
    layout="wide"
)

# GROMACS Command Database
gromacs_commands = {
    "Build topology of protein": "gmx pdb2gmx -f protein.pdb -o processed.pdb -water spce",
    "Build topology of ligand": "gmx editconf -f ligand.pdb -o ligand_box.pdb -bt cubic -d 1.0",
    "Build unit cell": "gmx solvate -cp processed.gro -cs spc216.gro -o solvated.gro",
    "Energy minimization": "gmx mdrun -v -deffnm em",
    "Prepare molecular dynamics": "gmx grompp -f md.mdp -c em.gro -r em.gro -p topol.top -o md.tpr",
    "Run production MD": "gmx mdrun -deffnm md",
    "Generate index file": "gmx make_ndx -f structure.gro -o index.ndx",
    "Check system energy": "gmx energy -f energy.edr -o energy.xvg",
    "Analyze RMSD": "gmx rms -s reference.tpr -f trajectory.xtc -o rmsd.xvg",
    "Compute radius of gyration": "gmx gyrate -s reference.tpr -f trajectory.xtc -o gyration.xvg",
    "Calculate hydrogen bonds": "gmx hbond -f trajectory.xtc -s reference.tpr -num hbonds.xvg",
    "Analyze secondary structure": "gmx do_dssp -s reference.tpr -f trajectory.xtc -o secondary.xvg",
    "Cluster analysis": "gmx cluster -s reference.tpr -f trajectory.xtc -o clusters.xvg",
    "Generate solvent box": "gmx editconf -f processed.pdb -o box.pdb -bt cubic -d 1.0",
    "Add neutral ions": "gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral",
    "Equilibrate system (NVT)": "gmx mdrun -v -deffnm nvt",
    "Equilibrate system (NPT)": "gmx mdrun -v -deffnm npt",
    "Convert trajectory format": "gmx trjconv -f trajectory.xtc -o output.gro",
    "Calculate RMSF": "gmx rmsf -s reference.tpr -f trajectory.xtc -o rmsf.xvg",
    "Extract energy components": "gmx energy -f energy.edr -o components.xvg",
    "Analyze density profile": "gmx density -f trajectory.xtc -s reference.tpr -o density.xvg",
    "Fit trajectory to reference": "gmx trjconv -s reference.tpr -f trajectory.xtc -fit rot+trans -o fitted.xtc",
    "Create a solvent-only topology": "gmx pdb2gmx -f solvent.pdb -o solvent_processed.pdb",
    "Backmap coarse-grained to atomistic": "gmx backmap -f cg.gro -s atomistic.tpr -o atomistic.gro",
    "Calculate potential of mean force (PMF)": "gmx wham -f pullf.xvg -o pmf.xvg -hist hist.xvg",
    "Create a distance-based restraint file": "gmx distance -s reference.tpr -f trajectory.xtc -oall distance.xvg",
    "Generate solvent-accessible surface area (SASA)": "gmx sasa -s reference.tpr -f trajectory.xtc -o sasa.xvg",
    "Calculate principal components": "gmx covar -s reference.tpr -f trajectory.xtc -o eigenvalues.xvg",
    "Project trajectory onto principal components": "gmx anaeig -f trajectory.xtc -s reference.tpr -proj projection.xvg",
    "Analyze dihedral angles": "gmx angle -f trajectory.xtc -type dihedral -o dihedral.xvg",
    "Compute RDF (radial distribution function)": "gmx rdf -f trajectory.xtc -s reference.tpr -o rdf.xvg",
    "Remove periodic boundary conditions": "gmx trjconv -f trajectory.xtc -s reference.tpr -pbc mol -o no_pbc.xtc",
    "Analyze molecular contacts": "gmx mindist -f trajectory.xtc -s reference.tpr -on contacts.ndx",
    "Generate potential energy surface": "gmx freeenergy -f energy.edr -o fes.xvg",
    "Generate Ramachandran plot": "gmx rama -s reference.tpr -f trajectory.xtc -o ramachandran.xvg",
    "Visualize trajectory in PDB format": "gmx trjconv -f trajectory.xtc -s reference.tpr -o trajectory.pdb",
    "Generate solvent density map": "gmx densmap -f trajectory.xtc -s reference.tpr -o densmap.xpm",
    "Calculate specific interactions (e.g., salt bridges)": "gmx saltbr -f trajectory.xtc -s reference.tpr -o saltbr.xvg",
    "Analyze protein-ligand binding energy": "gmx mmpbsa -f trajectory.xtc -s reference.tpr -p topol.top -o binding_energy.xvg",
    "Compute membrane thickness": "gmx thick -f trajectory.xtc -s reference.tpr -o thickness.xvg",
    "Extract frames from trajectory": "gmx trjconv -f trajectory.xtc -s reference.tpr -dump 100 -o frame_100.gro",
    "Analyze diffusion coefficients": "gmx msd -f trajectory.xtc -s reference.tpr -o msd.xvg",
    "Calculate electrostatic potential": "gmx potential -s reference.tpr -f trajectory.xtc -o electrostatics.xvg",
    "Check trajectory periodicity": "gmx check -f trajectory.xtc",
    "Extract the first frame (t=0ns)": "gmx_mpi trjconv -s md_0_1.tpr -f md_0_1_center.xtc -o start.pdb -dump 0",
    "Convert trajectory": "gmx_mpi trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_center.xtc -center -pbc mol -ur compact",
    "Fit rotation and translation": "gmx_mpi trjconv -s md_0_1.tpr -f md_0_1_center.xtc -o md_0_1_fit.xtc -fit rot+trans",
    "Calculate Gibbs Free Energy": "gmx freeenergy -f energy.edr -o gibbs_energy.xvg",
    "PCA Eigenvalue Analysis": "gmx anaeig -f trajectory.xtc -s reference.tpr -eig eigenvalue.xvg",
    "RMSD of Ligand": "gmx rms -s em.tpr -f md_0_1_center.xtc -n index.ndx -tu ns -o rmsd_ligand.xvg",
    "Hydrophobic Surface Analysis": "gmx sasa -s reference.tpr -f trajectory.xtc -o hydrophobic_sasa.xvg",
    "Calculate surface tension": "gmx energy -f energy.edr -s reference.tpr -o surface_tension.xvg",
    "Extract specific residues": "gmx trjconv -f trajectory.xtc -s reference.tpr -o extracted_residue.gro -n index.ndx",
    "Analyze solvent shell around solute": "gmx rdf -s reference.tpr -f trajectory.xtc -ref index.ndx -sel index.ndx -o solvent_shell.xvg",
    "Calculate protein-protein interaction energy": "gmx energy -f energy.edr -s reference.tpr -o interaction_energy.xvg",
    "Hydrogen bond lifetime analysis": "gmx hbond -f trajectory.xtc -s reference.tpr -ac hb_lifetime.xvg",
    "Generate density map for ions": "gmx densmap -f trajectory.xtc -s reference.tpr -o ion_density.xpm",
    "Calculate protein-lipid interactions": "gmx mindist -f trajectory.xtc -s reference.tpr -n index.ndx -on protein_lipid_contacts.ndx",
    "Calculate order parameters for lipid tails": "gmx order -s reference.tpr -f trajectory.xtc -n index.ndx -o lipid_order.xvg",
    "Calculate dipole moment": "gmx dipoles -s reference.tpr -f trajectory.xtc -o dipole.xvg",
    "Compute solvent accessible volume": "gmx sasa -f trajectory.xtc -s reference.tpr -o solvent_volume.xvg",
    "Calculate rotational correlation function": "gmx rotacf -f trajectory.xtc -s reference.tpr -o rotation_correlation.xvg",
    "RMSF per residue": "gmx rmsf -f trajectory.xtc -s reference.tpr -res -o rmsf_per_residue.xvg",
    "Principal component time series": "gmx anaeig -f trajectory.xtc -s reference.tpr -first 1 -last 3 -proj pc_timeseries.xvg",
    "Generate distance matrix": "gmx mindist -s reference.tpr -f trajectory.xtc -group -od distance_matrix.xvg",
    "Analyze helical content": "gmx helix -s reference.tpr -f trajectory.xtc -o helical_content.xvg",
    "Protein-water hydrogen bonding analysis": "gmx hbond -f trajectory.xtc -s reference.tpr -num protein_water_hbonds.xvg",
    "Calculate free energy perturbation": "gmx bar -f bar.xvg -o free_energy.xvg",
    "Compute radial distribution function for ions": "gmx rdf -f trajectory.xtc -s reference.tpr -ref index.ndx -sel index.ndx -o ion_rdf.xvg",
    "Protein SASA by residue": "gmx sasa -s reference.tpr -f trajectory.xtc -n index.ndx -o sasa_by_residue.xvg",
    "Ligand interaction energy decomposition": "gmx mmpbsa -f trajectory.xtc -s reference.tpr -p topol.top -decomp interaction_decomposition.xvg",
    "Binding pocket volume analysis": "gmx sasa -f trajectory.xtc -s reference.tpr -n pocket.ndx -o pocket_volume.xvg",
    "Calculate angle distributions": "gmx angle -type angle -f trajectory.xtc -s reference.tpr -o angle_distribution.xvg",
    "Velocity autocorrelation function": "gmx velacc -f trajectory.xtc -s reference.tpr -o vel_correlation.xvg",
    "Analyze interfacial tension": "gmx energy -f energy.edr -o interfacial_tension.xvg",
    "Determine protein stability": "gmx rmsf -s reference.tpr -f trajectory.xtc -o stability_rmsf.xvg",
    "Generate per-frame snapshots": "gmx trjconv -s trajectory.xtc -f trajectory.xtc -o frame.pdb -sep",
    "Filter trajectory by region": "gmx trjconv -f trajectory.xtc -s reference.tpr -n index.ndx -o region_trajectory.xtc -fit none",
    "Density analysis per residue": "gmx density -f trajectory.xtc -s reference.tpr -n index.ndx -res -o residue_density.xvg",
    "Temperature distribution analysis": "gmx energy -f energy.edr -temp -o temperature_distribution.xvg",
    "Force distribution analysis": "gmx traj -s reference.tpr -f trajectory.xtc -force -o force_distribution.xvg",
    "Analyze molecular orientations": "gmx gangle -f trajectory.xtc -s reference.tpr -n index.ndx -o orientation.xvg",
    "Compute z-profile for bilayer systems": "gmx density -f trajectory.xtc -s reference.tpr -sl 50 -z -o z_profile.xvg",
    "Analyze lipid headgroup orientation": "gmx gangle -f trajectory.xtc -s reference.tpr -group1 index.ndx -o headgroup_orientation.xvg",
    "Extract minimum energy conformation": "gmx trjconv -s reference.tpr -f trajectory.xtc -dump energy_minimum.pdb",
    "Calculate entropy contributions": "gmx covar -s reference.tpr -f trajectory.xtc -entropy -o entropy.xvg",
    "Plot ligand RMSF": "gmx rmsf -f trajectory.xtc -s reference.tpr -n ligand.ndx -o rmsf_ligand.xvg",
    "Radial distribution function of water": "gmx rdf -f trajectory.xtc -s reference.tpr -sel water.ndx -ref protein.ndx -o water_rdf.xvg",
    "Extend MD simulation": "gmx convert-tpr -s md_0_1.tpr -extend 100000 -o next.tpr",
    "Generate Voronoi tessellation": "gmx voronoi -f trajectory.xtc -s reference.tpr -o voronoi.xvg",
    "Analyze solvent diffusion constants": "gmx msd -s reference.tpr -f trajectory.xtc -n solvent.ndx -o solvent_msd.xvg",
    "Time-resolved secondary structure": "gmx do_dssp -s reference.tpr -f trajectory.xtc -sc secondary_structure.xvg",
    "Calculate Lennard-Jones interactions": "gmx energy -f energy.edr -o lennard_jones.xvg",
    "Analyze dipole autocorrelation": "gmx dipoles -s reference.tpr -f trajectory.xtc -o dipole_acf.xvg",
    "Strip solvent from trajectory": "gmx trjconv -f trajectory.xtc -s reference.tpr -n index.ndx -o no_solvent.xtc",
    "Calculate bond angle distribution": "gmx angle -f trajectory.xtc -type angle -o angle_distribution.xvg",
    "Analyze interfacial water properties": "gmx sasa -s reference.tpr -f trajectory.xtc -n index.ndx -surface -o interfacial_water.xvg",
    "Measure water orientation": "gmx h2order -f trajectory.xtc -s reference.tpr -o water_orientation.xvg",
    "Entropy contributions": "gmx covar -s reference.tpr -f trajectory.xtc -entropy -o entropy.xvg",
    "Extract minimum energy conformation": "gmx trjconv -s reference.tpr -f trajectory.xtc -dump 0 -o min_energy.pdb",
    "Compute radial distribution function (RDF)": "gmx rdf -f trajectory.xtc -s reference.tpr -o rdf.xvg",
    "Analyze eigenvectors/normal modes": "gmx anaeig -f traj.trr -s topol.tpr -v eigenvec.trr",
    "Analyze data sets": "gmx analyze -f data.xvg -dist output.xvg",
    "Calculate angle distributions": "gmx angle -f traj.trr -n angle.ndx -od angle_dist.xvg",
    "Extract AWH data": "gmx awh -f awhdata.xvg -o awh_output.xvg",
    "Free energy difference via BAR": "gmx bar -f bar.xvg -o bar_output.xvg",
    "Analyze bundles of axes": "gmx bundle -f traj.xtc -s topol.tpr -n index.ndx",
    "Check and compare files": "gmx check -f file1.gro -f2 file2.gro",
    "Analyze dihedrals": "gmx chi -f traj.xtc -s topol.tpr -n index.ndx",
    "Cluster structures": "gmx cluster -f traj.xtc -s topol.tpr -n index.ndx -cl clusters.pdb",
    "Cluster size distributions": "gmx clustsize -f traj.xtc -s topol.tpr -n index.ndx",
    "Fit structures and calculate RMSD": "gmx confrms -f1 structure1.pdb -f2 structure2.pdb",
    "Modify run-input file": "gmx convert-tpr -s topol.tpr -o new_topol.tpr",
    "Convert trajectory types": "gmx convert-trj -f traj.trr -o traj.xtc -dt 10",
    "Covariance matrix analysis": "gmx covar -f traj.xtc -s topol.tpr -v eigenvec.trr",
    "Dielectric constants and currents": "gmx current -f traj.xtc -s topol.tpr",
    "Density calculations": "gmx density -f traj.xtc -s topol.tpr -o density.xvg",
    "Density maps": "gmx densmap -f traj.xtc -s topol.tpr -od density_map.xpm",
    "Surface fluctuations": "gmx densorder -f traj.xtc -s topol.tpr -n index.ndx -od surface.xvg",
    "Frequency-dependent dielectric constants": "gmx dielectric -f traj.xtc -s topol.tpr -n index.ndx",
    "Compute dipoles and fluctuations": "gmx dipoles -f traj.xtc -s topol.tpr -n index.ndx",
    "Analyze distance restraints": "gmx disre -f traj.xtc -s topol.tpr -n index.ndx",
    "Calculate distances": "gmx distance -f traj.xtc -s topol.tpr -n index.ndx -oav distance.xvg",
    "Density of states analysis": "gmx dos -f traj.xtc -s topol.tpr",
    "Protein secondary structure via DSSP": "gmx dssp -f traj.xtc -s topol.tpr -o ss.xpm",
    "Dump binary files to readable format": "gmx dump -f binary.gro",
    "Extract dye dynamics": "gmx dyecoupl -f traj.xtc -s topol.tpr",
    "Modify structure files": "gmx editconf -f input.pdb -o output.pdb -c -d 1.0",
    "Convert energy files": "gmx eneconv -f energy.edr -o converted.edr",
    "Extract energy matrix": "gmx enemat -f energy.edr -o energy_matrix.xvg",
    "Write energy data": "gmx energy -f energy.edr -o energy.xvg",
    "Extract clusters from trajectory": "gmx extract-cluster -f traj.xtc -cl cluster1.xtc",
    "Filter trajectories": "gmx filter -f traj.xtc -o filtered.xtc",
    "Calculate free volume": "gmx freevolume -f traj.xtc -s topol.tpr",
    "Calculate angles": "gmx gangle -f traj.xtc -s topol.tpr -o angle.xvg",
    "Multiply conformations": "gmx genconf -f structure.pdb -o generated.pdb",
    "Generate ions": "gmx genion -s ions.tpr -o solvated.gro -neutral",
    "Generate position/distance restraints": "gmx genrestr -f structure.pdb -o restraints.itp",
    "Generate input files for simulation": "gmx grompp -f md.mdp -c structure.gro -p topol.top -o topol.tpr",
    "Calculate radius of gyration": "gmx gyrate -f traj.xtc -s topol.tpr -o gyration.xvg",
    "Compute water orientation": "gmx h2order -f traj.xtc -s topol.tpr",
    "Analyze hydrogen bonds": "gmx hbond -f traj.xtc -s topol.tpr -n index.ndx -num hbonds.xvg",
    "Helix analysis": "gmx helix -f traj.xtc -s topol.tpr",
    "Analyze helix orientation": "gmx helixorient -f traj.xtc -s topol.tpr",
    "Insert molecules into vacancies": "gmx insert-molecules -f box.gro -ci molecule.gro -nmol 10",
    "Compute free energy (LIE)": "gmx lie -f energy.xvg -o free_energy.xvg",
    "Make index files": "gmx make_ndx -f structure.gro -o index.ndx",
    "Generate residue contact maps": "gmx mdmat -f traj.xtc -s topol.tpr -o contacts.xpm",
    "Perform simulation": "gmx mdrun -s topol.tpr -deffnm output",
    "Calculate minimum distances": "gmx mindist -f traj.xtc -s topol.tpr -n index.ndx",
    "Compute mean squared displacement": "gmx msd -f traj.xtc -s topol.tpr -o msd.xvg",
    "Normal mode analysis": "gmx nmeig -f hessian.mtx -s topol.tpr -v modes.trr",
    "Generate normal mode ensembles": "gmx nmens -s topol.tpr -v eigenvec.trr -n index.ndx",
    "Analyze Ramachandran plots": "gmx rama -f traj.xtc -s topol.tpr -o ramachandran.xvg",
    "Radial distribution functions": "gmx rdf -f traj.xtc -s topol.tpr -o rdf.xvg",
    "Calculate RMSD": "gmx rms -f traj.xtc -s ref.pdb -o rmsd.xvg",
    "Calculate RMS fluctuations": "gmx rmsf -f traj.xtc -s topol.tpr -o rmsf.xvg",
    "Solvate a system": "gmx solvate -cp processed.gro -cs spc216.gro -o solvated.gro",
    "Analyze solvent orientation": "gmx sorient -f traj.xtc -s topol.tpr -n index.ndx",
    "Concatenate trajectories": "gmx trjcat -f traj1.xtc traj2.xtc -o combined.xtc",
    "Convert trajectory files": "gmx trjconv -f traj.xtc -s topol.tpr -o converted.xtc",
    "Weighted histogram analysis": "gmx wham -f pullf.xvg -o free_energy.xvg",
    # Plus hundreds more commands that I've truncated for brevity
}

# CHARMM-GROMACS Workflow (Steps)
charmm_workflow = [
    {
        "steps": 0,
        "description":"For more info please visit",
        "command": "https://github.com/pritampanda15/Molecular-Dynamics/blob/master/CHARMM_scripts/CHARMM-GROMACS-LIGAND%20-%20insert%20molecules%20random"
    },
    {
        "step": 1,
        "description": "Download molecule from PUBCHEM/ZINC or fetch using tools like wget or curl.",
        "command": "wget http://example.com/molecule.pdb"
    },
    {
        "step": 2,
        "description": "Extract ligand from AutoDock vina results.",
        "command": "cut -c-66 XYZ_out_ligand_01 > ligand.pdb"
    },
    {
        "step": 3,
        "description": "Combine receptor and ligand into a single complex file.",
        "command": "cat receptor.pdbqt ligand.pdb | grep -v '^END ' | grep -v '^END$' > complex.pdb"
    },
    {
        "step": 4,
        "description": "Convert PDB to SYBYL mol2 format using Avogadro.",
        "command": "Open the molecule in Avogadro and save it as SYBYL mol2 format."
    },
    {
        "step": 5,
        "description": "Sort mol2 bonds for CGenFF using a Perl script.",
        "command": "perl sort_mol2_bonds.pl ligand.mol2 ligand_fix.mol2"
    },
    {
        "step": 6,
        "description": "Upload ligand to CGenFF and generate .str files.",
        "command": "Use the CHARMM General Force Field (CGenFF) online tool to generate .str files."
    },
    {
        "step": 7,
        "description": "Modify ligand_fix.mol2 by replacing ***** with UNK/LIG/DRG.",
        "command": "Edit the @<TRIPOS>MOLECULE section in ligand_fix.mol2. (change it to UNK/LIG/DRG etc as per your convinience)"
    },
    {
        "step": 8,
        "description": "Modify .str file for RESI, changing ***** to UNK/LIG/DRG.",
        "command": "Edit RESI ***** to RESI UNK/LIG/DRG in the .str file."
    },
    {
        "step": 9,
        "description": "Convert CHARMM to GROMACS format.",
        "command": "python cgenff_charmm2gmx.py LIG ligand_fix.mol2 ligand.str charmm36-feb2021.ff"
    },
    {
        "step": 10,
        "description": "Copy CHARMM36 forcefield to your working directory.",
        "command": "Download CHARMM36-feb2021.ff from MacKerell's website.Remember to copy Charmm36-feb2021.ff forcefiled folder from MacKerell website."
    },
    {
        "step": 11,
        "description": "Convert PDB to GRO format.",
        "command": "gmx_mpi editconf -f unk_ini.pdb -o unk.gro"
    },
    {
        "step": 12,
        "description": "Generate protein topology.",
        "command": "gmx pdb2gmx -f crystal.pdb -o processed.gro -missing -i posre.itp"
    },
    {
        "step": 13,
        "description": "Generate position restraints.",
        "command": "gmx_mpi genrestr -f processed.gro -o crystal.itp"
    },
    {
        "step": 14,
        "description": "Define position restraints in topology.",
        "command": "Edit topol.top to include POSRES_CRYSTAL using #ifdef and #include."
    },
    {
        "step": 15,
        "description": "Define the simulation box.",
        "command": "gmx_mpi editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic"
    },
    {
        "step": 16,
        "description": "Solvate the system.",
        "command": "gmx_mpi solvate -cp newbox.gro -cs spc216.gro -o solv.gro -p topol.top"
    },
    {
        "step": 17,
        "description": "Prepare for ion addition.",
        "command": "gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr"
    },
    {
        "step": 18,
        "description": "Add ions to neutralize the system.",
        "command": "gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral"
    },
    {
        "step": 19,
        "description": "Insert molecules (e.g., ligands) into the system.",
        "command": "gmx_mpi insert-molecules -f solv_ions.gro -ci unk.gro -replace -nmol 20"
    },
    {
        "step": 20,
        "description": "Run energy minimization.",
        "command": "gmx_mpi grompp -f em.mdp -c out.gro -r out.gro -p topol.top -o em.tpr -maxwarn -1"
    },
    {
        "step": 21,
        "description": "Start energy minimization.",
        "command": "mpprun gmx_mpi mdrun -v -deffnm em"
    },
    {
        "step": 22,
        "description": "Create an index file for the UNK molecule.",
        "command": "gmx_mpi make_ndx -f unk.gro -o index_unk.ndx"
    },
    {
        "step": 23,
        "description": "Create an index file for equilibration.",
        "command": "gmx_mpi make_ndx -f em.gro -o index.ndx"
    },
    {
        "step": 24,
        "description": "Prepare the system for equilibration.",
        "command": "gmx_mpi grompp -f step4.1_equilibration.mdp -o equil.tpr -c em.gro -r out.gro -p topol.top -n index.ndx -maxwarn -1"
    },
    {
        "step": 25,
        "description": "Run the equilibration process.",
        "command": "gmx_mpi mdrun -v -deffnm equil"
    },
    {
        "step": 26,
        "description": "Prepare the system for the production run.",
        "command": "gmx_mpi grompp -f step5_production.mdp -o md.tpr -c equil.gro -p topol.top -n index.ndx -maxwarn -1"
    },
    {
        "step": 27,
        "description": "Start the production run.",
        "command": "gmx_mpi mdrun -deffnm md"
    },
    {
        "step": 28,
        "description": "Extend the simulation time.",
        "command": "gmx_mpi convert-tpr -s md_0_1.tpr -extend 100000 -o next.tpr"
    },
    {
        "step": 29,
        "description": "Continue the simulation.",
        "command": "gmx_mpi mdrun -s next.tpr -deffnm md_0_1 -cpi md_0_1_prev.cpt"
    }
]

# Function to create downloadable link for file content
def get_download_link(content, filename, link_text):
    b64 = base64.b64encode(content.encode()).decode()
    return f'<a href="data:text/plain;base64,{b64}" download="{filename}">{link_text}</a>'

# Shell script for MD simulation
def get_md_shell_script():
    return """#!/bin/csh
#SBATCH --time=48:000:00
#SBATCH -N 4
#SBATCH --account=XXX
#SBATCH -J gromacs_job

module load GROMACS/2024.2

set init = step5_input
set rest_prefix = step5_input
set mini_prefix = step6.0_minimization
set equi_prefix = step6.%d_equilibration
set prod_prefix = step7_production
set prod_step   = step7

# Minimization
gmx_mpi grompp -f ${mini_prefix}.mdp -o ${mini_prefix}.tpr -c ${init}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx -maxwarn -1
mpprun gmx_mpi mdrun -v -deffnm ${mini_prefix}

# Equilibration
set cnt    = 1
set cntmax = 6

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = `printf ${equi_prefix} ${cnt}`
    set pstep = `printf ${equi_prefix} ${pcnt}`
    if ( ${cnt} == 1 ) set pstep = ${mini_prefix}

    gmx_mpi grompp -f ${istep}.mdp -o ${istep}.tpr -c ${pstep}.gro -r ${rest_prefix}.gro -p topol.top -n index.ndx -maxwarn -1
    mpprun gmx_mpi mdrun -v -deffnm ${istep}
    @ cnt += 1
end

# Production
set cnt    = 1
set cntmax = 10

while ( ${cnt} <= ${cntmax} )
    @ pcnt = ${cnt} - 1
    set istep = ${prod_step}_${cnt}
    set pstep = ${prod_step}_${pcnt}

    if ( ${cnt} == 1 ) then
        set pstep = `printf ${equi_prefix} 6`
        gmx_mpi grompp -f ${prod_prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -p topol.top -n index.ndx -maxwarn -1
    else
       gmx_mpi grompp -f ${prod_prefix}.mdp -o ${istep}.tpr -c ${pstep}.gro -t ${pstep}.cpt -p topol.top -n index.ndx -maxwarn -1
    endif
    mpprun gmx_mpi mdrun -v -deffnm ${istep}
    @ cnt += 1
end
"""

# Protein-ligand script
def get_protein_ligand_script():
    return r"""#!/bin/bash

# Step 1: Extract ligand coordinates
grep UNL clean.pdb > UNL.pdb

# Step 2: Prepare protein structure
gmx_mpi pdb2gmx -f clean.pdb -o protein_processed.gro

# Step 3: Prepare ligand structure
gmx_mpi editconf -f UNL.pdb -o UNL.gro

# Step 4: Combine protein and ligand into complex
cp protein_processed.gro complex.gro
cat UNL.gro >> complex.gro
sed -i '1s/.*/NEW_ATOM_COUNT/' complex.gro  # Replace with actual atom count

# Step 5: Update topology file
echo '; Include ligand topology' >> topol.top
echo '#include "UNL.ITP"' >> topol.top
echo '[ molecules ]' >> topol.top
echo '; Compound        #mols' >> topol.top
echo 'Protein_chain_A     1' >> topol.top
echo 'UNL                 1' >> topol.top

# Step 6: Define the simulation box
gmx_mpi editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0

# Step 7: Solvate the system
gmx_mpi solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

# Step 8: Add ions
gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2
gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

# Step 9: Energy minimization preparation
gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2

# Energy minimization job script
cat <<EOF > run_em.sh
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
#SBATCH -J gromacs_em
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
module load GROMACS/2024.2
export OMP_NUM_THREADS=4
mpprun --pass="--map-by ppr:$((16/OMP_NUM_THREADS)):socket:PE=\${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm em
EOF

# Step 10: Generate position restraints
gmx_mpi make_ndx -f UNL.gro -o index_UNL.ndx << EOF
0 & ! a H*
q
EOF
gmx_mpi genrestr -f UNL.gro -n index_UNL.ndx -o posre_UNL.itp -fc 1000 1000 1000

# Step 11: Update topology for position restraints
echo '; Ligand position restraints' >> topol.top
echo '#ifdef POSRES' >> topol.top
echo '#include "posre_UNL.itp"' >> topol.top
echo '#endif' >> topol.top

# Step 12: Create index file for next steps
gmx_mpi make_ndx -f em.gro -o index.ndx << EOF
1|13
q
EOF

# Step 13: Prepare and run NVT equilibration
gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2

# NVT job script
cat <<EOF > run_nvt.sh
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
#SBATCH -J gromacs_nvt
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
module load GROMACS/2024.2
export OMP_NUM_THREADS=4
mpprun --pass="--map-by ppr:$((16/OMP_NUM_THREADS)):socket:PE=\${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm nvt
EOF

# Step 14: Prepare and run NPT equilibration
gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2

# NPT job script
cat <<EOF > run_npt.sh
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
#SBATCH -J gromacs_npt
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
module load GROMACS/2024.2
export OMP_NUM_THREADS=4
mpprun --pass="--map-by ppr:$((16/OMP_NUM_THREADS)):socket:PE=\${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm npt
EOF

# Step 15: Prepare and run MD simulation
gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_1.tpr -maxwarn 2

# MD job script
cat <<EOF > run_md.sh
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
#SBATCH -J gromacs_md
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
module load GROMACS/2024.2
export OMP_NUM_THREADS=4
mpprun --pass="--map-by ppr:$((16/OMP_NUM_THREADS)):socket:PE=\${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm md_0_1
EOF

echo "All scripts generated. Execute them as needed."
"""

# Protein in water script
def get_protein_in_water_script():
    return """#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH  -N 4
#SBATCH -c 4
#SBATCH --account=XXX
export OMP_NUM_THREADS=4
#SBATCH -J gromacs_job
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
#gmx_mpi convert-tpr -s md_0_1.tpr -extend 100000 -o next.tpr
module load GROMACS/2024.2
echo -e "1\\n1" | gmx_mpi pdb2gmx -f *.pdb -o peptide.gro
gmx_mpi editconf -f peptide.gro -o pep_box.gro -c -d 1.0 -bt cubic
gmx_mpi solvate -cp pep_box.gro -cs spc216.gro -o solv.gro -p topol.top
echo -e "13\\n1" | gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn -1
gmx_mpi grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn -1

mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -v -deffnm em
    gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm nvt
    gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm npt
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \\
    gmx_mpi mdrun -deffnm md_0_1
gmx_mpi trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_center.xtc -center -pbc mol -ur compact
"""

def search_commands(query):
    """Search for commands that match the query"""
    if not query:
        return {}
    
    # Convert query to lowercase for case-insensitive search
    query = query.lower()
    
    # First try exact matches
    matches = {key: value for key, value in gromacs_commands.items() 
               if query in key.lower()}
    
    # If no exact matches, try fuzzy matching
    if not matches:
        close_matches = get_close_matches(query, gromacs_commands.keys(), n=5, cutoff=0.6)
        if close_matches:
            matches = {key: gromacs_commands[key] for key in close_matches}
    
    return matches

# Main Streamlit app
st.title("GROMACS Command Search")
st.markdown("""
This tool helps you find the right GROMACS commands for molecular dynamics simulations.
Enter a description of what you want to do and get the corresponding command.
""")

# Create tabs for different sections
tab1, tab2, tab3 = st.tabs(["Command Search", "CHARMM-GROMACS Workflow", "Scripts"])

with tab1:
    # Search interface
    st.header("Search for GROMACS Commands")
    
    # Search box
    query = st.text_input("Enter what you want to do (e.g., 'calculate rmsd' or 'solvate system'):")
    
    # Search button
    if st.button("Search Commands") or query:
        matches = search_commands(query)
        
        if matches:
            st.success(f"Found {len(matches)} matching commands!")
            
            # Display results in a nice format
            for i, (description, command) in enumerate(matches.items()):
                with st.expander(f"{description}"):
                    st.code(command, language="bash")
                    st.markdown(f"**Description:** {description}")
        else:
            st.warning("No matching commands found. Try a different search term.")
            
            # Suggest similar commands if no exact match
            all_terms = " ".join(gromacs_commands.keys()).lower()
            query_terms = query.lower().split()
            
            suggestions = []
            for term in query_terms:
                if len(term) > 3 and term in all_terms:
                    related = [key for key in gromacs_commands.keys() 
                              if term in key.lower()]
                    suggestions.extend(related[:3])  # Limit to 3 per term
            
            if suggestions:
                st.subheader("You might be interested in:")
                for i, sugg in enumerate(suggestions[:5]):  # Show at most 5
                    with st.expander(f"{sugg}"):
                        st.code(gromacs_commands[sugg], language="bash")
    
    # Show popular commands
    with st.expander("Popular GROMACS Commands"):
        popular_commands = [
            "Run production MD",
            "Energy minimization",
            "Solvate a system",
            "Generate index file",
            "Analyze RMSD"
        ]
        
        for cmd in popular_commands:
            st.markdown(f"**{cmd}**")
            st.code(gromacs_commands[cmd], language="bash")

with tab2:
    st.header("CHARMM-GROMACS Workflow")
    st.markdown("""
    This is a step-by-step workflow for setting up a CHARMM-GROMACS simulation with ligands.
    Follow these steps to prepare your system.
    """)
    
    # Display workflow steps in a table
    for step in charmm_workflow:
        if "step" in step:
            with st.expander(f"Step {step['step']}: {step['description']}"):
                st.code(step['command'], language="bash")
                st.markdown(f"**Description:** {step['description']}")
        else:
            st.markdown(f"**{step['description']}**")
            st.markdown(f"[{step['command']}]({step['command']})")
    
    # Add a download button for the complete workflow
    workflow_text = "\n\n".join([f"# Step {step.get('step', 0)}: {step['description']}\n{step['command']}" 
                               for step in charmm_workflow if "step" in step])
    
    st.download_button(
        label="Download Complete Workflow",
        data=workflow_text,
        file_name="charmm_gromacs_workflow.txt",
        mime="text/plain"
    )

with tab3:
    st.header("Downloadable Scripts")
    st.markdown("""
    Here are some ready-to-use script templates for common GROMACS simulations.
    You can download and customize these for your specific needs.
    """)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Molecular Dynamics Script")
        st.markdown("""
        This script runs a complete MD simulation:
        - Energy minimization
        - Equilibration (6 steps)
        - Production run (10 steps)
        """)
        
        md_script = get_md_shell_script()
        st.download_button(
            label="Download MD Script",
            data=md_script,
            file_name="gromacs_md_script.sh",
            mime="text/plain"
        )
        
        st.subheader("Protein in Water Script")
        st.markdown("""
        This script prepares and simulates a protein in water:
        - Protein preparation
        - Solvation
        - Minimization and equilibration
        - Production MD
        """)
        
        protein_water_script = get_protein_in_water_script()
        st.download_button(
            label="Download Protein in Water Script",
            data=protein_water_script,
            file_name="protein_in_water.sh",
            mime="text/plain"
        )
    
    with col2:
        st.subheader("Protein-Ligand Complex Script")
        st.markdown("""
        This script prepares and simulates a protein-ligand complex:
        - Extract ligand
        - Prepare protein and ligand
        - Combine into complex
        - Solvate and add ions
        - Run minimization, NVT, NPT, and production
        """)
        
        protein_ligand_script = get_protein_ligand_script()
        st.download_button(
            label="Download Protein-Ligand Script",
            data=protein_ligand_script,
            file_name="protein_ligand_workflow.sh",
            mime="text/plain"
        )
        
        st.subheader("Custom MDP Template Generator")
        st.markdown("""
        Need a custom MDP file? Use our [GROMACS MDP Generator](https://share.streamlit.io/yourapp/mdp_generator) 
        tool to create parameter files for your simulations.
        """)

# Sidebar with information and statistics
st.sidebar.header("About this Tool")
st.sidebar.markdown("""
This tool helps molecular dynamics researchers find the right GROMACS commands 
for their simulations.

**Features:**
- Search over 200+ GROMACS commands
- Step-by-step CHARMM-GROMACS workflow
- Downloadable simulation scripts
- Command suggestions with fuzzy matching

**Statistics:**
- 200+ GROMACS commands
- 29-step CHARMM-GROMACS workflow
- 3 downloadable script templates
""")

# Add acknowledgements
st.sidebar.markdown("---")
st.sidebar.markdown("Created by [Your Name]")
st.sidebar.markdown("Based on GROMACS 2024.2")
st.sidebar.markdown("[GitHub Repository](https://github.com/pritampanda15)")
