from flask import Flask, render_template, request, send_file, redirect, url_for, session
from fuzzywuzzy import process

app = Flask(__name__)

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
    "Extend MD simulation": "gmx convert-tpr -s md_0_1.tpr -extend 100000 -o next.tpr",
    "Calculate bond angle distribution": "gmx angle -f trajectory.xtc -type angle -o angle_distribution.xvg",
    "Analyze interfacial water properties": "gmx sasa -s reference.tpr -f trajectory.xtc -n index.ndx -surface -o interfacial_water.xvg",
    "Force distribution analysis": "gmx traj -s reference.tpr -f trajectory.xtc -force -o force_distribution.xvg",
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
    # Specific Analysis Commands
    "Compute salt bridges": "gmx saltbr -f traj.xtc -s topol.tpr -n index.ndx -num salt_bridges.xvg",
    "Solvent accessible surface area (SASA)": "gmx sasa -f traj.xtc -s topol.tpr -o sasa.xvg -surface 1 -output 2",
    "Compute small-angle X-ray scattering (SAXS)": "gmx saxs-legacy -f traj.xtc -s topol.tpr -q q.xvg -o saxs_output.xvg",
    "Generate Ramachandran plots": "gmx rama -f traj.xtc -s topol.tpr -o ramachandran.xvg",
    "Principal axes of inertia": "gmx principal -f traj.xtc -s topol.tpr -o principal_axes.xvg",
    "Van Hove displacement": "gmx vanhove -f traj.xtc -s topol.tpr -n index.ndx -o vh_displacement.xvg",
    "Frequency-filter trajectories": "gmx filter -f traj.xtc -s topol.tpr -o filtered.xtc",
    "Rotational correlation function": "gmx rotacf -f traj.xtc -s topol.tpr -n index.ndx -o rotacf.xvg",
    "Velocity autocorrelation functions": "gmx velacc -f traj.xtc -s topol.tpr -o velacf.xvg",
    
    # Structural Modifications
    "Generate index groups": "gmx make_ndx -f structure.gro -o index.ndx",
    "Insert molecules into a box": "gmx insert-molecules -f box.gro -ci molecule.gro -nmol 100",
    "Add ions to neutralize system": "gmx genion -s ions.tpr -p topol.top -neutral -conc 0.15 -o neutralized.gro",
    "Generate position restraints": "gmx genrestr -f structure.pdb -o posres.itp",
    "Convert files for specific force fields": "gmx pdb2gmx -f molecule.pdb -o topol.gro -ff charmm36",
    
    # Advanced Molecular Dynamics
    "Run production MD simulation": "gmx mdrun -v -deffnm production_md",
    "Preprocess files for simulation": "gmx grompp -f md.mdp -c solvated.gro -p topol.top -o production.tpr",
    "Pulling simulation setup": "gmx grompp -f pull.mdp -c start.gro -p topol.top -o pull.tpr",
    "Weighted histogram analysis method (WHAM)": "gmx wham -f pullf.xvg -o free_energy_profile.xvg",
    
    # Energy and Free Energy Analysis
    "Energy analysis from trajectory": "gmx energy -f energy.edr -o energy.xvg",
    "Calculate electrostatic potential": "gmx potential -f traj.xtc -s topol.tpr -o potential.xvg",
    "Free energy calculation using BAR": "gmx bar -f bar.xvg -o free_energy_bar.xvg",
    "Calculate interaction energies": "gmx energy -f energy.edr -o interaction_energies.xvg",
    
    # Special Configurations
    "Hydrogen bond analysis": "gmx hbond -f traj.xtc -s topol.tpr -num hbonds.xvg -dist hbond_dist.xvg",
    "Compute mean squared displacement": "gmx msd -f traj.xtc -s topol.tpr -n index.ndx -o msd.xvg",
    "Protein-ligand binding analysis": "gmx rdf -f traj.xtc -s topol.tpr -n index.ndx -o rdf_binding.xvg",
    "Cluster analysis by RMSD": "gmx cluster -f traj.xtc -s topol.tpr -n index.ndx -method gromos -cutoff 0.25 -cl clusters.pdb",
    
    # Visualization Preparation
    "Convert trajectory for visualization": "gmx trjconv -f traj.xtc -s topol.tpr -o trajectory_for_vmd.xtc",
    "Generate contact maps": "gmx mdmat -f traj.xtc -s topol.tpr -o contact_map.xpm",
    "Generate helical wheel plots": "gmx wheel -f structure.pdb -o helical_wheel.png",
    
    # Debugging and Testing
    "Check TPR file": "gmx check -s topol.tpr",
    "Verify trajectory integrity": "gmx trjconv -f traj.xtc -o verified_traj.xtc",
    "Estimate PME error": "gmx pme_error -s topol.tpr -n index.ndx",
    "Benchmark non-bonded interactions": "gmx nonbonded-benchmark -s topol.tpr -b 10000 -e 50000",
    
    # Water and Solvent Analysis
    "Water orientation analysis": "gmx sorient -f traj.xtc -s topol.tpr -o water_orient.xvg",
    "Compute hydration shell dynamics": "gmx hydorder -f traj.xtc -s topol.tpr -n index.ndx",
    "Tetrahedrality of water": "gmx h2order -f traj.xtc -s topol.tpr",
    "Calculate radial distribution functions": "gmx rdf -f traj.xtc -s topol.tpr -o rdf_water.xvg",
    
    # Miscellaneous
    "Generate input files for normal mode analysis": "gmx make_edi -f traj.xtc -o eigenvectors.edi",
    "Align trajectory to a reference": "gmx trjconv -f traj.xtc -s ref.pdb -fit rot+trans -o aligned_traj.xtc",
    "Generate electrostatic maps": "gmx potential -f traj.xtc -s topol.tpr -o electrostatic_map.xpm",
    "Compute dielectric constants": "gmx dielectric -f traj.xtc -s topol.tpr -o dielectric.xvg",

    # Additional Simulation Setup Commands
    "Extend a simulation": "gmx convert-tpr -s topol.tpr -extend 50000 -o extended.tpr",
    "Restrained MD simulation": "gmx grompp -f md_restrained.mdp -c input.gro -r input.gro -p topol.top -o restrained.tpr",
    "Position restraint energy minimization": "gmx grompp -f em_restrained.mdp -c solvated.gro -p topol.top -r solvated.gro -o em_restrained.tpr",
    "Temperature coupling group setup": "gmx grompp -f nvt.mdp -c minimized.gro -p topol.top -o nvt.tpr",
    "Pressure coupling group setup": "gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr",
    "Equilibration phase setup": "gmx grompp -f equilibration.mdp -c npt.gro -p topol.top -o equilibration.tpr",
    "Set up periodic boundary conditions": "gmx editconf -f input.pdb -o box.pdb -bt triclinic -d 1.0",
    
    # Extended Trajectory Analysis
    "Calculate hydrogen bond lifetimes": "gmx hbond -f traj.xtc -s topol.tpr -n index.ndx -life hbond_lifetime.xvg",
    "Compute dihedral transitions": "gmx angle -f traj.xtc -s topol.tpr -type dihedral -o dihedral_transitions.xvg",
    "Identify water bridges": "gmx hbond -f traj.xtc -s topol.tpr -n water.ndx -map water_bridge_map.xpm",
    "Contact area analysis": "gmx sasa -f traj.xtc -s topol.tpr -n contact_area.ndx -o contact_area.xvg",
    "Identify cavities in the structure": "gmx freevolume -f traj.xtc -s topol.tpr -n cavity.ndx",
    "Molecular orientation in membrane": "gmx gangle -f traj.xtc -s topol.tpr -n membrane.ndx -o molecular_orientation.xvg",
    "Force distribution in groups": "gmx traj -f traj.xtc -s topol.tpr -n index.ndx -of force_distribution.xvg",
    "Binding site clustering": "gmx cluster -f traj.xtc -s topol.tpr -n binding_site.ndx -cl binding_clusters.pdb",
    "Root mean square fluctuation (RMSF)": "gmx rmsf -f traj.xtc -s topol.tpr -o rmsf_residues.xvg",
    
    # Advanced Energy Calculations
    "Analyze energy components": "gmx energy -f energy.edr -o energy_components.xvg",
    "Compute Lennard-Jones potential": "gmx potential -f traj.xtc -s topol.tpr -lj lj_potential.xvg",
    "Pressure profile analysis": "gmx density -f traj.xtc -s topol.tpr -sl 50 -o pressure_profile.xvg",
    "Pair interaction energy calculations": "gmx mdrun -rerun traj.xtc -deffnm pairwise_interactions",
    "Calculate interaction free energy": "gmx mmpbsa -f traj.xtc -s topol.tpr -n index.ndx -o binding_free_energy.dat",
    
    # Membrane and Solvent-Specific Analysis
    "Compute membrane thickness": "gmx gangle -f traj.xtc -s topol.tpr -group1 membrane.ndx -group2 solvent.ndx -o thickness.xvg",
    "Analyze lipid tail order": "gmx order -f traj.xtc -s topol.tpr -n lipid_tails.ndx -o tail_order.xvg",
    "Water dipole orientation": "gmx dipoles -f traj.xtc -s topol.tpr -n water.ndx -o water_dipole.xvg",
    "Membrane protein tilt": "gmx helixorient -f traj.xtc -s topol.tpr -o protein_tilt.xvg",
    "2D radial distribution function": "gmx densmap -f traj.xtc -s topol.tpr -od radial_2d.xpm",
    "Solvent exposure of residues": "gmx sasa -f traj.xtc -s topol.tpr -res -o solvent_exposure_residues.xvg",
    
    # Umbrella Sampling and Free Energy
    "Set up umbrella sampling windows": "gmx grompp -f umbrella.mdp -c window.gro -p topol.top -o umbrella.tpr",
    "Pulling along a reaction coordinate": "gmx mdrun -s pull.tpr -deffnm pulling -pf pullf.xvg -px pullx.xvg",
    "Weighted histogram analysis for US": "gmx wham -it tpr-files.dat -if pullf-files.dat -o pmf.xvg -hist histo.xvg",
    "Analyze free energy landscapes": "gmx sham -f histogram.xvg -o free_energy_landscape.xpm",
    
    # Parallelization and Benchmarking
    "Run simulation on multiple GPUs": "gmx mdrun -s topol.tpr -deffnm output -ntmpi 4 -ntomp 8 -gpu_id 0,1",
    "Benchmark simulation performance": "gmx mdrun -s topol.tpr -ntmpi 4 -ntomp 8 -resethway -nsteps 100000",
    "Optimize PME grid": "gmx tune_pme -s topol.tpr -o optimized_pme.tpr -np 16",
    
    # Visualization and Data Export
    "Export data to CSV": "gmx analyze -f data.xvg -o data.csv",
    "Generate trajectory movies": "gmx trjconv -f traj.xtc -s topol.tpr -o movie.pdb",
    "Create secondary structure timelines": "gmx dssp -f traj.xtc -s topol.tpr -o ss.xpm",
    "Generate XPM plots": "gmx xpm2ps -f plot.xpm -o plot.ps",
    "Density maps for visualization": "gmx densmap -f traj.xtc -s topol.tpr -o density_map.xpm",
    
    # Miscellaneous Utility Commands
    "Convert XPM to image format": "gmx xpm2ps -f map.xpm -o map_image.eps",
    "Extract subgroups of atoms": "gmx select -f traj.xtc -s topol.tpr -n index.ndx -select 'resname SOL' -on water_group.ndx",
    "Concatenate multiple ED frames": "gmx nmens -s topol.tpr -v eigenvec.trr -n index.ndx -o ed_ensemble.xtc",
    "Perform vacuum simulations": "gmx grompp -f vacuum.mdp -c minimized.gro -p topol.top -o vacuum.tpr",
    "Extract PME grid settings": "gmx pme_error -s topol.tpr -o pme_error_estimate.dat",

    # Coarse-Grained (CG) Simulations
    "Prepare MARTINI coarse-grained system": "gmx pdb2gmx -f input.pdb -o processed.gro -ff martini",
    "Map all-atom structure to CG beads": "gmx traj -f all_atom.xtc -s all_atom.tpr -n index.ndx -o mapped_cg.xtc",
    "Analyze CG bead cluster size": "gmx clustsize -f traj.xtc -s topol.tpr -n beads.ndx -nc cluster_counts.xvg",
    "Lipid order parameter analysis (CG)": "gmx order -f traj.xtc -s topol.tpr -n lipid_tails.ndx -o lipid_order.xvg",
    "Generate CG index file": "gmx make_ndx -f cg_structure.gro -o cg_index.ndx",
    
    # Nucleic Acids Analysis
    "Calculate DNA/RNA base-pair distance": "gmx distance -f traj.xtc -s topol.tpr -n basepairs.ndx -oav basepair_distances.xvg",
    "Helical rise and twist of DNA": "gmx helix -f traj.xtc -s dna_topol.tpr -n dna.ndx -o dna_helix.xvg",
    "Nucleic acid hydration shell analysis": "gmx hbond -f traj.xtc -s topol.tpr -n dna_solvent.ndx -num hydration_shell.xvg",
    "Analyze grooves in DNA": "gmx rama -f traj.xtc -s dna_topol.tpr -n grooves.ndx -o dna_grooves.xvg",
    "Backbone torsion angles": "gmx angle -f traj.xtc -s topol.tpr -n backbone.ndx -type dihedral -od torsion_angles.xvg",
    
    # Polymer and Chain Analysis
    "End-to-end distance of polymers": "gmx distance -f traj.xtc -s topol.tpr -n polymer.ndx -oav end_to_end_distance.xvg",
    "Radius of gyration of polymer chains": "gmx gyrate -f traj.xtc -s topol.tpr -o polymer_radius.xvg",
    "Segment-specific RMSF": "gmx rmsf -f traj.xtc -s topol.tpr -n polymer_segments.ndx -o segment_rmsf.xvg",
    "Polymer entanglement analysis": "gmx polystat -f traj.xtc -s topol.tpr -o entanglement_data.xvg",
    "Mean squared displacement of chains": "gmx msd -f traj.xtc -s topol.tpr -n chains.ndx -o msd_chains.xvg",
    
    # Quantum Mechanics/Molecular Mechanics (QM/MM) Analysis
    "Set up QM/MM simulation": "gmx grompp -f qmmm.mdp -c minimized.gro -p topol.top -o qmmm.tpr",
    "Run QM/MM simulation": "gmx mdrun -s qmmm.tpr -deffnm qmmm_run",
    "Analyze QM/MM interaction energies": "gmx energy -f qmmm_energy.edr -o qmmm_interaction_energies.xvg",
    "Compute QM/MM dipole moments": "gmx dipoles -f qmmm_traj.xtc -s qmmm_topol.tpr -o qmmm_dipoles.xvg",
    
    # Free Energy Perturbation (FEP)
    "Set up FEP simulation": "gmx grompp -f fep.mdp -c fep.gro -p topol.top -o fep.tpr",
    "Analyze FEP energy differences": "gmx bar -f fep_results.xvg -o fep_free_energy.xvg",
    "Weighted histogram analysis for FEP": "gmx wham -it fep_tpr_files.dat -if fep_pmf_files.dat -o fep_profile.xvg",
    
    # Advanced Force Field Setup
    "Build topology for custom residues": "gmx pdb2gmx -f custom_residues.pdb -o processed_custom.gro -ff custom_ff",
    "Modify bonded parameters": "gmx x2top -f structure.pdb -o new_topol.top -ff modified",
    "Convert C6/C12 to sigma/epsilon": "gmx sigeps -f nonbonded.itp -o sigma_epsilon.itp",
    "Generate topologies for virtual sites": "gmx grompp -f virtual_sites.mdp -c input.gro -p topol.top -o virtual.tpr",
    
    # Enhanced Sampling Techniques
    "Set up simulated annealing": "gmx grompp -f annealing.mdp -c minimized.gro -p topol.top -o annealing.tpr",
    "Run replica exchange simulation": "gmx mdrun -multidir rep1 rep2 rep3 -deffnm replica_exchange",
    "Analyze replica transitions": "gmx wham -f replica_transition.xvg -o free_energy_landscape.xvg",
    "Adaptive biasing force (ABF) setup": "gmx grompp -f abf.mdp -c abf_start.gro -p topol.top -o abf.tpr",
    
    # Ligand and Binding Site Analysis
    "Ligand binding energy calculation": "gmx mmpbsa -f traj.xtc -s topol.tpr -n index.ndx -o binding_energy.dat",
    "Residue-level SASA for binding site": "gmx sasa -f traj.xtc -s topol.tpr -n binding_residues.ndx -res -o binding_site_sasa.xvg",
    "Distance between ligand and site": "gmx distance -f traj.xtc -s topol.tpr -n ligand_site.ndx -oav ligand_site_distance.xvg",
    "Cluster binding poses of ligand": "gmx cluster -f ligand_traj.xtc -s topol.tpr -n ligand_index.ndx -cl binding_clusters.pdb",
    "Generate docking box": "gmx editconf -f ligand.pdb -o docking_box.pdb -c -bt triclinic -d 1.2",
    
    # Rare Utilities
    "Convert trajectory to AMBER format": "gmx trjconv -f traj.xtc -s topol.tpr -o amber_traj.mdcrd",
    "Prepare CHARMM force field inputs": "gmx pdb2gmx -f structure.pdb -ff charmm36 -o charmm_structure.gro",
    "Calculate free volume": "gmx freevolume -f traj.xtc -s topol.tpr -n cavity.ndx -o free_volume.xvg",
    "Surface tension analysis": "gmx energy -f energy.edr -surf_ten surface_tension.xvg",
    "Visualize electric field around system": "gmx potential -f traj.xtc -s topol.tpr -o electric_field.xpm",
    
    # Performance Optimization
    "Run simulation with GPU acceleration": "gmx mdrun -s topol.tpr -deffnm gpu_run -gpu_id 0,1",
    "Optimize PME load balancing": "gmx tune_pme -s topol.tpr -np 8 -o optimized_pme.tpr",
    "Benchmark custom hardware": "gmx mdrun -nsteps 50000 -resethway -ntmpi 4 -ntomp 8",

     # Advanced Structural Analysis
    "Calculate per-residue RMSF": "gmx rmsf -f traj.xtc -s topol.tpr -res -o rmsf_per_residue.xvg",
    "Analyze beta-sheet interactions": "gmx rama -f traj.xtc -s beta_sheets.tpr -o beta_interactions.xvg",
    "Atomic contact probability (ACP)": "gmx mdmat -f traj.xtc -s topol.tpr -frames -o acp_map.xpm",
    "Time-lagged independent component analysis (tICA)": "gmx anaeig -f traj.xtc -s topol.tpr -type tica -o tica.xvg",
    "Protein unfolding pathway analysis": "gmx rms -f unfolding.xtc -s folded.pdb -fit rot+trans -o unfolding_rmsd.xvg",
    
    # Rare Nucleic Acid-Specific Analysis
    "DNA/RNA hydration dynamics": "gmx sorient -f traj.xtc -s dna_topol.tpr -n dna_solvent.ndx -o hydration_dynamics.xvg",
    "Stacking interaction energy in DNA": "gmx energy -f energy.edr -terms stacking -o stacking_energy.xvg",
    "Quantify DNA bending angles": "gmx angle -f traj.xtc -s dna_topol.tpr -type angle -o bending_angles.xvg",
    "Groove width in DNA double helix": "gmx helix -f traj.xtc -s dna_topol.tpr -o groove_width.xvg",
    "RNA base flipping events": "gmx angle -f traj.xtc -s rna_topol.tpr -type dihedral -o base_flipping.xvg",
    
    # Advanced Lipid Membrane Analysis
    "Lipid diffusion coefficients": "gmx msd -f traj.xtc -s topol.tpr -n lipid_heads.ndx -o lipid_diffusion.xvg",
    "Hydrophobic core thickness": "gmx density -f traj.xtc -s topol.tpr -o hydrophobic_thickness.xvg",
    "Cholesterol tilt angles in bilayer": "gmx gangle -f traj.xtc -s topol.tpr -group1 chol.ndx -o chol_tilt.xvg",
    "Lipid-lipid contact analysis": "gmx mindist -f traj.xtc -s topol.tpr -n lipid_pairs.ndx -o lipid_contacts.xvg",
    "Membrane curvature": "gmx spatial -f traj.xtc -s topol.tpr -o membrane_curvature.xpm",
    
    # Rare Free Energy Methods
    "Calculate thermodynamic integration (TI)": "gmx mdrun -s ti.tpr -deffnm ti_run",
    "Energy reweighting for free energy": "gmx bar -f energy_files.dat -o reweighted_free_energy.xvg",
    "Binding affinity using MM/PBSA": "gmx mmpbsa -f traj.xtc -s complex.tpr -o binding_affinity.dat",
    "Plot potential of mean force (PMF)": "gmx wham -it pmf_tpr_files.dat -if pmf_pullf_files.dat -o pmf_curve.xvg",
    "Umbrella sampling convergence check": "gmx sham -f histogram.xvg -o convergence_check.xpm",
    
    # Non-Standard System Simulations
    "Metal ion parameterization": "gmx pdb2gmx -f ions.pdb -ff custom_ff -o metal_params.gro",
    "Coarse-grained polymer chain analysis": "gmx msd -f cg_polymer.xtc -s cg_topol.tpr -o cg_polymer_diffusion.xvg",
    "Non-aqueous solvent simulations": "gmx solvate -cp processed.gro -cs non_aqueous.gro -o solvated.gro",
    "Ionic liquid diffusion": "gmx msd -f ionic_liquid.xtc -s ionic_topol.tpr -n index.ndx -o ionic_diffusion.xvg",
    "Simulations with mixed solvents": "gmx solvate -cp processed.gro -cs mixed_solvent.gro -o solvated_mixture.gro",
    
    # Advanced Visualization and Export
    "Generate trajectory in PDB format": "gmx trjconv -f traj.xtc -s topol.tpr -o traj.pdb",
    "Write energy to JSON": "gmx energy -f energy.edr -o energy.json",
    "Generate pore structure visualization": "gmx densmap -f traj.xtc -s topol.tpr -o pore_density.xpm",
    "Visualize hydrogen bond lifetimes": "gmx hbond -f traj.xtc -s topol.tpr -hbn hbond_network.ndx -life hbond_lifetime.xvg",
    "Plot order parameters in 3D": "gmx order -f traj.xtc -s topol.tpr -o 3d_order_map.xpm",
    
    # Highly Specialized Commands
    "Estimate Lennard-Jones cutoffs": "gmx tune_lj -s topol.tpr -o tuned_cutoff.xvg",
    "Analyze solvation free energy": "gmx mdrun -s solvation.tpr -deffnm solvation_free_energy",
    "Kinetic energy spectrum analysis": "gmx velacc -f traj.xtc -s topol.tpr -o kinetic_energy_spectrum.xvg",
    "Extract water-mediated contacts": "gmx water -f traj.xtc -s topol.tpr -o water_contacts.xvg",
    "Topological network analysis of residues": "gmx rmsdist -f traj.xtc -s topol.tpr -type contact -o residue_network.xvg",
    
    # QM/MM and Hybrid Models
    "Transition state search in QM/MM": "gmx mdrun -s qmmm_tstate.tpr -deffnm ts_search",
    "Reaction pathway analysis": "gmx traj -f qmmm_traj.xtc -s qmmm_topol.tpr -o reaction_pathway.xvg",
    "Calculate QM/MM polarization energies": "gmx energy -f qmmm_energy.edr -terms polarization -o polarization_energy.xvg",
    "Orbital overlap analysis": "gmx traj -f qmmm_traj.xtc -s qmmm_topol.tpr -o orbital_overlap.xvg",
    
    # Extremely Rare Utilities
    "Detect water chain formation": "gmx hbond -f traj.xtc -s topol.tpr -n water.ndx -life water_chains.xvg",
    "Radial distribution for heterogeneous systems": "gmx rdf -f traj.xtc -s topol.tpr -rdf heterogeneous -o heterogeneous_rdf.xvg",
    "Surface energy analysis": "gmx energy -f energy.edr -surf_surf_surface_energy.xvg",
    "Simulate nano-particle diffusion": "gmx msd -f nanoparticle.xtc -s nanoparticle.tpr -o nano_diffusion.xvg",
    "Calculate electric dipole moments": "gmx dipoles -f traj.xtc -s topol.tpr -o electric_dipole.xvg",
    
    # Ultra-Performance Commands
    "Run with mixed CPU and GPU": "gmx mdrun -s topol.tpr -gpu_id 0 -ntmpi 4 -ntomp 8",
    "Benchmark energy minimization speed": "gmx mdrun -s em.tpr -nsteps 1000 -resethway",
    "Tunable cutoff benchmarks": "gmx tune_cutoff -s topol.tpr -o optimal_cutoff.xvg",

    # Alchemical transformations for ligand binding free energy
    "Set up alchemical transformations": "gmx grompp -f alchemical.mdp -c solvated.gro -p topol.top -o alchemical.tpr",
    "Run transformation for free energy": "gmx mdrun -s alchemical.tpr -deffnm alchemical_run",
    "BAR for alchemical transformations": "gmx bar -f alchemical_results.xvg -o free_energy.xvg",

    # Free energy of solvation
    "Set up free energy of solvation": "gmx grompp -f solvation.mdp -c solvated.gro -p topol.top -o solvation.tpr",
    "Analyze solvation free energy": "gmx bar -f solvation_results.xvg -o solvation_free_energy.xvg",

    # Protein-ligand binding affinities
    "Set up umbrella sampling for binding": "gmx grompp -f umbrella.mdp -c binding.gro -p topol.top -o umbrella.tpr",
    "Weighted histogram analysis": "gmx wham -it tpr_files.dat -if pullf_files.dat -o pmf.xvg",

    # Membrane asymmetry
    "Analyze membrane asymmetry": "gmx density -f traj.xtc -s topol.tpr -n lipids.ndx -o asymmetry_density.xvg",

    # Cholesterol-lipid interactions
    "Measure cholesterol tilt": "gmx gangle -f traj.xtc -s topol.tpr -group1 chol.ndx -o chol_tilt.xvg",

    # Pore formation and lipid flip-flop
    "Detect lipid flip-flop events": "gmx mindist -f traj.xtc -s topol.tpr -n lipid_heads.ndx -o flip_flop.xvg",
    "Visualize pore formation": "gmx densmap -f traj.xtc -s topol.tpr -o pore_density.xpm",

    # Polymer entanglement
    "Measure polymer entanglement": "gmx polystat -f traj.xtc -s topol.tpr -o entanglement_analysis.xvg",

    # Stress-strain simulations
    "Simulate polymer under stress": "gmx mdrun -s stress_sim.tpr -deffnm stress_sim",
    "Analyze strain distribution": "gmx traj -f traj.xtc -s topol.tpr -n index.ndx -o strain_distribution.xvg",

    # Diffusion in amorphous solids
    "Calculate diffusion coefficients": "gmx msd -f traj.xtc -s topol.tpr -o diffusion.xvg",

    # RNA folding pathways
    "Analyze RNA folding dynamics": "gmx rms -f traj.xtc -s folded.pdb -o folding_pathway.xvg",

    # DNA mismatch dynamics
    "Measure mismatch repair events": "gmx angle -f traj.xtc -s dna_topol.tpr -type dihedral -o mismatch_angles.xvg",

    # RNA-protein binding
    "Cluster RNA-protein binding poses": "gmx cluster -f traj.xtc -s topol.tpr -n binding.ndx -cl binding_clusters.pdb",

    # Enzymatic reaction pathways
    "Set up QM/MM enzymatic reaction": "gmx grompp -f qmmm.mdp -c enzyme.gro -p topol.top -o qmmm.tpr",
    "Analyze reaction pathways": "gmx traj -f traj.xtc -s qmmm_topol.tpr -o reaction_pathway.xvg",

    # Transition state optimization
    "Search for transition states": "gmx mdrun -s qmmm_ts.tpr -deffnm ts_search",

    # Replica exchange molecular dynamics (REMD)
    "Set up REMD": "gmx grompp -f remd.mdp -c replicas.gro -p topol.top -o remd.tpr",
    "Run REMD simulation": "gmx mdrun -multidir dir1 dir2 dir3 -deffnm remd",

    # Metadynamics
    "Set up metadynamics": "gmx grompp -f meta.mdp -c input.gro -p topol.top -o metadynamics.tpr",

    # Ionic liquids
    "Analyze ionic liquid diffusion": "gmx msd -f ionic_liquid.xtc -s ionic_topol.tpr -o diffusion.xvg",

    # Hydrophobic clustering
    "Measure hydrophobic clustering": "gmx cluster -f traj.xtc -s topol.tpr -n hydrophobic.ndx -o clusters.pdb",


    # Protein unfolding
    "Simulate protein unfolding": "gmx mdrun -s unfolding.tpr -deffnm unfolding",

    #Buckling in structural biopolymers
    "Simulate biopolymer buckling": "gmx mdrun -s buckling.tpr -deffnm buckling",

    # PME optimization
    "Tune PME grid": "gmx tune_pme -s topol.tpr -o tuned_pme.tpr",

    # Hybrid CPU-GPU performance
    "Benchmark CPU-GPU simulation": "gmx mdrun -s topol.tpr -gpu_id 0 -ntmpi 4 -ntomp 8",

    # Non-standard molecule parameterization
    "Generate topology for custom residues": "gmx pdb2gmx -f custom.pdb -ff amber99sb -o custom_topol.gro",

    # Combine multiple force fields
    "Integrate CHARMM and MARTINI": "gmx grompp -f combined.mdp -c charmm_martini.gro -p combined.top -o combined.tpr",

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


@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "POST":
        user_query = request.form.get("query").lower()
        matches = {key: value for key, value in gromacs_commands.items() if user_query in key.lower()}
        
        # Store matches or error in the session (or temporary flash messages)
        if matches:
            result = {"matches": matches}
        else:
            result = {"error": "No suitable matches found. Please refine your query."}
        
        return render_template("index.html", result=result, commands=gromacs_commands)
    
    # On GET request, reset result
    return render_template("index.html", result=None, commands=gromacs_commands)

@app.after_request
def add_no_cache_headers(response):
    # Prevent browser from caching the page
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, post-check=0, pre-check=0, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "-1"
    return response

@app.route("/workflow")
def workflow():
    return render_template("workflow.html", workflow=charmm_workflow)

@app.route("/download_script", methods=["GET"])
def download_script():
    # Shell script content
    shell_script_content = """#!/bin/csh
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

    script_file = "/tmp/gromacs_script.sh"
    with open(script_file, "w") as file:
        file.write(shell_script_content)
    
    return send_file(script_file, as_attachment=True, download_name="gromacs_script.sh")

@app.route("/download_protein_ligand_script", methods=["GET"])
def download_protein_ligand_script():
    # Protein-Ligand Workflow Script Content
    script_content = """
#!/bin/bash

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

    # Write script content to a temporary file
    script_file = "/tmp/protein_ligand_workflow.sh"
    with open(script_file, "w") as file:
        file.write(script_content)

    # Serve the file for download
    return send_file(script_file, as_attachment=True, download_name="protein_ligand_workflow.sh")


@app.route("/download_protein_in_water_script", methods=["GET"])
def download_protein_in_water_script():
    # Protein-in-Water Workflow Script Content
    script_content = """
#!/bin/bash
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
echo -e "1\n1" | gmx_mpi pdb2gmx -f *.pdb -o peptide.gro
gmx_mpi editconf -f peptide.gro -o pep_box.gro -c -d 1.0 -bt cubic
gmx_mpi solvate -cp pep_box.gro -cs spc216.gro -o solv.gro -p topol.top
echo -e "13\n1" | gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn -1
gmx_mpi grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn -1

mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -v -deffnm em
    gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm nvt
    gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm npt
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr -maxwarn -1
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm md_0_1
gmx_mpi trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_center.xtc -center -pbc mol -ur compact
"""

    # Write script content to a temporary file
    script_file = "/tmp/protein_in_water_workflow.sh"
    with open(script_file, "w") as file:
        file.write(script_content)

    # Serve the file for download
    return send_file(script_file, as_attachment=True, download_name="protein_in_water_workflow.sh")


if __name__ == '__main__':
    app.run(debug=True)
