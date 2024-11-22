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
    "Add ions": "gmx genion -s ions.tpr -o solvated_ions.gro -p topol.top -pname NA -nname CL -neutral",
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
    "Radial distribution function of water": "gmx rdf -f trajectory.xtc -s reference.tpr -sel water.ndx -ref protein.ndx -o water_rdf.xvg"
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

module load GROMACS/2024.2-nsc1-gcc-9.3.0-bare

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
# Protein-Ligand Complex: Series of Commands to execute

grep UNL clean.pdb > UNL.pdb
gmx_mpi pdb2gmx -f clean.pdb -o protein_processed.gro
Generate Ligand topology by extracting the coordinates from the complex using ATB or PRODRG if you prefer GROMOS96 forcefield
gmx_mpi editconf -f UNL.pdb -o UNL.gro
Copy the protein_processed.gro to complex.gro and copy the UNL.gro coordinates to the end of complex.gro. Make sure to update the number of atoms at the beginning.
Open topol.top file and add 
; Include ligand topology
#include "UNL.ITP"
[ molecules ]
; Compound        #mols
Protein_chain_A     1
UNL                 1

gmx_mpi editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
gmx_mpi solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 2 (Name of the ligand should correspond to the ITP file)
gmx_mpi genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral
gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 2


SNIC gromacs run file for energy minimization
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH  -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
export OMP_NUM_THREADS=4
#SBATCH -J gromacs_job
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
module load GROMACS/2024.2 
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm em
gmx_mpi make_ndx -f UNL.gro -o index_UNL.ndx
0 & ! a H*
q
gmx_mpi genrestr -f UNL.gro -n index_UNL.ndx -o posre_UNL.itp -fc 1000 1000 1000
3

Then add this to the topol.top file

; Ligand position restraints
#ifdef POSRES
#include "posre_UNL.itp"
#endif

gmx_mpi make_ndx -f em.gro -o index.ndx
1|13
q
gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2

For NVT, NPT, MD
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH  -N 4
#SBATCH -c 4
#SBATCH --account=XXXX
export OMP_NUM_THREADS=4
#SBATCH -J gromacs_job
#SBATCH --exclusive
#SBATCH --ntasks-per-node=32
#gmx_mpi convert-tpr -s md_0_1.tpr -extend 100000 -o next.tpr
module load GROMACS/2019.2-nsc1-gcc-7.3.0-bare 
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm nvt
    gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm npt
    gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md_0_1.tpr -maxwarn 2
mpprun --pass="--map-by
ppr:$((16/OMP_NUM_THREADS)):socket:PE=${OMP_NUM_THREADS}" \
    gmx_mpi mdrun -deffnm md_0_1
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
