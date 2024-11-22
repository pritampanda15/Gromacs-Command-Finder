from flask import Flask, render_template, request
from fuzzywuzzy import process

app = Flask(__name__)

# GROMACS Command Database
gromacs_commands = {
    "Build topology of protein": "gmx pdb2gmx -f protein.pdb -o processed.pdb -water spce",
    "Build topology of ligand": "gmx editconf -f ligand.pdb -o ligand_box.pdb -bt dodecahedron -d 1.0",
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
    "Check trajectory periodicity": "gmx check -f trajectory.xtc"
    # CHARMM-GROMACS Workflow (Steps)
    }
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
    result = None
    if request.method == "POST":
        user_query = request.form.get("query")
        best_match, confidence = process.extractOne(user_query, gromacs_commands.keys())
        if confidence > 60:  # Confidence threshold
            result = {"task": best_match, "command": gromacs_commands[best_match]}
        else:
            result = {"task": "No suitable match found.", "command": "Please refine your query."}
    return render_template("index.html", result=result, commands=gromacs_commands)

@app.route("/workflow")
def workflow():
    return render_template("workflow.html", workflow=charmm_workflow)



if __name__ == '__main__':
    app.run(debug=True)
