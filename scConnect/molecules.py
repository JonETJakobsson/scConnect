# Used to add molecules manually. Check that ligands are named
# correcly as in the Ligands table
# Current version expect mouse genes, and automatically find orthogonal 
# genes when a species specific library is constructed. These files can be edited manually if needed


def get_molecules():
    molecules = [{
        "ligand": "L-glutamic acid",
        "type": "molecule",
        "synthesis": None,
        "transport": ["Slc17a6", "Slc17a7", "Slc17a8"],
        "reuptake": ["Slc1a1", "Slc1a2", "Slc1a6"],
        "excluded": None
    }, {
        "ligand": "GABA",
        "type": "molecule",
        "synthesis": ["Gad1", "Gad2"],
        "transport": ["Slc32a1"],
        "reuptake": ["Slc6a1"],
        "excluded": None
    }, {
        "ligand": "glycine",
        "type": "molecule",
        "synthesis": None,
        "transport": ["Slc32a1"],
        "reuptake": ["Slc6a9"],
        "excluded": None
    }, {
        "ligand": "dopamine",
        "type": "molecule",
        "synthesis": ["Th", "Ddc"],
        "transport": ["Slc18a2", "Slc18a3"],
        "reuptake": ["Slc6a3"],
        "excluded": ["Dbh"]
    }, {
        "ligand": "(-)-noradrenaline",
        "type": "molecule",
        "synthesis": ["Th", "Ddc", "Dbh"],
        "transport": ["Slc18a2", "Slc18a3"],
        "reuptake": ["Slc6a3"],
        "excluded": ["Pnmt"]
    }, {
        "ligand": "(-)-adrenaline",
        "type": "molecule",
        "synthesis": ["Th", "Ddc", "Dbh", "Pnmt"],
        "transport": ["Slc18a2"],
        "reuptake": None,
        "excluded": None
    }, {
        "ligand": "5-hydroxytryptamine",
        "type": "molecule",
        "synthesis": ["Tph2", "Ddc"],
        "transport": ["Slc18a2"],
        "reuptake": ["Slc6a4"],
        "excluded": None
    }, {
        "ligand": "acetylcholine",
        "type": "molecule",
        "synthesis": ["Chat"],
        "transport": ["Slc18a3"],
        "reuptake": ["Ache", "Slc5a7"],
        "excluded": None
    }]
    return molecules
