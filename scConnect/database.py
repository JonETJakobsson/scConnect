
version = "2020-5"  # this is currently imported in other modules, so set verions here
organism = "mmusculus"

# Clean original GTP tables-----------------------------------------------------------------------------------------


def import_GTP_tables():
    """Imports all GTP tables into a dictionary of tables."""

    import pandas as pd
    import pkg_resources
    print(f"fetching GTF version {version}")
    # pkg_resourses can be used to access the local files in the package, independent of where the package is installed.
    i_path = pkg_resources.resource_filename(
        __name__, f"data/GTP_tables/{version}/interactions.csv")
    l_path = pkg_resources.resource_filename(
        __name__, f"data/GTP_tables/{version}/ligands.csv")
    p_path = pkg_resources.resource_filename(
        __name__, f"data/GTP_tables/{version}/peptides.csv")
    t_path = pkg_resources.resource_filename(
        __name__, f"data/GTP_tables/{version}/targets_and_families.csv")

    # imports full tables from Guidtopharmacology database
    interactions = pd.read_csv(i_path)
    interactions.index.name = "interaction_id"
    ligands = pd.read_csv(l_path, index_col=0)
    peptides = pd.read_csv(p_path, index_col=0)
    targets = pd.read_csv(t_path, index_col=3)

    return {"interactions": interactions, "ligands": ligands, "peptides": peptides, "targets": targets}


# Cleans tha names of targets and ligands in the target and ligands tables
def remove_html_tags(text, verbouse=True):
    """Remove html tags from a string"""

    import html
    import re

    clean = re.compile('<.*?>')
    if type(text) is str:
        try:
            clean_text = html.unescape(re.sub(clean, '', text))
            return clean_text
        except:
            if verbouse:
                print("could not clean ", text)
            return text


def clean_table(table):
    """Clean Html from any table and return the cleaned table."""

    return table.applymap(
        lambda text: remove_html_tags(text, verbouse=False)
    )


# Use this function to clean any new database tables.
def clean_GTP_tables():
    """Cleans all GTP tables from html and saves them as as new csv files."""

    import pkg_resources

    data = import_GTP_tables()
    for key, table in data.items():
        data[key] = clean_table(table)
        data[key].to_csv(pkg_resources.resource_filename(
            __name__, (f"data/GTP_tables_clean/{version}/{key}.csv")
        ))


# Input and Output----------------------------------------------------------------------------------------------
def get_data(data):
    """loads the data as a pandas DataFrame. 

    available tables are:

    - ligands
    - targets
    - interactions
    - peptides

    returns: pandas DataFrame
    """
    import pkg_resources
    import pandas as pd
    path = pkg_resources.resource_filename(
        __name__, (f"data/GTP_tables_clean/{version}/{data}.csv"))
    df = pd.read_csv(path)
    return df


# Annotation support ------------------------------------------------------------------------------------
    # Searches for orthogonal genes given an organism for which the gene exist
    # and a target organism for which to search.
def find_orth_gene(gene, organism, target):
    """Find orthogonal gene via Gprofiler

    returns a list of gene(s)"""
    import pandas as pd
    from gprofiler import GProfiler
    gp = GProfiler()

    if organism == target:  # do not search if original gene is known
        genes = list(set([gene, ]))
    else:
        results = pd.DataFrame(
            gp.orth(query=gene, organism=organism, target=target))
        results.dropna(subset=["name"], axis=0)
        genes = [gene for gene in results.name if gene != "N/A"]

    return genes


def get_peptide_ligands(organism=organism, save=True, verbouse=False):
    """Find genes encoding peptidergic ligands using gprofiler orthogonal gene search and
    the interactions table from GTP

    Saves updated peptide ligands to data/Gene_annotation/{organism}/peptide_ligands.csv
    Returns: pandas DataFrame
    """

    from scConnect.database import get_data
    import pandas as pd
    import numpy as np
    import pkg_resources
    import sys

    interactions = get_data("interactions")

    # Clean the dataframe removing duplicates of ligands from the same species
    dropped_dup = interactions.drop_duplicates(
        subset=["ligand", "ligand_species"])

    # All ligands need to have a gene and species attribute
    dropped_na = dropped_dup.dropna(
        axis=0, subset=["ligand_species", "ligand_gene_symbol"])

    # Subset dataframe to ligands, species and gene symbols
    gene_ligands = dropped_na[[
        "ligand", "ligand_species", "ligand_gene_symbol"]]

    # splits gene information from several species into induvidual records
    new_gene_list = []
    for row in gene_ligands.iterrows():
        ligand = row[1].ligand
        species = row[1].ligand_species.split("|")
        genes = row[1].ligand_gene_symbol.split("|")

        # Checks if two genes were reffered to.
        if len(species) < len(genes):
            genes = "|".join(genes)
            species_list = [ligand, str(species[0]), genes]
            new_gene_list.append(species_list)

        else:
            for i in np.arange(len(species)):
                species_list = [ligand, species[i], genes[i]]
                new_gene_list.append(species_list)

    gene_ligands_df = pd.DataFrame(new_gene_list, columns=[
        "ligand",
        "ligand_species",
        "ligand_gene_symbol"], dtype=str
    )
    # pivot gene symbols for each ligand and species.
    gene_pivot = gene_ligands_df.pivot(
        index="ligand", columns="ligand_species", values="ligand_gene_symbol")

    inferred_genes = list()

    # Add organism names. Key should match column names from GTP,
    # and values should match species names from gprofiler.
    # Note that the order of the dictionary matters (Py3.5>),
    # and if one gene is found ,the search complete.
    organism_dict = {
        "Mouse": "mmusculus",
        "Rat": "rnorvegicus",
        "Human": "hsapiens"
    }

    inferred_genes = list()
    for ligand in gene_pivot.index:
        genes = gene_pivot.loc[ligand]
        for name in organism_dict:
            if str(genes[name]) != "nan":
                gene = find_orth_gene(
                    gene=genes[name], organism=organism_dict[name], target=organism)
            else:
                continue

            if len(gene) > 0:
                comment = f"inferred from {name}"
                break
            else:
                gene = "nan"
                comment = "could not find gene"

        result = [ligand, gene, comment]

        inferred_genes.append(result)

        # add some output estimating % done
        percent = len(inferred_genes) / len(gene_pivot) * 100
        if verbouse:
            sys.stdout.write("\r" + str(np.around(percent, decimals=2)
                                        ) + "% done. \n Adding ligand: " + ligand)

    inferred_genes = pd.DataFrame(inferred_genes, columns=[
                                  "ligand", "gene", "comment"])

    inferred_genes = inferred_genes[inferred_genes.gene != "nan"]

    ligands = pd.DataFrame({"ligand": inferred_genes.ligand,
                            "ligand_type": "peptide",
                            "preprogene": inferred_genes.gene,
                            "synthesis": None,
                            "transport": None,
                            "reuptake": None,
                            "excluded": None,
                            "comment": inferred_genes.comment}
                           )

    if save == True:
        ligands.to_csv(pkg_resources.resource_filename(
            __name__, (f"data/Gene_annotation/{version}/{organism}/peptide_ligands.csv")))

    return ligands


def get_molecule_ligands(organism=organism, save=True):
    '''
    Creates a file with the molecule ligands

    Saves updated molecule ligands to data/Gene_annotation/molecule_ligands.csv
    returns: pandas DataFrame

    '''
    from scConnect.molecules import get_molecules
    import pkg_resources
    import pandas as pd

    mol = get_molecules()

    molecules = pd.DataFrame({
        "ligand": [molecule["ligand"] for molecule in mol],
        "ligand_type": [molecule["type"] for molecule in mol],
        "preprogene": [None for molecule in mol],
        "synthesis": [molecule["synthesis"] for molecule in mol],
        "transport": [molecule["transport"] for molecule in mol],
        "reuptake": [molecule["reuptake"] for molecule in mol],
        "excluded": [molecule["excluded"] for molecule in mol],
        "comment": ["Manually added genes" for molecule in mol]
    })

    drop_ligand = []  # detect if a critical gene category is not detected, then remove the ligand
    if organism != "mmusculus":
        types = ["synthesis", "transport", "reuptake", "excluded"]
        for ty in types:
            inferred_genes = {}
            for i, genes in enumerate(molecules[ty]):
                if genes != None:
                    inferred_genes[i] = list()
                    for gene in genes:
                        inferred_gene = find_orth_gene(
                            gene, organism="mmusculus", target=organism)
                        for gene in inferred_gene:
                            inferred_genes[i].append(gene)
                    if not len(inferred_genes[i]) > 0:
                        drop_ligand.append(i)
                else:
                    inferred_genes[i] = None
            molecules[ty] = inferred_genes.values()

    molecules.drop(index=list(set(drop_ligand)), inplace=True)
    molecules.reindex()

    if save is True:
        molecules.to_csv(pkg_resources.resource_filename(
            __name__, (f"data/Gene_annotation/{version}/{organism}/molecule_ligands.csv")))

    return molecules


def merge_ligand_lists(organism=organism, save=True):
    '''Returns a list with peptides and molecule ligands merged, 
    ready to be used with the package tools.

    saves updated ligands list to data/Gene_annotation/ligands.csv
    returns: pandas DataFrame
    '''
    import pandas as pd
    import pkg_resources

    molecules = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/molecule_ligands.csv")), index_col=0)

    peptides = pd.read_csv(pkg_resources.resource_filename(
        __name__, (f"data/Gene_annotation/{version}/{organism}/peptide_ligands.csv")), index_col=0)

    ligands = peptides.append(molecules, ignore_index=True, sort=False)

    if save is True:
        ligands.to_csv(pkg_resources.resource_filename(
            __name__, (f"data/Gene_annotation/{version}/{organism}/ligands.csv")))
    return ligands


def get_receptors(organism=organism, receptor_types=["gpcr", "lgic"], save=True):
    '''
    returns a list of all receptors with known genes. 
    Finds ortholog genes for human, then mouse then rat.

    Many receptor types are available, defaults to only GPCRs and LGIC:
    enzyme                1259
    transporter            510
    gpcr                   416
    catalytic_receptor     307
    other_protein          228
    vgic                   145
    lgic                    86
    other_ic                54
    nhr                     49

    Saves updated receptor list to data/Gene_annotation/receptors.csv
    Returns: pandas DataFrame
    '''
    from scConnect.database import get_data
    import pandas as pd
    import pkg_resources

    receptors = get_data("targets")

    # Add all receptors to a list with correponding gene.
    # When two genes are found, these are added as a list.

    organism_dict = {  # decides in what order to look for orthological genes
        "HGNC symbol": "hsapiens",
        "MGI symbol": "mmusculus",
        "RGD symbol": "rnorvegicus"
    }

    receptor_list = list()
    for i, receptor in receptors.iterrows():
        if receptor["Type"] in receptor_types:  # select only receptors of specified types

            name = receptor["Target name"]
            family = receptor["Family name"]
            receptor_type = receptor["Type"]
            for species in organism_dict:
                if str(receptor[species]) != "nan":
                    genes = find_orth_gene(
                        gene=receptor[species], organism=organism_dict[species], target=organism)
                    break
                else:
                    genes = []
            if len(genes) > 0:
                receptor_list.append((name, family, receptor_type, genes))

    receptor_genes = pd.DataFrame(receptor_list, columns=[
                                  "receptor", "family", "type", "gene"])

    if save is True:
        receptor_genes.to_csv(pkg_resources.resource_filename(
            __name__, (f"data/Gene_annotation/{version}/{organism}/receptors.csv")))
    return receptor_genes


# build interaction table
def get_interactions(save=True):
    """build interaction table from GTP interaction df with only relevant data.

    Saves updated interaction table to data/Gene_annotation/interactions.csv.
    Returns: pandas DataFrame
    """
    import pkg_resources

    # Load data
    interactions = get_data("interactions")
    # remove unused receptors
    interactions = interactions.dropna(axis=0, subset=["target"])
    # Make a multiindex for quick indexing when finding interactions
    interactions.index = [interactions.ligand, interactions.target]

    # Keep useful columns and remove nan values for targets
    interactions = interactions[['endogenous',
                                 'ligand_species', 'target_species', 'pubmed_id', 'action']]

    # use ";" to seperate values due to "," being used in ligand names
    if save is True:
        interactions.to_csv(pkg_resources.resource_filename(
            __name__, (f"data/Gene_annotation/{version}/interactions.csv")), sep=";")

    return interactions


def setup_database(
    organism=organism,
    receptor_types=[
        "enzyme", "transporter", "gpcr", "catalytic_receptor",
        "other_protein", "vgic", "lgic", "other_ic", "nhr"]):
    '''Use this function to recalculate all
     ligands and receptor gene lists based on
     GTP tables located under data/GTP_tables.
     Manually create a folder for new species
     before running this function under /data/Gene_annotation.

     i.e. /data/Gene_annotation/mmusculus -- for mouse


     .. note:: 
        Download new `Guide to phalmacology data`_ and save the tables to data/GTP_tables.

    .. _Guide to phalmacology data: https://www.guidetopharmacology.org/download.jsp

    Many receptor types are available, defaults to only GPCRs and LGIC:
    enzyme                1259
    transporter            510
    gpcr                   416
    catalytic_receptor     307
    other_protein          228
    vgic                   145
    lgic                    86
    other_ic                54
    nhr                     49
    '''
    import scConnect as cn

    # Clean the new GTP tables from HTML
    cn.database.clean_GTP_tables()

    # Compute ligand genes from known genes and from manually added ligands
    # found in molecules.py
    print("getting peptide ligands...")
    pep = cn.database.get_peptide_ligands(organism=organism)
    print(f'found {pep.shape[0]} peptide ligands')

    print("getting molecule ligands...")
    mol = cn.database.get_molecule_ligands(organism=organism)
    print(f'found {mol.shape[0]} molecular ligands')

    print("merging ligands...")
    cn.database.merge_ligand_lists(organism=organism)

    print("getting receptors...")
    rec = cn.database.get_receptors(
        organism=organism, receptor_types=receptor_types)
    print(f'found {rec.shape[0]} receptors')
    print(rec.type.value_counts())

    print("getting interactions...")
    cn.database.get_interactions()

    print("Database build complete")
