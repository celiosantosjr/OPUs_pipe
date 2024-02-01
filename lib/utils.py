import os
import gzip
import lzma

emapper_ofile_fmt=['originalname', 'seed eggNOG ortholog',
                   'seed ortholog evalue', 'seed ortholog score',
                   'Predicted taxonomic group', 'Predicted protein name',
                   'Gene Ontology terms', 'EC number', 'KEGG_ko',
                   'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction',
                   'KEGG_rclass', 'BRITE', 'KEGG_TC', 'CAZy',
                   'BiGG Reaction', 'tax_scope', 'eggNOG OGs',
                   'bestOG', 'COG Functional Category',
                   'eggNOG free text description']


def cr(infile, ofile):
    with open(infile, 'rb') as input_file, lzma.open(ofile, 'wb') as output_file:
        for data in input_file:
            output_file.write(data)


def highest_value_with_percentage(d):
    total_sum = sum(d.values())
    highest_key = max(d, key=d.get)
    highest_value = d[highest_key]
    percentage = (highest_value / total_sum) * 100
    return f"{highest_key} ({percentage:.2f}%)"


def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    else:
        print(f"{directory_path} already exists.")


def corruption_test(infile):
    try:
        if infile.endswith(".gz"):
            with gzip.GzipFile(infile, 'rb') as gzfile:
                for _ in gzfile:  # Try iterating through the entire file
                    pass
                return False
        elif infile.endswith(".xz"):
            with lzma.open(infile, "rb") as xzfile:
                for _ in xzfile:  # Try iterating through the entire file
                    pass
                return False
    except (OSError, EOFError):
        return True


def mark_bad(seq, minsize=10):
    xokay, nox, comp = [False, False, False]
    amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    leftover = set(seq.upper()) - set(amino_acids)
    if leftover:
        if (len(leftover) == 1) and ('X' in leftover):
            if len(seq.upper().split('X')[0]) >= minsize:
                xokay = True
    else:
        nox = True
    if len(seq) >= minsize:
        comp = True
    if comp and nox:
        return True
    elif xokay:
        return True
    else:
        return False
