import requests
import csv


# ============================================================
# 1. UniProt API query - NEGATIVE DATASET
# ============================================================

base_url = "https://rest.uniprot.org/uniprotkb/search"

query = (
    "(reviewed:true) AND "
    "(taxonomy_id:2759) AND "
    "(length:[40 TO *]) AND "
    "(existence:1) AND "
    "(fragment:false) AND "
    "NOT (ft_signal:*) AND "
    "("
        "(cc_scl_term_exp:SL-0091) OR "   # Cytosol
        "(cc_scl_term_exp:SL-0191) OR "   # Nucleus
        "(cc_scl_term_exp:SL-0173) OR "   # Mitochondrion
        "(cc_scl_term_exp:SL-0209) OR "   # Plastid
        "(cc_scl_term_exp:SL-0204) OR "   # Peroxisome
        "(cc_scl_term_exp:SL-0039)"       # Cell membrane
    ")"
)


# ============================================================
# 2. Retrieve ALL results using UniProt pagination
# ============================================================

all_results = []

url = base_url

params = {
    "query": query,
    "format": "json",
    "size": 500
}

page = 1

while True:

    print(f"Retrieving page {page}...")

    response = requests.get(
        url,
        params=params
    )

    response.raise_for_status()

    data = response.json()

    results = data.get("results", [])

    all_results.extend(results)

    print(f"  Retrieved this page: {len(results)}")
    print(f"  Total retrieved so far: {len(all_results)}")


    # ========================================================
    # Get next page from HTTP Link header
    # ========================================================

    link_header = response.headers.get("Link")

    next_url = None

    if link_header:

        links = link_header.split(",")

        for link in links:

            if 'rel="next"' in link:

                next_url = (
                    link
                    .split(";")[0]
                    .strip()
                    .strip("<>")
                )

                break


    # ========================================================
    # Stop when there is no next page
    # ========================================================

    if next_url is None:
        break


    # ========================================================
    # Continue to next page
    # ========================================================

    url = next_url

    # Next URL already contains query parameters
    params = None

    page += 1


# ============================================================
# 3. Overall number retrieved
# ============================================================

print()
print("==============================================")
print("TOTAL PROTEINS RETRIEVED:", len(all_results))
print("==============================================")


# ============================================================
# 4. Prepare negative proteins
# ============================================================

negative_proteins = []

for entry in all_results:

    accession = entry["primaryAccession"]

    organism = entry["organism"]["scientificName"]

    protein_length = entry["sequence"]["length"]

    sequence = entry["sequence"]["value"]


    # --------------------------------------------------------
    # Determine eukaryotic kingdom
    # --------------------------------------------------------

    lineage = entry["organism"]["lineage"]

    if "Metazoa" in lineage:

        kingdom = "Metazoa"

    elif "Fungi" in lineage:

        kingdom = "Fungi"

    elif "Viridiplantae" in lineage:

        kingdom = "Plants"

    else:

        kingdom = "Other"


    # --------------------------------------------------------
    # Safety check:
    # make sure there is NO signal peptide
    # --------------------------------------------------------

    has_signal_peptide = False

    for feature in entry.get("features", []):

        if feature["type"] == "Signal":

            has_signal_peptide = True

            break


    # Skip proteins with a signal peptide
    if has_signal_peptide:

        continue


    # --------------------------------------------------------
    # Store negative protein
    # --------------------------------------------------------

    negative_proteins.append({

        "accession": accession,

        "organism": organism,

        "kingdom": kingdom,

        "length": protein_length,

        "sequence": sequence

    })


# ============================================================
# 5. Number of final negative proteins
# ============================================================

print()
print(
    "NEGATIVE PROTEINS AFTER FILTERING:",
    len(negative_proteins)
)


# ============================================================
# 6. Save TSV
# ============================================================

tsv_file = "../data/raw/negative.tsv"

with open(
    tsv_file,
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.writer(
        f,
        delimiter="\t"
    )

    writer.writerow([
        "accession",
        "organism",
        "kingdom",
        "protein_length"
    ])

    for protein in negative_proteins:

        writer.writerow([
            protein["accession"],
            protein["organism"],
            protein["kingdom"],
            protein["length"]
        ])


# ============================================================
# 7. Save FASTA
# ============================================================

fasta_file = "../data/raw/negative.fasta"

with open(
    fasta_file,
    "w",
    encoding="utf-8"
) as f:

    for protein in negative_proteins:

        f.write(
            f">{protein['accession']}|"
            f"{protein['organism']}|"
            f"{protein['kingdom']}\n"
        )

        sequence = protein["sequence"]

        for i in range(
            0,
            len(sequence),
            60
        ):

            f.write(
                sequence[i:i + 60] + "\n"
            )


# ============================================================
# 8. Final summary
# ============================================================

print()
print("==============================================")
print("FINAL SUMMARY")
print("==============================================")

print(
    "Total UniProt proteins retrieved:",
    len(all_results)
)

print(
    "Final negative proteins:",
    len(negative_proteins)
)

print()
print("Files created:")

print(tsv_file)
print(fasta_file)
