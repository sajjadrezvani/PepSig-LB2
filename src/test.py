import requests
import csv


# ============================================================
# 1. UniProt API query
# ============================================================

base_url = "https://rest.uniprot.org/uniprotkb/search"

query = (
    "(reviewed:true) AND "
    "(taxonomy_id:2759) AND "
    "(length:[40 TO *]) AND "
    "(existence:1) AND "
    "(fragment:false) AND "
    "(ft_signal_exp:*)"
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

        # Example:
        # <https://rest.uniprot.org/...&cursor=...>; rel="next"

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

    # The next URL already contains the query parameters
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
# 4. Filter positive proteins
# ============================================================

positive_proteins = []

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
    # Find signal peptide
    # --------------------------------------------------------

    signal_peptide = None

    for feature in entry.get("features", []):

        if feature["type"] == "Signal":

            location = feature.get("location", {})

            start = location.get("start", {}).get("value")

            end = location.get("end", {}).get("value")

            if start is not None and end is not None:

                signal_peptide = (start, end)

                break


    # --------------------------------------------------------
    # Skip if no usable signal peptide
    # --------------------------------------------------------

    if signal_peptide is None:

        continue


    signal_start, cleavage_position = signal_peptide


    # --------------------------------------------------------
    # Signal peptide length
    # --------------------------------------------------------

    signal_length = (
        cleavage_position - signal_start + 1
    )


    # --------------------------------------------------------
    # Project requirement:
    # signal peptide length > 13
    # --------------------------------------------------------

    if signal_length <= 13:

        continue


    # --------------------------------------------------------
    # Store protein
    # --------------------------------------------------------

    positive_proteins.append({

        "accession": accession,

        "organism": organism,

        "kingdom": kingdom,

        "length": protein_length,

        "cleavage_position": cleavage_position,

        "sequence": sequence

    })


# ============================================================
# 5. Number of positive proteins
# ============================================================

print()
print(
    "POSITIVE PROTEINS AFTER FILTERING:",
    len(positive_proteins)
)


# ============================================================
# 6. Save TSV
# ============================================================

tsv_file = "../data/raw/positive.tsv"

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
        "protein_length",
        "cleavage_position"
    ])

    for protein in positive_proteins:

        writer.writerow([
            protein["accession"],
            protein["organism"],
            protein["kingdom"],
            protein["length"],
            protein["cleavage_position"]
        ])


# ============================================================
# 7. Save FASTA
# ============================================================

fasta_file = "../data/raw/positive.fasta"

with open(
    fasta_file,
    "w",
    encoding="utf-8"
) as f:

    for protein in positive_proteins:

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
    "Positive proteins after filtering:",
    len(positive_proteins)
)

print()
print("Files created:")

print(tsv_file)
print(fasta_file)
print()
print("==============================================")
print("FINAL NUMBER OF RESULTS")
print("==============================================")
print("Before filtering:", len(all_results))
print("After filtering:", len(positive_proteins))
