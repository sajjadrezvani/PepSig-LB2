import requests
import csv
from pathlib import Path


# Project paths
project_root = Path(__file__).resolve().parent.parent
raw_dir = project_root / "data" / "raw"
raw_dir.mkdir(parents=True, exist_ok=True)


# Positive dataset: reviewed, eukaryotic, >=40 aa, protein-level evidence,
# not a fragment, and experimentally supported signal peptide
query_pos = (
    "(reviewed:true) AND "
    "(taxonomy_id:2759) AND "
    "(length:[40 TO *]) AND "
    "(existence:1) AND "
    "(fragment:false) AND "
    "(ft_signal_exp:*)"
)

# Negative dataset: same general criteria, no signal peptide at any evidence
# level, and experimental localization to the specified cellular compartments
query_neg = (
    "(reviewed:true) AND "
    "(taxonomy_id:2759) AND "
    "(length:[40 TO *]) AND "
    "(existence:1) AND "
    "(fragment:false) AND "
    "NOT (ft_signal:*) AND "
    "("
    "(cc_scl_term_exp:SL-0091) OR "
    "(cc_scl_term_exp:SL-0191) OR "
    "(cc_scl_term_exp:SL-0173) OR "
    "(cc_scl_term_exp:SL-0209) OR "
    "(cc_scl_term_exp:SL-0204) OR "
    "(cc_scl_term_exp:SL-0039)"
    ")"
)


def retrieve_all_uniprot(query):
    """Retrieve all UniProt entries matching a query using pagination."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    all_results = []
    url = base_url
    params = {"query": query, "format": "json", "size": 500}
    page = 1

    while True:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        results = data.get("results", [])
        all_results.extend(results)
        print(f"Page {page}: +{len(results)} | Total: {len(all_results)}")

        # UniProt provides the next page URL in the HTTP Link header
        link_header = response.headers.get("Link")
        next_url = None

        if link_header:
            for link in link_header.split(","):
                if 'rel="next"' in link:
                    next_url = link.split(";")[0].strip().strip("<>")
                    break

        if next_url is None:
            break

        url = next_url
        params = None
        page += 1

    print(f"Total retrieved: {len(all_results)}")
    return all_results


def get_kingdom(entry):
    """Assign the project kingdom category from the UniProt lineage."""
    lineage = entry["organism"]["lineage"]

    if "Metazoa" in lineage:
        return "Metazoa"
    if "Fungi" in lineage:
        return "Fungi"
    if "Viridiplantae" in lineage:
        return "Plants"
    return "Other"


def process_positive(results):
    """Apply the two positive-specific filters and collect required fields."""
    positive_proteins = []

    for entry in results:
        accession = entry["primaryAccession"]
        organism = entry["organism"]["scientificName"]
        protein_length = entry["sequence"]["length"]
        kingdom = get_kingdom(entry)

        # Find the signal peptide feature
        signal_peptide = None

        for feature in entry.get("features", []):
            if feature["type"] != "Signal":
                continue

            location = feature.get("location", {})
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")

            # A numerical end position means the cleavage position is known
            if start is not None and end is not None:
                signal_peptide = (start, end)
                break

        # Exclude proteins with an unknown cleavage site
        if signal_peptide is None:
            continue

        signal_start, cleavage_position = signal_peptide
        signal_length = cleavage_position - signal_start + 1

        # Project requirement: signal peptide must be longer than 13 aa
        if signal_length <= 13:
            continue

        positive_proteins.append({
            "accession": accession,
            "organism": organism,
            "kingdom": kingdom,
            "length": protein_length,
            "cleavage_position": cleavage_position
        })

    return positive_proteins


def process_negative(results):
    """Collect negative proteins and record TM helix status in first 90 aa."""
    negative_proteins = []

    for entry in results:
        accession = entry["primaryAccession"]
        organism = entry["organism"]["scientificName"]
        protein_length = entry["sequence"]["length"]
        kingdom = get_kingdom(entry)

        # Safety check: exclude any entry containing a signal peptide
        has_signal_peptide = any(
            feature["type"] == "Signal"
            for feature in entry.get("features", [])
        )

        if has_signal_peptide:
            continue

        # Check whether a transmembrane helix starts within the first 90 residues
        tm_first_90 = False

        for feature in entry.get("features", []):
            if feature["type"] != "Transmembrane":
                continue

            location = feature.get("location", {})
            start = location.get("start", {}).get("value")

            if start is not None and start <= 90:
                tm_first_90 = True
                break

        negative_proteins.append({
            "accession": accession,
            "organism": organism,
            "kingdom": kingdom,
            "length": protein_length,
            "tm_first_90": tm_first_90
        })

    return negative_proteins


def save_positive_tsv(proteins, output_file):
    """Save the required positive dataset information as TSV."""
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "accession",
            "organism",
            "kingdom",
            "protein_length",
            "cleavage_position"
        ])

        for protein in proteins:
            writer.writerow([
                protein["accession"],
                protein["organism"],
                protein["kingdom"],
                protein["length"],
                protein["cleavage_position"]
            ])


def save_negative_tsv(proteins, output_file):
    """Save the required negative dataset information as TSV."""
    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow([
            "accession",
            "organism",
            "kingdom",
            "protein_length",
            "tm_helix_first_90"
        ])

        for protein in proteins:
            writer.writerow([
                protein["accession"],
                protein["organism"],
                protein["kingdom"],
                protein["length"],
                protein["tm_first_90"]
            ])


def download_fasta(accessions, output_file, batch_size=50):
    """Download FASTA sequences in small batches to avoid large API requests."""
    with open(output_file, "w") as f:
        for i in range(0, len(accessions), batch_size):
            batch = accessions[i:i + batch_size]
            query = " OR ".join(f"accession:{acc}" for acc in batch)

            response = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": query,
                    "format": "fasta",
                    "size": batch_size
                },
                timeout=120
            )
            response.raise_for_status()
            f.write(response.text)

            print(f"FASTA: {min(i + batch_size, len(accessions))}/{len(accessions)}")

# Retrieve all positive and negative candidates
print("Retrieving positive data...")
positive_results = retrieve_all_uniprot(query_pos)

print("Retrieving negative data...")
negative_results = retrieve_all_uniprot(query_neg)


# Apply project-specific processing
positive_proteins = process_positive(positive_results)
negative_proteins = process_negative(negative_results)


# Get accessions of the final selected proteins
positive_accessions = [p["accession"] for p in positive_proteins]
negative_accessions = [p["accession"] for p in negative_proteins]


# Save TSV files
positive_tsv = raw_dir / "positive.tsv"
negative_tsv = raw_dir / "negative.tsv"

save_positive_tsv(positive_proteins, positive_tsv)
save_negative_tsv(negative_proteins, negative_tsv)


# Download FASTA files directly from UniProt
positive_fasta = raw_dir / "positive.fasta"
negative_fasta = raw_dir / "negative.fasta"

download_fasta(positive_accessions, positive_fasta)
download_fasta(negative_accessions, negative_fasta)


# Final summary
print(
    f"\nPositive: {len(positive_results)} retrieved → "
    f"{len(positive_proteins)} final"
)
print(
    f"Negative: {len(negative_results)} retrieved → "
    f"{len(negative_proteins)} final"
)
print(f"Files saved to: {raw_dir}")


