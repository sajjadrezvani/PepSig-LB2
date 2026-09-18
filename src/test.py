import requests
import json

accession = "P01165"

url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
data = requests.get(url).json()

print(json.dumps(data["features"], indent=2))