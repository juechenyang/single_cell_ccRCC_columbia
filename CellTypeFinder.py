import json
import requests
from pandas.io.json import json_normalize
import pandas as pd

def get_cell_type(genes_str):
    #database example: 'KEGG_2015'
    database = "Azimuth_Cell_Types_2021"
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }
    try:
        response = requests.post(ENRICHR_URL, files=payload, timeout=5)
    except:
        print('connection error')
        return 'error'

    if not response.ok:
        return 'error'
    data = json.loads(response.text)
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'
    user_list_id = data['userListId']
    gene_set_library = database
    response = requests.get(
        ENRICHR_URL + query_string % (user_list_id, gene_set_library), timeout=10
    )
    if not response.ok:
        return 'error'
    data = json.loads(response.text)
    printed_data = pd.DataFrame(data=data[database],
                                columns=['Rank', 'Term name', 'P-value', 'Z-score', 'Combined score',
                                         'Overlapping genes', 'Adjusted p-value', 'Old p-value',
                                         'Old adjusted p-value'])
    return printed_data
