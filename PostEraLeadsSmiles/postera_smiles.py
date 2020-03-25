# -*- coding: utf-8 -*-
# Module: postera_smiles.py
# Tested on python version: 3.6.7

def extract_postera_smiles(url=None, save_as=None):
    """
    Extract the lead compounds submission id and smile
    code from PostEra.ai COVID lead submisson page.
    Input:
    url: "https://covid.postera.ai/covid/submissions" if None
    save_as: valid csv filename or None
    Output:
    A list of (id, smiles) tuples if save_as is None, else
    a csv file.
    """
    import requests
    from bs4 import BeautifulSoup as btfsoup

    # PostEra Leads Submissions site:
    if url is None:
        url = "https://covid.postera.ai/covid/submissions"

    r = requests.get(url)
    if r.status_code == 200: 
        soup = btfsoup(r.content, 'html.parser')
    else:
        print(f'Oops? Request status: {r.status_code}')
        return None
    
    def get_list(soup):
        # The div(card h-100) container has all the info needed: 
        # submission id: from the first a.href
        # smiles code: from the img.alt
        data = []
        for card in soup.find_all('div', class_='card h-100'):
            submit_id = card.find_next('a').get('href').split('/')[1]
            smiles = card.find_next('img').get('alt')
            data.append((submit_id, smiles))
        return data
    
    id_smiles = get_list(soup)
    
    if save_as is None:
        return id_smiles
    else:
        with open(save_as, 'w+') as fh:
            fh.write('submission_id, submission_smiles\n')
            for i, v in enumerate(id_smiles):
                fh.write(f'{v[0]}, {v[1]}\n')
        print(f'Submission id and smiles saved in {save_as}.')  
        
        
if '__name__' == '__main__':
    import sys
    
    fname = sys.argv[1:].strip()
    extract_postera_smiles(save_as=fname)