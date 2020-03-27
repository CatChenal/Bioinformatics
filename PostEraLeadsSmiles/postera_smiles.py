# -*- coding: utf-8 -*-
# Module: postera_smiles.py
# Tested on python version: 3.6.7
__doc__="""
# Module: postera_smiles.py
The main function `get_submissions_dict()` performs two extractions:
* First extraction on https://covid.postera.ai/covid/submissions:
  The submitter's details page address from the anchor.href of the h5 (class="mb-1") tag.
* Second extraction on the the submitter's details page:
  All the alt.string from the img (class="card-img") tags.

Submission page entry example on https://covid.postera.ai/covid/submissions:
```
<div class="card h-100">
    <div class="card-header">
      <h5 class="mb-1"><a href="submissions/fedd1f79-9ed6-4487-a801-e6f14abf8e11">PET-SGC-fed</a></h5>
    </div>
    <img src="/synthesize/CC1C=CN=CC=1NC(=O)CCCCN1CCN(C(=O)COC2C=CC(C)=CC=2)CC1" class="card-img" alt="CC1C=CN=CC=1NC(=O)CCCCN1CCN(C(=O)COC2C=CC(C)=CC=2)CC1">
    <div class="card-body">
      <a style="width: 100%" href="/covid/submissions/fedd1f79-9ed6-4487-a801-e6f14abf8e11" class="mb-2 btn btn-primary">View</a>
    </div>
</div>
```
Details page example (here with single entry): 
```
<div id="smiles_list" class="row row-cols-1 row-cols-md-4">  
    <div class="col col-md-3 mb-4">
      <div class="card h-100">
        <img src="/synthesize/CC1C=CN=CC=1NC(=O)CCCCN1CCN(C(=O)COC2C=CC(C)=CC=2)CC1" class="card-img" alt="CC1C=CN=CC=1NC(=O)CCCCN1CCN(C(=O)COC2C=CC(C)=CC=2)CC1">
        <div class="card-body text-left">
          <p><strong>PET-SGC-fed-1</strong></p>
          <p>CC1C=CN=CC=1NC(=O)CCCCN1CCN(C(=O)COC2C=CC(C)=CC=2)CC1</p>
        </div>
      </div>
    </div>
</div>
```
Call examples:
--------------
## 1. Command line: Optional argument is a folder path. 
The json filename has pattern <2020_12_31_23_59>_postera_smiles.json.
Command:
```
python postera_smiles.py ./data/intermediate
```
Output:
```
Getting PostEra.ai COVID submissions SMILES codes (Note: slow process!)...
Submissions SMILES json file saved as:
    data/intermediate/2020_03_27_15_38_postera_smiles.json
```

## 2. Using the postera_smiles.py module:
```
import postera_smiles as psmiles
from pathlib import Path

submitters_dict = psmiles.get_submissions_dict()

fname = Path.cwd().joinpath('data', 'intermediate', psmiles.get_stamped_filename())
psmiles.save_as_json(fname, submitters_dict)

# to load:

jdict = psmiles.load_json(fname)
list(jdict.keys())[:5]
```
"""
import requests
from bs4 import BeautifulSoup as btfsoup, SoupStrainer
import json
import datetime

POSTERA_COVID_URL = "https://covid.postera.ai/covid/"
FNAME_SUFFIX = 'postera_smiles.json'


def get_stamped_filename(dto=datetime.datetime.today(),
                         frmt='%Y_%m_%d_%H_%M',
                         suffix=FNAME_SUFFIX):
    """
    Create csv filename with default pattern today().strftime(frmt)_suffix.
    """
    fname = F"{dto.strftime(frmt)}_{suffix}"
    return fname


def save_as_json(fname, obj):
    with open(fname, 'w') as fh:
        fh.write(json.dumps(obj))


def load_json(fname):
    with open(fname) as fh:
        jdict = json.load(fh)
    return jdict


def get_soup(url, only_parse=None):
    r = requests.get(url)
    if r.status_code == 200: 
        return btfsoup(r.content, 'html.parser', parse_only=only_parse)
    else:
        print(F'Failed request? url :: status:\n\t{url} :: {r.status_code}')
        return None


def get_details_list(details_url):
    """
    Get the smiles codes on PostEra.ai COVID submissions details page, using img.alt.
    Submission details page: https://covid.postera.ai/covid/submissions/<hash>.
    """
    molec_img = SoupStrainer('img', attrs={'class':'card-img'})
    img_soup = get_soup(details_url, only_parse=molec_img)
    
    if img_soup is not None:
        smiles = [img.get('alt') for img in img_soup] 
    else:
        smiles = []
    return smiles


def get_submissions_dict():
    """
    Return a dict of submitter_id:[submitted_molecules.smiles].
    
    The h5 header of div(card h-100) has the anchor to the details page.
    Submission details url: from h5.a.href;
    """
    submitters = {}
    
    domain_url = POSTERA_COVID_URL
    submission_url = domain_url + 'submissions'
    
    h5_tags = SoupStrainer('h5', attrs={'class':'mb-1'})
    h5_soup = get_soup(submission_url, only_parse=h5_tags)
    
    if h5_soup is not None:
        for h5 in h5_soup:
            a = h5.find_next('a')
            details_url = domain_url + a.get('href')
            submitters[a.string] = get_details_list(details_url)
    return submitters


if '__name__' == '__main__':
    import sys
    from pathlib import Path
    
    if len(sys.argv[1:]):
        data_path = sys.argv[1:].strip()
    else:
        data_path = Path.cwd()
        
    print('Getting PostEra.ai COVID submissions SMILES codes (Note: slow process!)...')
    submitters_dict = get_submissions_dict()
    
    fname = data_path.joinpath(get_stamped_filename())
    save_as_json(fname, submitters_dict)
        
    print(F'Submission SMILES json file saved as:\n\t{fname}')
