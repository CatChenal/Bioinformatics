# PostEraLeadsSmiles
### Update: Obsolete project (2020-04-04)
The data files of the lead compounds submitted to the [COVID Moonshot](https://discuss.postera.ai/t/round-2-is-finished-now-a-call-for-help/816/11) are now hosted on Github: https://github.com/mc-robinson/COVID_moonshot_submissions.

---
---
### Python: 3.6.7
---
Update:
2020-03-27: Refactored main function to output a json file: submitter_id:[smiles of each entry]
---
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
