# PostEraLeadsSmiles
### Module: postera_smiles.py
### Python: 3.6.7
---
Update:
2020-03-26: Fix write mode, 'w+' -> 'w'.

On-going review due to design flaw:
* 1. The extracting function is likely missing additional entries past the first one due to the fact that it is not visiting the submission details page. Unlike what I initially thought, it is not sufficient.
* 2. Need to come up with an incremental process of the SMILES extraction; right now, all listing (on the submmission apge) are extracted, not just the new one => Need a way to determine new ones. 
---
## Function `extract_postera_smiles()`: parses https://covid.postera.ai/covid/submissions in order to retrieve the submission id and its [SMILES](https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system) code.
### Note: If the `save_as` argument is None (default), the function returns a list of (id, smiles) tuples.
NOTE: The following description is erroneous (the function retrieves only the first of possible multiple submissions for a submitter): will be updated with a new verion of the function.

The div(card h-100) container on the submission page has all the info needed[1]: 
* submission id: from the first a.href
* smiles code: from the img.alt
```
<div class="card h-100">
    <div class="card-header">
      <h5 class="mb-1"><a href="submissions/fd8d85a5-db4a-4ecf-93d1-416eee50b961">JAN-GHE-fd8</a></h5>
    </div>
    <img src="/synthesize/C1(C=C(C%23N)C=NC=1NC1CC1)CNCCNC(C)=O" class="card-img" alt="C1(C=C(C#N)C=NC=1NC1CC1)CNCCNC(C)=O">
    <div class="card-body">
      <a style="width: 100%" href="/covid/submissions/fd8d85a5-db4a-4ecf-93d1-416eee50b961" class="mb-2 btn btn-primary">View</a>

    </div>
</div>
```
[1] Given the html structure as at 03-25-2020.

## Call:

### As a function:
```
import postera_smiles

postera_smiles.extract_postera_smiles(save_as='data/postera_smiles.csv')
```

### At the comman line:
```
python postera_smiles.py data/postera_smiles.csv
```
---