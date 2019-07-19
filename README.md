# DrugMineR

`DrugMineR` is a package to retrieve information for a given drug, based on information found on [KEGG](https://www.kegg.jp), [PubChem](https://pubchem.ncbi.nlm.nih.gov), and [DrugBank](https://www.drugbank.ca), making use of `REST` APIs and data pre-processed. 

## Data pre-processing
Included with the package are data objects extracted from KEGG and DrugBank. These can be easily updated with the provided functions `kegg_processing` and `parse_drugbank`, though this is not recommended (a few tweaks internally need to happen.) For completeness, the process for extracting these data is described below.

### Processing KEGG data
The KEGG data processing is done with the `kegg_processing` function, which queries the `PUG REST` service of KEGG using the `httr` package.  Briefly, this can be done as:

```{r}
kegg_data <- kegg_processing()
```

The `kegg_data` object contains two lists, each containing diverse information from KEGG, such as drug/compound name, pathways, and so on.  Please refer to the data documentation in the package by doing:

```{r}
?KEGG
```

or

```{r}
?DRUG
```

### Parsing DrugBank
DrugMineR also comes with data from DrugBank, which has been processed using the `dbparser` package. From DrugMineR, you can parse the DrugBank data by doing:

```{r} 
drugbank_data <- parse_drugbank(xml_file)
```

where `xml_file` is the full data base download from DrugBank (you need this file before running the parser). The output is a set of tibbles that extract relevant information from the XML file. See what information is captured by doing:

```{r}
?DRUGBANK
```


## Running DrugMineR
The best way to use the DrugMineR package is to make use of its `compound_query` wrapper. As an example, if we want to extract all of the known information on aspirin, we can just call the wrapper as:

```{r}
compound_query("aspirin")
```

which will return six tibbles:

 - `general`, for general information about a drug. Tibble includes:
   - drug name
   - KEGG ID
   - PubChem CID
   - DrugBank ID
   - CAS registration number
   - Chemical formula
   - Molecular weight
   - ChEBI ID
   - Canonical SMILES
   - InChIKey
 - `targets`, for known targets of a given drug;
 - `pathways`, for known pathways that are affected by the drug;
 - `diseases`, for diseases/syndromes that the drug *can* be used against (not a comprehensive list)
 - `drug_uses`, for indications for a given drug (e.g., analgesic, anti-neoplastic)
 - `drug_groups`, for KEGG mappings for drugs into larger classes (e.g., calcium blocker or SSRI)

The `compound_query` function will call the `kegg_query` function first, to map a given compound string to information present in KEGG (looking in the compound, drug, and synonyms tables.) When a compound in *not* found in KEGG, it will then be searched directly on PubChem, using the `name2cid` and the `property_extractor` functions, in sequence. What these functions do is to 1) translate the compound string into a PubChem CID (in `name2cid`) and 2) get all of the information on the compound by querying PubChem through its `PUG REST` API. Of note, we also query the [Chemical Translation Service](http://cts.fiehnlab.ucdavis.edu)  at UC Davis to extract the CAS number for a given PubChem CID through their `REST` services. In the event that no PubChem ID is found for a given query, then NAs are returned.


## Future developments
I am always open to adding features and functionalities to `DrugMineR` and would be excited about possible collaborations to develop the package further.  Please [contact me](mailto:diogo.camacho@wyss.harvard.edu) if you want to chat. 