# Version 0.1.2
### Minor changes
 - Fixed dependency issue where packages were not being loaded upon starting `drugminer`. (for future reference: there seems to be an issue with Imports versus Depends in the DESCRIPTION file, where the former does not seem to work that well while the latter always does. see more [here](https://stackoverflow.com/questions/37568884/r-package-does-not-load-dependencies).)

# Version 0.1.1
## Major changes
 - Was missing the `pubchem2cas` function that is called on the `compound_query` wrapper.
 
## Minor changes
 - Bug fixes on `compound_query` function
 - Fixed dependency issue where packages were not being loaded upon starting `drugminer`

# Version 0.0.2
## Major changes
  - Added the `retrieve_data`, `parse_drugbank`, and `kegg_processing` functions, for completeness. These will allow for the generation of the data frames that are used internally.

## Minor changes
  - Made changes in documentation of all functions to include examples, as well as fixing minor typos.


## Bug Fixes
  - Fixed a couple of bugs on function linking in the documentation.

# Version 0.0.1
First release version of the DrugMineR package