common_wrangler
==============================
[//]: # (Badges)
[![Travis Build Status](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/common_wrangler.png)](https://travis-ci.org/REPLACE_WITH_OWNER_ACCOUNT/common_wrangler)
[![AppVeyor Build status](https://ci.appveyor.com/api/projects/status/REPLACE_WITH_APPVEYOR_LINK/branch/master?svg=true)](https://ci.appveyor.com/project/REPLACE_WITH_OWNER_ACCOUNT/common_wrangler/branch/master)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/common_wrangler/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/common_wrangler/branch/master)

Scripts for file and data manipulation, as well as common scripts used in this project as well as gaussian_wrangler and md_wrangler

### Copyright

Copyright (c) 2019, hmayes

This package contains a common.py that is used in other packages from this developer, such as:
https://github.com/team-mayes/gaussian_wrangler

There are also several scripts that are useful in many scenarios:

*fill_tpl*: set up so that a user can provide a configuration file with values for any parameters, and a 
corresponding template file that has parameter names in brackets. For example, the configuration file can 
have `project_name = common_wrangler` and the template would have `{project_name}` as a field to fill in.

*rename_files*: meant to be used in any terminal to rename any type of file.

#### Acknowledgements
 
Project based on the 
[Computational Chemistry Python Cookiecutter](https://github.com/choderalab/cookiecutter-python-comp-chem)
