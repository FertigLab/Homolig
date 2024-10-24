# installation

- OPTIONAL consider using [pyproject.toml](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/) instead of setup.py?
- DONE - python 3.11 in docker fixing installation that fails with 
  ```
  building 'Bio.cpairwise2' extension
      clang -fno-strict-overflow -Wsign-compare -Wunreachable-code -fno-common -dynamic -DNDEBUG -g -O3 -Wall -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk -I/Users/dlvovs1/envs/homolig/include -I/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/include/python3.13 -c Bio/cpairwise2module.c -o build/temp.macosx-14.0-arm64-cpython-313/Bio/cpairwise2module.o
      Bio/cpairwise2module.c:60:22: error: call to undeclared function 'PyEval_CallObject'; ISO C99 and later do not support implicit function declarations [-Wimplicit-function-declaration]
          if(!(py_result = PyEval_CallObject(py_match_fn, py_arglist)))
                           ^
      Bio/cpairwise2module.c:60:20: error: incompatible integer to pointer conversion assigning to 'PyObject *' (aka 'struct _object *') from 'int' [-Wint-conversion]
          if(!(py_result = PyEval_CallObject(py_match_fn, py_arglist)))
                         ^ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      2 errors generated.
      error: command '/usr/bin/clang' failed with exit code 1
      [end of output]
  
  note: This error originates from a subprocess, and is likely not a problem with pip.
  ERROR: Failed building wheel for biopython
  building 'Bio.cpairwise2' extension
      clang -fno-strict-overflow -Wsign-compare -Wunreachable-code -fno-common -dynamic -DNDEBUG -g -O3 -Wall -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk -I/Users/dlvovs1/envs/homolig/include -I/opt/homebrew/opt/python@3.13/Frameworks/Python.framework/Versions/3.13/include/python3.13 -c Bio/cpairwise2module.c -o build/temp.macosx-14.0-arm64-cpython-313/Bio/cpairwise2module.o
      Bio/cpairwise2module.c:60:22: error: call to undeclared function 'PyEval_CallObject'; ISO C99 and later do not support implicit function declarations [-Wimplicit-function-declaration]
          if(!(py_result = PyEval_CallObject(py_match_fn, py_arglist)))
                           ^
      Bio/cpairwise2module.c:60:20: error: incompatible integer to pointer conversion assigning to 'PyObject *' (aka 'struct _object *') from 'int' [-Wint-conversion]
          if(!(py_result = PyEval_CallObject(py_match_fn, py_arglist)))
                         ^ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      2 errors generated.
      error: command '/usr/bin/clang' failed with exit code 1
      [end of output]
  
  note: This error originates from a subprocess, and is likely not a problem with pip.
  ERROR: Failed building wheel for biopython
  Failed to build homolig biopython
  ERROR: ERROR: Failed to build installable wheels for some pyproject.toml based projects (homolig, biopython)
  ```
- remove venv, add requirements.txt (I think its's cleaner to gen this from inside the docker)
- scanpy and scikit-learn are missing from dependencies -> add to setup.py

# running
- fix path that causes 
```
root@39dd7c8e5b11:/homolig# bash testcode-python-clustering-module.sh 
python3: can't open file '/homolig/./homolig/homolig_format_checker.py ': [Errno 2] No such file or directory
testcode-python-clustering-module.sh: line 24: -i: command not found
```
- fix anndata futurewarnings
```
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_csv from `anndata` is deprecated. Import anndata.io.read_csv instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_excel from `anndata` is deprecated. Import anndata.io.read_excel instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_hdf from `anndata` is deprecated. Import anndata.io.read_hdf instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_loom from `anndata` is deprecated. Import anndata.io.read_loom instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_mtx from `anndata` is deprecated. Import anndata.io.read_mtx instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_text from `anndata` is deprecated. Import anndata.io.read_text instead.
  warnings.warn(msg, FutureWarning)
/usr/local/lib/python3.11/site-packages/anndata-0.11.0rc3-py3.11.egg/anndata/utils.py:429: FutureWarning: Importing read_umi_tools from `anndata` is deprecated. Import anndata.io.read_umi_tools instead.
  warnings.warn(msg, FutureWarning)
```

# documentation

- updating [metatada](https://github.com/FertigLab/Homolig/blob/acc46726d922f11981de9e5944a6c1fd4e54c15d/setup.py#L27) - description is dummy, github link returns 404

- remove step "Activate a pre-made virtual environment with all homolig dependencies pre-installed", venv is very much bound to developer's machine

# tests

- add tests - no `tests/` directory found

# misc
 - cleaning misc files in the `/homolig`:
 ```
 (homolig) dlvovs1@Dmitrijss-MacBook-Pro Homolig % ls homolig | grep cluster
clusterHomolig.py
clusterHomolig.v2.py
clusterHomolig.v3.py
clusterHomolig.v4.py
```
- OPTIONAL simplifying, there is R, cpp, and python in this tool, are all of them really needed?