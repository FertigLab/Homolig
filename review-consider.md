# installation

- consider using [pyproject.toml](https://packaging.python.org/en/latest/guides/writing-pyproject-toml/) instead of setup.py?
- fixing installation that fails with 
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

# documentation

- updating [metatada](https://github.com/FertigLab/Homolig/blob/acc46726d922f11981de9e5944a6c1fd4e54c15d/setup.py#L27) - description is dummy, github link returns 404

- stating "Activate a pre-made virtual environment with all homolig dependencies pre-installed" recommended? or removing at all

# tests

- adding tests as no `tests/` directory found 

# misc
 - cleaning misc files in the `/homolig`:
 ```
 (homolig) dlvovs1@Dmitrijss-MacBook-Pro Homolig % ls homolig | grep cluster
clusterHomolig.py
clusterHomolig.v2.py
clusterHomolig.v3.py
clusterHomolig.v4.py
```
- simplifying, there is R, cpp, and python in this tool, are all of them really needed?