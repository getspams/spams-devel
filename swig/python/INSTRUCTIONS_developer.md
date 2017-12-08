# INSTRUCTIONS for developers

This files describes the different bash scripts that can be used to test and
package the Python SPAMS library.


### Test package building and installation

```bash
./package_build_test_XX.sh
./package_install_test_XX.sh
```

### Test run

```bash
./package_run_test_XX.sh
```

### Create sources

```bash
./make_XX_package.sh
```

### uploard on PyPI

```bash
./make_pypi_package.#!/bin/sh

twine upload -r <repos> dist/*.tar.gz
```

`<respos>` can be any PyPI repository. A simple trick is to define a
'~/.pypirc' file in your home with the following content:

```
[distutils]
index-servers=
        testpypi

[testpypi]
repository: https://test.pypi.org/legacy/
username: <your_username>
password: <your_password>
```
