#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)  # developmental hack, remove later!
import introSpect

__name__, __package__, invoked_directly = introSpect.cmdSupport(
    __name__, __package__, __file__
)
hint = introSpect.hint

latex_container = "docker://blang/latex:ctanfull"

environment = """name: survivor
channels:
  - brown-data-science
  - anaconda
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=0_gnu
  - attrs=19.3.0=py_0
  - autograd=1.3=py_0
  - autograd-gamma=0.4.1=py_0
  - backcall=0.1.0=py_0
  - beautifulsoup4=4.8.2=py38_0
  - bleach=3.1.1=py_0
  - ca-certificates=2020.1.1=0
  - cachetools=3.1.1=py_0
  - cairo=1.16.0=hcf35c78_1003
  - certifi=2020.4.5.1=py38_0
  - cffi=1.14.0=py38hd463f26_0
  - chardet=3.0.4=py38_1003
  - cryptography=2.8=py38h72c5cf5_1
  - curl=7.68.0=hf8cf82a_0
  - cycler=0.10.0=py_2
  - decorator=4.4.2=py_0
  - defusedxml=0.6.0=py_0
  - entrypoints=0.3=py38_1000
  - expat=2.2.9=he1b5a44_2
  - fontconfig=2.13.1=h86ecdb6_1001
  - freetype=2.10.0=he983fc9_1
  - future=0.18.2=py38_0
  - gettext=0.19.8.1=hc5be6a0_1002
  - git=2.26.0=pl526hf241897_0
  - glib=2.58.3=py38h6f030ca_1002
  - gmp=6.2.0=he1b5a44_2
  - google-api-core=1.16.0=py38_1
  - google-api-python-client=1.8.0=pyh8c360ce_0
  - google-auth=1.11.2=py_0
  - google-auth-httplib2=0.0.3=py_3
  - googleapis-common-protos=1.51.0=py38_1
  - graphite2=1.3.13=he1b5a44_1001
  - harfbuzz=2.4.0=h9f30f68_3
  - httplib2=0.17.0=py38_0
  - icu=64.2=he1b5a44_1
  - idna=2.9=py_1
  - importlib_metadata=1.5.0=py38_0
  - ipykernel=5.1.4=py38h5ca1d4c_0
  - ipython=7.13.0=py38h5ca1d4c_0
  - ipython_genutils=0.2.0=py_1
  - jedi=0.16.0=py38_0
  - jinja2=2.11.1=py_0
  - jpeg=9c=h14c3975_1001
  - json5=0.9.0=py_0
  - jsonschema=3.2.0=py38_0
  - jupyter_client=6.0.0=py_0
  - jupyter_core=4.6.3=py38_0
  - jupyterlab=2.0.0=py_1
  - jupyterlab_server=1.0.6=py_0
  - jupytext=1.4.2=pyh9f0ad1d_0
  - kiwisolver=1.1.0=py38hc9558a2_0
  - krb5=1.16.4=h2fd8d38_0
  - ld_impl_linux-64=2.33.1=h53a641e_8
  - libblas=3.8.0=14_openblas
  - libcblas=3.8.0=14_openblas
  - libcurl=7.68.0=hda55be3_0
  - libedit=3.1.20170329=hf8c457e_1001
  - libffi=3.2.1=he1b5a44_1006
  - libgcc-ng=9.2.0=h24d8f2e_2
  - libgfortran-ng=7.3.0=hdf63c60_5
  - libgomp=9.2.0=h24d8f2e_2
  - libiconv=1.15=h516909a_1006
  - liblapack=3.8.0=14_openblas
  - libopenblas=0.3.7=h5ec1e0e_6
  - libpng=1.6.37=hed695b0_0
  - libprotobuf=3.11.4=h8b12597_0
  - libsodium=1.0.17=h516909a_0
  - libssh2=1.8.2=h22169c7_2
  - libstdcxx-ng=9.2.0=hdf63c60_2
  - libtiff=4.1.0=hc7e4089_6
  - libuuid=2.32.1=h14c3975_1000
  - libwebp-base=1.1.0=h516909a_3
  - libxcb=1.13=h14c3975_1002
  - libxml2=2.9.10=hee79883_0
  - lifelines=0.23.9=py_0
  - lz4-c=1.8.3=he1b5a44_1001
  - markupsafe=1.1.1=py38h516909a_0
  - matplotlib-base=3.1.3=py38h250f245_0
  - mistune=0.8.4=py38h516909a_1000
  - mpfr=4.0.2=he80fd80_0
  - nbconvert=5.6.1=py38_0
  - nbformat=5.0.4=py_0
  - ncurses=6.1=hf484d3e_1002
  - notebook=6.0.3=py38_0
  - numpy=1.18.1=py38h95a1406_0
  - oauth2client=4.1.3=py_0
  - openjpeg=2.3.1=h981e76c_3
  - openssl=1.1.1g=h7b6447c_0
  - pandas=1.0.1=py38hb3f55d8_0
  - pandoc=2.9.2=0
  - pandocfilters=1.4.2=py_1
  - parso=0.6.1=py_0
  - patsy=0.5.1=py_0
  - pcre=8.44=he1b5a44_0
  - perl=5.26.2=h516909a_1006
  - pexpect=4.8.0=py38_0
  - pickleshare=0.7.5=py38_1000
  - pip=20.0.2=py_2
  - pixman=0.38.0=h516909a_1003
  - pkg-config=0.29.2=h516909a_1006
  - poppler=0.65.0=h14e79db_0
  - poppler-data=0.4.9=1
  - predictmd-texlive=20170524=0
  - prometheus_client=0.7.1=py_0
  - prompt_toolkit=3.0.3=py_0
  - protobuf=3.11.4=py38he1b5a44_0
  - pthread-stubs=0.4=h14c3975_1001
  - ptyprocess=0.6.0=py_1001
  - pyasn1=0.4.8=py_0
  - pyasn1-modules=0.2.7=py_0
  - pycparser=2.20=py_0
  - pydrive=1.3.1=py_1
  - pygments=2.5.2=py_0
  - pyopenssl=19.1.0=py_1
  - pyparsing=2.4.6=py_0
  - pyrsistent=0.15.7=py38h516909a_0
  - pysocks=1.7.1=py38_0
  - python=3.8.2=h9d8adfe_1_cpython
  - python-dateutil=2.8.1=py_0
  - pytz=2019.3=py_0
  - pyyaml=5.1.2=py38h516909a_1
  - pyzmq=19.0.0=py38h1768529_0
  - readline=8.0=hf8c457e_0
  - requests=2.23.0=pyh8c360ce_2
  - rsa=4.0=py_0
  - scipy=1.4.1=py38h921218d_0
  - seaborn=0.10.0=py_1
  - send2trash=1.5.0=py_0
  - setuptools=45.2.0=py38_0
  - simplejson=3.17.0=py38h516909a_0
  - six=1.14.0=py38_0
  - soupsieve=2.0.1=py_0
  - sqlite=3.30.1=hcee41ef_0
  - statsmodels=0.11.1=py38h516909a_0
  - terminado=0.8.3=py38_0
  - testpath=0.4.4=py_0
  - texlive-core=20180414=pl526h89d1741_1
  - tikzmagic=0.1.1=py_0
  - tk=8.6.10=hed695b0_0
  - tornado=6.0.3=py38h516909a_4
  - traitlets=4.3.3=py38_0
  - uritemplate=3.0.1=py_0
  - urllib3=1.25.7=py38_0
  - wcwidth=0.1.8=py_0
  - webencodings=0.5.1=py_1
  - wheel=0.34.2=py_1
  - xorg-kbproto=1.0.7=h14c3975_1002
  - xorg-libice=1.0.10=h516909a_0
  - xorg-libsm=1.2.3=h84519dc_1000
  - xorg-libx11=1.6.9=h516909a_0
  - xorg-libxau=1.0.9=h14c3975_0
  - xorg-libxdmcp=1.1.3=h516909a_0
  - xorg-libxext=1.3.4=h516909a_0
  - xorg-libxrender=0.9.10=h516909a_1002
  - xorg-renderproto=0.11.1=h14c3975_1002
  - xorg-xextproto=7.3.0=h14c3975_1002
  - xorg-xproto=7.0.31=h14c3975_1007
  - xz=5.2.4=h14c3975_1001
  - yaml=0.2.2=h516909a_1
  - zeromq=4.3.2=he1b5a44_2
  - zipp=3.0.0=py_0
  - zlib=1.2.11=h516909a_1006
  - zstd=1.4.4=h3b9ef0a_2
  - pip:
    - xenapython==1.0.10
"""


def recipe(*, verbose: bool = False,) -> dict:

    """
    Create a yaml file with all dependencies of the package that can be used as template
    for conda.

    Parameters
    ----------
    verbose
        Flag to turn on messages.
    
    Returns
    -------
    The yaml of dependencies.
    
    """

    ### Just retrieve the settings
    hint(verbose, environment)

    return environment


def main():

    """
    Define what should be executed when invoked from command line.
    """
    modified_kws = {
        "verbose": (
            0,
            "-v",
            "--verbose",
            {"dest": "verbose", "default": False, "action": "store_true"},
        ),
        "outFile": (
            1,
            "-o",
            "--outFile",
            {
                "dest": "outFile",
                "help": "Location where results should be saved. If not specified, STDOUT will be used.",
            },
        ),
    }
    mainFunction = introSpect.commandLines.cmdConnect(recipe, modified_kws)
    mainFunction.eval()
    mainFunction.save()
    return


__doc__ = recipe.__doc__
if invoked_directly:
    main()
