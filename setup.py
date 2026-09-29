from setuptools import setup, Extension

module = Extension(
    'symnmfmodule',
    sources=[
        'symnmfmodule.c',  
        'symnmf.c',  
        'matrix.c'      
    ]
)

setup(
    name='symnmfmodule',
    version='1.0',
    description='Python C-extension for SymNMF',
    ext_modules=[module]
)
