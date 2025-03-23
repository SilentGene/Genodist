from setuptools import setup, find_packages

setup(
    name="genodist",
    version="0.2.0",
    description="Genomic Distance Analysis Tool",
    license='GPL-3.0',
    author="Heyu Lin",
    author_email="heyu.lin@qut.edu.au",
    packages=find_packages(),  # Automatically find all packages in the project
    py_modules=["genodist"],
    install_requires=[
        "ruamel.yaml>=0.15.99",
        "snakemake>=8.0.0",
    ],
    entry_points={ 
        'console_scripts': [
            'genodist=genodist.genodist:main',  # 命令名称 = 包名.模块名:函数名
        ],
    },
    classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
    python_requires=">=3.8",
)
