# -*- coding: utf-8 -*-

"""
Pip setup for `pipeline' package
"""

import os
from subprocess import check_call
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install

CURDIR = os.getcwd() # this gets the src directory of python package
USERDIR = os.path.expanduser("~")

class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        config_script = os.path.join(self.script_dir,
                                     "pipeline_conf-riboprof.sh")
        develop.run(self)
        check_call(["cp", config_script, USERDIR])


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        config_script = os.path.join(self.install_scripts,
                                     "pipeline_conf-riboprof.sh")
        install.run(self)
        check_call(["cp", config_script, CURDIR])


setup(
    name='pipeline',
    version='0.1.0',
    description='Tools for analysing HTSeq data',
#   url='NA',
    author='Bulak Arpat, René Dreos',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=find_packages(),
    scripts=[
        'scripts/meta_pipeline.sh',
        'scripts/pipeline_bwt2_single.sh',
        'scripts/pipeline_common.sh',
        'scripts/subpipeline.sh',
        'scripts/pipeline_conf-riboprof.sh',
        'scripts/pipeline_conf-riboprof_umi.sh',
        'scripts/pipeline_whitelist.sh',
        'scripts/pipeline_setconf.sh'
    ],
    entry_points={
        'console_scripts': [
            'make_sample_db = pipeline.make_sample_db:main',
            'concat_map_logs = pipeline.concat_map_logs:main',
            'consume = pipeline.consume:main',
            'filter_sam = pipeline.filter_sam:main',
            'filterUmiFromSam = pipeline.filterUmiFromSam:main',
            'split_barcode = pipeline.split_barcode:main'
        ]},
    install_requires=[
        'xopen',
        'gatlab-tools-accessories==0.1.0'
    ],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand},
    zip_safe=False
)
