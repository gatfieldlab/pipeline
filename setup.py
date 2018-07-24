# -*- coding: utf-8 -*-

"""
Pip setup for process_reads package
"""

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        # PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION
        print("Do you see me?")
        develop.run(self)


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        # PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION
        install.run(self)

setup(
    name='pipeline',
    version='0.1.0',
    description='Tools for analysing RPF-Seq data',
#   url='NA',
    author='Bulak Arpat',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=find_packages(),
    scripts=['scripts/meta_pipeline.sh',
             'scripts/pipeline_bwt2_single.sh',
             'scripts/pipeline_common.sh',
             'scripts/subpipeline.sh'],
    entry_points={
        'console_scripts': ['make_sample_db = pipeline.make_sample_db:main',
                            'concat_map_logs = pipeline.concat_map_logs:main']},
    install_requires=[],
    dependency_links=[],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand,
    },
    zip_safe=False
)
