# -*- coding: utf-8 -*-

"""
Pip setup for process_reads package
"""

import os
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install
from distutils.command.build_py import build_py
from subprocess import check_call


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        develop.run(self)
        check_call(["echo", "'{}'".format(self.install_dir), ">test.txt"])

class BuildConfiguresScriptDir(build_py):
    def run(self):
        build_py.run(self)
        if self.dry_run:
            return
        target = os.path.join(self.build_lib, 'scripts', 'pipeline_conf-riboprof.sh')
        check_call(["echo", "'{}'".format(target), ">test.txt"])

class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        check_call(["touch", "SomeNewFile"])
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
        'build_py': BuildConfiguresScriptDir
    },
    zip_safe=False
)
