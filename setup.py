# -*- coding: utf-8 -*-

"""
Pip setup for process_reads package
"""

import os
from subprocess import check_call
from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install


class PostDevelopCommand(develop):
    """Post-installation for development mode."""
    def run(self):
        config_script = os.path.join(self.script_dir, "pipeline_conf-riboprof.sh")
        develop.run(self)
        check_call(["cp", config_script, "."])


class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        config_script = os.path.join(self.script_dir, "pipeline_conf-riboprof.sh")
        install.run(self)
        check_call(["cp", config_script, "."])


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
             'scripts/subpipeline.sh',
             'scripts/pipeline_conf-riboprof.sh'],
    entry_points={
        'console_scripts': ['make_sample_db = pipeline.make_sample_db:main',
                            'concat_map_logs = pipeline.concat_map_logs:main']},
    install_requires=[],
    dependency_links=[],
    cmdclass={
        'develop': PostDevelopCommand,
        'install': PostInstallCommand},
    zip_safe=False
)
