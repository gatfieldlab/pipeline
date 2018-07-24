# -*- coding: utf-8 -*-

"""
Pip setup for process_reads package
"""

from setuptools import setup, find_packages
from setuptools.command.develop import develop
from setuptools.command.install import install



def _post_install(my_setup):
    def _post_actions():
        print("THIS WORKED")
    _post_actions()
    return my_setup

_post_install(
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
        zip_safe=False
    )
)
