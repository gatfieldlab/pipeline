from setuptools import setup

setup(
    name='gatlab-pipeline',
    version='0.3.0',
    description='Tools for analysing HTSeq data',
    author='Bulak Arpat, René Dreos',
    author_email='bulak.arpat@gmail.com',
    license='GPLv3',
    packages=['gatlab.pipeline'],
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
    zip_safe=False
)
