from setuptools import setup, find_packages

setup(name='decontam_pipeline',
      version='0.1',
      description="""The pipeline used for decontamination of the megadaph 
                     Illumina nuclear assemblies"""
      url='http://github.com/fennerm/megadaph/decontam_pipeline/',
      author='Fenner Macrae',
      author_email='fennermacrae@gmail.com',
      license='MIT',
      packages=find_packages(exclude=["*.test", "*.test.*", "test.*", 
                                      "test"]),
      zip_safe=False)
