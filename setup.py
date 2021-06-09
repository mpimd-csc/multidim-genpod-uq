from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(name='gen_pod_uq',
      version='0.0.1',
      description='Core module for genpod for UQ',
      # license="GPLv3",
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Jan Heiland',
      author_email='jnhlnd@gmail.com',
      install_requires=['numpy', 'scipy',
                        'dolfin-navier-scipy>=1.0.0',
                        'sadptprj-riclyap-adi>=1.0.0',
                        'multidim-galerkin-pod>=1.0.2',
                        ],
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          "Operating System :: OS Independent",
          ]
      )
