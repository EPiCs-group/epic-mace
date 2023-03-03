import setuptools

with open('README.md', 'r', encoding = 'utf-8') as inpf:
    long_description = inpf.read()

setuptools.setup(
      name = 'epic_mace',
      version = '0.5.0',
      description = 'Python library and command-line tool for generation of 3D coordinates for complexes of d-/f-elements',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      license = 'GPLv3+',
      classifiers = [
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Chemistry',
      ],
      keywords = 'stereomers 3D embedding geometry complex',
      url = 'http://github.com/EPiCs-group/mace',
      author = 'Ivan Yu. Chernyshov',
      author_email = 'ivan.chernyshoff@gmail.com',
      packages = setuptools.find_packages(),
      entry_points = {
          'console_scripts': [
              'epic-mace = mace.__main__:main',
              'epic-mace-quickstart = mace._cli_quickstart:main'
          ]
      },
      install_requires = [
          'numpy',
          'pyyaml'
      ],
      python_requires = '>=3.7'
)
