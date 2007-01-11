
from archimedes import use_merlin
use_merlin()

from merlin import setup
#import os

setup(
    name = 'Specfem3DBasin', 
    version = '1.4',
    url = 'http://www.gps.caltech.edu/~jtromp/research/downloads.html',
    author = 'Dimitri Komatitsch and Jeroen Tromp',
    author_email = 'jtromp AT caltech.edu',
    packages = [ 'Specfem3DBasin' ],
    
    install_requires = [
    'cig >= 1.0dev-r4449, < 2.0a, == dev',
    'pythia[mpi] >= 0.8.1.3, < 0.8.2a',
    ],
    
    dependency_links = [
    'http://geodynamics.org/svn/cig/cs/framework/trunk#egg=cig-dev',
    'http://geodynamics.org/svn/cig/cs/pythia/trunk#egg=pythia-0.8.1.3', # temporary
    ],

    #interpreter = os.path.join(os.getcwd(), "pyspecfem3D"),
    entry_points = {
    'console_scripts': [
    'xspecfem3D = Specfem3DBasin.Specfem:main',
    #'xcreate_movie_AVS_DX = Specfem3DBasin.Specfem:create_movie_AVS_DX',
    ],
    },
)
