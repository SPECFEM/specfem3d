from distutils.core import setup
import glob

setup(name='GEOCUBIT',
      version='3',
      description='CUBIT plugin',
      author='emanuele casarotti',
      author_email='emanuele.casarotti@ingv.it',
      license='GNU General Public License v3 or later',
      url='',
      packages=['','geocubitlib'],
      scripts=['GEOCUBIT.py']
)

